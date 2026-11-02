"""Independent analog Monte Carlo reference for layered-slab contact sensors."""

from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Any
import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class LayerOptics:
    """Monochromatic optical properties for one homogeneous slab layer."""

    thickness_mm: float
    absorption_per_mm: float
    reduced_scattering_per_mm: float
    anisotropy: float

    @property
    def scattering_per_mm(self) -> float:
        """Return the total scattering coefficient in inverse millimeters."""
        return self.reduced_scattering_per_mm / (1 - self.anisotropy)


def _scatter_directions(
    *,
    directions: NDArray[np.float64],
    anisotropy: NDArray[np.float64],
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Sample HG scattering directions independently of MomentRT kernels."""
    uniform = rng.random(len(directions))
    cosine = np.empty_like(uniform)
    isotropic = np.abs(anisotropy) < 1e-12
    cosine[isotropic] = 2 * uniform[isotropic] - 1
    g = anisotropy[~isotropic]
    term = (1 - g**2) / (1 - g + 2 * g * uniform[~isotropic])
    cosine[~isotropic] = np.clip((1 + g**2 - term**2) / (2 * g), -1, 1)
    sine = np.sqrt(np.maximum(0, 1 - cosine**2))
    azimuth = 2 * np.pi * rng.random(len(directions))
    reference = np.zeros_like(directions)
    reference[:, 2] = 1
    near_pole = np.abs(directions[:, 2]) > 0.9
    reference[near_pole] = [1, 0, 0]
    first_axis = np.cross(directions, reference)
    first_axis /= np.linalg.norm(first_axis, axis=1)[:, None]
    second_axis = np.cross(directions, first_axis)
    scattered = (
        cosine[:, None] * directions
        + sine[:, None] * np.cos(azimuth)[:, None] * first_axis
        + sine[:, None] * np.sin(azimuth)[:, None] * second_axis
    )
    return scattered / np.linalg.norm(scattered, axis=1)[:, None]


def simulate_layered_slab_reference(
    *,
    layers: tuple[LayerOptics, ...],
    sample_width_mm: float,
    sample_height_mm: float,
    source_radius_mm: float,
    detector_offsets_mm: tuple[float, ...],
    detector_radius_mm: float,
    photon_count: int,
    seed: int,
    optical_power_w: float = 1e-3,
    batch_size: int = 20_000,
) -> dict[str, Any]:
    """Estimate diffuse detector powers by boundary-injected photon transport.

    Photons enter uniformly over a circular source disk with a cosine-weighted
    inward direction, scatter from Henyey-Greenstein phase functions, and lose
    weight continuously through absorption. Every sample face is index matched
    and transmitting; there is no Fresnel reflection. The output is outward
    power crossing each detector aperture, multiplied by the same explicit
    Lambertian detector acceptance fraction used by MomentRT (one for the full
    hemisphere). The independent trajectories provide standard errors.
    """
    if not layers or any(
        layer.thickness_mm <= 0
        or layer.absorption_per_mm < 0
        or layer.reduced_scattering_per_mm <= 0
        or not 0 <= layer.anisotropy < 1
        for layer in layers
    ):
        raise ValueError(
            "Layer thicknesses and scattering must be positive; optics valid."
        )
    if sample_height_mm <= 0 or sample_width_mm <= 0:
        raise ValueError("Sample dimensions must be positive.")
    if sum(layer.thickness_mm for layer in layers) != sample_height_mm:
        raise ValueError("Layer thicknesses must sum exactly to sample_height_mm.")
    if photon_count < 2 or batch_size < 1:
        raise ValueError("photon_count must be at least two and batch_size positive.")
    if source_radius_mm <= 0 or detector_radius_mm <= 0:
        raise ValueError("Source and detector radii must be positive.")

    rng = np.random.default_rng(seed)
    interfaces = np.cumsum([layer.thickness_mm for layer in layers])
    absorption = np.array([layer.absorption_per_mm for layer in layers])
    scattering = np.array([layer.scattering_per_mm for layer in layers])
    anisotropy = np.array([layer.anisotropy for layer in layers])
    detector_offsets = np.asarray(detector_offsets_mm, dtype=float)
    tallied_sum = np.zeros(len(detector_offsets))
    tallied_square_sum = np.zeros(len(detector_offsets))
    processed = 0
    boundary_epsilon = min(sample_width_mm, sample_height_mm) * 1e-12

    while processed < photon_count:
        count = min(batch_size, photon_count - processed)
        positions = np.zeros((count, 3))
        source_radius = source_radius_mm * np.sqrt(rng.random(count))
        source_azimuth = 2 * np.pi * rng.random(count)
        positions[:, 0] = sample_width_mm / 2 + source_radius * np.cos(source_azimuth)
        positions[:, 1] = sample_width_mm / 2 + source_radius * np.sin(source_azimuth)
        positions[:, 2] = boundary_epsilon
        cosine_in = np.sqrt(rng.random(count))
        azimuth = 2 * np.pi * rng.random(count)
        directions = np.column_stack(
            (
                np.sqrt(1 - cosine_in**2) * np.cos(azimuth),
                np.sqrt(1 - cosine_in**2) * np.sin(azimuth),
                cosine_in,
            )
        )
        weights = np.ones(count)
        active = np.ones(count, dtype=bool)
        scores = np.zeros((len(detector_offsets), count))
        events = 0

        while np.any(active):
            events += 1
            if events > 100_000:
                raise RuntimeError("Monte Carlo packet exceeded the event limit.")
            indices = np.flatnonzero(active)
            position = positions[indices]
            direction = directions[indices]
            layer_index = np.searchsorted(interfaces, position[:, 2], side="right")
            layer_index = np.minimum(layer_index, len(layers) - 1)
            total = absorption[layer_index] + scattering[layer_index]
            free_path = -np.log(np.maximum(rng.random(len(indices)), 1e-15)) / total

            boundary_distance = np.full((len(indices), 3), np.inf)
            for axis in (0, 1):
                positive = direction[:, axis] > 0
                plane = np.where(positive, sample_width_mm, 0.0)
                boundary_distance[:, axis] = (plane - position[:, axis]) / direction[
                    :, axis
                ]

            positive_z = direction[:, 2] > 0
            next_interface = np.empty(len(indices))
            next_interface[positive_z] = interfaces[
                np.minimum(layer_index[positive_z] + 1, len(interfaces) - 1)
            ]
            next_interface[~positive_z] = np.where(
                layer_index[~positive_z] == 0,
                0.0,
                interfaces[layer_index[~positive_z] - 1],
            )
            boundary_distance[:, 2] = (next_interface - position[:, 2]) / direction[
                :, 2
            ]
            boundary_distance = np.maximum(boundary_distance, 0)
            nearest_axis = np.argmin(boundary_distance, axis=1)
            distance_to_boundary = boundary_distance[
                np.arange(len(indices)), nearest_axis
            ]
            scatters = free_path < distance_to_boundary
            travel = np.where(scatters, free_path, distance_to_boundary)
            position += direction * travel[:, None]
            positions[indices] = position
            weights[indices] *= np.exp(-absorption[layer_index] * travel)

            scatter_indices = indices[scatters]
            if len(scatter_indices):
                directions[scatter_indices] = _scatter_directions(
                    directions=direction[scatters],
                    anisotropy=anisotropy[layer_index[scatters]],
                    rng=rng,
                )
                roulette = weights[scatter_indices] < 1e-3
                if np.any(roulette):
                    roulette_indices = scatter_indices[roulette]
                    survive = (
                        rng.random(len(roulette_indices))
                        < weights[roulette_indices] / 1e-3
                    )
                    weights[roulette_indices[survive]] = 1e-3
                    active[roulette_indices[~survive]] = False

            boundary_indices = indices[~scatters]
            if len(boundary_indices):
                axes = nearest_axis[~scatters]
                outward = np.zeros(len(boundary_indices), dtype=bool)
                z_axis = axes == 2
                z_positions = position[~scatters, 2]
                outward[z_axis] = z_positions[z_axis] <= boundary_epsilon * 2
                inward = z_axis & (
                    z_positions >= sample_height_mm - boundary_epsilon * 2
                )
                escaped = outward | inward | (axes != 2)

                top_escape_local = np.flatnonzero(outward)
                if len(top_escape_local):
                    escaped_positions = position[~scatters][top_escape_local]
                    dx = (
                        escaped_positions[:, 0, None]
                        - sample_width_mm / 2
                        - detector_offsets
                    )
                    dy = escaped_positions[:, 1, None] - sample_width_mm / 2
                    hit = dx**2 + dy**2 <= detector_radius_mm**2
                    packet_scores = (
                        weights[boundary_indices[top_escape_local], None] * hit
                    )
                    scores[:, boundary_indices[top_escape_local]] = packet_scores.T

                active[boundary_indices[escaped]] = False
                internal = np.flatnonzero(z_axis & ~escaped)
                if len(internal):
                    crossed = boundary_indices[internal]
                    positions[crossed, 2] += directions[crossed, 2] * boundary_epsilon

        tallied_sum += scores.sum(axis=1)
        tallied_square_sum += np.square(scores).sum(axis=1)
        processed += count

    mean_fraction = tallied_sum / photon_count
    variance = np.maximum(
        0,
        (tallied_square_sum - photon_count * mean_fraction**2) / (photon_count - 1),
    )
    standard_error_fraction = np.sqrt(variance / photon_count)
    power = optical_power_w * mean_fraction
    standard_error = optical_power_w * standard_error_fraction

    return {
        "method": "analog photon transport with continuous absorption weighting",
        "photon_count": photon_count,
        "seed": seed,
        "source_power_w": optical_power_w,
        "detector_offsets_mm": detector_offsets_mm,
        "detector_power_w": power.tolist(),
        "standard_error_w": standard_error.tolist(),
        "relative_standard_error": np.divide(
            standard_error,
            power,
            out=np.full_like(power, np.inf),
            where=power > 0,
        ).tolist(),
        "layers": [asdict(layer) for layer in layers],
    }
