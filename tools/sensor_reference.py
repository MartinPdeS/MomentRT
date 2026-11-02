"""Independent weighted Monte Carlo reference for layered-slab contact sensors."""

from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Any
import numpy as np
from numpy.typing import NDArray

REFERENCE_REVISION = 2


@dataclass(frozen=True)
class LayerOptics:
    """Monochromatic optical properties for one homogeneous slab layer.

    Parameters
    ----------
    thickness_mm : float
        Positive layer thickness in millimeters.
    absorption_per_mm : float
        Nonnegative absorption coefficient in inverse millimeters.
    reduced_scattering_per_mm : float
        Nonnegative reduced scattering coefficient in inverse millimeters.
        Zero permits independent transparent and purely absorbing limit checks.
    anisotropy : float
        Henyey-Greenstein anisotropy in [0, 1).
    """

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
    weight continuously through absorption. Free paths use the total scattering
    coefficient, not absorption plus scattering. Paths stop at every interface;
    a new scattering distance in the next layer uses exponential memorylessness.
    Every sample face is index matched
    and transmitting; there is no Fresnel reflection. The output is outward
    power crossing each detector aperture, multiplied by the same explicit
    Lambertian detector acceptance fraction used by MomentRT (one for the full
    hemisphere). The independent trajectories provide standard errors, per-face
    escape and unscattered powers, and absorbed-power and balance diagnostics.

    All lengths and inverse-length coefficients in this standalone numerical
    reference use millimeters. Returned optical powers and standard errors use
    watts. Zero scattering is supported for analytical limit checks.

    Parameters
    ----------
    layers : tuple of LayerOptics
        Ordered homogeneous layers, starting at z=0, spanning the sample height.
    sample_width_mm : float
        Positive width of both lateral sample dimensions in millimeters.
    sample_height_mm : float
        Positive sample depth in millimeters, equal to the layer thickness sum.
    source_radius_mm : float
        Radius of the centered top-face source disk in millimeters.
    detector_offsets_mm : tuple of float
        Top-face aperture x offsets from the source center in millimeters.
        Empty collects only whole-face and absorption diagnostics.
    detector_radius_mm : float
        Common aperture radius in millimeters.
    photon_count : int
        Number of independent histories, at least two.
    seed : int
        Seed passed to NumPy's random generator.
    optical_power_w : float, optional
        Positive incoming diffuse optical power in watts. Default is 0.001.
    batch_size : int, optional
        Maximum histories held in memory per batch. Default is 20000.

    Returns
    -------
    dict of str to object
        Revision, inputs, detector and face powers, and standard errors in watts.
        The balance standard error uses per-history correlated residuals, so it
        also accounts for the unbiased stochastic weight changes from roulette.
    """
    if not layers or any(
        not np.all(
            np.isfinite(
                [
                    layer.thickness_mm,
                    layer.absorption_per_mm,
                    layer.reduced_scattering_per_mm,
                    layer.anisotropy,
                ]
            )
        )
        or layer.thickness_mm <= 0
        or layer.absorption_per_mm < 0
        or layer.reduced_scattering_per_mm < 0
        or not 0 <= layer.anisotropy < 1
        for layer in layers
    ):
        raise ValueError(
            "Layer thicknesses must be positive; coefficients finite and nonnegative; anisotropy in [0, 1)."
        )
    if (
        not np.all(np.isfinite([sample_height_mm, sample_width_mm]))
        or sample_height_mm <= 0
        or sample_width_mm <= 0
    ):
        raise ValueError("Sample dimensions must be positive.")
    if not np.isclose(
        sum(layer.thickness_mm for layer in layers),
        sample_height_mm,
        rtol=8 * np.finfo(float).eps,
        atol=0,
    ):
        raise ValueError("Layer thicknesses must sum to sample_height_mm.")
    if (
        any(
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            for value in (photon_count, batch_size)
        )
        or photon_count < 2
        or batch_size < 1
    ):
        raise ValueError("photon_count must be at least two and batch_size positive.")
    if (
        not np.all(np.isfinite([source_radius_mm, detector_radius_mm]))
        or source_radius_mm <= 0
        or detector_radius_mm <= 0
    ):
        raise ValueError("Source and detector radii must be positive.")
    if not np.isfinite(optical_power_w) or optical_power_w <= 0:
        raise ValueError("optical_power_w must be finite and positive.")
    detector_offsets = np.asarray(detector_offsets_mm, dtype=float)
    if detector_offsets.ndim != 1 or not np.all(np.isfinite(detector_offsets)):
        raise ValueError("Detector offsets must be a finite vector.")
    if source_radius_mm > sample_width_mm / 2 or np.any(
        np.abs(detector_offsets) + detector_radius_mm > sample_width_mm / 2
    ):
        raise ValueError("Source and detector disks must fit on the sample face.")

    rng = np.random.default_rng(seed)
    interfaces = np.cumsum([layer.thickness_mm for layer in layers])
    interfaces[-1] = sample_height_mm
    absorption = np.array([layer.absorption_per_mm for layer in layers])
    scattering = np.array([layer.scattering_per_mm for layer in layers])
    anisotropy = np.array([layer.anisotropy for layer in layers])
    face_names = ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
    detector_count = len(detector_offsets)
    # Scores include detectors, six faces, six unscattered face powers,
    # absorbed power, and the correlated energy-balance residual.
    tallied_sum = np.zeros(detector_count + 14)
    tallied_square_sum = np.zeros(detector_count + 14)
    processed = 0

    while processed < photon_count:
        count = min(batch_size, photon_count - processed)
        positions = np.zeros((count, 3))
        source_radius = source_radius_mm * np.sqrt(rng.random(count))
        source_azimuth = 2 * np.pi * rng.random(count)
        positions[:, 0] = sample_width_mm / 2 + source_radius * np.cos(source_azimuth)
        positions[:, 1] = sample_width_mm / 2 + source_radius * np.sin(source_azimuth)
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
        layer_indices = np.zeros(count, dtype=np.intp)
        scattering_counts = np.zeros(count, dtype=np.intp)
        scores = np.zeros((len(detector_offsets), count))
        face_scores = np.zeros((6, count))
        unscattered_scores = np.zeros((6, count))
        absorbed_scores = np.zeros(count)
        events = 0

        while np.any(active):
            events += 1
            if events > 100_000:
                raise RuntimeError("Monte Carlo packet exceeded the event limit.")
            indices = np.flatnonzero(active)
            position = positions[indices]
            direction = directions[indices]
            layer_index = layer_indices[indices]
            # Continuous absorption weighting samples scattering events at mu_s.
            # Sampling at mu_a + mu_s and always scattering biases transport.
            free_path = np.divide(
                rng.exponential(size=len(indices)),
                scattering[layer_index],
                out=np.full(len(indices), np.inf),
                where=scattering[layer_index] > 0,
            )

            boundary_distance = np.full((len(indices), 3), np.inf)
            for axis in (0, 1):
                positive = direction[:, axis] > 0
                plane = np.where(positive, sample_width_mm, 0.0)
                np.divide(
                    plane - position[:, axis],
                    direction[:, axis],
                    out=boundary_distance[:, axis],
                    where=direction[:, axis] != 0,
                )

            positive_z = direction[:, 2] > 0
            next_interface = np.empty(len(indices))
            next_interface[positive_z] = interfaces[layer_index[positive_z]]
            next_interface[~positive_z] = np.where(
                layer_index[~positive_z] == 0,
                0.0,
                interfaces[layer_index[~positive_z] - 1],
            )
            np.divide(
                next_interface - position[:, 2],
                direction[:, 2],
                out=boundary_distance[:, 2],
                where=direction[:, 2] != 0,
            )
            boundary_distance = np.maximum(boundary_distance, 0)
            nearest_axis = np.argmin(boundary_distance, axis=1)
            distance_to_boundary = boundary_distance[
                np.arange(len(indices)), nearest_axis
            ]
            scatters = free_path < distance_to_boundary
            travel = np.where(scatters, free_path, distance_to_boundary)
            position += direction * travel[:, None]
            positions[indices] = position
            absorption_depth = absorption[layer_index] * travel
            absorbed_scores[indices] += weights[indices] * -np.expm1(-absorption_depth)
            weights[indices] *= np.exp(-absorption_depth)

            scatter_indices = indices[scatters]
            if len(scatter_indices):
                scattering_counts[scatter_indices] += 1
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
                z_axis = axes == 2
                outward = (
                    z_axis & ~positive_z[~scatters] & (layer_index[~scatters] == 0)
                )
                inward = (
                    z_axis
                    & positive_z[~scatters]
                    & (layer_index[~scatters] == len(layers) - 1)
                )
                escaped = outward | inward | (axes != 2)

                escaped_indices = boundary_indices[escaped]
                escaped_axes = axes[escaped]
                escaped_faces = 2 * escaped_axes + (
                    directions[escaped_indices, escaped_axes] > 0
                ).astype(np.intp)
                face_scores[escaped_faces, escaped_indices] = weights[escaped_indices]
                unscattered = scattering_counts[escaped_indices] == 0
                unscattered_scores[
                    escaped_faces[unscattered], escaped_indices[unscattered]
                ] = weights[escaped_indices[unscattered]]

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
                    # Track the side of an exact interface explicitly, avoiding
                    # geometric nudges that omit absorption or skip thin layers.
                    positions[crossed, 2] = next_interface[~scatters][internal]
                    layer_indices[crossed] += np.where(
                        directions[crossed, 2] > 0, 1, -1
                    )

        residual_scores = 1 - absorbed_scores - face_scores.sum(axis=0)
        all_scores = np.vstack(
            (scores, face_scores, unscattered_scores, absorbed_scores, residual_scores)
        )
        tallied_sum += all_scores.sum(axis=1)
        tallied_square_sum += np.square(all_scores).sum(axis=1)
        processed += count

    mean_fraction = tallied_sum / photon_count
    variance = np.maximum(
        0,
        (tallied_square_sum - photon_count * mean_fraction**2) / (photon_count - 1),
    )
    standard_error_fraction = np.sqrt(variance / photon_count)
    power = optical_power_w * mean_fraction
    standard_error = optical_power_w * standard_error_fraction
    detector_power = power[:detector_count]
    detector_standard_error = standard_error[:detector_count]

    return {
        "method": "scattering-only free paths with continuous absorption weighting",
        "reference_revision": REFERENCE_REVISION,
        "photon_count": photon_count,
        "seed": seed,
        "source_power_w": optical_power_w,
        "detector_offsets_mm": detector_offsets_mm,
        "detector_power_w": detector_power.tolist(),
        "standard_error_w": detector_standard_error.tolist(),
        "relative_standard_error": np.divide(
            detector_standard_error,
            detector_power,
            out=np.full_like(detector_power, np.inf),
            where=detector_power > 0,
        ).tolist(),
        "escaping_power_by_face_w": dict(
            zip(face_names, power[detector_count : detector_count + 6].tolist())
        ),
        "escaping_standard_error_by_face_w": dict(
            zip(
                face_names, standard_error[detector_count : detector_count + 6].tolist()
            )
        ),
        "unscattered_power_by_face_w": dict(
            zip(face_names, power[detector_count + 6 : detector_count + 12].tolist())
        ),
        "unscattered_standard_error_by_face_w": dict(
            zip(
                face_names,
                standard_error[detector_count + 6 : detector_count + 12].tolist(),
            )
        ),
        "absorbed_power_w": float(power[-2]),
        "absorbed_standard_error_w": float(standard_error[-2]),
        "energy_balance_residual_w": float(power[-1]),
        "energy_balance_standard_error_w": float(standard_error[-1]),
        "layers": [asdict(layer) for layer in layers],
    }
