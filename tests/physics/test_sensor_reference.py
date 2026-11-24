"""Analytical and independently discretized RTE checks for the photon reference."""

from __future__ import annotations
from typing import Any
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm
from scipy.special import expn
from tools.sensor_reference import (
    LayerOptics,
    _scatter_directions,
    simulate_layered_slab_reference,
)


def plane_slab(*, layers: tuple[LayerOptics, ...], seed: int = 1729) -> dict[str, Any]:
    """Make lateral boundaries negligible for plane-slab analytical comparisons."""
    return simulate_layered_slab_reference(
        layers=layers,
        sample_width_mm=1e6,
        sample_height_mm=sum(layer.thickness_mm for layer in layers),
        source_radius_mm=0.1,
        detector_offsets_mm=(),
        detector_radius_mm=0.1,
        photon_count=50_000,
        seed=seed,
        optical_power_w=1.0,
    )


@pytest.mark.parametrize("top_thickness_mm", [0.25, 1e-10])
def test_absorbing_layers_match_cosine_weighted_beer_lambert(top_thickness_mm):
    """Every layer contributes its own optical depth, including very thin ones."""
    layers = (
        LayerOptics(
            thickness_mm=top_thickness_mm,
            absorption_per_mm=0.75 / top_thickness_mm,
            reduced_scattering_per_mm=0,
            anisotropy=0,
        ),
        LayerOptics(
            thickness_mm=1 - top_thickness_mm,
            absorption_per_mm=0.1 / (1 - top_thickness_mm),
            reduced_scattering_per_mm=0,
            anisotropy=0,
        ),
    )
    result = plane_slab(layers=layers)
    # Diffuse incidence has p(mu)=2mu. Integrating exp(-tau/mu) gives
    # 2 E_3(tau), independent of how the optical depth is split into layers.
    expected = 2 * expn(3, 0.85)
    transmission = result["escaping_power_by_face_w"]["z_max"]
    standard_error = result["escaping_standard_error_by_face_w"]["z_max"]
    assert transmission == pytest.approx(expected, abs=5 * standard_error)
    assert result["unscattered_power_by_face_w"]["z_max"] == transmission
    assert result["escaping_power_by_face_w"]["z_min"] == 0
    assert result["absorbed_power_w"] == pytest.approx(1 - transmission, abs=1e-12)
    assert abs(result["energy_balance_residual_w"]) < 1e-12


def test_ballistic_power_distinguishes_scattering_from_extinction_sampling():
    """Unscattered transmission includes scattering survival and absorption once."""
    result = plane_slab(
        layers=(
            LayerOptics(
                thickness_mm=1,
                absorption_per_mm=0.7,
                reduced_scattering_per_mm=0.3,
                anisotropy=0,
            ),
        ),
    )
    expected = 2 * expn(3, 1.0)
    measured = result["unscattered_power_by_face_w"]["z_max"]
    standard_error = result["unscattered_standard_error_by_face_w"]["z_max"]
    assert measured == pytest.approx(expected, abs=5 * standard_error)
    biased_extinction_sampling = 2 * expn(3, 1.7)
    assert abs(measured - biased_extinction_sampling) > 10 * standard_error


def isotropic_rte(*, layers: tuple[LayerOptics, ...], order: int) -> np.ndarray:
    """Solve the independent plane RTE by angular quadrature and matrix exponentials.

    Vacuum boundaries have unit inward radiance at z=0 and zero inward
    radiance at the bottom. Gauss-Legendre ordinates discretize
    mu dI/dz = -(mu_a+mu_s) I + mu_s/2 integral(I dmu).
    Small optical thickness keeps the transfer matrices well conditioned.
    """
    direction_cosines, angular_weights = leggauss(order)
    transfer = np.eye(order)
    for layer in layers:
        assert layer.anisotropy == 0
        scattering = layer.reduced_scattering_per_mm
        extinction = layer.absorption_per_mm + scattering
        generator = (
            0.5 * scattering * angular_weights[None, :] - extinction * np.eye(order)
        ) / direction_cosines[:, None]
        transfer = expm(generator * layer.thickness_mm) @ transfer
    negative = np.flatnonzero(direction_cosines < 0)
    positive = np.flatnonzero(direction_cosines > 0)
    incoming = np.ones(len(positive))
    outgoing_top = np.linalg.solve(
        transfer[np.ix_(negative, negative)],
        -transfer[np.ix_(negative, positive)] @ incoming,
    )
    bottom = transfer @ np.concatenate((outgoing_top, incoming))
    incoming_flux = np.sum(angular_weights[positive] * direction_cosines[positive])
    reflected = (
        -np.sum(angular_weights[negative] * direction_cosines[negative] * outgoing_top)
        / incoming_flux
    )
    transmitted = (
        np.sum(
            angular_weights[positive] * direction_cosines[positive] * bottom[positive]
        )
        / incoming_flux
    )
    return np.array([reflected, transmitted])


def test_layered_photon_transport_matches_independent_rte():
    """Compare reflected and transmitted powers with a refined angular RTE solve."""
    layers = (
        LayerOptics(
            thickness_mm=0.3,
            absorption_per_mm=0.2,
            reduced_scattering_per_mm=0.3,
            anisotropy=0,
        ),
        LayerOptics(
            thickness_mm=0.7,
            absorption_per_mm=0.02,
            reduced_scattering_per_mm=0.15,
            anisotropy=0,
        ),
    )
    coarse = isotropic_rte(layers=layers, order=32)
    fine = isotropic_rte(layers=layers, order=64)
    quadrature_change = np.abs(fine - coarse)
    assert np.all(quadrature_change < 5e-4)
    result = plane_slab(layers=layers)
    for index, face in enumerate(("z_min", "z_max")):
        standard_error = result["escaping_standard_error_by_face_w"][face]
        assert result["escaping_power_by_face_w"][face] == pytest.approx(
            fine[index], abs=5 * standard_error + quadrature_change[index]
        )


def test_inserting_identical_material_interfaces_preserves_transport():
    """Artificial layer boundaries do not change homogeneous reflected power."""
    properties = {
        "absorption_per_mm": 0.1,
        "reduced_scattering_per_mm": 0.5,
        "anisotropy": 0.6,
    }
    homogeneous = plane_slab(layers=(LayerOptics(thickness_mm=1, **properties),))
    partitioned = plane_slab(
        layers=tuple(
            LayerOptics(thickness_mm=thickness, **properties)
            for thickness in (0.1, 0.3, 0.6)
        ),
        seed=1730,
    )
    for face in ("z_min", "z_max"):
        uncertainty = np.hypot(
            homogeneous["escaping_standard_error_by_face_w"][face],
            partitioned["escaping_standard_error_by_face_w"][face],
        )
        assert homogeneous["escaping_power_by_face_w"][face] == pytest.approx(
            partitioned["escaping_power_by_face_w"][face], abs=5 * uncertainty
        )


@pytest.mark.parametrize("anisotropy", [0, 0.8])
def test_hg_scattering_preserves_norm_and_known_moments(anisotropy):
    """HG angular moments have P1=g and P2=g squared for any initial direction."""
    count = 50_000
    direction = np.array([1, 2, 3]) / np.sqrt(14)
    scattered = _scatter_directions(
        directions=np.tile(direction, (count, 1)),
        anisotropy=np.full(count, anisotropy, dtype=float),
        rng=np.random.default_rng(1729),
    )
    np.testing.assert_allclose(np.linalg.norm(scattered, axis=1), 1, atol=1e-14)
    cosine = scattered @ direction
    for values, expected in (
        (cosine, anisotropy),
        ((3 * cosine**2 - 1) / 2, anisotropy**2),
    ):
        standard_error = np.std(values, ddof=1) / np.sqrt(count)
        assert np.mean(values) == pytest.approx(expected, abs=5 * standard_error)


def test_nonabsorbing_sample_accounts_for_all_six_exit_faces():
    """No absorption or roulette means all launched weight exits through one face."""
    result = simulate_layered_slab_reference(
        layers=(
            LayerOptics(
                thickness_mm=1,
                absorption_per_mm=0,
                reduced_scattering_per_mm=1,
                anisotropy=0,
            ),
        ),
        sample_width_mm=2,
        sample_height_mm=1,
        source_radius_mm=0.1,
        detector_offsets_mm=(0.4,),
        detector_radius_mm=0.2,
        photon_count=10_000,
        seed=1729,
        optical_power_w=1,
    )
    assert all(power > 0 for power in result["escaping_power_by_face_w"].values())
    assert sum(result["escaping_power_by_face_w"].values()) == pytest.approx(
        1, abs=1e-12
    )
    assert result["absorbed_power_w"] == 0
    assert abs(result["energy_balance_residual_w"]) < 1e-12
    assert result["detector_power_w"][0] <= result["escaping_power_by_face_w"]["z_min"]
