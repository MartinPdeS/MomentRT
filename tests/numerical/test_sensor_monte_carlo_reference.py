"""Checks for the independent layered-slab photon reference."""

import numpy as np
from tools.sensor_reference import LayerOptics, simulate_layered_slab_reference


def test_layered_slab_reference_is_seeded_and_reports_uncertainty() -> None:
    """Repeated seeded runs agree and return finite power uncertainties."""
    layers = (
        LayerOptics(
            thickness_mm=0.5,
            absorption_per_mm=0.1,
            reduced_scattering_per_mm=1.0,
            anisotropy=0.8,
        ),
        LayerOptics(
            thickness_mm=0.5,
            absorption_per_mm=0.05,
            reduced_scattering_per_mm=0.8,
            anisotropy=0.7,
        ),
    )
    arguments = {
        "layers": layers,
        "sample_width_mm": 4.0,
        "sample_height_mm": 1.0,
        "source_radius_mm": 0.2,
        "detector_offsets_mm": (0.8, 1.5),
        "detector_radius_mm": 0.2,
        "photon_count": 5_000,
        "seed": 1729,
    }

    first = simulate_layered_slab_reference(**arguments)
    second = simulate_layered_slab_reference(**arguments)

    np.testing.assert_array_equal(first["detector_power_w"], second["detector_power_w"])
    np.testing.assert_array_equal(first["standard_error_w"], second["standard_error_w"])
    assert np.all(np.asarray(first["detector_power_w"]) > 0)
    assert np.all(np.asarray(first["standard_error_w"]) > 0)
    assert np.all(np.isfinite(first["relative_standard_error"]))
