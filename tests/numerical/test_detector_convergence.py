"""Convergence criteria, area aliases, unresolved levels, and exact slab signals."""

import numpy as np
import pytest
import MomentRT
from MomentRT import (
    CartesianMesh,
    DetectorConvergenceResult,
    DetectorResult,
    Diffusion,
    Material,
    ReflectingBoundary,
    RobinBoundary,
    Simulation,
    SimulationResult,
    StructuredMedium,
    SurfaceDetector,
    SurfaceSource,
    UniformSurfaceSource,
    ureg,
)


def configuration():
    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.array([2, 2, 4]) * ureg.mm,
    )

    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    sample.add_background(material=material)

    source = UniformSurfaceSource(
        face="z_min",
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    boundaries = {
        face: ReflectingBoundary()
        for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
    }

    boundaries["z_min"] = RobinBoundary(effective_reflectance=0.3)

    boundaries["z_max"] = RobinBoundary(effective_reflectance=0.3)

    detector = SurfaceDetector(
        face="z_max",
        position=np.array([1, 1, 4]) * ureg.mm,
        radius=0.7 * ureg.mm,
        collection_axis=[0, 0, 1],
        acceptance_angle=90 * ureg.degree,
    )

    return Simulation(
        medium=sample,
        mesh=CartesianMesh.from_sample(sample=sample, shape=(5, 5, 9)),
        source=source,
        solver=Diffusion(boundaries=boundaries),
        detectors={"transmission": detector},
    )


def synthetic_report(*, powers, area_ratios, power_atol=None):
    simulation = configuration()

    detector = simulation.detectors["transmission"]

    geometric_area = np.pi * detector.radius**2

    meshes = [
        CartesianMesh.from_sample(sample=simulation.medium, shape=(count, count, count))
        for count in [3, 5, 9, 17, 33, 65][: len(powers)]
    ]

    results = [
        SimulationResult(
            fluence=None,
            diagnostics={"source_quadrature_area": 4 * ureg.mm**2},
            detectors={
                "transmission": DetectorResult(
                    collected_power=power * ureg.watt,
                    sampled_area=ratio * geometric_area,
                    geometric_area=geometric_area,
                    angular_fraction=1,
                    wavelength=800 * ureg.nm,
                    minimum_escaping_irradiance=0 * ureg.W / ureg.m**2,
                    detector=detector,
                ),
            },
        )
        for power, ratio in zip(powers, area_ratios)
    ]

    return DetectorConvergenceResult(
        meshes=meshes,
        results=results,
        detectors=simulation.detectors,
        source_geometric_area=4 * ureg.mm**2,
        wavelength=800 * ureg.nm,
        level_errors=[{} for _ in powers],
        power_rtol=0.01,
        area_rtol=0.02,
        power_atol=0 * ureg.W if power_atol is None else power_atol,
        required_consecutive=2,
    )


def test_recommendation_needs_consecutive_power_and_geometric_area_checks():
    report = synthetic_report(
        powers=[0.0015, 0.00101, 0.001005, 0.001004],
        area_ratios=[0.9, 0.99, 1.001, 1.002],
    )

    assert report.recommended_shape == (17, 17, 17)

    np.testing.assert_array_equal(report.stable_refinements, [0, 0, 1, 2])

    assert report.power_relative_change["transmission"][2] == pytest.approx(
        0.000005 / 0.001005
    )

    assert np.isnan(report.power_relative_change["transmission"][0])

    assert report.records[-1]["detectors"]["transmission"]["collected_power"].to(
        "mW"
    ).magnitude == pytest.approx(1.004)

    unresolved_area = synthetic_report(
        powers=[0.001, 0.001, 0.001],
        area_ratios=[0.5, 0.5, 0.5],
    )

    assert not unresolved_area.converged

    assert unresolved_area.recommended_mesh is None


def test_zero_signals_and_absolute_power_tolerance_do_not_divide_by_zero():
    zero = synthetic_report(powers=[0, 0, 0], area_ratios=[1, 1, 1])

    assert zero.converged

    np.testing.assert_array_equal(zero.power_relative_change["transmission"][1:], 0)

    absolute = synthetic_report(
        powers=[1e-9, 0, 0],
        area_ratios=[1, 1, 1],
        power_atol=2e-9 * ureg.W,
    )

    assert absolute.converged

    assert np.isinf(absolute.power_relative_change["transmission"][1])


def test_late_failure_or_accumulated_drift_invalidates_a_coarse_recommendation():
    late_failure = synthetic_report(powers=[1, 1, 1, 2], area_ratios=[1, 1, 1, 1])

    assert not late_failure.converged

    drift = synthetic_report(
        powers=[1, 1.009, 1.018, 1.027, 1.036, 1.045],
        area_ratios=[1, 1, 1, 1, 1, 1],
    )

    assert drift.recommended_index == 4

    assert drift.power_relative_to_finest["transmission"][2] > drift.power_rtol


def test_real_study_preserves_configuration_and_matches_linear_slab():
    simulation = configuration()

    original_mesh = simulation.mesh

    original_source = simulation.source

    study = simulation.convergence_study(
        mesh_shapes=[(5, 5, 9), (9, 9, 17), (13, 13, 25)],
        power_rtol=0.15,
        area_rtol=0.15,
    )

    beta = 0.7 / 2.6

    current = (2 / 1.3) * 250 / (2 + beta * 0.004 / (1 / 3000))

    np.testing.assert_allclose(
        study.detector_power["transmission"].to("W").magnitude,
        current * study.detector_sampled_area["transmission"].to("m^2").magnitude,
        rtol=1e-11,
    )

    np.testing.assert_allclose(study.source_sampled_area.to("mm^2").magnitude, 4)

    assert study.converged

    summary = study.results[-1].summary(convergence=study)

    assert "Supplied convergence study: converged" in summary

    assert "Verify that the study applies" in summary

    assert simulation.mesh is original_mesh

    assert simulation.source is original_source

    assert not any(study.level_errors)


def test_unresolved_source_level_is_recorded_and_finer_mesh_still_runs():
    simulation = configuration()

    simulation.source = SurfaceSource(
        face="z_min",
        position=np.array([0.75, 0.75, 0]) * ureg.mm,
        radius=0.03 * ureg.mm,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    study = simulation.convergence_study(
        mesh_shapes=[(3, 3, 5), (17, 17, 33)],
    )

    assert study.results[0] is None

    assert "source" in study.level_errors[0]

    assert study.source_sampled_area[0].magnitude == 0

    assert study.results[1] is not None

    assert not study.converged


def test_unresolved_detector_does_not_abort_study():
    simulation = configuration()

    simulation.detectors["tiny"] = SurfaceDetector(
        face="z_max",
        position=np.array([0.75, 0.75, 4]) * ureg.mm,
        radius=0.03 * ureg.mm,
        collection_axis=[0, 0, 1],
        acceptance_angle=90 * ureg.degree,
    )

    study = simulation.convergence_study(mesh_shapes=[(3, 3, 5), (5, 5, 9)])

    assert all("detector:tiny" in errors for errors in study.level_errors)

    assert all(result is not None for result in study.results)

    np.testing.assert_array_equal(study.detector_sampled_area["tiny"].magnitude, 0)

    assert not study.converged


def test_invalid_interface_is_not_hidden_by_an_unresolved_source():
    simulation = configuration()

    other = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0.01 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.5,
    )

    simulation.medium.add_layer(
        name="unsupported index jump", thickness=1 * ureg.mm, material=other
    )

    simulation.source = SurfaceSource(
        face="z_min",
        position=np.array([0.75, 0.75, 0]) * ureg.mm,
        radius=1e-8 * ureg.mm,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    with pytest.raises(ValueError, match="refractive-index"):
        simulation.convergence_study(mesh_shapes=[(3, 3, 5), (5, 5, 9)])


@pytest.mark.parametrize(
    "arguments",
    [
        {"mesh_shapes": [(5, 5, 5)]},
        {"mesh_shapes": [(5, 5, 5), (3, 9, 9)]},
        {"mesh_shapes": [(5, 5, 5), (5, 5, 5)]},
        {"mesh_shapes": [(3, 3, 3), (5, 5, 5)], "power_rtol": -1},
        {"mesh_shapes": [(3, 3, 3), (5, 5, 5)], "required_consecutive": True},
        {"mesh_shapes": [(3, 3, 3), (5, 5, 5)], "power_atol": 1},
    ],
)
def test_invalid_study_inputs_fail_before_solving(arguments):
    with pytest.raises((ValueError, TypeError)):
        configuration().convergence_study(**arguments)
