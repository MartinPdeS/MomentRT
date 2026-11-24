"""Face placement, preflight, and reporting with real finite solutions."""

import numpy as np
import pytest
import MomentRT


def configuration():
    ureg = MomentRT.ureg

    material = MomentRT.Material(
        absorption=0.02 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
        anisotropy=0.8,
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )

    sample = MomentRT.StructuredMedium(
        origin=np.array([-2, -3, 1]) * ureg.mm,
        size=np.array([4, 6, 3]) * ureg.mm,
    )

    sample.add_background(material=material)

    source = MomentRT.SurfaceSource.on_face(
        sample=sample,
        face="z_min",
        coordinates=np.array([-1, 0]) * ureg.mm,
        radius=0.7 * ureg.mm,
        wavelength=850 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    detector = MomentRT.SurfaceDetector.on_face(
        sample=sample,
        face="z_min",
        coordinates=np.array([1, 0]) * ureg.mm,
        radius=0.7 * ureg.mm,
        acceptance_angle=90 * ureg.degree,
        orientation="outward_normal",
    )

    boundary = MomentRT.RobinBoundary(effective_reflectance=0)

    boundaries = MomentRT.Boundaries.uniform(condition=boundary)

    solver = MomentRT.SP3(boundaries=boundaries)

    mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=(9, 9, 7))

    return MomentRT.Simulation(
        medium=sample,
        source=source,
        solver=solver,
        mesh=mesh,
        detectors={"contact": detector},
    )


@pytest.mark.parametrize(
    "face,coordinates,position,axis",
    [
        ("x_min", [0, 2], [-2, 0, 2], [-1, 0, 0]),
        ("x_max", [0, 2], [2, 0, 2], [1, 0, 0]),
        ("y_min", [0, 2], [0, -3, 2], [0, -1, 0]),
        ("y_max", [0, 2], [0, 3, 2], [0, 1, 0]),
        ("z_min", [0, 0], [0, 0, 1], [0, 0, -1]),
        ("z_max", [0.25, -0.75], [0.25, -0.75, 4], [0, 0, 1]),
    ],
)
def test_face_coordinates_and_outward_axis(face, coordinates, position, axis):
    simulation = configuration()

    detector = MomentRT.SurfaceDetector.on_face(
        sample=simulation.medium,
        face=face,
        coordinates=np.array(coordinates) * MomentRT.ureg.mm,
        radius=0.2 * MomentRT.ureg.mm,
        acceptance_angle=45 * MomentRT.ureg.degree,
        orientation="outward_normal",
    )

    np.testing.assert_allclose(detector.position.to("mm").magnitude, position)

    np.testing.assert_array_equal(detector.collection_axis, axis)


def test_boundaries_override_is_a_snapshot():
    condition = MomentRT.RobinBoundary(effective_reflectance=0)

    boundaries = MomentRT.Boundaries.uniform(condition=condition)

    overridden = boundaries.with_overrides(
        conditions={"z_max": MomentRT.ReflectingBoundary()}
    )

    condition.effective_reflectance = 0.5

    assert boundaries["z_min"].effective_reflectance == 0

    assert isinstance(boundaries["z_max"], MomentRT.RobinBoundary)

    assert isinstance(overridden["z_max"], MomentRT.ReflectingBoundary)

    with pytest.raises(ValueError, match="Unknown"):
        boundaries.with_overrides(conditions={"wrong": condition})


def test_preflight_footprint_areas_match_solution_and_does_not_solve(monkeypatch):
    simulation = configuration()

    solver = simulation.solver

    original = solver.solve

    def forbidden(**kwargs):
        raise AssertionError("Preflight must not solve the PDE")

    monkeypatch.setattr(solver, "solve", forbidden)

    report = simulation.validate()

    assert report.valid

    assert "does not establish mesh convergence" in report.summary()

    monkeypatch.setattr(solver, "solve", original)

    result = simulation.run()

    assert report.source_sampled_area.to("m^2").magnitude == pytest.approx(
        result.source_sampled_area.to("m^2").magnitude
    )

    assert report.detector_sampled_areas["contact"].to(
        "m^2"
    ).magnitude == pytest.approx(
        result.detectors["contact"].sampled_area.to("m^2").magnitude
    )

    assert "Detector contact" in result.summary()

    assert "Residual power" in result.summary()


def test_preflight_collects_independent_errors():
    simulation = configuration()

    simulation.medium.background.anisotropy = None

    simulation.solver.boundaries["z_min"] = MomentRT.ReflectingBoundary()

    simulation.detectors["contact"].position = np.array([20, 0, 1]) * MomentRT.ureg.mm

    report = simulation.validate()

    assert not report.valid

    assert any("anisotropy" in error for error in report.errors)

    assert any("illuminated face" in error for error in report.errors)

    assert any("Detector 'contact'" in error for error in report.errors)


def test_preflight_reports_spectral_coverage_and_unresolved_detector():
    simulation = configuration()

    material = MomentRT.SpectralMaterial(
        wavelengths=np.array([800, 900]) * MomentRT.ureg.nm,
        absorption=np.array([0.02, 0.03]) / MomentRT.ureg.mm,
        reduced_scattering=np.array([1.0, 1.1]) / MomentRT.ureg.mm,
        refractive_index=1.4,
        anisotropy=[0.7, 0.8],
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )

    simulation.medium.add_background(material=material)

    simulation.detectors["contact"].radius = 1e-8 * MomentRT.ureg.mm

    report = simulation.validate(wavelengths=np.array([750, 850]) * MomentRT.ureg.nm)

    assert not report.valid

    assert any("outside" in error for error in report.errors)

    assert any("not sampled" in error for error in report.errors)

    assert report.detector_sampled_areas["contact"].magnitude == 0

    assert isinstance(
        material.at(wavelength=850 * MomentRT.ureg.nm).phase_function,
        MomentRT.phase_function.HenyeyGreenstein,
    )


def test_face_placement_rejects_footprints_outside_sample():
    simulation = configuration()

    with pytest.raises(ValueError, match="fit"):
        MomentRT.SurfaceSource.on_face(
            sample=simulation.medium,
            face="z_min",
            coordinates=np.array([2, 0]) * MomentRT.ureg.mm,
            radius=0.5 * MomentRT.ureg.mm,
            wavelength=850 * MomentRT.ureg.nm,
            optical_power=1 * MomentRT.ureg.mW,
        )


def test_material_phase_function_is_visible_and_checked():
    simulation = configuration()

    assert isinstance(
        simulation.medium.background.phase_function,
        MomentRT.phase_function.HenyeyGreenstein,
    )

    with pytest.raises(ValueError, match="phase_function"):
        MomentRT.Material(
            absorption=1 / MomentRT.ureg.m,
            reduced_scattering=1 / MomentRT.ureg.mm,
            refractive_index=1.4,
            phase_function="unknown",
        )


def test_spectral_detector_measurements_and_result_summary():
    simulation = configuration()

    spectrum = simulation.sweep(
        wavelengths=np.array([850, 800]) * MomentRT.ureg.nm,
    )

    measurements = spectrum.detectors["contact"]

    assert [
        measurement.wavelength.to("nm").magnitude for measurement in measurements
    ] == [
        850,
        800,
    ]

    for index, measurement in enumerate(measurements):
        assert measurement is spectrum.results[index].detectors["contact"]

        assert measurement.collected_power.to("W").magnitude == pytest.approx(
            spectrum.detector_power["contact"][index].to("W").magnitude
        )

        assert measurement.sampled_area.to("mm^2").magnitude > 0

        assert measurement.geometric_area.to("mm^2").magnitude == pytest.approx(
            np.pi * 0.7**2
        )

    summary = spectrum.results[0].summary()

    assert "Model: steady-state tetrahedral SP3 FEM" in summary

    assert "Wavelength: 850" in summary

    assert "Mesh shape: (9, 9, 7)" in summary

    assert "convergence is not established" in summary


def test_configuration_errors_identify_affected_face_and_bounds():
    simulation = configuration()

    with pytest.raises(ValueError) as error:
        MomentRT.SurfaceDetector.full_hemisphere_on_face(
            sample=simulation.medium,
            face="z_min",
            coordinates=np.array([2, 0]) * MomentRT.ureg.mm,
            radius=0.5 * MomentRT.ureg.mm,
        )

    message = str(error.value)

    assert "z_min" in message and "radius=" in message

    assert "x=[" in message and "y=[" in message

    boundaries = MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0),
    ).with_overrides(conditions={"x_max": MomentRT.ReflectingBoundary()})

    with pytest.raises(ValueError, match="x_max"):
        MomentRT.P3(boundaries=boundaries)


def test_refractive_index_errors_identify_materials():
    simulation = configuration()

    other = MomentRT.Material(
        absorption=0.02 / MomentRT.ureg.mm,
        reduced_scattering=1 / MomentRT.ureg.mm,
        refractive_index=1.5,
        anisotropy=0.8,
        phase_function=MomentRT.HenyeyGreenstein(),
    )

    simulation.medium.add_layer(
        name="upper layer",
        thickness=0.5 * MomentRT.ureg.mm,
        material=other,
    )

    report = simulation.validate()

    assert any(
        "background" in error and "upper layer" in error for error in report.errors
    )

    with pytest.raises(ValueError, match="upper layer"):
        simulation.run()
