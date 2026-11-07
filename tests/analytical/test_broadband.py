"""Independent slab powers, electrical quadrature, and spectral API checks."""

import copy
import numpy as np
import pytest
import MomentRT


def spectrum():
    return MomentRT.SourceSpectrum(
        wavelengths=np.array([700, 750, 900]) * MomentRT.ureg.nm,
        spectral_power=np.array([0, 2, 1]) * MomentRT.ureg.uW / MomentRT.ureg.nm,
    )


def response():
    return MomentRT.DetectorResponsivity(
        wavelengths=np.array([700, 750, 900]) * MomentRT.ureg.nm,
        responsivity=np.array([0.2, 0.3, 0.6]) * MomentRT.ureg.A / MomentRT.ureg.W,
    )


def simulation(*, solver_type=MomentRT.Diffusion, slab=False, responsivity=True):
    unit = MomentRT.ureg
    material = MomentRT.Material(
        phase_function=MomentRT.HenyeyGreenstein(),
        absorption=(0 if slab else 0.03) / unit.mm,
        reduced_scattering=1 / unit.mm,
        refractive_index=1.4,
        anisotropy=0.8,
    )
    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * unit.mm,
        size=np.array([2, 2, 4]) * unit.mm,
    )
    sample.add_background(material=material)
    source = MomentRT.UniformSurfaceSource(face="z_min", spectrum=spectrum())
    detector = MomentRT.SurfaceDetector.full_hemisphere_on_face(
        sample=sample,
        face="z_max",
        coordinates=np.array([1, 1]) * unit.mm,
        radius=0.7 * unit.mm,
        responsivity=response() if responsivity else None,
    )
    boundaries = {
        face: (
            MomentRT.ReflectingBoundary()
            if slab and face[0] != "z"
            else MomentRT.RobinBoundary(effective_reflectance=0)
        )
        for face in MomentRT.boundary.FACES
    }
    return MomentRT.Simulation(
        medium=sample,
        source=source,
        mesh=MomentRT.TetrahedralMesh.from_sample(
            sample=sample,
            resolution=0.7 * unit.mm,
        ),
        solver=solver_type(boundaries=boundaries),
        detectors={"sensor": detector},
    )


def test_broadband_power_and_current_against_exact_slab():
    config = simulation(slab=True)
    assert config.validate().valid
    result = config.run()
    assert isinstance(result, MomentRT.BroadbandSimulationResult)
    measurement = result.detectors["sensor"]
    reference = result.reference_results.detectors["sensor"][0]
    # Exact nonabsorbing P1 slab: D = 1/3 mm; the configured partial-current source factor
    # is two, giving T = 4D/(L+4D) = 1/4.
    fraction = 0.25 * reference.sampled_area.to("mm^2").magnitude / 4
    np.testing.assert_allclose(
        measurement.spectral_power.to("uW/nm").magnitude,
        np.array([0, 2, 1]) * fraction,
        rtol=1e-10,
        atol=1e-14,
    )
    # Hand-evaluated nonuniform trapezoids, including the zero-density endpoint.
    assert result.source_spectrum.optical_power.to("uW").magnitude == pytest.approx(275)
    assert measurement.collected_power.to("uW").magnitude == pytest.approx(
        275 * fraction, rel=1e-10
    )
    assert measurement.photocurrent.to("uA").magnitude == pytest.approx(
        105 * fraction, rel=1e-10
    )
    assert result.diagnostics["geometry_preparations"] == 1
    assert result.reference_optical_power == 1 * MomentRT.ureg.W
    assert config.source.spectrum is not None
    assert config.source.optical_power.to("uW").magnitude == pytest.approx(275)


@pytest.mark.parametrize("solver_type", [MomentRT.Diffusion, MomentRT.SP3, MomentRT.P3])
def test_broadband_matches_fresh_monochromatic_solves(solver_type):
    config = simulation(solver_type=solver_type, responsivity=False)
    broadband = config.run()
    powers = []
    for wavelength in config.source.spectrum.wavelengths:
        monochromatic = copy.copy(config)
        monochromatic.source = MomentRT.UniformSurfaceSource(
            face="z_min",
            wavelength=wavelength,
            optical_power=2 * MomentRT.ureg.mW,
        )
        powers.append(
            monochromatic.run().detector_power["sensor"].to("W").magnitude / 0.002
        )
    np.testing.assert_allclose(
        broadband.detectors["sensor"].spectral_power.to("uW/nm").magnitude,
        np.array(powers) * [0, 2, 1],
        rtol=1e-11,
        atol=1e-14,
    )
    assert broadband.detectors["sensor"].photocurrent is None


def test_responsivity_coverage_fails_before_preparation(monkeypatch):
    config = simulation()
    config.detectors["sensor"].responsivity = MomentRT.DetectorResponsivity(
        wavelengths=np.array([720, 900]) * MomentRT.ureg.nm,
        responsivity=np.ones(2) * MomentRT.ureg.A / MomentRT.ureg.W,
    )
    assert not config.validate().valid

    def unexpected(**kwargs):
        raise AssertionError("Must reject response coverage before FEM preparation")

    monkeypatch.setattr(config.solver, "_prepare", unexpected)
    with pytest.raises(ValueError, match="outside"):
        config.run()


def test_units_negative_outputs_and_spectral_grid_convergence():
    unit = MomentRT.ureg
    currents = []
    for count in (3, 5, 9):
        wavelength_nm = np.linspace(700, 900, count)
        measured = MomentRT.BroadbandDetectorResult(
            wavelengths=wavelength_nm / 1000 * unit.um,
            spectral_power=(wavelength_nm - 700) / 200 * unit.uW / unit.nm,
            responsivity=MomentRT.DetectorResponsivity(
                wavelengths=np.array([700, 900]) * unit.nm,
                responsivity=np.array([0, 1]) * unit.A / unit.W,
            ),
        )
        # Integral of the linear density is exactly 100 uW. The current
        # integrand is quadratic with exact integral 200/3 uA.
        assert measured.collected_power.to("uW").magnitude == pytest.approx(100)
        currents.append(measured.photocurrent.to("uA").magnitude)
    errors = np.array(currents) - 200 / 3
    np.testing.assert_allclose(errors[:-1] / errors[1:], 4, rtol=1e-10)
    negative = MomentRT.BroadbandDetectorResult(
        wavelengths=np.array([700, 900]) * unit.nm,
        spectral_power=np.array([-1, -1]) * unit.uW / unit.nm,
    )
    assert negative.collected_power.to("uW").magnitude == pytest.approx(-200)


def test_monochromatic_response_and_snapshot_preservation():
    config = simulation()
    config.source = MomentRT.UniformSurfaceSource(
        face="z_min",
        wavelength=750 * MomentRT.ureg.nm,
        optical_power=1 * MomentRT.ureg.mW,
    )
    result = config.run()
    measured = result.detectors["sensor"]
    assert (measured.photocurrent / measured.collected_power).to(
        "A/W"
    ).magnitude == pytest.approx(0.3)
    config.detectors["sensor"].responsivity.responsivity[:] = (
        0 * MomentRT.ureg.A / MomentRT.ureg.W
    )
    assert measured.detector.responsivity.at(
        wavelength=750 * MomentRT.ureg.nm
    ).magnitude == pytest.approx(0.3)


def test_broadband_snapshots_and_emission_mutation_guards():
    config = simulation()
    result = config.run()
    stored = result.detectors["sensor"].spectral_power.copy()
    config.source.spectrum.spectral_power[:] *= 2
    np.testing.assert_array_equal(
        result.source_spectrum.spectral_power.to("uW/nm").magnitude, [0, 2, 1]
    )
    np.testing.assert_array_equal(
        result.detectors["sensor"].spectral_power.magnitude, stored.magnitude
    )
    with pytest.raises(ValueError, match="spectrum"):
        config.source.optical_power = 1 * MomentRT.ureg.W
    with pytest.raises(ValueError, match="spectrum"):
        config.source.wavelength = 800 * MomentRT.ureg.nm


def test_broadband_geometry_plot_keeps_source_and_detector():
    config = simulation()
    figure = config.plot_system(show=False)
    assert any(trace.name == "Full-face illumination" for trace in figure.data)
    assert any(trace.name == "Surface detector" for trace in figure.data)


def test_emission_exclusivity_factories_and_invalid_tables():
    unit = MomentRT.ureg
    config = simulation()
    with pytest.raises(ValueError, match="exclusively"):
        MomentRT.UniformSurfaceSource(
            face="z_min", spectrum=spectrum(), wavelength=800 * unit.nm
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        MomentRT.SourceSpectrum(
            wavelengths=np.array([800, 700]) * unit.nm,
            spectral_power=np.ones(2) * unit.W / unit.nm,
        )
    with pytest.raises(ValueError, match="nonnegative"):
        MomentRT.SourceSpectrum(
            wavelengths=np.array([700, 800]) * unit.nm,
            spectral_power=np.array([1, -1]) * unit.W / unit.nm,
        )
    with pytest.raises(TypeError, match="Pint quantity"):
        MomentRT.SourceSpectrum(
            wavelengths=np.array([700, 800]) * unit.nm, spectral_power=np.ones(2)
        )
    disk = MomentRT.SurfaceSource.on_face(
        sample=config.medium,
        face="z_min",
        coordinates=np.array([1, 1]) * unit.mm,
        radius=0.4 * unit.mm,
        spectrum=spectrum(),
    )
    assert disk.spectrum is not None
    with pytest.raises(ValueError, match="monochromatic members"):
        MomentRT.SurfaceSourceGroup(sources=[disk])
    with pytest.raises(ValueError, match="no single wavelength"):
        _ = disk.wavelength
    with pytest.raises(ValueError, match="use run"):
        config.sweep(wavelengths=np.array([700, 800]) * unit.nm)
    with pytest.raises(ValueError, match="monochromatic source"):
        config.convergence_study(mesh_resolutions=[1 * unit.mm, 0.7 * unit.mm])


def test_face_acceptance_shortcut_and_legacy_factory():
    config = simulation()
    unit = MomentRT.ureg
    parameters = {
        "sample": config.medium,
        "face": "z_max",
        "coordinates": np.array([1, 1]) * unit.mm,
        "radius": 0.7 * unit.mm,
    }
    detector = MomentRT.SurfaceDetector.on_face(
        **parameters,
        acceptance="full_hemisphere",
    )
    legacy = MomentRT.SurfaceDetector.full_hemisphere_on_face(**parameters)
    assert detector.angular_fraction == pytest.approx(1)
    np.testing.assert_array_equal(detector.collection_axis, legacy.collection_axis)
    assert detector.acceptance_angle == legacy.acceptance_angle
    with pytest.raises(ValueError, match="exclusively"):
        MomentRT.SurfaceDetector.on_face(
            **parameters,
            acceptance="full_hemisphere",
            acceptance_angle=45 * unit.degree,
        )
    with pytest.raises(ValueError, match="Supply acceptance"):
        MomentRT.SurfaceDetector.on_face(**parameters)
    with pytest.raises(ValueError, match="full_hemisphere"):
        MomentRT.SurfaceDetector.on_face(**parameters, acceptance="automatic")


def test_measurements_share_area_interface_and_snapshot():
    config = simulation()
    broadband = config.run().detectors["sensor"]
    assert broadband.sampled_area == broadband.reference_measurements[0].sampled_area
    assert broadband.geometric_area == np.pi * (0.7 * MomentRT.ureg.mm) ** 2
    assert broadband.angular_fraction == pytest.approx(1)
    config.detectors["sensor"].radius *= 2
    assert broadband.detector.radius == 0.7 * MomentRT.ureg.mm
    monochromatic = broadband.reference_measurements[0]
    assert monochromatic.sampled_area == broadband.sampled_area
    assert monochromatic.geometric_area == broadband.geometric_area


def test_source_resampling_preserves_knots_units_and_power():
    source = spectrum()
    fine = source.resample(
        wavelengths=np.array([0.7, 0.725, 0.75, 0.825, 0.9]) * MomentRT.ureg.um,
    )
    np.testing.assert_allclose(
        fine.spectral_power.to("uW/nm").magnitude, [0, 1, 2, 1.5, 1]
    )
    assert fine.optical_power.to("uW").magnitude == pytest.approx(275)
    fine.spectral_power[1] *= 2
    np.testing.assert_array_equal(
        source.spectral_power.to("uW/nm").magnitude, [0, 2, 1]
    )
    with pytest.raises(ValueError, match="outside"):
        source.resample(wavelengths=np.array([650, 900]) * MomentRT.ureg.nm)
    source_figure = source.plot(wavelength_unit="um", power_unit="uW/nm", show=False)
    np.testing.assert_allclose(source_figure.data[0].x, [0.7, 0.75, 0.9])
    response_figure = response().plot(show=False)
    np.testing.assert_allclose(response_figure.data[0].y, [0.2, 0.3, 0.6])


def test_spectral_convergence_against_integrated_slab_response():
    config = simulation(slab=True)
    original = config.source.spectrum.wavelengths.copy()
    study = config.spectral_convergence(
        wavelength_grids=[
            np.linspace(700, 900, count) * MomentRT.ureg.nm for count in (5, 9, 17)
        ],
        power_rtol=1e-10,
        photocurrent_rtol=0.01,
    )
    assert study.converged
    assert study.recommended_index == 2
    assert len(study.recommended_wavelengths) == 17
    np.testing.assert_array_equal(study.sample_counts, [5, 9, 17])
    np.testing.assert_allclose(study.power_relative_change["sensor"][1:], 0, atol=1e-12)
    # Integrate the product of piecewise linear source and response exactly:
    # the two intervals contribute 40/3 and 195/2 microamperes before transport.
    aperture_fraction = (
        study.results[-1].detectors["sensor"].sampled_area.to("mm^2").magnitude / 4
    )
    exact_current = (40 / 3 + 195 / 2) * 0.25 * aperture_fraction
    errors = abs(
        study.detector_photocurrent["sensor"].to("uA").magnitude - exact_current
    )
    np.testing.assert_allclose(errors[:-1] / errors[1:], 4, rtol=1e-8)
    np.testing.assert_array_equal(
        config.source.spectrum.wavelengths.magnitude, original.magnitude
    )
    assert study.records[-1]["meets_tolerances"]
    assert "Converged: True" in study.summary()
    assert len(study.plot(show=False).data) == 2
    assert config.mesh_convergence.__func__ is config.convergence_study.__func__


def test_spectral_grids_keep_source_material_and_response_knots():
    config = simulation(slab=True)
    unit = MomentRT.ureg
    config.medium.add_background(
        material=MomentRT.SpectralMaterial(
            phase_function=MomentRT.HenyeyGreenstein(),
            wavelengths=np.array([700, 730, 900]) * unit.nm,
            absorption=np.zeros(3) / unit.mm,
            reduced_scattering=np.ones(3) / unit.mm,
            refractive_index=1.4,
        ),
    )
    config.detectors["sensor"].responsivity = MomentRT.DetectorResponsivity(
        wavelengths=np.array([700, 710, 900]) * unit.nm,
        responsivity=np.array([0.1, 0.8, 0.5]) * unit.A / unit.W,
    )
    grids = [np.linspace(700, 900, count) * unit.nm for count in (5, 9, 17)]
    study = config.spectral_convergence(wavelength_grids=grids)
    for grid, result in zip(study.wavelength_grids, study.results):
        for knot in (710, 730, 750):
            assert np.any(np.isclose(grid.to("nm").magnitude, knot))
        assert result.source_spectrum.optical_power.to("uW").magnitude == pytest.approx(
            275
        )
    plain = config.spectral_convergence(
        wavelength_grids=grids, include_table_knots=False
    )
    assert len(plain.wavelength_grids[0]) == 5
    with pytest.raises(ValueError, match="source interval"):
        config.spectral_convergence(
            wavelength_grids=[np.array([710, 900]) * unit.nm] * 2
        )
    with pytest.raises(ValueError, match="nested"):
        config.spectral_convergence(
            wavelength_grids=[np.array([700, 800, 900]) * unit.nm] * 2
        )


def test_failed_final_spectral_refinement_cannot_recommend(monkeypatch):
    config = simulation(slab=True)
    original_run = MomentRT.Simulation.run

    def fail_fine(self):
        if (
            self.source.spectrum is not None
            and len(self.source.spectrum.wavelengths) >= 33
        ):
            raise RuntimeError("Injected fine-level failure")
        return original_run(self)

    monkeypatch.setattr(MomentRT.Simulation, "run", fail_fine)
    study = config.spectral_convergence(
        wavelength_grids=[
            np.linspace(700, 900, count) * MomentRT.ureg.nm for count in (5, 9, 17, 33)
        ],
        required_consecutive=1,
        photocurrent_rtol=0.1,
    )
    assert np.any(study.meets_tolerances[:-1])
    assert study.results[-1] is None
    assert np.isnan(study.detector_power["sensor"].magnitude[-1])
    assert not study.converged
    assert study.recommended_wavelengths is None
    assert "Injected fine-level failure" in study.records[-1]["errors"]["simulation"]


def test_negative_final_integrated_output_is_retained(monkeypatch):
    config = simulation(slab=True)
    original_run = MomentRT.Simulation.run

    def negative_fine(self):
        result = original_run(self)
        if (
            self.source.spectrum is not None
            and len(self.source.spectrum.wavelengths) >= 17
        ):
            result.detectors["sensor"].collected_power *= -1
        return result

    monkeypatch.setattr(MomentRT.Simulation, "run", negative_fine)
    study = config.spectral_convergence(
        wavelength_grids=[
            np.linspace(700, 900, count) * MomentRT.ureg.nm for count in (5, 9, 17)
        ],
        required_consecutive=1,
        photocurrent_rtol=0.1,
    )
    assert study.detector_power["sensor"][-1].magnitude < 0
    assert not study.converged
    assert "Negative integrated" in study.level_errors[-1]["detector:sensor"]
