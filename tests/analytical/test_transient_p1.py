"""Time-dependent P1 checks against a frequency-domain slab solution."""

import numpy as np
import pytest
import MomentRT


def configuration(
    *, step_ps=2.0, duration_ps=1000.0, layered=False, detector=False, reflectance=0
):
    """Build a sample with uniform illumination and reflecting lateral faces."""
    unit = MomentRT.ureg
    material = MomentRT.Material(
        absorption=0.02 / unit.mm,
        reduced_scattering=1 / unit.mm,
        refractive_index=1.4,
        phase_function=MomentRT.Isotropic(),
    )
    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * unit.mm,
        size=np.array([2, 2, 2]) * unit.mm,
    )
    sample.add_background(material=material)
    if layered:
        sample.add_layer(
            name="Upper",
            thickness=1 * unit.mm,
            material=MomentRT.Material(
                absorption=0.04 / unit.mm,
                reduced_scattering=1.5 / unit.mm,
                refractive_index=1.4,
                phase_function=MomentRT.Isotropic(),
            ),
        )
    boundaries = {
        face: MomentRT.ReflectingBoundary() for face in MomentRT.boundary.FACES
    }
    boundaries["z_min"] = MomentRT.RobinBoundary(effective_reflectance=reflectance)
    boundaries["z_max"] = MomentRT.RobinBoundary(effective_reflectance=reflectance)
    source = MomentRT.UniformSurfaceSource(
        face="z_min",
        wavelength=850 * unit.nm,
        optical_power=1 * unit.mW,
    )
    mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=(3, 3, 13))
    times = np.arange(round(duration_ps / step_ps) + 1) * step_ps * unit.ps
    profile = 1 - np.cos(2 * np.pi * times.to("ps").magnitude / 100)
    solver = MomentRT.TransientP1(
        boundaries=boundaries,
        times=times,
        power_profile=profile,
        energy_speed=2e8 * unit.m / unit.s,
    )
    detectors = {}
    if detector:
        detectors = {
            "out": MomentRT.SurfaceDetector.full_hemisphere_on_face(
                sample=sample,
                face="z_max",
                coordinates=np.array([1, 1]) * unit.mm,
                radius=0.7 * unit.mm,
            )
        }
    return MomentRT.Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=solver,
        detectors=detectors,
    )


def exact_slab_transmission(*, angular_frequency):
    """Return independent complex P1 fluence at the far face of a slab."""
    speed = 2e8
    absorption = 20.0
    transport = 1020.0
    thickness = 0.002
    beta = 0.5
    incoming_irradiance = 0.001 / 0.002**2
    diffusion = 1 / (3 * (transport + 1j * angular_frequency / speed))
    wave_number = np.sqrt((absorption + 1j * angular_frequency / speed) / diffusion)
    attenuation = diffusion * wave_number
    denominator = 2 * beta * np.cosh(wave_number * thickness) + (
        attenuation + beta**2 / attenuation
    ) * np.sinh(wave_number * thickness)
    return 2 * incoming_irradiance / denominator


def test_periodic_slab_against_independent_complex_p1_solution():
    """Late transient modulation agrees with an independent slab calculation."""
    simulation = configuration()
    result = simulation.run()
    omega = 2 * np.pi / (100e-12)
    steady = exact_slab_transmission(angular_frequency=0).real
    harmonic = exact_slab_transmission(angular_frequency=omega)
    for index in (-1, -13, -26, -38):
        time = result.times[index].to("s").magnitude
        expected = steady - np.real(harmonic * np.exp(1j * omega * time))
        calculated = (
            result.frames[index].fluence.values.to("W/m^2").magnitude[:, :, -1].mean()
        )
        assert calculated == pytest.approx(expected, rel=0.06)


@pytest.mark.parametrize("reflectance", [0, 0.4])
def test_layered_pulse_preserves_energy_units_and_detector_outputs(reflectance):
    """Layered time steps satisfy midpoint conservation and retain signed signals."""
    simulation = configuration(
        step_ps=2,
        duration_ps=80,
        layered=True,
        detector=True,
        reflectance=reflectance,
    )
    assert simulation.validate().valid
    result = simulation.run()
    assert isinstance(result, MomentRT.TransientResult)
    assert result.times.shape == (41,)
    assert len(result.frames) == 41
    assert result.detector_power["out"].shape == (41,)
    assert result.detector_power["out"].check("[power]")
    assert result.stored_energy.check("[energy]")
    np.testing.assert_allclose(
        result.net_outward_power.to("W").magnitude,
        (result.escaping_power - result.incoming_power).to("W").magnitude,
    )
    assert set(result.frames[-1].outgoing_power_by_face) == set(MomentRT.boundary.FACES)
    assert result.frames[0].fluence.values.max().magnitude == 0
    assert np.max(np.abs(result.interval_energy_residual.to("W").magnitude)) < 1e-12
    assert result.frames[-1].flux.values.shape == (
        len(result.frames[-1].flux.values),
        3,
    )
    assert result.frames[-1].fluence.metadata["time"] == result.times[-1]
    with pytest.raises(ValueError, match="does not support wavelength sweeps"):
        simulation.sweep(wavelengths=np.array([800, 850]) * MomentRT.ureg.nm)


def test_time_refinement_reduces_layered_solution_error():
    """Midpoint time refinement converges on a fixed heterogeneous mesh."""
    values = []
    for step in (4, 2, 1, 0.5):
        result = configuration(step_ps=step, duration_ps=200, layered=True).run()
        values.append(result.frames[-1].fluence.values.to("W/m^2").magnitude[1, 1, -1])
    coarse_error = abs(values[0] - values[-1])
    medium_error = abs(values[1] - values[-1])
    fine_error = abs(values[2] - values[-1])
    assert coarse_error > medium_error > fine_error
    assert coarse_error / medium_error > 3


def test_layered_continuous_illumination_approaches_steady_p1():
    """The time-dependent current relaxes to the steady layered FEM solution."""
    simulation = configuration(step_ps=4, duration_ps=800, layered=True)
    simulation.solver.power_profile[:] = 1
    simulation.solver.power_profile[0] = 0
    transient = simulation.run()
    steady = MomentRT.Simulation(
        medium=simulation.medium,
        source=simulation.source,
        mesh=simulation.mesh,
        solver=MomentRT.P1(boundaries=simulation.solver.boundaries),
    ).run()
    np.testing.assert_allclose(
        transient.frames[-1].fluence.values.to("W/m^2").magnitude,
        steady.fluence.values.to("W/m^2").magnitude,
        rtol=0.001,
    )


def test_transient_input_validation():
    """Time grid, source pulse, and energy speed must be explicit and valid."""
    unit = MomentRT.ureg
    boundary = MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0)
    )
    arguments = dict(
        boundaries=boundary,
        times=np.array([0, 1, 2]) * unit.ps,
        power_profile=[0, 1, 0],
        energy_speed=2e8 * unit.m / unit.s,
    )
    for replacement in (
        {"times": np.array([0, 1, 1]) * unit.ps},
        {"times": np.array([1, 2, 3]) * unit.ps},
        {"power_profile": [1, 1, 0]},
        {"power_profile": [0, -1, 0]},
        {"energy_speed": 0 * unit.m / unit.s},
    ):
        with pytest.raises((TypeError, ValueError)):
            MomentRT.TransientP1(**(arguments | replacement))


def test_transient_broadband_source_is_rejected_by_preflight_and_run():
    """Preflight reports the same monochromatic limitation as the solver."""
    unit = MomentRT.ureg
    simulation = configuration(step_ps=2, duration_ps=10)
    source = MomentRT.UniformSurfaceSource(
        face="z_min",
        spectrum=MomentRT.SourceSpectrum(
            wavelengths=np.array([800, 900]) * unit.nm,
            spectral_power=np.array([1, 1]) * unit.uW / unit.nm,
        ),
    )
    simulation = MomentRT.Simulation(
        medium=simulation.medium,
        source=source,
        mesh=simulation.mesh,
        solver=simulation.solver,
    )
    report = simulation.validate()
    assert any("monochromatic source" in error for error in report.errors)
    with pytest.raises(ValueError, match="monochromatic source"):
        simulation.run()
