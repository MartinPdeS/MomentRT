"""Independent constant-current slab and Lambertian projected-solid-angle checks."""

import numpy as np
import pytest
import MomentRT
from MomentRT import (
    CartesianMesh,
    Diffusion,
    Material,
    ReflectingBoundary,
    RobinBoundary,
    Simulation,
    SpectralMaterial,
    StructuredMedium,
    SurfaceDetector,
    UniformSurfaceSource,
    ureg,
)


def detector(*, face, axis, angle):
    return SurfaceDetector(
        face=face,
        position=np.array([1, 1, 0 if face == "z_min" else 4]) * ureg.mm,
        radius=0.7 * ureg.mm,
        collection_axis=axis,
        acceptance_angle=angle * ureg.degree,
    )


def test_projected_acceptance_for_aligned_and_tilted_cones():
    for half_angle in (0, 15, 30, 60, 90):
        aperture = detector(face="z_max", axis=[0, 0, 1], angle=half_angle)

        assert aperture.angular_fraction == pytest.approx(
            np.sin(np.deg2rad(half_angle)) ** 2, abs=1e-12
        )

    aperture = detector(face="z_min", axis=[1, 0, -1], angle=30)

    assert aperture.angular_fraction == pytest.approx(
        np.cos(np.pi / 4) * 0.25, abs=1e-12
    )

    # For two intersecting hemispheres the cosine-weighted fraction is (1+cos tilt)/2.
    aperture = detector(face="z_max", axis=[1, 0, 1], angle=90)

    assert aperture.angular_fraction == pytest.approx(
        (1 + np.cos(np.pi / 4)) / 2, abs=1e-10
    )


def test_contact_collection_in_exact_nonabsorbing_slab():
    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.array([2, 2, 4]) * ureg.mm,
    )

    sample.add_background(material=material)

    mesh = CartesianMesh.from_sample(sample=sample, shape=(9, 9, 17))

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

    transmitted = detector(face="z_max", axis=[0, 0, 1], angle=30)

    reflected = detector(face="z_min", axis=[0, 0, -1], angle=90)

    simulation = Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=Diffusion(boundaries=boundaries),
        detectors={"transmitted": transmitted, "reflected": reflected},
    )

    result = simulation.run()

    # Independent solution: Phi is linear, and Jz is constant throughout the slab.
    beta = 0.7 / 2.6

    diffusion_length = 1 / 3000

    incoming = 0.001 / (0.002**2)

    current = (2 / 1.3) * incoming / (2 + beta * 0.004 / diffusion_length)

    measurement = result.detectors["transmitted"]

    expected_transmission = (
        current * measurement.sampled_area.to("m^2").magnitude * 0.25
    )

    assert measurement.collected_power.to("W").magnitude == pytest.approx(
        expected_transmission, rel=1e-11
    )

    reflection = result.detectors["reflected"]

    expected_reflection = (incoming - current) * reflection.sampled_area.to(
        "m^2"
    ).magnitude

    assert reflection.collected_power.to("W").magnitude == pytest.approx(
        expected_reflection, rel=1e-11
    )

    assert measurement.sampled_area.to("mm^2").magnitude == pytest.approx(
        np.pi * 0.7**2, rel=0.08
    )

    for face, field in result.boundary_fields.items():
        assert field.escaping_power.to("W").magnitude == pytest.approx(
            result.outgoing_power_by_face[face].to("W").magnitude, abs=1e-15
        )

    assert result.detector_power["transmitted"] == measurement.collected_power

    # Reading a solved boundary field remains independent of later source edits.
    source.optical_power = 9 * ureg.mW

    repeated = reflected.measure(result=result)

    assert repeated.collected_power == reflection.collected_power

    tiny = SurfaceDetector(
        face="z_max",
        position=np.array([1, 1, 4]) * ureg.mm,
        radius=1e-8 * ureg.mm,
        collection_axis=[0, 0, 1],
        acceptance_angle=90 * ureg.degree,
    )

    with pytest.raises(ValueError, match="resolve the detector"):
        tiny.measure(result=result)

    outside = SurfaceDetector(
        face="z_max",
        position=np.array([0.1, 1, 4]) * ureg.mm,
        radius=0.7 * ureg.mm,
        collection_axis=[0, 0, 1],
        acceptance_angle=90 * ureg.degree,
    )

    with pytest.raises(ValueError, match="entire detector"):
        outside.measure(result=result)


@pytest.mark.parametrize(
    "axis,angle",
    [([0, 0, 0], 30), ([0, 0, -1], 30), ([1, 0, 0], 30), ([0, 0, 1], 91)],
)
def test_invalid_collection_configuration(axis, angle):
    with pytest.raises(ValueError):
        detector(face="z_max", axis=axis, angle=angle)


def test_acceptance_requires_an_explicit_quantity():
    with pytest.raises(TypeError, match="carry units"):
        SurfaceDetector(
            face="z_max",
            position=np.ones(3) * ureg.mm,
            radius=0.5 * ureg.mm,
            collection_axis=[0, 0, 1],
            acceptance_angle=30,
        )


def test_spectral_slab_detector_signals_match_constant_current_solutions():
    material = SpectralMaterial(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        wavelengths=np.array([700, 900]) * ureg.nm,
        absorption=np.zeros(2) / ureg.mm,
        reduced_scattering=np.array([0.5, 1.5]) / ureg.mm,
        refractive_index=1.4,
    )

    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.array([2, 2, 4]) * ureg.mm,
    )

    sample.add_background(material=material)

    mesh = CartesianMesh.from_sample(sample=sample, shape=(5, 5, 9))

    source = UniformSurfaceSource(
        face="z_min",
        wavelength=700 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    boundaries = {
        face: ReflectingBoundary()
        for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
    }

    boundaries["z_min"] = RobinBoundary(effective_reflectance=0.3)

    boundaries["z_max"] = RobinBoundary(effective_reflectance=0.3)

    aperture = detector(face="z_max", axis=[0, 0, 1], angle=90)

    simulation = Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=Diffusion(boundaries=boundaries),
        detectors={"transmission": aperture},
    )

    spectrum = simulation.sweep(wavelengths=np.array([700, 800, 900]) * ureg.nm)

    beta = 0.7 / 2.6

    measured = spectrum.detector_power["transmission"].to("W").magnitude

    expected = []

    for result, scattering in zip(spectrum.results, [500, 1000, 1500]):
        diffusion_length = 1 / (3 * scattering)

        current = (2 / 1.3) * 250 / (2 + beta * 0.004 / diffusion_length)

        expected.append(
            current * result.detectors["transmission"].sampled_area.to("m^2").magnitude
        )

    np.testing.assert_allclose(measured, expected, rtol=1e-11)

    assert source.wavelength == 700 * ureg.nm

    assert sample.background is material

    assert all(
        result.detectors["transmission"].detector is not aperture
        for result in spectrum.results
    )


@pytest.mark.parametrize("face", MomentRT.boundary.FACES)
def test_full_hemisphere_detector_on_each_face(face):
    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.full(3, 4) * ureg.mm,
    )

    aperture = SurfaceDetector.full_hemisphere_on_face(
        sample=sample,
        face=face,
        coordinates=np.array([2, 2]) * ureg.mm,
        radius=0.5 * ureg.mm,
    )

    axis = "xyz".index(face[0])

    expected = np.zeros(3)

    expected[axis] = -1 if face.endswith("min") else 1

    np.testing.assert_array_equal(aperture.collection_axis, expected)

    assert aperture.position[axis].to("mm").magnitude == (
        0 if face.endswith("min") else 4
    )

    assert aperture.angular_fraction == pytest.approx(1)

    assert aperture.acceptance_angle.to("degree").magnitude == 90
