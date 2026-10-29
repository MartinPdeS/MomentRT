"""Tabulated optics, independent Green functions, and nonmutating sweeps."""

import numpy as np
import pytest

import MomentRT
from MomentRT import (
    CartesianMesh,
    HomogeneousMedium,
    InfiniteDiffusion,
    Material,
    PointSource,
    Simulation,
    SpectralMaterial,
    SpectralSimulationResult,
    StructuredMedium,
    ureg,
)


def table():
    return SpectralMaterial(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        wavelengths=np.array([700, 800, 900]) * ureg.nm,
        absorption=np.array([10, 20, 40]) / ureg.m,
        reduced_scattering=np.array([1.2, 1.0, 0.8]) / ureg.mm,
        refractive_index=np.array([1.42, 1.41, 1.4]),
    )


def test_spectral_interpolation_and_coverage():
    spectral = table()

    material = spectral.at(wavelength=0.75 * ureg.micrometer)

    assert material.absorption.to("1/m").magnitude == pytest.approx(15)

    assert material.reduced_scattering.to("1/mm").magnitude == pytest.approx(1.1)

    assert material.refractive_index == pytest.approx(1.415)

    assert spectral.at(wavelength=900 * ureg.nm).absorption == 40 / ureg.m

    for wavelength in (699, 901):
        with pytest.raises(ValueError, match="outside"):
            spectral.at(wavelength=wavelength * ureg.nm)

    with pytest.raises(TypeError, match="carry units"):
        spectral.at(wavelength=750)

    with pytest.raises(ValueError, match="shape"):
        spectral.at(wavelength=np.array([750, 800]) * ureg.nm)


@pytest.mark.parametrize(
    "wavelengths,absorption,scattering,refractive_index",
    [
        ([800, 700], [1, 2], [1, 2], 1.4),
        ([700, 700], [1, 2], [1, 2], 1.4),
        ([700, 800], [-1, 2], [1, 2], 1.4),
        ([700, 800], [1, 2], [0, 2], 1.4),
        ([700, 800], [1, 2], [1, 2], [1.4]),
        ([700, 800], [1, 2], [1, 2], True),
    ],
)
def test_invalid_spectral_tables(wavelengths, absorption, scattering, refractive_index):
    with pytest.raises(ValueError):
        SpectralMaterial(
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
            wavelengths=np.array(wavelengths) * ureg.nm,
            absorption=np.array(absorption) / ureg.mm,
            reduced_scattering=np.array(scattering) / ureg.mm,
            refractive_index=refractive_index,
        )


def test_sweep_matches_green_functions_and_preserves_configuration():
    spectral = table()

    medium = HomogeneousMedium(material=spectral)

    source = PointSource(
        position=np.zeros(3) * ureg.mm,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    mesh = CartesianMesh(
        shape=(4, 4, 4),
        spacing=np.ones(3) * ureg.mm,
        origin=np.full(3, 0.5) * ureg.mm,
    )

    simulation = Simulation(
        medium=medium, source=source, mesh=mesh, solver=InfiniteDiffusion()
    )

    wavelengths = np.array([850, 750, 800]) * ureg.nm

    powers = np.array([1, 2, 3]) * ureg.mW

    spectrum = simulation.sweep(wavelengths=wavelengths, optical_power=powers)

    radii = np.linalg.norm(mesh.points.to("m").magnitude, axis=-1)

    # Explicit values independently interpolated from the supplied table.
    for result, absorption, scattering, power in zip(
        spectrum.results, [30, 15, 20], [900, 1100, 1000], [0.001, 0.002, 0.003]
    ):
        diffusion_length = 1 / (3 * (absorption + scattering))

        exact = (
            power
            * np.exp(-np.sqrt(absorption / diffusion_length) * radii)
            / (4 * np.pi * diffusion_length * radii)
        )

        np.testing.assert_allclose(
            result.fluence.values.to("W/m^2").magnitude, exact, rtol=1e-12
        )

    np.testing.assert_array_equal(spectrum.wavelengths.magnitude, [850, 750, 800])

    assert source.wavelength == 800 * ureg.nm

    assert source.optical_power == 1 * ureg.mW

    assert simulation.mesh is mesh

    assert medium.material is spectral

    assert spectrum.detector_power == {}

    baseline = simulation.sweep(wavelengths=np.array([800, 800]) * ureg.nm)

    np.testing.assert_array_equal(
        baseline.results[0].fluence.values.magnitude,
        baseline.results[1].fluence.values.magnitude,
    )


def test_spectral_sample_snapshots_and_optical_sections():
    spectral = table()

    fixed = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=1 / ureg.m,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.ones(3) * 4 * ureg.mm,
    )

    sample.add_background(material=spectral)

    sample.add_layer(name="fixed", thickness=1 * ureg.mm, material=fixed)

    resolved = sample.at(wavelength=750 * ureg.nm)

    assert resolved.background.absorption.to("1/m").magnitude == pytest.approx(15)

    assert resolved.structures[0].material is fixed

    assert sample.background is spectral

    resolved.structures[0].thickness *= 2

    assert sample.structures[0].thickness == 1 * ureg.mm

    assert resolved.structures[0].thickness == 2 * ureg.mm

    section = sample.section(
        normal="z", position=2 * ureg.mm, resolution=(4, 4), wavelength=800 * ureg.nm
    )

    np.testing.assert_allclose(section.fields["absorption"].to("1/m").magnitude, 20)

    with pytest.raises(ValueError, match="Supply wavelength"):
        sample.section(normal="z", position=2 * ureg.mm)


def test_sweep_rejects_invalid_vectors_and_power_units():
    medium = HomogeneousMedium(material=table())

    source = PointSource(
        position=np.zeros(3) * ureg.mm,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    mesh = CartesianMesh(
        shape=(2, 2, 2),
        spacing=np.ones(3) * ureg.mm,
        origin=np.full(3, 0.5) * ureg.mm,
    )

    simulation = Simulation(
        medium=medium,
        source=source,
        mesh=mesh,
        solver=InfiniteDiffusion(),
    )

    with pytest.raises(ValueError, match="nonempty vector"):
        simulation.sweep(wavelengths=800 * ureg.nm)

    with pytest.raises(ValueError, match="nonempty vector"):
        simulation.sweep(wavelengths=np.array([]) * ureg.nm)

    with pytest.raises(ValueError, match="shaped like"):
        simulation.sweep(
            wavelengths=np.array([700, 800]) * ureg.nm,
            optical_power=np.ones(3) * ureg.mW,
        )

    with pytest.raises(ValueError, match="outside"):
        simulation.sweep(wavelengths=np.array([800, 950]) * ureg.nm)


def test_spectral_result_rejects_empty_and_mismatched_data():
    with pytest.raises(ValueError, match="nonempty vector"):
        SpectralSimulationResult(
            wavelengths=np.array([]) * ureg.nm,
            results=[],
        )

    with pytest.raises(ValueError, match="one SimulationResult"):
        SpectralSimulationResult(
            wavelengths=np.array([800]) * ureg.nm,
            results=[],
        )
