"""Independent normalization, moment, and solver checks for phase functions."""

import numpy as np
import pytest

import MomentRT


@pytest.mark.parametrize(
    "phase,anisotropy,expected",
    [
        (MomentRT.Isotropic(), None, [1, 0, 0, 0]),
        (MomentRT.Rayleigh(), None, [1, 0, 0.1, 0]),
        (MomentRT.HenyeyGreenstein(), 0.7, [1, 0.7, 0.49, 0.343]),
    ],
)
def test_density_and_moments_match_independent_spherical_integration(
    phase, anisotropy, expected
):
    mu, weights = np.polynomial.legendre.leggauss(160)

    density = phase.evaluate(cosine=mu, anisotropy=anisotropy)

    polynomials = np.array(
        [np.ones_like(mu), mu, (3 * mu**2 - 1) / 2, (5 * mu**3 - 3 * mu) / 2]
    )

    integrated = 2 * np.pi * np.sum(polynomials * density * weights, axis=1)

    np.testing.assert_allclose(integrated, expected, atol=1e-12)

    np.testing.assert_allclose(
        phase.moments(anisotropy=anisotropy), expected, atol=1e-12
    )


def test_tabulated_phase_integrates_piecewise_linear_density():
    mu = np.array([-1, -0.4, 0.25, 1])

    phase = MomentRT.phase_function.Tabulated(cosine=mu, values=1 + 0.6 * mu)

    np.testing.assert_allclose(phase.moments(), [1, 0.2, 0, 0], atol=1e-14)

    np.testing.assert_allclose(
        phase.evaluate(cosine=[-1, 0, 1]), np.array([0.4, 1, 1.6]) / (4 * np.pi)
    )

    with pytest.raises(ValueError, match="first moment"):
        phase.moments(anisotropy=0.5)


@pytest.mark.parametrize(
    "cosine,values",
    [
        ([-0.9, 1], [1, 1]),
        ([-1, 0, 0, 1], [1, 1, 1, 1]),
        ([-1, 1], [0, 0]),
        ([-1, 1], [1, -1]),
    ],
)
def test_tabulated_phase_rejects_invalid_data(cosine, values):
    with pytest.raises(ValueError):
        MomentRT.phase_function.Tabulated(cosine=cosine, values=values)


def simulation(*, phase, solver_type, anisotropy=None):
    ureg = MomentRT.ureg

    material = MomentRT.Material(
        absorption=0.1 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
        anisotropy=anisotropy,
        phase_function=phase,
    )

    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * ureg.mm, size=np.ones(3) * 2 * ureg.mm
    )

    sample.add_background(material=material)

    mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=(3, 3, 3))

    source = MomentRT.UniformSurfaceSource(
        face="z_min", wavelength=850 * ureg.nm, optical_power=1 * ureg.mW
    )

    boundary = MomentRT.RobinBoundary(effective_reflectance=0)

    solver = solver_type(boundaries=MomentRT.Boundaries.uniform(condition=boundary))

    return MomentRT.Simulation(medium=sample, mesh=mesh, source=source, solver=solver)


@pytest.mark.parametrize("solver_type", [MomentRT.SP3, MomentRT.P3])
def test_solver_uses_higher_phase_moments(solver_type):
    isotropic = simulation(phase=MomentRT.Isotropic(), solver_type=solver_type)

    hg_zero = simulation(
        phase=MomentRT.HenyeyGreenstein(), anisotropy=0, solver_type=solver_type
    )

    rayleigh = simulation(phase=MomentRT.Rayleigh(), solver_type=solver_type)

    assert isotropic.validate().valid

    isotropic_field = isotropic.run().fluence.values.magnitude

    hg_field = hg_zero.run().fluence.values.magnitude

    rayleigh_field = rayleigh.run().fluence.values.magnitude

    np.testing.assert_allclose(isotropic_field, hg_field, rtol=1e-12)

    assert not np.allclose(rayleigh_field, isotropic_field, rtol=1e-5)


def test_p1_needs_only_reduced_scattering():
    isotropic = simulation(phase=MomentRT.Isotropic(), solver_type=MomentRT.Diffusion)

    rayleigh = simulation(phase=MomentRT.Rayleigh(), solver_type=MomentRT.Diffusion)

    np.testing.assert_allclose(
        isotropic.run().fluence.values.magnitude,
        rayleigh.run().fluence.values.magnitude,
    )


@pytest.mark.parametrize("solver_type", [MomentRT.SP3, MomentRT.P3])
def test_tabulated_phase_drives_solvers_without_explicit_anisotropy(solver_type):
    phase = MomentRT.phase_function.Tabulated(cosine=[-1, 0, 1], values=[0.4, 1, 1.6])

    tabulated = simulation(phase=phase, solver_type=solver_type)

    hg = simulation(
        phase=MomentRT.HenyeyGreenstein(), anisotropy=0.2, solver_type=solver_type
    )

    assert tabulated.validate().valid

    tabulated_result = tabulated.run()

    hg_result = hg.run()

    assert not np.allclose(
        tabulated_result.fluence.values.magnitude,
        hg_result.fluence.values.magnitude,
        rtol=1e-5,
    )

    assert (
        abs(tabulated_result.energy_balance["residual_power"].to("W").magnitude) < 1e-12
    )


def test_spectral_phase_preserves_consistency():
    ureg = MomentRT.ureg

    phase = MomentRT.Rayleigh()

    material = MomentRT.SpectralMaterial(
        wavelengths=np.array([800, 900]) * ureg.nm,
        absorption=np.array([0.1, 0.1]) / ureg.mm,
        reduced_scattering=np.array([1, 1]) / ureg.mm,
        refractive_index=1.4,
        phase_function=phase,
    )

    resolved = material.at(wavelength=850 * ureg.nm)

    np.testing.assert_allclose(resolved._phase_moments(), [1, 0, 0.1, 0])

    model = simulation(phase=phase, solver_type=MomentRT.SP3)

    model.medium.add_background(material=material)

    spectrum = model.sweep(wavelengths=np.array([800, 850, 900]) * ureg.nm)

    assert spectrum.diagnostics["geometry_preparations"] == 1

    np.testing.assert_allclose(
        spectrum.results[0].fluence.values.magnitude,
        spectrum.results[-1].fluence.values.magnitude,
        rtol=1e-12,
    )

    with pytest.raises(ValueError, match="anisotropy"):
        MomentRT.Material(
            absorption=0.1 / ureg.mm,
            reduced_scattering=1 / ureg.mm,
            refractive_index=1.4,
            anisotropy=0.8,
            phase_function=phase,
        )


def test_solver_rejects_unrealizable_custom_moments():
    class InvalidMoments(MomentRT.PhaseFunction):
        def evaluate(self, *, cosine, anisotropy=None):
            return np.ones_like(cosine, dtype=float) / (4 * np.pi)

        def moments(self, *, anisotropy=None):
            return np.array([1, 0, -1, 0])

    model = simulation(phase=InvalidMoments(), solver_type=MomentRT.SP3)

    assert not model.validate().valid

    with pytest.raises(ValueError, match="nonnegative scattering distribution"):
        model.run()


def test_material_requires_explicit_phase_function_object():
    arguments = {
        "absorption": 0.1 / MomentRT.ureg.mm,
        "reduced_scattering": 1 / MomentRT.ureg.mm,
        "refractive_index": 1.4,
    }

    with pytest.raises(TypeError, match="phase_function"):
        MomentRT.Material(**arguments)

    with pytest.raises(ValueError, match="PhaseFunction object"):
        MomentRT.Material(**arguments, phase_function="henyey_greenstein")

    material = MomentRT.Material(
        **arguments,
        anisotropy=0.7,
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )

    np.testing.assert_allclose(
        material.phase_function.moments(anisotropy=material.anisotropy),
        [1, 0.7, 0.49, 0.343],
    )


def test_spectral_material_requires_explicit_phase_function_object():
    with pytest.raises(TypeError, match="phase_function"):
        MomentRT.SpectralMaterial(
            wavelengths=np.array([800, 900]) * MomentRT.ureg.nm,
            absorption=np.array([0.1, 0.1]) / MomentRT.ureg.mm,
            reduced_scattering=np.array([1, 1]) / MomentRT.ureg.mm,
            refractive_index=1.4,
        )
