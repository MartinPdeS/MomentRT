"""Independent matrix-exponential slab reference for SP3."""

import numpy as np
import pytest
from scipy.linalg import expm

import MomentRT


def configuration(*, count, anisotropy=0.7, layered=False, layer_anisotropy=0.4):
    ureg = MomentRT.ureg

    material = MomentRT.Material(
        absorption=0.02 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
        anisotropy=anisotropy,
    )

    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.array([2, 2, 4]) * ureg.mm,
    )

    sample.add_background(material=material)

    if layered:
        layer_material = MomentRT.Material(
            absorption=0.04 / ureg.mm,
            reduced_scattering=2 / ureg.mm,
            refractive_index=1.4,
            anisotropy=layer_anisotropy,
        )

        sample.add_layer(name="upper", thickness=2 * ureg.mm, material=layer_material)

    boundary = MomentRT.ReflectingBoundary()

    boundaries = {face: boundary for face in MomentRT.boundary.FACES}

    boundaries["z_min"] = MomentRT.RobinBoundary(effective_reflectance=0)

    boundaries["z_max"] = MomentRT.RobinBoundary(effective_reflectance=0)

    source = MomentRT.UniformSurfaceSource(
        face="z_min", wavelength=970 * ureg.nm, optical_power=1 * ureg.mW
    )

    detector = MomentRT.SurfaceDetector(
        face="z_max",
        position=np.array([1, 1, 4]) * ureg.mm,
        radius=0.8 * ureg.mm,
        collection_axis=np.array([0, 0, 1]),
        acceptance_angle=90 * ureg.degree,
    )

    mesh = MomentRT.CartesianMesh.from_sample(
        sample=sample, shape=(count, count, 2 * count - 1)
    )

    solver = MomentRT.SP3(boundaries=boundaries)

    return MomentRT.Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=solver,
        detectors={"sensor": detector},
    )


def exact_slab(*, depths, anisotropy, layered=False, layer_anisotropy=0.4):
    """Solve first-order moment ODEs and both Marshak boundary conditions."""

    def generator(*, mu_a, mu_s1, anisotropy):
        mu_1 = mu_a + mu_s1

        mu_2 = mu_a + mu_s1 * (1 + anisotropy)

        mu_3 = mu_a + mu_s1 * (1 + anisotropy + anisotropy**2)

        return np.array(
            [
                [0, 0, -3 * mu_1, 0],
                [0, 0, 0, -7 * mu_3],
                [-mu_a, 2 * mu_a / 3, 0, 0],
                [2 * mu_a / 3, -(4 * mu_a + 5 * mu_2) / 9, 0, 0],
            ]
        )

    lower = generator(mu_a=20.0, mu_s1=1000.0, anisotropy=anisotropy)

    upper = (
        generator(mu_a=40.0, mu_s1=2000.0, anisotropy=layer_anisotropy)
        if layered
        else lower
    )

    def transfer(*, depth):
        if depth <= 0.002:
            return expm(upper * depth)

        return expm(lower * (depth - 0.002)) @ expm(upper * 0.002)

    incoming = 0.001 / 0.002**2

    bottom = np.array([[-0.5, 0.125, 1, 0], [0.125, -7 / 24, 0, 1]])

    # At z=0, j = -B[u,v] + [2I,-I/2]; at z=L, j = B[u,v].
    top = np.array([[0.5, -0.125, 1, 0], [-0.125, 7 / 24, 0, 1]])

    constraints = np.vstack((top, bottom @ transfer(depth=0.004)))

    initial = np.linalg.solve(
        constraints, np.array([2 * incoming, -incoming / 2, 0, 0])
    )

    states = np.array([transfer(depth=z) @ initial for z in depths])

    return states[:, 0] - 2 * states[:, 1] / 3, states[-1, 2]


@pytest.mark.parametrize("layered", [False, True])
def test_sp3_slab_converges_to_independent_reference(layered):
    errors = []

    for count in (3, 5, 9):
        simulation = configuration(count=count, layered=layered)

        result = simulation.run()

        depths = simulation.mesh.points.to("m").magnitude[0, 0, :, 2]

        expected, outgoing = exact_slab(depths=depths, anisotropy=0.7, layered=layered)

        actual = result.fluence.values.to("W/m^2").magnitude[0, 0, :]

        errors.append(np.linalg.norm(actual - expected) / np.linalg.norm(expected))

        assert result.outgoing_power_by_face["z_max"].to(
            "W"
        ).magnitude == pytest.approx(outgoing * 0.002**2, rel=0.08)

        assert abs(result.energy_balance["residual_power"].to("W").magnitude) < 1e-12

        assert result.detectors["sensor"].collected_power.magnitude > 0

    assert errors[2] < errors[1] < errors[0]

    assert errors[-1] < 0.01


def test_sp3_rejects_unavailable_boundary_reflection_model():
    simulation = configuration(count=3)

    simulation.solver.boundaries["z_max"] = MomentRT.RobinBoundary(
        effective_reflectance=0.2
    )

    with pytest.raises(ValueError, match="zero-reflectance"):
        MomentRT.SP3(boundaries=simulation.solver.boundaries)

    with pytest.raises(ValueError, match="zero-reflectance"):
        simulation.run()


def test_sp3_sweep_reuses_geometry():
    simulation = configuration(count=3)

    spectrum = simulation.sweep(wavelengths=np.array([800, 900]) * MomentRT.ureg.nm)

    assert spectrum.diagnostics["geometry_preparations"] == 1

    assert spectrum.results[1].diagnostics["preparation_reused"] is True

    assert spectrum.detector_power["sensor"].shape == (2,)


@pytest.mark.parametrize(
    "anisotropy",
    [-1, 1, float("nan"), True, np.bool_(False), 0.5 * MomentRT.ureg.dimensionless],
)
def test_sp3_requires_physical_anisotropy(anisotropy):
    with pytest.raises(ValueError, match="anisotropy"):
        MomentRT.Material(
            absorption=0.02 / MomentRT.ureg.mm,
            reduced_scattering=1 / MomentRT.ureg.mm,
            refractive_index=1.4,
            anisotropy=anisotropy,
        )


def test_sp3_requires_anisotropy_for_each_region_and_p1_does_not():
    simulation = configuration(count=3, layered=True)

    simulation.medium.structures[0].material.anisotropy = None

    with pytest.raises(ValueError, match="anisotropy on every material"):
        simulation.run()

    simulation.solver = MomentRT.Diffusion(boundaries=simulation.solver.boundaries)

    assert simulation.run().fluence.values.magnitude.size > 0


def test_sp3_has_no_solver_anisotropy_argument():
    simulation = configuration(count=3)

    with pytest.raises(TypeError, match="anisotropy"):
        MomentRT.SP3(boundaries=simulation.solver.boundaries, anisotropy=0.7)


def test_sp3_spectral_anisotropy_interpolates_and_changes_solution():
    ureg = MomentRT.ureg

    material = MomentRT.SpectralMaterial(
        wavelengths=np.array([800, 1000]) * ureg.nm,
        absorption=np.array([0.02, 0.02]) / ureg.mm,
        reduced_scattering=np.array([1.0, 1.0]) / ureg.mm,
        refractive_index=1.4,
        anisotropy=[0.2, 0.8],
    )

    assert material.at(wavelength=900 * ureg.nm).anisotropy == pytest.approx(0.5)

    simulation = configuration(count=5)

    simulation.medium.add_background(material=material)

    wavelengths = np.array([800, 900, 1000]) * ureg.nm

    spectrum = simulation.sweep(wavelengths=wavelengths)

    for index, anisotropy in enumerate([0.2, 0.5, 0.8]):
        reference = configuration(count=5, anisotropy=anisotropy).run()

        np.testing.assert_allclose(
            spectrum.results[index].fluence.values.magnitude,
            reference.fluence.values.magnitude,
            rtol=1e-12,
        )

    assert not np.allclose(
        spectrum.results[0].fluence.values.magnitude,
        spectrum.results[-1].fluence.values.magnitude,
    )


@pytest.mark.parametrize(
    "anisotropy",
    [
        [0.2],
        [0.1, 1],
        [0.1, float("nan")],
        [True, False],
        [0.1, 0.2] * MomentRT.ureg.dimensionless,
    ],
)
def test_spectral_anisotropy_rejects_invalid_tables(anisotropy):
    with pytest.raises(ValueError, match="anisotropy"):
        MomentRT.SpectralMaterial(
            wavelengths=np.array([800, 1000]) * MomentRT.ureg.nm,
            absorption=np.array([0.02, 0.02]) / MomentRT.ureg.mm,
            reduced_scattering=np.array([1, 1]) / MomentRT.ureg.mm,
            refractive_index=1.4,
            anisotropy=anisotropy,
        )


def test_detector_convergence_accepts_sp3():
    simulation = configuration(count=3)

    study = simulation.convergence_study(
        mesh_shapes=[(3, 3, 5), (5, 5, 9), (7, 7, 13)],
        area_rtol=0.8,
        power_rtol=0.8,
        required_consecutive=1,
    )

    assert len(study.results) == 3

    assert study.detector_power["sensor"].shape == (3,)
