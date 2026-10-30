"""Independent angular and radiative-transfer references for full P3."""

import numpy as np
import pytest
from scipy.linalg import expm
import MomentRT


def cube(*, count, source_face="z_min", absorption=1.0, scattering=1e-7):
    ureg = MomentRT.ureg

    material = MomentRT.Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=absorption / ureg.mm,
        reduced_scattering=scattering / ureg.m,
        refractive_index=1.4,
        anisotropy=0.0,
    )

    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.ones(3) * 2 * ureg.mm,
    )

    sample.add_background(material=material)

    mesh = MomentRT.CartesianMesh.from_sample(
        sample=sample, shape=(count, count, count)
    )

    source = MomentRT.UniformSurfaceSource(
        face=source_face,
        wavelength=850 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    boundary = MomentRT.RobinBoundary(effective_reflectance=0)

    boundaries = MomentRT.Boundaries.uniform(condition=boundary)

    solver = MomentRT.P3(boundaries=boundaries)

    return MomentRT.Simulation(medium=sample, mesh=mesh, source=source, solver=solver)


def test_planar_angular_reduction_has_exact_p3_coefficients():
    """Compare all planar diffusion and boundary coefficients with P3 algebra."""
    streaming, _, boundary, source = MomentRT.transport.p3._angular_operators()

    odd_sigma = np.array(
        [1 if degree == 1 else 2 for degree, _ in MomentRT.transport.p3._ODD]
    )

    diffusion = streaming[2].T @ np.diag(1 / odd_sigma) @ streaming[2]

    # Azimuthally symmetric l=0,2 moments, transformed to the conventional
    # P3 composite variables u = phi_0 + 2 phi_2 and v = 3 phi_2.
    indices = [0, 3]

    transform = np.array([[1, -2 / 3], [0, np.sqrt(5) / 3]]) / np.sqrt(4 * np.pi)

    planar_diffusion = (
        4 * np.pi * transform.T @ diffusion[np.ix_(indices, indices)] @ transform
    )

    planar_boundary = (
        4
        * np.pi
        * transform.T
        @ boundary["z_min"][np.ix_(indices, indices)]
        @ transform
    )

    planar_source = 4 * np.pi * transform.T @ source["z_min"][indices]

    np.testing.assert_allclose(planar_diffusion, np.diag([1 / 3, 1 / 14]), atol=1e-14)

    np.testing.assert_allclose(
        planar_boundary, [[1 / 2, -1 / 8], [-1 / 8, 7 / 24]], atol=1e-14
    )

    np.testing.assert_allclose(planar_source, [1, -1 / 4], atol=1e-14)


def test_planar_p3_fem_against_independent_transfer_solution():
    """Compare the reduced angular operator to an analytical homogeneous slab."""
    streaming, _, boundary, source = MomentRT.transport.p3._angular_operators()

    mu_a, mu_s1, anisotropy = 20.0, 1000.0, 0.7

    mu_1 = mu_a + mu_s1

    mu_2 = mu_a + mu_s1 * (1 + anisotropy)

    mu_3 = mu_a + mu_s1 * (1 + anisotropy + anisotropy**2)

    odd_sigma = np.array(
        [mu_1 if degree == 1 else mu_3 for degree, _ in MomentRT.transport.p3._ODD]
    )

    indices = [0, 3]

    transform = np.array([[1, -2 / 3], [0, np.sqrt(5) / 3]]) / np.sqrt(4 * np.pi)

    angular_diffusion = streaming[2].T @ np.diag(1 / odd_sigma) @ streaming[2]

    diffusion = (
        4
        * np.pi
        * transform.T
        @ angular_diffusion[np.ix_(indices, indices)]
        @ transform
    )

    reaction = 4 * np.pi * transform.T @ np.diag([mu_a, mu_2]) @ transform

    face = (
        4
        * np.pi
        * transform.T
        @ boundary["z_min"][np.ix_(indices, indices)]
        @ transform
    )

    drive = 4 * np.pi * transform.T @ source["z_min"][indices]

    generator = np.array(
        [
            [0, 0, -3 * mu_1, 0],
            [0, 0, 0, -7 * mu_3],
            [-mu_a, 2 * mu_a / 3, 0, 0],
            [2 * mu_a / 3, -(4 * mu_a + 5 * mu_2) / 9, 0, 0],
        ]
    )

    upper = np.array([[0.5, -0.125, 1, 0], [-0.125, 7 / 24, 0, 1]])

    lower = np.array([[-0.5, 0.125, 1, 0], [0.125, -7 / 24, 0, 1]])

    incoming = 250.0

    state = np.linalg.solve(
        np.vstack((upper, lower @ expm(generator * 0.004))),
        [2 * incoming, -incoming / 2, 0, 0],
    )

    errors = []

    for count in (5, 9, 17):
        depth = np.linspace(0, 0.004, count)

        spacing = depth[1] - depth[0]

        stiffness = np.zeros((count, count))

        mass = np.zeros((count, count))

        for index in range(count - 1):
            stiffness[index : index + 2, index : index + 2] += (
                np.array([[1, -1], [-1, 1]]) / spacing
            )

            mass[index : index + 2, index : index + 2] += (
                spacing * np.array([[2, 1], [1, 2]]) / 6
            )

        endpoints = np.zeros((count, count))

        endpoints[0, 0] = endpoints[-1, -1] = 1

        matrix = (
            np.kron(diffusion, stiffness)
            + np.kron(reaction, mass)
            + np.kron(face, endpoints)
        )

        rhs = np.zeros((2, count))

        rhs[:, 0] = 2 * incoming * drive

        values = np.linalg.solve(matrix, rhs.ravel()).reshape(2, count)

        actual = values[0] - 2 * values[1] / 3

        states = np.array([expm(generator * z) @ state for z in depth])

        expected = states[:, 0] - 2 * states[:, 1] / 3

        errors.append(np.linalg.norm(actual - expected) / np.linalg.norm(expected))

    assert errors[2] < errors[1] < errors[0]

    assert errors[-1] < 0.01


def test_3d_p3_against_independent_ballistic_ray_integral():
    """Compare cube-center fluence with the near-zero-scattering RTE limit."""
    # Integrate Lambertian radiance backwards from the cube center. A ray
    # reaches the illuminated upper face first when its z cosine exceeds its
    # absolute x and y cosines. Optical depth to that face is 1 / mu.
    mu, weights = np.polynomial.legendre.leggauss(200)

    mu = (mu[:, None] + 1) / 2

    weights = weights[:, None] / 2

    azimuth = (2 * np.pi * np.arange(800) / 800)[None, :]

    transverse = np.sqrt(1 - mu**2)

    selected = (mu >= np.abs(transverse * np.cos(azimuth))) & (
        mu >= np.abs(transverse * np.sin(azimuth))
    )

    incoming_irradiance = 0.001 / 0.002**2

    reference = (
        incoming_irradiance
        / np.pi
        * np.sum(weights * (2 * np.pi / 800) * np.exp(-1 / mu) * selected)
    )

    errors = []

    for count in (3, 5, 7):
        result = cube(count=count).run()

        center = result.fluence.values.to("W/m^2").magnitude[
            count // 2, count // 2, count // 2
        ]

        errors.append(abs(center - reference) / reference)

        assert abs(result.energy_balance["residual_power"].to("W").magnitude) < 1e-12

        assert len(result.moments) == 6

    assert errors[2] < errors[1] < errors[0]

    # This tolerance is a model-comparison bound for this specific ballistic
    # cube; it is not a universal accuracy claim for P3.
    assert errors[-1] < 0.13


def test_p3_rotational_symmetry_and_spectral_reuse():
    z_simulation = cube(count=5, source_face="z_min", absorption=0.1, scattering=1000)

    x_simulation = cube(count=5, source_face="x_min", absorption=0.1, scattering=1000)

    z_result = z_simulation.run()

    x_result = x_simulation.run()

    np.testing.assert_allclose(
        z_result.fluence.values.magnitude.transpose(2, 1, 0),
        x_result.fluence.values.magnitude,
        rtol=1e-11,
        atol=1e-10,
    )

    spectrum = z_simulation.sweep(wavelengths=np.array([800, 900]) * MomentRT.ureg.nm)

    assert spectrum.diagnostics["geometry_preparations"] == 1

    assert spectrum.results[1].diagnostics["preparation_reused"]

    np.testing.assert_allclose(
        spectrum.results[0].fluence.values.magnitude, z_result.fluence.values.magnitude
    )


def test_p3_rejects_unsupported_boundaries_and_missing_anisotropy():
    simulation = cube(count=3)

    bad_boundaries = dict(simulation.solver.boundaries)

    bad_boundaries["x_min"] = MomentRT.ReflectingBoundary()

    with pytest.raises(ValueError, match="zero-reflectance"):
        MomentRT.P3(boundaries=bad_boundaries)

    simulation.medium.background.anisotropy = None

    assert not simulation.validate().valid

    with pytest.raises(ValueError, match="anisotropy"):
        simulation.run()


def test_p3_detector_integrates_boundary_quadrature():
    simulation = cube(count=5, absorption=0.1, scattering=1000)

    detector = MomentRT.SurfaceDetector.on_face(
        sample=simulation.medium,
        face="z_max",
        coordinates=np.array([1, 1]) * MomentRT.ureg.mm,
        radius=0.75 * MomentRT.ureg.mm,
        acceptance_angle=90 * MomentRT.ureg.degree,
        orientation="outward_normal",
    )

    simulation.detectors["sensor"] = detector

    report = simulation.validate()

    assert report.valid

    result = simulation.run()

    field = result.boundary_fields["z_max"]

    displacement = field.positions.to("mm").magnitude[..., :2] - [1, 1]

    selected = np.linalg.norm(displacement, axis=-1) <= 0.75

    reference = (
        field.escaping_irradiance[selected] * field.area_weights[selected]
    ).sum()

    np.testing.assert_allclose(
        result.detector_power["sensor"].to("W").magnitude,
        reference.to("W").magnitude,
        rtol=1e-13,
    )

    assert result.detectors["sensor"].sampled_area.to("m^2").magnitude == pytest.approx(
        report.detector_sampled_areas["sensor"].to("m^2").magnitude
    )
