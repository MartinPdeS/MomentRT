"""Independent one-dimensional slab reference for the 3D FEM solution."""

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
    StructuredMedium,
    UniformSurfaceSource,
    ureg,
)


def slab(*, count, reflectance, layered=False, absorption=0.02):
    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=absorption / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm, size=np.array([2, 2, 4]) * ureg.mm
    )

    sample.add_background(material=material)

    if layered:
        upper_material = Material(
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
            absorption=0.04 / ureg.mm,
            reduced_scattering=2 / ureg.mm,
            refractive_index=1.4,
        )

        sample.add_layer(name="upper", thickness=2 * ureg.mm, material=upper_material)

    boundaries = {
        face: ReflectingBoundary()
        for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
    }

    boundaries["z_min"] = RobinBoundary(effective_reflectance=reflectance)

    boundaries["z_max"] = RobinBoundary(effective_reflectance=reflectance)

    source = UniformSurfaceSource(
        face="z_min", wavelength=970 * ureg.nm, optical_power=1 * ureg.mW
    )

    mesh = CartesianMesh.from_sample(sample=sample, shape=(count, count, 2 * count - 1))

    solver = Diffusion(boundaries=boundaries)

    simulation = Simulation(medium=sample, source=source, mesh=mesh, solver=solver)

    return simulation, simulation.run()


def transfer(*, distance, absorption, scattering):
    # Independent exact transfer of [fluence, z-directed current] in a slab.
    diffusion_length = 1 / (3 * (absorption + scattering))

    attenuation = np.sqrt(absorption / diffusion_length)

    hyperbolic_cosine = np.cosh(attenuation * distance)

    hyperbolic_sine = np.sinh(attenuation * distance)

    return np.array(
        [
            [hyperbolic_cosine, -hyperbolic_sine / (diffusion_length * attenuation)],
            [-diffusion_length * attenuation * hyperbolic_sine, hyperbolic_cosine],
        ]
    )


def exact_slab(*, z, reflectance, layered):
    layers = [(0.002, 40, 2000), (0.002, 20, 1000)] if layered else [(0.004, 20, 1000)]

    total = np.eye(2)

    for thickness, absorption, scattering in layers:
        total = (
            transfer(distance=thickness, absorption=absorption, scattering=scattering)
            @ total
        )

    beta = (1 - reflectance) / (2 * (1 + reflectance))

    incoming_drive = 2 / (1 + reflectance) * 0.001 / (0.002**2)

    combination = total[1, 1] - beta * total[0, 1]

    initial_fluence = (
        -combination
        * incoming_drive
        / (total[1, 0] - beta * total[0, 0] - beta * combination)
    )

    state = np.array([initial_fluence, incoming_drive - beta * initial_fluence])

    reference = np.empty_like(z)

    start = 0

    for thickness, absorption, scattering in layers:
        selected = (z >= start) & (z <= start + thickness)

        for index in np.flatnonzero(selected):
            propagated = (
                transfer(
                    distance=z[index] - start,
                    absorption=absorption,
                    scattering=scattering,
                )
                @ state
            )

            reference[index] = propagated[0]

        state = (
            transfer(distance=thickness, absorption=absorption, scattering=scattering)
            @ state
        )

        start += thickness

    return reference


@pytest.mark.parametrize("reflectance", [0, 0.5])
@pytest.mark.parametrize("layered", [False, True])
def test_slab_convergence_and_energy_balance(reflectance, layered):
    errors = []

    for count in (3, 5, 9):
        simulation, result = slab(count=count, reflectance=reflectance, layered=layered)

        z = simulation.mesh.points.to("m").magnitude[..., 2].ravel()

        reference = exact_slab(z=z, reflectance=reflectance, layered=layered)

        error = np.linalg.norm(
            result.fluence.values.magnitude.ravel() - reference
        ) / np.linalg.norm(reference)

        errors.append(error)

        assert abs(result.energy_balance["residual_power"].to("W").magnitude) < 1e-13

        assert result.energy_balance["absorbed_power"].magnitude > 0

        assert result.energy_balance["escaping_power"].magnitude > 0

        assert result.diagnostics["linear_relative_residual"] < 1e-10

        assert result.flux.values.shape[1] == 3

        assert result.flux.values.to("W/m**2").magnitude[:, 2].mean() > 0

    assert errors[1] < errors[0]

    assert errors[2] < errors[1] * 0.6

    assert errors[-1] < 0.003


def test_nonabsorbing_slab_has_exact_linear_fluence():
    simulation, result = slab(count=4, reflectance=0.3, absorption=0)

    z = simulation.mesh.points.to("m").magnitude[..., 2]

    beta = 0.7 / 2.6

    current = (2 / 1.3) * 250 / (2 + beta * 0.004 / (1 / 3000))

    exact = current / beta + current * (0.004 - z) / (1 / 3000)

    np.testing.assert_allclose(result.fluence.values.magnitude, exact, rtol=1e-12)

    assert result.energy_balance["absorbed_power"].magnitude == 0

    assert result.energy_balance["escaping_power"].to("mW").magnitude == pytest.approx(
        1
    )


@pytest.mark.parametrize("layered", [False, True])
def test_unstructured_slab_converges_to_transfer_solution(layered):
    simulation, _ = slab(count=3, reflectance=0.5, layered=layered)
    errors = []

    for resolution in (1.0, 0.5, 0.25):
        simulation.mesh = MomentRT.TetrahedralMesh.from_sample(
            sample=simulation.medium,
            resolution=resolution * ureg.mm,
        )
        result = simulation.run()
        depths = simulation.mesh.points.to("m").magnitude[:, 2]
        reference = exact_slab(z=depths, reflectance=0.5, layered=layered)
        errors.append(
            np.linalg.norm(result.fluence.values.magnitude - reference)
            / np.linalg.norm(reference)
        )
        assert abs(result.energy_balance["residual_power"].to("W").magnitude) < 1e-12

    assert errors[-1] < errors[1] < errors[0]
    assert errors[-1] < 0.01
