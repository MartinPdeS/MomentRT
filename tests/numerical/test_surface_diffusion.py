import runpy
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
    SurfaceSource,
    ureg,
)

FACES = ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")


def configuration(*, face="z_min"):
    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0.02 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm, size=np.full(3, 4) * ureg.mm
    )

    sample.add_background(material=material)

    mesh = CartesianMesh.from_sample(sample=sample, shape=(7, 7, 7))

    position = np.full(3, 2.0)

    position["xyz".index(face[0])] = 0 if face.endswith("min") else 4

    source = SurfaceSource(
        face=face,
        position=position * ureg.mm,
        radius=0.8 * ureg.mm,
        wavelength=970 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    boundary = RobinBoundary(effective_reflectance=0.4)

    solver = Diffusion(boundaries={name: boundary for name in FACES})

    return sample, mesh, source, solver


@pytest.mark.parametrize("face", FACES)
def test_disk_power_on_every_face(face):
    sample, mesh, source, solver = configuration(face=face)

    simulation = Simulation(medium=sample, mesh=mesh, source=source, solver=solver)

    result = simulation.run()

    assert abs(result.energy_balance["residual_power"].to("W").magnitude) < 1e-13

    assert result.diagnostics["source_quadrature_area"].to(
        "mm**2"
    ).magnitude == pytest.approx(np.pi * 0.8**2, rel=0.25)

    assert result.fluence.values.min().magnitude > 0

    assert sum(result.outgoing_power_by_face.values()).to(
        "W"
    ).magnitude == pytest.approx(
        result.energy_balance["escaping_power"].to("W").magnitude
    )


def test_configuration_validation():
    sample, mesh, source, solver = configuration()

    with pytest.raises(ValueError, match="boundary faces"):
        Diffusion(boundaries={"z_min": RobinBoundary(effective_reflectance=0)})

    with pytest.raises(ValueError, match="plain number"):
        RobinBoundary(effective_reflectance=1)

    with pytest.raises(TypeError):
        SurfaceSource(
            face="z_min",
            position=np.zeros(3) * ureg.mm,
            radius=1,
            wavelength=970 * ureg.nm,
            optical_power=1 * ureg.mW,
        )

    source.position = np.ones(3) * ureg.mm

    with pytest.raises(ValueError, match="position"):
        solver.solve(medium=sample, mesh=mesh, source=source)

    source.position = np.zeros(3) * ureg.mm

    with pytest.raises(ValueError, match="disk"):
        solver.solve(medium=sample, mesh=mesh, source=source)

    sample, mesh, source, solver = configuration()

    source.radius = 1e-9 * ureg.mm

    with pytest.raises(ValueError, match="resolve"):
        solver.solve(medium=sample, mesh=mesh, source=source)

    sample, mesh, source, solver = configuration()

    mesh.origin = np.ones(3) * ureg.mm

    with pytest.raises(ValueError, match="mesh"):
        solver.solve(medium=sample, mesh=mesh, source=source)

    sample, mesh, source, solver = configuration()

    solver.boundaries["z_min"] = ReflectingBoundary()

    with pytest.raises(ValueError, match="illuminated"):
        solver.solve(medium=sample, mesh=mesh, source=source)

    sample, mesh, source, solver = configuration()

    other = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0.02 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.5,
    )

    sample.add_layer(name="index jump", material=other, thickness=1 * ureg.mm)

    with pytest.raises(ValueError, match="refractive-index"):
        solver.solve(medium=sample, mesh=mesh, source=source)


def test_surface_example_and_flux_figures():
    pytest.importorskip("plotly")

    example = runpy.run_path("docs/examples/getting_started/surface_diffusion.py")

    simulation = example["simulation"]

    figure = simulation.plot_system(show=False)

    assert "Illumination spot" in [trace.name for trace in figure.data]

    # Use a smaller grid for this integration check; the analytical tests cover convergence.
    simulation.mesh = CartesianMesh.from_sample(
        sample=simulation.medium, shape=(9, 9, 9)
    )

    result = simulation.run()

    fluence = result.fluence.plot(show=False, scale="log")

    flux = result.flux.plot(show=False, max_vectors=40)

    assert fluence.data[0].type == "volume"

    assert flux.data[0].type == "cone"

    assert len(flux.data[0].x) == 40


def test_source_power_scales_fields_and_all_power_outputs():
    sample, mesh, source, solver = configuration()

    first = solver.solve(medium=sample, mesh=mesh, source=source)

    source.optical_power = 2 * ureg.mW

    second = solver.solve(medium=sample, mesh=mesh, source=source)

    np.testing.assert_allclose(
        second.fluence.values.magnitude, 2 * first.fluence.values.magnitude, rtol=1e-12
    )

    np.testing.assert_allclose(
        second.flux.values.magnitude, 2 * first.flux.values.magnitude, atol=1e-10
    )

    assert second.energy_balance["absorbed_power"].to("W").magnitude == pytest.approx(
        2 * first.energy_balance["absorbed_power"].to("W").magnitude
    )
