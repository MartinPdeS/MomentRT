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
    ureg,
)


def configuration():
    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0.02 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    medium = HomogeneousMedium(material=material)

    mesh = CartesianMesh(
        shape=(20, 20, 20),
        spacing=np.ones(3) * ureg.mm,
        origin=np.full(3, -9.5) * ureg.mm,
    )

    source = PointSource(
        position=np.zeros(3) * ureg.mm,
        wavelength=970 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    solver = InfiniteDiffusion()

    return Simulation(medium=medium, source=source, mesh=mesh, solver=solver)


def test_diffusion_energy_and_radial_symmetry():
    simulation = configuration()

    result = simulation.run()

    values = result.fluence.values.to("W/m**2").magnitude

    np.testing.assert_allclose(values, values[::-1, :, :])

    np.testing.assert_allclose(values, values.transpose(2, 1, 0))

    # Independent spherical energy balance: outward power plus absorbed
    # power inside radius r must equal source power.
    diffusion_length = 1 / (3 * 1020)

    attenuation = np.sqrt(20 / diffusion_length)

    radius = np.linalg.norm(simulation.mesh.points.to("m").magnitude, axis=-1)

    outward_power = (
        4 * np.pi * diffusion_length * radius * values * (1 + attenuation * radius)
    )

    absorbed_power = 0.001 * (
        1 - (1 + attenuation * radius) * np.exp(-attenuation * radius)
    )

    np.testing.assert_allclose(outward_power + absorbed_power, 0.001, rtol=1e-12)


def test_explicit_units_and_singularity():
    with pytest.raises(TypeError):
        Material(
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
            absorption=0.02,
            reduced_scattering=1 / ureg.mm,
            refractive_index=1.4,
        )

    simulation = configuration()

    simulation.source.position = simulation.mesh.points[0, 0, 0]

    with pytest.raises(ValueError, match="singular"):
        simulation.run()


@pytest.mark.parametrize("view", ["slices", "volume", "isosurfaces"])
def test_plotly_fields(view, tmp_path):

    simulation = configuration()

    field = simulation.run().fluence

    original = field.values.magnitude.copy()

    figure = field.plot(view=view, scale="log", show=False)

    assert (
        figure.data[0].type
        == {"slices": "surface", "volume": "volume", "isosurfaces": "isosurface"}[view]
    )

    assert figure.layout.scene.zaxis.autorange == "reversed"

    np.testing.assert_array_equal(field.values.magnitude, original)

    if view == "slices":
        plane = original[original.shape[0] // 2, :, :]

        np.testing.assert_allclose(figure.data[0].surfacecolor, np.log10(plane))

        np.testing.assert_allclose(figure.data[0].customdata, plane)
    else:
        np.testing.assert_allclose(figure.data[0].value, np.log10(original).ravel())

        density_slider = next(
            slider
            for slider in figure.layout.sliders
            if slider.currentvalue.prefix == "Iso-surface density: "
        )

        assert density_slider.active == {"volume": 4, "isosurfaces": 3}[view]

        density_step = next(step for step in density_slider.steps if step.label == "24")

        assert density_step.method == "restyle"

        assert density_step.args[0] == {"surface.count": 24}

        assert list(density_step.args[1]) == [0]

        if view == "volume":
            assert len(figure.layout.sliders) == 2

    figure.write_html(tmp_path / f"{view}.html")

    assert (tmp_path / f"{view}.html").stat().st_size > 1000


def test_plotly_system():

    simulation = configuration()

    figure = simulation.plot_system(show=False, length_unit="m")

    assert [trace.type for trace in figure.data] == ["mesh3d", "scatter3d", "scatter3d"]

    assert figure.layout.scene.xaxis.title.text == "x [m]"

    np.testing.assert_allclose(figure.data[0].x, [-0.0095] * 4 + [0.0095] * 4)
