import ast
import builtins
from pathlib import Path

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
        Material(absorption=0.02, reduced_scattering=1 / ureg.mm, refractive_index=1.4)

    simulation = configuration()

    simulation.source.position = simulation.mesh.points[0, 0, 0]

    with pytest.raises(ValueError, match="singular"):
        simulation.run()


def test_vtk_order_and_export(tmp_path):
    pv = pytest.importorskip("pyvista")

    simulation = configuration()

    field = simulation.run().fluence

    grid = field.to_pyvista()

    np.testing.assert_allclose(
        grid.points, field.mesh.points.to("mm").magnitude.reshape(-1, 3, order="F")
    )

    target = tmp_path / "fluence.vti"

    field.export(path=target)

    restored = pv.read(target)

    np.testing.assert_allclose(
        restored[field.name], field.values.magnitude.ravel(order="F")
    )

    assert restored.field_data["coordinate_units"][0] == "mm"


@pytest.mark.parametrize("view", ["slices", "volume", "isosurfaces"])
def test_plotly_fields(view, tmp_path):
    pytest.importorskip("plotly")

    simulation = configuration()

    field = simulation.run().fluence

    original = field.values.magnitude.copy()

    figure = field.plot(view=view, scale="log", show=False)

    assert (
        figure.data[0].type
        == {"slices": "surface", "volume": "volume", "isosurfaces": "isosurface"}[view]
    )

    np.testing.assert_array_equal(field.values.magnitude, original)

    if view == "slices":
        plane = original[original.shape[0] // 2, :, :]

        np.testing.assert_allclose(figure.data[0].surfacecolor, np.log10(plane))

        np.testing.assert_allclose(figure.data[0].customdata, plane)
    else:
        np.testing.assert_allclose(figure.data[0].value, np.log10(original).ravel())

    figure.write_html(tmp_path / f"{view}.html")

    assert (tmp_path / f"{view}.html").stat().st_size > 1000


def test_plotly_system():
    pytest.importorskip("plotly")

    simulation = configuration()

    figure = simulation.plot_system(show=False, length_unit="m")

    assert [trace.type for trace in figure.data] == ["mesh3d", "scatter3d", "scatter3d"]

    assert figure.layout.scene.xaxis.title.text == "x [m]"

    np.testing.assert_allclose(figure.data[0].x, [-0.0095] * 4 + [0.0095] * 4)


def test_optional_visualization(monkeypatch):
    original = builtins.__import__

    def without_plotly(name, *args, **kwargs):
        if name.startswith("plotly"):
            raise ImportError("blocked for test")

        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_plotly)

    simulation = configuration()

    assert simulation.run().fluence.values.shape == (20, 20, 20)

    with pytest.raises(ImportError, match="viz"):
        simulation.plot_system(show=False)


def test_no_python_annotations():
    for path in Path(MomentRT.__file__).parent.glob("*.py"):
        tree = ast.parse(path.read_text())

        for node in ast.walk(tree):
            assert not isinstance(node, ast.AnnAssign)

            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                assert node.returns is None

            if isinstance(node, ast.arg):
                assert node.annotation is None
