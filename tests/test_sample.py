import runpy

import numpy as np
import pytest

import MomentRT
from MomentRT import (
    Box,
    CartesianMesh,
    InfiniteDiffusion,
    Layer,
    Material,
    Sphere,
    StructuredMedium,
    ureg,
)


def material(*, absorption):
    return Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=absorption / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )


def sample_configuration():
    background = material(absorption=0.01)

    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.full(3, 4) * ureg.mm,
    )

    sample.add_background(material=background)

    layer = Layer(
        name="layer",
        depth=0 * ureg.mm,
        thickness=2 * ureg.mm,
        material=material(absorption=0.02),
    )

    sphere = Sphere(
        name="sphere",
        center=np.ones(3) * ureg.mm,
        radius=0.5 * ureg.mm,
        material=material(absorption=0.1),
    )

    sample.add_structures(layer, sphere)

    return sample


def test_material_precedence_and_background_replacement():
    sample = sample_configuration()

    points = np.array([[1, 1, 1], [3, 3, 1], [3, 3, 3], [2, 2, 2]]) * ureg.mm

    np.testing.assert_array_equal(sample.material_labels(points=points), [2, 1, 0, 1])

    replacement = material(absorption=0.04)

    assert sample.add_background(material=replacement) is None

    assert len(sample.structures) == 2

    assert sample.background is replacement

    np.testing.assert_array_equal(
        sample.material_labels(points=points.to("m")), [2, 1, 0, 1]
    )


def test_optical_fields_are_independent_snapshots():
    sample = sample_configuration()

    mesh = CartesianMesh(
        shape=(5, 5, 5),
        spacing=np.ones(3) * ureg.mm,
        origin=np.zeros(3) * ureg.mm,
    )

    fields = sample.optical_fields(mesh=mesh)

    assert fields["absorption"].values.to("1/mm").magnitude[1, 1, 1] == pytest.approx(
        0.1
    )

    assert fields["absorption"].values.to("1/mm").magnitude[3, 3, 3] == pytest.approx(
        0.01
    )

    sample.add_background(material=material(absorption=0.9))

    assert fields["absorption"].values.to("1/mm").magnitude[3, 3, 3] == pytest.approx(
        0.01
    )


def test_geometry_validation_and_atomic_addition():
    sample = sample_configuration()

    box = Box(
        name="box",
        origin=np.ones(3) * ureg.mm,
        size=np.ones(3) * ureg.mm,
        material=material(absorption=0.3),
    )

    outside = Sphere(
        name="outside",
        center=np.zeros(3) * ureg.mm,
        radius=1 * ureg.mm,
        material=material(absorption=0.3),
    )

    with pytest.raises(ValueError, match="outside"):
        sample.add_structures(box, outside)

    assert len(sample.structures) == 2

    sample.add_structures(box)

    np.testing.assert_array_equal(
        sample.material_labels(points=np.array([[1.5, 1.5, 1.5]]) * ureg.mm), [3]
    )

    with pytest.raises(ValueError, match="Duplicate"):
        sample.add_structures(box)

    with pytest.raises(TypeError):
        Sphere(
            name="bad",
            center=np.zeros(3),
            radius=1 * ureg.mm,
            material=sample.background,
        )

    with pytest.raises(ValueError, match="positive"):
        Layer(
            name="bad",
            depth=0 * ureg.mm,
            thickness=-1 * ureg.mm,
            material=sample.background,
        )

    with pytest.raises(ValueError, match="inside"):
        sample.material_labels(points=np.array([9, 9, 9]) * ureg.mm)


def test_unconfigured_sample_and_unsupported_solver():
    sample = StructuredMedium(origin=np.zeros(3) * ureg.mm, size=np.ones(3) * ureg.mm)

    with pytest.raises(ValueError, match="background"):
        sample.plot(show=False)

    solver = InfiniteDiffusion()

    with pytest.raises(TypeError, match="boundary-aware"):
        solver.solve(medium=sample, source=None, mesh=None)


def test_example_and_geometry_plot():
    pytest.importorskip("plotly")

    example = runpy.run_path("docs/examples/finite_samples/structured_sample.py")

    sample = example["sample"]

    figure = sample.plot(show=False, length_unit="um")

    assert [trace.name for trace in figure.data] == [
        "Background / sample extent",
        "Epidermis",
        "Papillary layer",
        "Reticular layer",
        "Absorbing inclusion",
        "Scattering inclusion",
        "Rectangular heterogeneity",
    ]

    assert max(figure.data[1].z) == pytest.approx(200)

    assert figure.layout.scene.xaxis.title.text == "x [um]"

    with pytest.raises(ValueError, match="opacity"):
        sample.plot(show=False, opacity=2)


def test_layer_stacking_mixed_units_and_atomic_failure():
    sample = StructuredMedium(
        origin=np.array([0, 0, -2]) * ureg.mm,
        size=np.array([4, 4, 4]) * ureg.mm,
    )

    background = material(absorption=0.01)

    sample.add_background(material=background)

    assert (
        sample.add_layer(name="first", thickness=500 * ureg.um, material=background)
        is None
    )

    sample.add_layer(name="second", thickness=1.5 * ureg.mm, material=background)

    assert sample.structures[0].depth.to("mm").magnitude == pytest.approx(-2)

    assert sample.structures[1].depth.to("mm").magnitude == pytest.approx(-1.5)

    explicit = Layer(
        name="explicit", depth=0 * ureg.mm, thickness=1 * ureg.mm, material=background
    )

    sample.add_structures(explicit)

    sample.add_layer(name="last", thickness=1 * ureg.mm, material=background)

    assert sample.structures[-1].depth.to("mm").magnitude == pytest.approx(1)

    with pytest.raises(ValueError, match="outside"):
        sample.add_layer(name="overflow", thickness=1 * ureg.mm, material=background)

    with pytest.raises(ValueError, match="Duplicate"):
        sample.add_layer(name="last", thickness=1 * ureg.mm, material=background)

    assert len(sample.structures) == 4


@pytest.mark.parametrize(
    "normal, axes", [("x", ("y", "z")), ("y", ("x", "z")), ("z", ("x", "y"))]
)
def test_sections_resolve_material_precedence(normal, axes):
    sample = sample_configuration()

    section = sample.section(normal=normal, position=1 * ureg.mm, resolution=(9, 7))

    assert section.axes == axes

    assert section.labels.shape == (7, 9)

    horizontal, vertical = np.meshgrid(
        section.horizontal.to("mm").magnitude, section.vertical.to("mm").magnitude
    )

    points = np.ones((7, 9, 3))

    points[..., "xyz".index(axes[0])] = horizontal

    points[..., "xyz".index(axes[1])] = vertical

    expected = sample.material_labels(points=points * ureg.mm)

    np.testing.assert_array_equal(section.labels, expected)

    assert 2 in section.labels

    np.testing.assert_allclose(
        section.fields["absorption"].to("1/mm").magnitude,
        np.array([0.01, 0.02, 0.1])[expected],
    )

    sample.background.absorption = 0.9 / ureg.mm

    np.testing.assert_allclose(
        section.fields["absorption"].to("1/mm").magnitude,
        np.array([0.01, 0.02, 0.1])[expected],
    )


def test_section_plots_and_validation():
    pytest.importorskip("plotly")

    sample = sample_configuration()

    section = sample.section(normal="y", position=1 * ureg.mm, resolution=(9, 7))

    figure = section.plot(show=False)

    np.testing.assert_array_equal(figure.data[0].z, section.labels)

    assert figure.layout.yaxis.range == (4, 0)

    property_figure = section.plot(
        quantity="absorption", value_unit="1/mm", scale="log", show=False
    )

    np.testing.assert_allclose(
        property_figure.data[0].z,
        np.log10(section.fields["absorption"].to("1/mm").magnitude),
    )

    assert property_figure.data[0].zsmooth is False

    with pytest.raises(TypeError):
        sample.section(normal="y", position=1)

    with pytest.raises(ValueError, match="normal"):
        sample.section(normal="bad", position=1 * ureg.mm)

    with pytest.raises(ValueError, match="inside"):
        sample.section(normal="y", position=-1 * ureg.mm)

    with pytest.raises(ValueError, match="resolution"):
        sample.section(normal="y", position=1 * ureg.mm, resolution=(1, 7))

    with pytest.raises(ValueError, match="linear"):
        section.plot(scale="log", show=False)

    with pytest.raises(ValueError, match="quantity"):
        section.plot(quantity="unknown", show=False)

    sample.background.absorption = 0 / ureg.mm

    zero_section = sample.section(normal="y", position=3 * ureg.mm)

    with pytest.raises(ValueError, match="positive"):
        zero_section.plot(quantity="absorption", scale="log", show=False)
