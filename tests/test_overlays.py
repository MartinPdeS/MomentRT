"""Spatial overlay coordinates, grouping, and shared source/detector hooks."""

import numpy as np
import pytest
import MomentRT
from MomentRT import (
    CartesianMesh,
    HomogeneousMedium,
    InfiniteDiffusion,
    Material,
    PointSource,
    ScalarField,
    Simulation,
    SurfaceDetector,
    SurfaceSource,
    UniformSurfaceSource,
    ureg,
)

pytest.importorskip("plotly")


@pytest.fixture
def field():
    mesh = CartesianMesh(
        shape=(3, 5, 3),
        spacing=np.array([1, 1, 2]) * ureg.mm,
        origin=np.array([-1, -2, 0]) * ureg.mm,
    )

    return ScalarField(
        mesh=mesh,
        values=np.ones(mesh.shape) * ureg.W / ureg.m**2,
        name="Fluence",
        metadata={
            "sample_origin": np.array([-1, -2, 0]) * ureg.mm,
            "sample_size": np.array([2, 4, 4]) * ureg.mm,
        },
    )


def test_named_overlays_convert_units_and_preserve_field_and_opacity_control(field):
    source = SurfaceSource(
        name="Illumination",
        face="z_min",
        position=np.array([0, 0, 0]) * ureg.mm,
        radius=0.4 * ureg.mm,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    detector = SurfaceDetector(
        name="Collection",
        face="z_max",
        position=np.array([0, 0, 4]) * ureg.mm,
        radius=0.5 * ureg.mm,
        collection_axis=[1, 0, 1],
        acceptance_angle=30 * ureg.degree,
    )

    original = field.values.copy()

    figure = field.plot(
        overlays=[source, detector],
        length_unit="um",
        show=False,
    )

    np.testing.assert_allclose(figure.data[1].z, 0)

    assert max(figure.data[1].x) == pytest.approx(400)

    np.testing.assert_allclose(figure.data[5].z, 4000)

    assert max(figure.data[5].x) == pytest.approx(500)

    np.testing.assert_allclose(
        np.array([figure.data[7].x[1], figure.data[7].y[1], figure.data[7].z[1]]),
        [1000 / np.sqrt(2), 0, 4000 + 1000 / np.sqrt(2)],
    )

    np.testing.assert_allclose(figure.data[3].z, [-800, 0])

    assert figure.data[2].type == "mesh3d"

    assert figure.data[6].type == "mesh3d"

    assert figure.data[4].type == "cone"

    assert figure.data[8].type == "cone"

    assert figure.data[1].name == "Illumination"

    assert figure.data[5].name == "Collection"

    assert figure.data[1].legendgroup == figure.data[2].legendgroup

    assert figure.data[5].legendgroup == figure.data[7].legendgroup

    assert figure.data[1].legendgroup != figure.data[5].legendgroup

    assert figure.layout.legend.orientation == "h"

    assert all(list(step.args[1]) == [0] for step in figure.layout.sliders[0].steps)

    np.testing.assert_array_equal(field.values.magnitude, original.magnitude)


@pytest.mark.parametrize(
    "face,axis,coordinate",
    [
        ("x_min", 0, -1),
        ("x_max", 0, 1),
        ("y_min", 1, -2),
        ("y_max", 1, 2),
        ("z_min", 2, 0),
        ("z_max", 2, 4),
    ],
)
def test_full_face_overlay_uses_physical_bounds(field, face, axis, coordinate):
    source = UniformSurfaceSource(
        face=face,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    figure = field.plot(
        view="isosurfaces",
        overlays=[source],
        show=False,
    )

    surface = figure.data[1]

    np.testing.assert_allclose(getattr(surface, "xyz"[axis]), coordinate)

    for other_axis, (lower, upper) in enumerate([(-1, 1), (-2, 2), (0, 4)]):
        if other_axis != axis:
            coordinates = np.asarray(getattr(surface, "xyz"[other_axis]))

            assert coordinates.min() == pytest.approx(lower)

            assert coordinates.max() == pytest.approx(upper)


def test_full_face_overlay_rejects_observation_window_without_sample_bounds(field):
    field.metadata.clear()

    source = UniformSurfaceSource(
        face="z_min",
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    with pytest.raises(ValueError, match="finite-sample bounds"):
        field.plot(overlays=[source], show=False)


def test_layer_overlay_does_not_mutate_sample_bounds(field):
    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0.02 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    layer = MomentRT.Layer(
        name="Dermis",
        depth=0.5 * ureg.mm,
        thickness=1 * ureg.mm,
        material=material,
    )

    original_origin = field.metadata["sample_origin"].copy()

    original_size = field.metadata["sample_size"].copy()

    field.plot(overlays=[layer], show=False)

    np.testing.assert_array_equal(
        field.metadata["sample_origin"].magnitude,
        original_origin.magnitude,
    )

    np.testing.assert_array_equal(
        field.metadata["sample_size"].magnitude,
        original_size.magnitude,
    )


def test_point_source_hook_is_shared_with_system_plot(field):
    source = PointSource(
        position=np.array([0.1, 0.2, 0.3]) * ureg.mm,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0.02 / ureg.mm,
        reduced_scattering=1 / ureg.mm,
        refractive_index=1.4,
    )

    simulation = Simulation(
        medium=HomogeneousMedium(material=material),
        source=source,
        mesh=field.mesh,
        solver=InfiniteDiffusion(),
    )

    standalone = field.plot(overlays=[source], length_unit="um", show=False)

    system = simulation.plot_system(length_unit="um", show=False)

    for name, expected in zip("xyz", [100, 200, 300]):
        assert getattr(standalone.data[-1], name)[0] == pytest.approx(expected)

        assert getattr(system.data[-1], name)[0] == pytest.approx(expected)


@pytest.mark.parametrize("overlays", ["source", [object()], 1])
def test_invalid_overlays_raise_actionable_error(field, overlays):
    with pytest.raises(TypeError, match="overlays|overlay"):
        field.plot(overlays=overlays, show=False)
