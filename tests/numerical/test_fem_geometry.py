import numpy as np
import pytest
from MomentRT import CartesianMesh, UniformSurfaceSource, fem, ureg
from MomentRT.fem import (
    FiniteElementGeometry,
    boundary_triangles,
    element_geometry,
    tetrahedralize,
    triangle_areas,
)


def test_prepared_source_load_uses_face_quadrature(monkeypatch):
    """Uniform unit-square illumination gives the exact triangle nodal loads."""
    mesh = CartesianMesh(
        shape=(2, 2, 2),
        spacing=np.ones(3) * ureg.meter,
        origin=np.zeros(3) * ureg.meter,
    )
    geometry = FiniteElementGeometry(mesh=mesh)
    source = UniformSurfaceSource(
        face="z_min",
        wavelength=850 * ureg.nm,
        optical_power=2 * ureg.mW,
    )
    monkeypatch.setattr(
        fem,
        "triangle_quadrature",
        lambda **kwargs: pytest.fail("Quadrature was recomputed."),
    )

    irradiance, sampled_area = geometry.source_irradiance(source=source)
    load, load_area = geometry.source_load(source=source)

    np.testing.assert_allclose(irradiance, 0.002)
    np.testing.assert_allclose(sampled_area, 1)
    np.testing.assert_allclose(load_area, 1)
    np.testing.assert_allclose(load, np.full((2, 3), 0.002 / 6))


def test_region_matrices_integrate_constant_and_linear_fields():
    """Region mass and stiffness reproduce integrals on a unit cube."""
    mesh = CartesianMesh(
        shape=(2, 2, 2),
        spacing=np.ones(3) * ureg.meter,
        origin=np.zeros(3) * ureg.meter,
    )
    geometry = FiniteElementGeometry(mesh=mesh)
    matrices = geometry.region_matrices(labels=np.array([0, 0, 0, 1, 1, 1]))
    ones = np.ones(len(geometry.points))
    x = geometry.points[:, 0]

    for stiffness, mass in matrices.values():
        np.testing.assert_allclose(ones @ mass @ ones, 0.5)
        np.testing.assert_allclose(x @ stiffness @ x, 0.5)
        np.testing.assert_allclose(stiffness @ ones, 0, atol=1e-14)


def test_tetrahedral_partition_and_linear_gradients():
    mesh = CartesianMesh(
        shape=(3, 4, 5),
        spacing=np.array([1, 2, 3]) * ureg.mm,
        origin=np.array([7, -4, 2]) * ureg.mm,
    )

    points = mesh.points.to("m").magnitude.reshape(-1, 3)

    tetrahedra = tetrahedralize(mesh=mesh)

    volumes, gradients = element_geometry(points=points, tetrahedra=tetrahedra)

    assert len(tetrahedra) == 6 * 2 * 3 * 4

    np.testing.assert_allclose(volumes.sum(), 0.002 * 0.006 * 0.012)

    coefficients = np.array([2, -3, 4])

    linear_values = 10 + points @ coefficients

    recovered = np.einsum("ti,tij->tj", linear_values[tetrahedra], gradients)

    np.testing.assert_allclose(
        recovered, np.broadcast_to(coefficients, recovered.shape), atol=1e-11
    )

    faces = boundary_triangles(tetrahedra=tetrahedra, shape=mesh.shape)

    expected = {"x": 0.006 * 0.012, "y": 0.002 * 0.012, "z": 0.002 * 0.006}

    for face, triangles in faces.items():
        np.testing.assert_allclose(
            triangle_areas(points=points, triangles=triangles).sum(), expected[face[0]]
        )
