import numpy as np
from MomentRT import CartesianMesh, ureg
from MomentRT.fem import (
    boundary_triangles,
    element_geometry,
    tetrahedralize,
    triangle_areas,
)


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
