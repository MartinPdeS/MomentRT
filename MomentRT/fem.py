"""Linear tetrahedral FEM kernels on a conforming Cartesian subdivision."""

from itertools import permutations

import numpy as np


def tetrahedralize(*, mesh):
    """Split each voxel into six tetrahedra sharing its lower/upper diagonal."""
    shape = np.asarray(mesh.shape)

    cells = np.stack(
        np.meshgrid(*(np.arange(n - 1) for n in shape), indexing="ij"), axis=-1
    ).reshape(-1, 3)

    tetrahedra = []

    for order in permutations(range(3)):
        offsets = np.zeros((4, 3), dtype=int)

        for step, axis in enumerate(order, start=1):
            offsets[step] = offsets[step - 1]

            offsets[step, axis] += 1

        indices = cells[:, None, :] + offsets

        tetrahedra.append(np.ravel_multi_index(indices.transpose(2, 0, 1), shape))

    return np.concatenate(tetrahedra)


def element_geometry(*, points, tetrahedra):
    vertices = points[tetrahedra]

    edges = vertices[:, 1:] - vertices[:, :1]

    volumes = np.abs(np.linalg.det(edges)) / 6

    gradients = np.empty((len(tetrahedra), 4, 3))

    gradients[:, 1:] = np.linalg.inv(edges).transpose(0, 2, 1)

    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)

    return volumes, gradients


def boundary_triangles(*, tetrahedra, shape):
    candidates = np.concatenate(
        [tetrahedra[:, face] for face in ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])]
    )

    unique, counts = np.unique(np.sort(candidates, axis=1), axis=0, return_counts=True)

    triangles = unique[counts == 1]

    indices = np.stack(np.unravel_index(triangles, shape), axis=-1)

    faces = {}

    for axis, name in enumerate("xyz"):
        faces[f"{name}_min"] = triangles[np.all(indices[..., axis] == 0, axis=1)]

        faces[f"{name}_max"] = triangles[
            np.all(indices[..., axis] == shape[axis] - 1, axis=1)
        ]

    return faces


def triangle_areas(*, points, triangles):
    vertices = points[triangles]

    return (
        np.linalg.norm(
            np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]),
            axis=1,
        )
        / 2
    )


def triangle_quadrature(*, points, triangles):
    """Return positions, barycentric basis values, and area weights in SI."""
    barycentric = np.array(
        [[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]]
    )

    quadrature = np.einsum("qa,tad->tqd", barycentric, points[triangles])

    areas = triangle_areas(points=points, triangles=triangles)

    weights = np.broadcast_to(areas[:, None] / 3, quadrature.shape[:2])

    return quadrature, barycentric, weights


def source_irradiance(*, quadrature, weights, source):
    """Evaluate the incoming diffuse profile with its sampled-area normalization."""
    from .source import SurfaceSource

    selected = np.ones(quadrature.shape[:2], dtype=bool)

    if isinstance(source, SurfaceSource):
        axis = "xyz".index(source.face[0])

        in_plane = [i for i in range(3) if i != axis]

        displacement = (
            quadrature[..., in_plane] - source.position.to("m").magnitude[in_plane]
        )

        selected = (
            np.linalg.norm(displacement, axis=-1) <= source.radius.to("m").magnitude
        )

    sampled_area = np.sum(weights * selected)

    if sampled_area == 0:
        raise ValueError("The mesh does not resolve the source footprint; refine it.")

    irradiance = (source.optical_power.to("W").magnitude / sampled_area) * selected

    return irradiance, sampled_area


def source_load(*, points, triangles, source):
    """Three-point triangle quadrature for the disk footprint and nodal load."""
    quadrature, barycentric, weights = triangle_quadrature(
        points=points, triangles=triangles
    )

    irradiance, sampled_area = source_irradiance(
        quadrature=quadrature, weights=weights, source=source
    )

    local_load = np.einsum("tq,qa->ta", weights * irradiance, barycentric)

    return local_load, sampled_area
