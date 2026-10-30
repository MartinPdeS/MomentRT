"""Linear tetrahedral FEM kernels for structured and unstructured meshes."""

from __future__ import annotations
from itertools import permutations
from typing import cast
import numpy as np
from numpy.typing import NDArray
from scipy import sparse
import MomentRT


def tetrahedralize(*, mesh: MomentRT.mesh.CartesianMesh) -> NDArray[np.int64]:
    """Split each voxel into six tetrahedra sharing its lower/upper diagonal.

    Parameters
    ----------
    mesh : CartesianMesh
        Cartesian nodes used for sampling or finite-element discretization.

    Returns
    -------
    ndarray of int, shape (T, 4)
        Tetrahedral connectivity into mesh.points reshaped to (N, 3),
        with T = 6 * (nx - 1) * (ny - 1) * (nz - 1).
    """
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


def element_geometry(
    *, points: NDArray[np.float64], tetrahedra: NDArray[np.int64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Calculate tetrahedral volumes and linear shape-function gradients.

    Parameters
    ----------
    points : ndarray of float, shape (N, 3)
        Flattened nodal coordinates in meters; this numerical kernel
        accepts SI magnitudes rather than quantities.
    tetrahedra : ndarray of int, shape (T, 4)
        Tetrahedral node indices into the flattened Cartesian point array.

    Returns
    -------
    tuple of ndarray
        Element volumes in cubic meters, shape (T,), and linear basis
        gradients in inverse meters, shape (T, 4, 3).
    """

    vertices = points[tetrahedra]

    edges = vertices[:, 1:] - vertices[:, :1]

    volumes = np.abs(np.linalg.det(edges)) / 6

    gradients = np.empty((len(tetrahedra), 4, 3))

    gradients[:, 1:] = np.linalg.inv(edges).transpose(0, 2, 1)

    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)

    return volumes, gradients


def boundary_triangles(
    *, tetrahedra: NDArray[np.int64], shape: tuple[int, int, int]
) -> dict[MomentRT.boundary.Face, NDArray[np.int64]]:
    """Extract exterior tetrahedron faces and group them by sample boundary.

    Parameters
    ----------
    tetrahedra : ndarray of int, shape (T, 4)
        Tetrahedral node indices into the flattened Cartesian point array.
    shape : tuple of int
        Node counts (nx, ny, nz), each at least two.

    Returns
    -------
    dict of str to ndarray
        Exterior triangle connectivity, shape (F, 3) per sample face.
    """

    candidates = np.concatenate(
        [tetrahedra[:, face] for face in ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])]
    )

    unique, counts = np.unique(np.sort(candidates, axis=1), axis=0, return_counts=True)

    triangles = unique[counts == 1]

    indices = np.stack(np.unravel_index(triangles, shape), axis=-1)

    faces: dict[MomentRT.boundary.Face, NDArray[np.int64]] = {}

    for axis, name in enumerate("xyz"):
        faces[cast(MomentRT.boundary.Face, f"{name}_min")] = triangles[
            np.all(indices[..., axis] == 0, axis=1)
        ]

        faces[cast(MomentRT.boundary.Face, f"{name}_max")] = triangles[
            np.all(indices[..., axis] == shape[axis] - 1, axis=1)
        ]

    return faces


def triangle_areas(
    *, points: NDArray[np.float64], triangles: NDArray[np.int64]
) -> NDArray[np.float64]:
    """Calculate boundary triangle areas from SI nodal coordinates.

    Parameters
    ----------
    points : ndarray of float, shape (N, 3)
        Flattened nodal coordinates in meters; this numerical kernel
        accepts SI magnitudes rather than quantities.
    triangles : ndarray of int, shape (F, 3)
        Boundary triangle node indices.

    Returns
    -------
    ndarray of float, shape (F,)
        Triangle areas in square meters.
    """

    vertices = points[triangles]

    return (
        np.linalg.norm(
            np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]),
            axis=1,
        )
        / 2
    )


def triangle_quadrature(
    *, points: NDArray[np.float64], triangles: NDArray[np.int64]
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Return positions, barycentric basis values, and area weights in SI.

    Parameters
    ----------
    points : ndarray of float, shape (N, 3)
        Flattened nodal coordinates in meters; this numerical kernel
        accepts SI magnitudes rather than quantities.
    triangles : ndarray of int, shape (F, 3)
        Boundary triangle node indices.

    Returns
    -------
    tuple of ndarray
        Quadrature positions (F, 3, 3) in meters, barycentric basis values
        (3, 3), and area weights (F, 3) in square meters. The three-point
        rule integrates quadratic polynomials on each triangle exactly.
    """
    barycentric = np.array(
        [[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]]
    )

    quadrature = np.einsum("qa,tad->tqd", barycentric, points[triangles])

    areas = triangle_areas(points=points, triangles=triangles)

    weights = np.broadcast_to(areas[:, None] / 3, quadrature.shape[:2])

    return quadrature, barycentric, weights


def source_irradiance(
    *,
    quadrature: NDArray[np.float64],
    weights: NDArray[np.float64],
    source: MomentRT.source.UniformSurfaceSource,
) -> tuple[NDArray[np.float64], float]:
    """Evaluate the incoming diffuse profile with its sampled-area normalization.

    Parameters
    ----------
    quadrature : ndarray of float, shape (F, 3, 3)
        Boundary quadrature positions in meters.
    weights : ndarray of float, shape (F, 3)
        Boundary quadrature area weights in square meters.
    source : UniformSurfaceSource or SurfaceSource
        Incoming diffuse illumination on a Robin face; optical power is
        already coupled into the sample.

    Returns
    -------
    tuple of ndarray and float
        Incoming irradiance (F, 3) in watts per square meter and sampled
        footprint area in square meters. Its weighted integral equals the
        configured incoming power.
    """

    selected = np.ones(quadrature.shape[:2], dtype=bool)

    if isinstance(source, MomentRT.source.SurfaceSourceGroup):
        irradiance = np.zeros(quadrature.shape[:2])
        sampled_area = 0.0

        for member in source.sources:
            member_irradiance, member_area = source_irradiance(
                quadrature=quadrature,
                weights=weights,
                source=member,
            )

            irradiance += member_irradiance

            sampled_area += member_area

        return irradiance, float(sampled_area)

    if isinstance(source, MomentRT.source.SurfaceSource):
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
        raise MomentRT.source.SourceResolutionError(
            "The mesh does not resolve the source footprint; refine it."
        )

    irradiance = (source.optical_power.to("W").magnitude / sampled_area) * selected

    return irradiance, sampled_area


def source_load(
    *,
    points: NDArray[np.float64],
    triangles: NDArray[np.int64],
    source: MomentRT.source.UniformSurfaceSource,
) -> tuple[NDArray[np.float64], float]:
    """Three-point triangle quadrature for the disk footprint and nodal load.

    Parameters
    ----------
    points : ndarray of float, shape (N, 3)
        Flattened nodal coordinates in meters; this numerical kernel
        accepts SI magnitudes rather than quantities.
    triangles : ndarray of int, shape (F, 3)
        Boundary triangle node indices.
    source : UniformSurfaceSource or SurfaceSource
        Incoming diffuse illumination on a Robin face; optical power is
        already coupled into the sample.

    Returns
    -------
    tuple of ndarray and float
        Triangle-local nodal loads (F, 3) in watts and sampled source area
        in square meters; the Robin source factor has not yet been applied.
    """
    quadrature, barycentric, weights = triangle_quadrature(
        points=points, triangles=triangles
    )

    irradiance, sampled_area = source_irradiance(
        quadrature=quadrature, weights=weights, source=source
    )

    local_load = np.einsum("tq,qa->ta", weights * irradiance, barycentric)

    return local_load, sampled_area


class _MeshGeometry:
    """Numerical geometry prepared once for a fixed finite tetrahedral mesh."""

    def __init__(
        self, *, mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh
    ) -> None:
        """Initialize the _MeshGeometry configuration.

        Parameters
        ----------
        mesh : CartesianMesh or TetrahedralMesh
            Nodes and tetrahedra used for finite-element discretization.
        """

        self.points = mesh.points.to("m").magnitude.reshape(-1, 3)

        self.tetrahedra = (
            mesh.tetrahedra
            if isinstance(mesh, MomentRT.mesh.TetrahedralMesh)
            else tetrahedralize(mesh=mesh)
        )

        self.volumes, self.gradients = element_geometry(
            points=self.points,
            tetrahedra=self.tetrahedra,
        )

        self.centers = self.points[self.tetrahedra].mean(axis=1)

        self.faces = (
            mesh.faces
            if isinstance(mesh, MomentRT.mesh.TetrahedralMesh)
            else boundary_triangles(tetrahedra=self.tetrahedra, shape=mesh.shape)
        )

        self.quadrature = {
            face: triangle_quadrature(points=self.points, triangles=triangles)
            for face, triangles in self.faces.items()
        }


def _region_matrices(
    *, geometry: _MeshGeometry, labels: NDArray[np.int64]
) -> dict[int, tuple[sparse.csr_matrix, sparse.csr_matrix]]:
    """Assemble unit diffusion and absorption matrices for each occupied region.

    Parameters
    ----------
    geometry : _MeshGeometry
        Prepared nodal coordinates, elements, and boundary quadrature in SI.
    labels : ndarray of int, shape (T,)
        Resolved material label at each tetrahedral centroid.

    Returns
    -------
    dict of int to tuple of scipy.sparse.csr_matrix
        Unit-coefficient stiffness and mass matrices for each occupied region.
        Multiply by diffusion length (m) and absorption (1/m), respectively.
    """

    stiffness = geometry.volumes[:, None, None] * np.einsum(
        "tik,tjk->tij", geometry.gradients, geometry.gradients
    )

    mass = geometry.volumes[:, None, None] * (np.ones((4, 4)) + np.eye(4)) / 20

    shape = (len(geometry.points), len(geometry.points))

    matrices = {}

    for label in np.unique(labels):
        selected = labels == label

        tetrahedra = geometry.tetrahedra[selected]

        local_shape = (len(tetrahedra), 4, 4)

        rows = np.broadcast_to(tetrahedra[:, :, None], local_shape).ravel()

        columns = np.broadcast_to(tetrahedra[:, None, :], local_shape).ravel()

        matrices[int(label)] = (
            sparse.coo_matrix(
                (stiffness[selected].ravel(), (rows, columns)), shape=shape
            ).tocsr(),
            sparse.coo_matrix(
                (mass[selected].ravel(), (rows, columns)), shape=shape
            ).tocsr(),
        )

    return matrices
