"""Prepared finite-element geometry and face-based calculations."""

from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from scipy import sparse
import MomentRT


class FiniteElementGeometry:
    """Prepare numerical geometry for a fixed finite tetrahedral mesh.

    Parameters
    ----------
    mesh : CartesianMesh or TetrahedralMesh
        Mesh providing nodes and tetrahedra. Cartesian cells are subdivided
        into six tetrahedra. Coordinates and quadrature are stored in SI units.

    Notes
    -----
    This numerical snapshot is prepared per solver call and must be rebuilt
    after changing the mesh. Element centers and gradients are not smoothed
    across material interfaces.
    """

    def __init__(
        self, *, mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh
    ) -> None:
        """Initialize finite-element geometry from a mesh.

        Parameters
        ----------
        mesh : CartesianMesh or TetrahedralMesh
            Nodes and tetrahedra used for finite-element discretization.
        """

        self.points = mesh.points.to("m").magnitude.reshape(-1, 3)

        self.tetrahedra = (
            mesh.tetrahedra
            if isinstance(mesh, MomentRT.mesh.TetrahedralMesh)
            else MomentRT.fem.tetrahedralize(mesh=mesh)
        )

        self.volumes, self.gradients = MomentRT.fem.element_geometry(
            points=self.points,
            tetrahedra=self.tetrahedra,
        )

        self.centers = self.points[self.tetrahedra].mean(axis=1)

        self.faces = (
            mesh.faces
            if isinstance(mesh, MomentRT.mesh.TetrahedralMesh)
            else MomentRT.fem.boundary_triangles(
                tetrahedra=self.tetrahedra,
                shape=mesh.shape,
            )
        )

        self.quadrature = {
            face: MomentRT.fem.triangle_quadrature(
                points=self.points,
                triangles=triangles,
            )
            for face, triangles in self.faces.items()
        }

    def source_irradiance(
        self, *, source: MomentRT.source.UniformSurfaceSource
    ) -> tuple[NDArray[np.float64], float]:
        """Evaluate incoming irradiance on the source's prepared face.

        Parameters
        ----------
        source : UniformSurfaceSource, SurfaceSource, or SurfaceSourceGroup
            Incoming diffuse illumination. Its optical power is already
            coupled into the sample.

        Returns
        -------
        tuple of ndarray and float
            Incoming irradiance in W/m², shape (F, 3), at the face quadrature
            points, and sampled source area in m². Group area is the sum of
            member areas, which can count overlapping footprints more than once.

        Raises
        ------
        SourceResolutionError
            If a source footprint contains no face quadrature points.
        """
        quadrature, _, weights = self.quadrature[source.face]

        return MomentRT.fem.source_irradiance(
            quadrature=quadrature,
            weights=weights,
            source=source,
        )

    def source_load(
        self, *, source: MomentRT.source.UniformSurfaceSource
    ) -> tuple[NDArray[np.float64], float]:
        """Integrate incoming irradiance against the prepared face basis.

        Parameters
        ----------
        source : UniformSurfaceSource, SurfaceSource, or SurfaceSourceGroup
            Incoming diffuse illumination already coupled into the sample.

        Returns
        -------
        tuple of ndarray and float
            Triangle-local nodal loads in watts, shape (F, 3), and sampled
            source area in m². The Robin source factor is not applied.

        Notes
        -----
        The precomputed three-point rule normalizes sampled incoming power;
        footprint-area agreement alone does not establish mesh convergence.
        """
        _, barycentric, weights = self.quadrature[source.face]

        irradiance, sampled_area = self.source_irradiance(source=source)

        local_load = np.einsum("tq,qa->ta", weights * irradiance, barycentric)

        return local_load, sampled_area

    def region_matrices(
        self, *, labels: NDArray[np.int64]
    ) -> dict[int, tuple[sparse.csr_matrix, sparse.csr_matrix]]:
        """Assemble unit diffusion and absorption matrices for each occupied region.

        Parameters
        ----------
        labels : ndarray of int, shape (T,)
            Resolved material label at each tetrahedral centroid.

        Returns
        -------
        dict of int to tuple of scipy.sparse.csr_matrix
            Unit-coefficient stiffness and mass matrices for each occupied region.
            Multiply by diffusion length (m) and absorption (1/m), respectively.
        """

        stiffness = self.volumes[:, None, None] * np.einsum(
            "tik,tjk->tij", self.gradients, self.gradients
        )

        mass = self.volumes[:, None, None] * (np.ones((4, 4)) + np.eye(4)) / 20

        shape = (len(self.points), len(self.points))

        matrices = {}

        for label in np.unique(labels):
            selected = labels == label

            tetrahedra = self.tetrahedra[selected]

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
