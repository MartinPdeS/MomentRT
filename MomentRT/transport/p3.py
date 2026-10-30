"""Steady-state full three-dimensional P3 transport using even-parity FEM."""

from __future__ import annotations
from collections.abc import Mapping
from math import factorial
import numpy as np
from scipy import sparse, special
from scipy.sparse import linalg as sparse_linalg
import MomentRT


class P3:
    """Full steady-state three-dimensional P3 transport on a finite sample.

    Parameters
    ----------
    boundaries : mapping of str to RobinBoundary
        Six matched-index faces with zero effective reflectance. Reflecting and
        index-mismatched boundaries are not implemented by this P3 solver.

    Notes
    -----
    Six even harmonics are solved with linear tetrahedral FEM; ten odd
    harmonics are reconstructed elementwise. Materials provide normalized
    scattering phase-function moments. Illumination is isotropic over
    the incoming hemisphere. P3 is an angular approximation and requires
    validation for the intended configuration.
    """

    model_notes: tuple[str, ...] = (
        "P3 is a steady-state full three-dimensional angular-moment approximation through degree 3, not a solution of the full radiative transfer equation.",
        "P3 assumes isotropic incoming diffuse radiance and supports only zero-reflectance Robin faces with equal refractive index throughout; reflecting faces and index jumps are unsupported.",
        "Detector angular acceptance uses a Lambertian approximation; validate the selected model against an independent reference for the intended use.",
    )

    _EVEN = [
        (degree, order) for degree in (0, 2) for order in range(-degree, degree + 1)
    ]
    _ODD = [
        (degree, order) for degree in (1, 3) for order in range(-degree, degree + 1)
    ]

    @staticmethod
    def _harmonics(
        *, directions: np.ndarray, modes: list[tuple[int, int]]
    ) -> np.ndarray:
        """Evaluate orthonormal real harmonics in the documented ordering."""
        z = np.clip(directions[:, 2], -1, 1)

        phi = np.arctan2(directions[:, 1], directions[:, 0])

        values = []

        for degree, order in modes:
            absolute = abs(order)

            scale = np.sqrt(
                (2 * degree + 1)
                / (4 * np.pi)
                * factorial(degree - absolute)
                / factorial(degree + absolute)
            )

            legendre = special.lpmv(absolute, degree, z)

            if order < 0:
                values.append(np.sqrt(2) * scale * legendre * np.sin(absolute * phi))
            elif order > 0:
                values.append(np.sqrt(2) * scale * legendre * np.cos(absolute * phi))
            else:
                values.append(scale * legendre)

        return np.column_stack(values)

    @staticmethod
    def _angular_quadrature(
        *, normal: np.ndarray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Integrate the sphere or the outward hemisphere about an axis-aligned normal."""
        nodes, weights = np.polynomial.legendre.leggauss(8)

        if normal is not None:
            nodes = (nodes + 1) / 2

            weights = weights / 2

        azimuth = 2 * np.pi * np.arange(16) / 16

        mu, angle = np.meshgrid(nodes, azimuth, indexing="ij")

        transverse = np.sqrt(1 - mu**2)

        local = np.column_stack(
            (
                transverse.ravel() * np.cos(angle.ravel()),
                transverse.ravel() * np.sin(angle.ravel()),
                mu.ravel(),
            )
        )

        if normal is not None:
            tangent = np.eye(3)[(int(np.argmax(np.abs(normal))) + 1) % 3]

            other = np.cross(normal, tangent)

            local = (
                local[:, 0, None] * tangent
                + local[:, 1, None] * other
                + local[:, 2, None] * normal
            )

        return local, np.repeat(weights, len(azimuth)) * (2 * np.pi / len(azimuth))

    @classmethod
    def _angular_operators(
        cls,
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        dict[MomentRT.boundary.Face, np.ndarray],
        dict[MomentRT.boundary.Face, np.ndarray],
    ]:
        """Build exact P3 streaming and matched-index half-range matrices."""
        directions, weights = cls._angular_quadrature()

        even = cls._harmonics(directions=directions, modes=cls._EVEN)

        odd = cls._harmonics(directions=directions, modes=cls._ODD)

        streaming = np.einsum("q,qo,qi,qe->ioe", weights, odd, directions, even)

        current = np.einsum("q,qo,qi->io", weights, odd, directions)

        boundary = {}

        illumination = {}

        for face in MomentRT.boundary.FACES:
            normal = np.zeros(3)

            normal["xyz".index(face[0])] = -1 if face.endswith("min") else 1

            outward, face_weights = cls._angular_quadrature(normal=normal)

            face_even = cls._harmonics(directions=outward, modes=cls._EVEN)

            weighted_mu = face_weights * (outward @ normal)

            boundary[face] = 2 * np.einsum(
                "q,qa,qb->ab", weighted_mu, face_even, face_even
            )

            illumination[face] = np.einsum("q,qa->a", weighted_mu, face_even) / np.pi

        return streaming, current, boundary, illumination

    def __init__(
        self,
        *,
        boundaries: Mapping[
            MomentRT.boundary.Face,
            MomentRT.boundary.RobinBoundary | MomentRT.boundary.ReflectingBoundary,
        ],
    ) -> None:
        """Configure matched-index P3 boundaries.

        Parameters
        ----------
        boundaries : mapping of str to RobinBoundary
            Six zero-reflectance Robin faces.
        """
        validated = MomentRT.transport.diffusion.Diffusion(boundaries=boundaries)

        unsupported = [
            face
            for face, boundary in validated.boundaries.items()
            if not isinstance(boundary, MomentRT.boundary.RobinBoundary)
            or boundary.effective_reflectance != 0
        ]

        if unsupported:
            raise ValueError(
                "P3 supports only zero-reflectance Robin boundaries on all faces. "
                f"Unsupported faces: {', '.join(unsupported)}."
            )

        self.boundaries = validated.boundaries

    def _prepare(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> _PreparedP3:
        validated = P3(boundaries=self.boundaries)

        return _PreparedP3(
            medium=medium, source=source, mesh=mesh, boundaries=validated.boundaries
        )

    def solve(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> MomentRT.results.SimulationResult:
        """Solve the P3 moment equations at the source wavelength.

        Parameters
        ----------
        medium : StructuredMedium
            Finite sample with explicitly configured material anisotropies.
        source : UniformSurfaceSource or SurfaceSource
            Incoming isotropic diffuse source on a sample face.
        mesh : CartesianMesh or TetrahedralMesh
            Nodes spanning the sample.

        Returns
        -------
        SimulationResult
            Fluence, element-centered flux, detector boundary fields, and six
            named nodal even-moment coefficients.
        """
        prepared = self._prepare(medium=medium, source=source, mesh=mesh)

        return prepared.solve(medium=medium, source=source, mesh=mesh)


class _PreparedP3(MomentRT.transport.diffusion._PreparedDiffusion):
    """Prepare physical geometry once for one solve or wavelength sweep."""

    def solve(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> MomentRT.results.SimulationResult:
        materials = MomentRT.transport.diffusion._finite_materials_at(
            medium=self.medium, wavelength=source.wavelength
        )

        names = ["background", *(shape.name for shape in self.medium.structures)]

        missing = [
            name
            for name, material in zip(names, materials)
            if material.anisotropy is None
            and isinstance(
                material.phase_function, MomentRT.phase_function.HenyeyGreenstein
            )
        ]

        if missing:
            raise ValueError(
                "P3 requires anisotropy on every material using Henyey–Greenstein; missing for: "
                + ", ".join(missing)
            )

        MomentRT.units.validate_units(
            value=source.optical_power,
            unit="W",
            name="optical_power",
            shape=(),
            positive=True,
        )

        absorption = np.array(
            [material.absorption.to("1/m").magnitude for material in materials]
        )

        scattering = np.array(
            [material.reduced_scattering.to("1/m").magnitude for material in materials]
        )

        factors = np.array([material._transport_factors() for material in materials])

        streaming, current, face_matrices, face_source = (
            MomentRT.transport.p3.P3._angular_operators()
        )

        geometry = self.geometry

        count = len(geometry.points)
        elements = geometry.tetrahedra
        element_count = len(elements)
        local_shape = (element_count, 4, 4)
        rows = np.broadcast_to(elements[:, :, None], local_shape).ravel()
        columns = np.broadcast_to(elements[:, None, :], local_shape).ravel()
        mass_template = (np.ones((4, 4)) + np.eye(4)) / 20

        blocks = [
            [sparse.csr_matrix((count, count)) for _ in MomentRT.transport.p3.P3._EVEN]
            for _ in MomentRT.transport.p3.P3._EVEN
        ]

        region_operators = {}

        for label in self.region_matrices:
            selection = self.labels == label

            odd_sigma = np.array(
                [
                    absorption[label] + scattering[label] * factors[label, degree]
                    for degree, _ in MomentRT.transport.p3.P3._ODD
                ]
            )

            diffusion = np.einsum(
                "ioe,o,jof->eijf", streaming, 1 / odd_sigma, streaming
            )

            even_sigma = np.array(
                [
                    absorption[label]
                    if degree == 0
                    else absorption[label] + scattering[label] * factors[label, 2]
                    for degree, _ in MomentRT.transport.p3.P3._EVEN
                ]
            )

            region_operators[label] = (selection, diffusion, even_sigma)

        for a in range(len(MomentRT.transport.p3.P3._EVEN)):
            for b in range(len(MomentRT.transport.p3.P3._EVEN)):
                local = np.zeros(local_shape)

                for label, (
                    selection,
                    diffusion,
                    even_sigma,
                ) in region_operators.items():
                    gradients = geometry.gradients[selection]
                    volumes = geometry.volumes[selection]

                    local[selection] = volumes[:, None, None] * np.einsum(
                        "tai,ij,tbj->tab",
                        gradients,
                        diffusion[a, :, :, b],
                        gradients,
                    )

                    if a == b:
                        local[selection] += (
                            even_sigma[a] * volumes[:, None, None] * mass_template
                        )

                blocks[a][b] = sparse.coo_matrix(
                    (local.ravel(), (rows, columns)), shape=(count, count)
                ).tocsr()

        for face, triangles in geometry.faces.items():
            areas = geometry.quadrature[face][2].sum(axis=1)

            face_mass = areas[:, None, None] * (np.ones((3, 3)) + np.eye(3)) / 12

            rows = np.broadcast_to(triangles[:, :, None], face_mass.shape).ravel()

            columns = np.broadcast_to(triangles[:, None, :], face_mass.shape).ravel()

            matrix = sparse.coo_matrix(
                (face_mass.ravel(), (rows, columns)), shape=(count, count)
            ).tocsr()

            for a in range(len(MomentRT.transport.p3.P3._EVEN)):
                for b in range(len(MomentRT.transport.p3.P3._EVEN)):
                    blocks[a][b] += face_matrices[face][a, b] * matrix

        operator = sparse.bmat(blocks, format="csr")

        power = source.optical_power.to("W").magnitude

        load = np.concatenate(
            [
                self.unit_load * power * coefficient
                for coefficient in face_source[self.source.face]
            ]
        )

        solution = sparse_linalg.spsolve(operator, load)

        if not np.all(np.isfinite(solution)):
            raise RuntimeError("The P3 solve returned nonfinite moments.")

        residual = np.linalg.norm(operator @ solution - load) / np.linalg.norm(load)

        if residual > 1e-8:
            raise RuntimeError(f"P3 linear residual is too large: {residual:g}.")

        even = solution.reshape(len(MomentRT.transport.p3.P3._EVEN), count)

        scalar = np.sqrt(4 * np.pi) * even[0]

        gradients = np.einsum(
            "tek,tki->tei",
            even[:, geometry.tetrahedra].transpose(1, 0, 2),
            geometry.gradients,
        )

        odd_degrees = [degree for degree, _ in MomentRT.transport.p3.P3._ODD]

        odd_sigma = (
            absorption[self.labels, None]
            + scattering[self.labels, None] * factors[self.labels][:, odd_degrees]
        )

        odd = -np.einsum("ioe,tei->to", streaming, gradients) / odd_sigma

        flux = np.einsum("io,to->ti", current, odd)

        absorbed = np.sum(
            absorption[self.labels]
            * geometry.volumes
            * scalar[geometry.tetrahedra].mean(axis=1)
        )

        unit = MomentRT.units.ureg.watt / MomentRT.units.ureg.meter**2

        boundary_fields: dict[str, MomentRT.fields.BoundaryField] = {}

        escaping = {}

        for face, triangles in geometry.faces.items():
            quadrature, barycentric, weights = geometry.quadrature[face]

            even_face = np.einsum("qa,eta->etq", barycentric, even[:, triangles])

            incoming = self.unit_irradiance * power if face == self.source.face else 0

            net = np.sqrt(4 * np.pi) * (
                np.einsum("a,atq->tq", face_matrices[face][0], even_face)
                - 2 * face_source[face][0] * incoming
            )

            outgoing = net + incoming

            escaping[face] = float(np.sum(outgoing * weights))

            boundary_fields[face] = MomentRT.fields.BoundaryField(
                face=face,
                positions=quadrature * MomentRT.units.ureg.meter,
                area_weights=weights * MomentRT.units.ureg.meter**2,
                escaping_irradiance=outgoing * unit,
                sample_origin=self.medium.origin,
                sample_size=self.medium.size,
                wavelength=source.wavelength,
            )

        result_mesh = self.mesh.copy()

        metadata = {
            "transport": "steady-state full P3 even-parity FEM",
            "wavelength": source.wavelength.copy(),
            "sample_origin": self.medium.origin.copy(),
            "sample_size": self.medium.size.copy(),
            "coefficient_sampling": "element centroid",
            "source_face": self.source.face,
        }

        moments = {
            f"l{degree}_m{order}": MomentRT.fields.ScalarField(
                mesh=result_mesh,
                values=even[index].reshape(self.mesh.shape) * unit,
                name=f"P3 harmonic l={degree}, m={order}",
                metadata=metadata,
            )
            for index, (degree, order) in enumerate(MomentRT.transport.p3.P3._EVEN)
        }

        reused = self.solve_count > 0

        self.solve_count += 1

        energy_balance, outgoing_power_by_face = MomentRT.results._energy_accounting(
            incoming_power=source.optical_power,
            absorbed_watts=absorbed,
            escaping_watts=escaping,
        )

        return MomentRT.results.SimulationResult(
            fluence=MomentRT.fields.ScalarField(
                mesh=result_mesh,
                values=scalar.reshape(self.mesh.shape) * unit,
                name="Fluence rate",
                metadata=metadata,
            ),
            flux=MomentRT.fields.VectorField(
                positions=geometry.centers * MomentRT.units.ureg.meter,
                values=flux * unit,
                name="P3 flux",
            ),
            moments=moments,
            energy_balance=energy_balance,
            outgoing_power_by_face=outgoing_power_by_face,
            diagnostics={
                "linear_relative_residual": residual,
                "source_quadrature_area": self.source_sampled_area
                * MomentRT.units.ureg.meter**2,
                "source_geometric_area": self.source_geometric_area
                * MomentRT.units.ureg.meter**2,
                "tetrahedron_count": len(geometry.tetrahedra),
                "minimum_fluence": scalar.min() * unit,
                "preparation_reused": reused,
            },
            boundary_fields=boundary_fields,
        )
