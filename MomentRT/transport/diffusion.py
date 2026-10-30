"""Analytical steady-state diffusion Green function in an infinite medium."""

from __future__ import annotations
import copy
from collections.abc import Mapping
import numpy as np
from pint import Quantity
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg
import MomentRT


class InfiniteDiffusion:
    """Evaluate Phi = P exp(-sqrt(mu_a/D) r) / (4 pi D r).

    No boundary, Fresnel interface, or finite-source correction is applied.
    The observation window does not truncate the physical medium.

    Notes
    -----
    The Green function is Phi(r) = P exp(-sqrt(mu_a / D) r) / (4 pi D r).
    No physical outer boundary, Fresnel interface, or finite-source correction
    is applied. The observation window does not truncate the medium.
    """

    model_notes: tuple[str, ...] = (
        "InfiniteDiffusion is an analytical steady-state P1 Green function for an embedded isotropic point source in a homogeneous, unbounded medium.",
        "The point-source singularity is excluded; boundaries, finite source footprints, and refractive interfaces are not modeled.",
        "P1 diffusion is least reliable close to the source and where scattering does not dominate absorption.",
    )

    def solve(
        self,
        *,
        medium: MomentRT.medium.HomogeneousMedium,
        source: MomentRT.source.PointSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> MomentRT.fields.ScalarField:
        """Solve the transport equation for the supplied physical configuration.

        Parameters
        ----------
        medium : HomogeneousMedium
            Infinite homogeneous material resolved at the source wavelength.
        source : PointSource
            Embedded isotropic point source; no grid node may coincide with it.
        mesh : CartesianMesh
            Cartesian nodes used for sampling or finite-element discretization.

        Returns
        -------
        ScalarField
            Analytical infinite-medium fluence rate in watts per square meter.

        Raises
        ------
        TypeError
            If the medium or source is incompatible with the infinite model.
        ValueError
            If a grid node coincides with the singular point source.
        """

        if not isinstance(medium, MomentRT.medium.HomogeneousMedium):
            raise TypeError(
                "InfiniteDiffusion supports only HomogeneousMedium; finite structured samples require a boundary-aware solver."
            )

        if not isinstance(source, MomentRT.source.PointSource):
            raise TypeError(
                "InfiniteDiffusion requires an embedded PointSource, not surface illumination."
            )

        if not isinstance(mesh, MomentRT.mesh.CartesianMesh):
            raise TypeError("InfiniteDiffusion requires a Cartesian observation mesh.")

        medium = medium.at(wavelength=source.wavelength)

        absorption = medium.material.absorption.to("1/m").magnitude

        scattering = medium.material.reduced_scattering.to("1/m").magnitude

        diffusion_length = 1 / (3 * (absorption + scattering))

        displacement = mesh.points.to("m").magnitude - source.position.to("m").magnitude

        distance = np.linalg.norm(displacement, axis=-1)

        if np.any(distance == 0):
            raise ValueError(
                "The point source is singular; choose a grid that excludes its position."
            )

        values = (
            source.optical_power.to("W").magnitude
            * np.exp(-np.sqrt(absorption / diffusion_length) * distance)
            / (4 * np.pi * diffusion_length * distance)
        )

        return MomentRT.fields.ScalarField(
            mesh=mesh,
            values=values * MomentRT.units.ureg.watt / MomentRT.units.ureg.meter**2,
            name="Fluence rate",
            metadata={
                "transport": "infinite-medium steady-state P1 diffusion",
                "wavelength": source.wavelength,
                "source_position": source.position,
            },
        )


def _finite_materials_at(
    *, medium: MomentRT.medium.StructuredMedium, wavelength: Quantity
) -> list[MomentRT.material.Material]:
    """Resolve finite material coefficients and the supported interface convention.

    Parameters
    ----------
    medium : StructuredMedium
        Finite sample with a background and no internal refractive-index jumps.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.

    Returns
    -------
    list of Material
        Resolved background followed by structure materials; internal index
        jumps are rejected.
    """
    if medium.background is None:
        raise ValueError("Configure the sample background before solving.")

    materials = [
        medium.background.at(wavelength=wavelength),
        *(shape.material.at(wavelength=wavelength) for shape in medium.structures),
    ]

    if any(
        material.refractive_index != materials[0].refractive_index
        for material in materials
    ):
        names = ["background", *(shape.name for shape in medium.structures)]

        indices = ", ".join(
            f"{name!r}={material.refractive_index:g}"
            for name, material in zip(names, materials)
        )

        raise ValueError(
            "Internal refractive-index jumps are not supported by this "
            f"continuous-fluence FEM solver at {wavelength:~P}: {indices}."
        )

    return materials


class Diffusion:
    """Steady-state P1 diffusion with linear tetrahedral finite elements.

    Geometry preparation is scoped to a solve or wavelength sweep. No persistent
    cache survives changes to a sample, source, mesh, or boundary configuration.

    Parameters
    ----------
    boundaries : mapping of str to RobinBoundary or ReflectingBoundary
        Conditions for exactly all six sample faces. Reflecting faces cannot
        carry illumination.

    Notes
    -----
    The steady-state P1 equation is -div(D grad(Phi)) + mu_a Phi = 0,
    with J = -D grad(Phi). Linear tetrahedra use centroid-sampled material
    coefficients. Fluence is continuous across interfaces; internal
    refractive-index jumps are unsupported. Footprint and unaligned interface
    errors require mesh refinement. Geometry preparation is scoped to a
    solve or sweep, with no persistent cache.
    """

    model_notes: tuple[str, ...] = (
        "Diffusion is a steady-state P1 approximation on a finite sample, using linear tetrahedra and centroid-sampled material coefficients.",
        "Fluence is continuous across internal interfaces; refractive-index jumps are unsupported, while absorption and reduced scattering may vary.",
        "P1 diffusion is least reliable close to illumination and in weakly scattering regions; mesh refinement does not remove this model error.",
        "Robin effective_reflectance is supplied explicitly and is not inferred from refractive index.",
    )

    def __init__(
        self,
        *,
        boundaries: Mapping[
            MomentRT.boundary.Face,
            MomentRT.boundary.RobinBoundary | MomentRT.boundary.ReflectingBoundary,
        ],
    ) -> None:
        """Initialize the Diffusion configuration.

        Parameters
        ----------
        boundaries : mapping of str to RobinBoundary or ReflectingBoundary
            Conditions for exactly all six sample faces. Reflecting faces cannot
            carry illumination.
        """

        if set(boundaries) != set(MomentRT.boundary.FACES):
            raise ValueError(
                f"Configure exactly these boundary faces: {MomentRT.boundary.FACES}."
            )

        if any(
            not isinstance(
                condition,
                (MomentRT.boundary.RobinBoundary, MomentRT.boundary.ReflectingBoundary),
            )
            for condition in boundaries.values()
        ):
            raise TypeError("Each face requires RobinBoundary or ReflectingBoundary.")

        self.boundaries = dict(boundaries)

    def _prepare(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> _PreparedDiffusion:
        """Prepare frozen finite-element geometry and source sampling.

        Parameters
        ----------
        medium : StructuredMedium
            Finite sample with a background and no internal refractive-index jumps.
        source : UniformSurfaceSource or SurfaceSource
            Incoming diffuse illumination on a Robin face; optical power is
            already coupled into the sample.
        mesh : CartesianMesh or TetrahedralMesh
            Nodal coordinates for sampling or finite-element discretization.

        Returns
        -------
        _PreparedDiffusion
            Frozen geometry and source footprint prepared for this solve or sweep.
        """

        return _PreparedDiffusion(
            medium=medium,
            source=source,
            mesh=mesh,
            boundaries=self.boundaries,
        )

    def solve(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> MomentRT.results.SimulationResult:
        """Solve the transport equation for the supplied physical configuration.

        Parameters
        ----------
        medium : StructuredMedium
            Finite sample with a background and no internal refractive-index jumps.
        source : UniformSurfaceSource or SurfaceSource
            Incoming diffuse illumination on a Robin face; optical power is
            already coupled into the sample.
        mesh : CartesianMesh or TetrahedralMesh
            Nodal coordinates for sampling or finite-element discretization.

        Returns
        -------
        SimulationResult
            Finite-element fluence, element flux, boundary fields, and power balance.

        Raises
        ------
        TypeError
            If the medium or source is incompatible with the finite model.
        ValueError
            If the mesh does not span the sample, interfaces have refractive-index
            jumps, or a source does not fit on an illuminated Robin face.
        SourceResolutionError
            If no boundary quadrature point resolves the illumination footprint.
        """

        problem = self._prepare(medium=medium, source=source, mesh=mesh)

        return problem.solve(medium=medium, source=source, mesh=mesh)


class _PreparedDiffusion:
    """One frozen instrument geometry with wavelength-dependent bulk operators."""

    def __init__(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
        boundaries: Mapping[
            MomentRT.boundary.Face,
            MomentRT.boundary.RobinBoundary | MomentRT.boundary.ReflectingBoundary,
        ],
    ) -> None:
        """Initialize the _PreparedDiffusion configuration.

        Parameters
        ----------
        medium : StructuredMedium
            Finite sample with a background and no internal refractive-index jumps.
        source : UniformSurfaceSource or SurfaceSource
            Incoming diffuse illumination on a Robin face; optical power is
            already coupled into the sample.
        mesh : CartesianMesh or TetrahedralMesh
            Nodal coordinates for sampling or finite-element discretization.
        boundaries : mapping of str to RobinBoundary or ReflectingBoundary
            Conditions for exactly all six sample faces. Reflecting faces cannot
            carry illumination.
        """

        if not isinstance(medium, MomentRT.medium.StructuredMedium):
            raise TypeError("Diffusion requires a finite StructuredMedium.")

        if not isinstance(source, MomentRT.source.UniformSurfaceSource):
            raise TypeError("Diffusion requires SurfaceSource or UniformSurfaceSource.")

        if medium.background is None:
            raise ValueError("Configure the sample background before solving.")

        if isinstance(mesh, MomentRT.mesh.TetrahedralMesh):
            mesh._validate_sample(sample=medium)

        MomentRT.units.validate_units(
            value=source.optical_power,
            unit="W",
            name="optical_power",
            shape=(),
            positive=True,
        )

        self.medium = copy.deepcopy(medium)

        self.source = copy.deepcopy(source)

        self.boundaries = copy.deepcopy(boundaries)

        self.mesh = mesh.copy()

        lower = self.medium.origin.to("m").magnitude

        upper = (self.medium.origin + self.medium.size).to("m").magnitude

        mesh_lower, mesh_upper = self.mesh.bounds.to("m").magnitude

        if not np.allclose(
            mesh_lower, lower, rtol=1e-12, atol=1e-15
        ) or not np.allclose(mesh_upper, upper, rtol=1e-12, atol=1e-15):
            raise ValueError(
                "The mesh must span exactly the finite sample; use a mesh from_sample()."
            )

        condition = self.boundaries[self.source.face]

        if not isinstance(condition, MomentRT.boundary.RobinBoundary):
            raise ValueError(
                f"The illuminated face {self.source.face!r} requires a RobinBoundary."
            )

        axis = "xyz".index(self.source.face[0])

        in_plane = [index for index in range(3) if index != axis]

        if isinstance(self.source, MomentRT.source.SurfaceSourceGroup):
            for source in self.source.sources:
                position = source.position.to("m").magnitude

                face_coordinate = (
                    lower[axis] if source.face.endswith("min") else upper[axis]
                )

                if not np.isclose(
                    position[axis], face_coordinate, rtol=1e-12, atol=1e-15
                ):
                    raise ValueError("Source position must lie on its configured face.")

                MomentRT.boundary._validate_face_disk(
                    sample=self.medium,
                    face=source.face,
                    position=source.position,
                    radius=source.radius,
                )

            self.source_geometric_area = sum(
                np.pi * source.radius.to("m").magnitude ** 2
                for source in self.source.sources
            )
        elif isinstance(self.source, MomentRT.source.SurfaceSource):
            position = self.source.position.to("m").magnitude

            face_coordinate = (
                lower[axis] if self.source.face.endswith("min") else upper[axis]
            )

            if not np.isclose(position[axis], face_coordinate, rtol=1e-12, atol=1e-15):
                raise ValueError("Source position must lie on its configured face.")

            radius = self.source.radius.to("m").magnitude

            if np.any(position[in_plane] - radius < lower[in_plane]) or np.any(
                position[in_plane] + radius > upper[in_plane]
            ):
                raise ValueError(
                    MomentRT.boundary._disk_bounds_message(
                        face=self.source.face,
                        position=self.source.position,
                        radius=self.source.radius,
                        lower=self.medium.origin,
                        upper=self.medium.origin + self.medium.size,
                        label="illumination disk",
                    )
                )

            self.source_geometric_area = np.pi * radius**2
        else:
            self.source_geometric_area = float(np.prod((upper - lower)[in_plane]))

        self.geometry = MomentRT.fem._MeshGeometry(mesh=self.mesh)

        self.labels = self.medium.material_labels(
            points=self.geometry.centers * MomentRT.units.ureg.meter
        )

        self.region_matrices = MomentRT.fem._region_matrices(
            geometry=self.geometry, labels=self.labels
        )

        node_count = len(self.geometry.points)

        rows = []

        columns = []

        data = []

        self.nodal_leakage = {}

        for face, triangles in self.geometry.faces.items():
            _, _, weights = self.geometry.quadrature[face]

            areas = weights.sum(axis=1)

            beta = self.boundaries[face].leakage_coefficient

            boundary_mass = (
                beta * areas[:, None, None] * (np.ones((3, 3)) + np.eye(3)) / 12
            )

            rows.append(
                np.broadcast_to(triangles[:, :, None], boundary_mass.shape).ravel()
            )

            columns.append(
                np.broadcast_to(triangles[:, None, :], boundary_mass.shape).ravel()
            )

            data.append(boundary_mass.ravel())

            leakage = np.zeros(node_count)

            np.add.at(leakage, triangles.ravel(), np.repeat(beta * areas / 3, 3))

            self.nodal_leakage[face] = leakage

        self.boundary_matrix = sparse.coo_matrix(
            (np.concatenate(data), (np.concatenate(rows), np.concatenate(columns))),
            shape=(node_count, node_count),
        ).tocsr()

        quadrature, barycentric, weights = self.geometry.quadrature[self.source.face]

        incoming, self.source_sampled_area = MomentRT.fem.source_irradiance(
            quadrature=quadrature,
            weights=weights,
            source=self.source,
        )

        self.unit_irradiance = incoming / self.source.optical_power.to("W").magnitude

        local_load = np.einsum("tq,qa->ta", weights * self.unit_irradiance, barycentric)

        self.unit_load = np.zeros(node_count)

        np.add.at(
            self.unit_load,
            self.geometry.faces[self.source.face].ravel(),
            local_load.ravel() * condition.source_factor,
        )

        self.solve_count = 0

    def solve(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> MomentRT.results.SimulationResult:
        """Solve the transport equation for the supplied physical configuration.

        Parameters
        ----------
        medium : StructuredMedium
            Finite sample with a background and no internal refractive-index jumps.
        source : UniformSurfaceSource or SurfaceSource
            Incoming diffuse illumination on a Robin face; optical power is
            already coupled into the sample.
        mesh : CartesianMesh or TetrahedralMesh
            Nodal coordinates for sampling or finite-element discretization.

        Returns
        -------
        SimulationResult
            Finite-element fluence, element flux, boundary fields, and power balance.
        """

        MomentRT.units.validate_units(
            value=source.optical_power,
            unit="W",
            name="optical_power",
            shape=(),
            positive=True,
        )

        materials = _finite_materials_at(
            medium=self.medium, wavelength=source.wavelength
        )

        absorption_by_region = np.array(
            [material.absorption.to("1/m").magnitude for material in materials]
        )

        scattering_by_region = np.array(
            [material.reduced_scattering.to("1/m").magnitude for material in materials]
        )

        diffusion_by_region = 1 / (3 * (absorption_by_region + scattering_by_region))

        matrix = self.boundary_matrix.copy()

        for label, (stiffness, mass) in self.region_matrices.items():
            matrix += (
                diffusion_by_region[label] * stiffness
                + absorption_by_region[label] * mass
            )

        power = source.optical_power.to("W").magnitude

        load = self.unit_load * power

        values = sparse_linalg.spsolve(matrix, load)

        if not np.all(np.isfinite(values)):
            raise RuntimeError("The diffusion solve returned nonfinite fluence.")

        relative_residual = np.linalg.norm(matrix @ values - load) / np.linalg.norm(
            load
        )

        if relative_residual > 1e-8:
            raise RuntimeError(
                f"Diffusion linear residual is too large: {relative_residual:g}."
            )

        geometry = self.geometry

        absorption = absorption_by_region[self.labels]

        coefficient = diffusion_by_region[self.labels]

        absorbed = np.sum(
            absorption * geometry.volumes * values[geometry.tetrahedra].mean(axis=1)
        )

        escaping = {
            face: float(weights @ values)
            for face, weights in self.nodal_leakage.items()
        }

        escaping[self.source.face] -= (
            2 * self.boundaries[self.source.face].leakage_coefficient * power
        )

        boundary_fields: dict[str, MomentRT.fields.BoundaryField] = {}

        for face, triangles in geometry.faces.items():
            quadrature, barycentric, weights = geometry.quadrature[face]

            face_fluence = np.einsum("qa,ta->tq", barycentric, values[triangles])

            incoming = self.unit_irradiance * power if face == self.source.face else 0

            beta = self.boundaries[face].leakage_coefficient

            boundary_fields[face] = MomentRT.fields.BoundaryField(
                face=face,
                positions=quadrature * MomentRT.units.ureg.meter,
                area_weights=weights * MomentRT.units.ureg.meter**2,
                escaping_irradiance=beta
                * (face_fluence - 2 * incoming)
                * MomentRT.units.ureg.watt
                / MomentRT.units.ureg.meter**2,
                sample_origin=self.medium.origin,
                sample_size=self.medium.size,
                wavelength=source.wavelength,
            )

        flux_values = -coefficient[:, None] * np.einsum(
            "ti,tij->tj", values[geometry.tetrahedra], geometry.gradients
        )

        fluence = MomentRT.fields.ScalarField(
            mesh=self.mesh.copy(),
            values=values.reshape(self.mesh.shape)
            * MomentRT.units.ureg.watt
            / MomentRT.units.ureg.meter**2,
            name="Fluence rate",
            metadata={
                "transport": "steady-state tetrahedral P1 FEM diffusion",
                "wavelength": source.wavelength.copy(),
                "coefficient_sampling": "element centroid",
                "source_face": self.source.face,
                "sample_origin": self.medium.origin.copy(),
                "sample_size": self.medium.size.copy(),
            },
        )

        flux = MomentRT.fields.VectorField(
            positions=geometry.centers * MomentRT.units.ureg.meter,
            values=flux_values
            * MomentRT.units.ureg.watt
            / MomentRT.units.ureg.meter**2,
            name="Diffusive flux",
        )

        reused = self.solve_count > 0

        self.solve_count += 1

        energy_balance, outgoing_power_by_face = MomentRT.results._energy_accounting(
            incoming_power=source.optical_power,
            absorbed_watts=absorbed,
            escaping_watts=escaping,
        )

        return MomentRT.results.SimulationResult(
            fluence=fluence,
            flux=flux,
            energy_balance=energy_balance,
            outgoing_power_by_face=outgoing_power_by_face,
            diagnostics={
                "linear_relative_residual": relative_residual,
                "source_quadrature_area": self.source_sampled_area
                * MomentRT.units.ureg.meter**2,
                "source_geometric_area": self.source_geometric_area
                * MomentRT.units.ureg.meter**2,
                "tetrahedron_count": len(geometry.tetrahedra),
                "minimum_fluence": values.min()
                * MomentRT.units.ureg.watt
                / MomentRT.units.ureg.meter**2,
                "preparation_reused": reused,
            },
            boundary_fields=boundary_fields,
        )
