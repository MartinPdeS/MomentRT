"""Analytical steady-state diffusion Green function in an infinite medium."""

import numpy as np

from .fields import ScalarField
from .medium import HomogeneousMedium
from .source import PointSource
from .units import ureg


class InfiniteDiffusion:
    """Evaluate Phi = P exp(-sqrt(mu_a/D) r) / (4 pi D r).

    No boundary, Fresnel interface, or finite-source correction is applied.
    The observation window does not truncate the physical medium.
    """

    def solve(self, *, medium, source, mesh):
        if not isinstance(medium, HomogeneousMedium):
            raise TypeError(
                "InfiniteDiffusion supports only HomogeneousMedium; finite structured samples require a boundary-aware solver."
            )

        if not isinstance(source, PointSource):
            raise TypeError(
                "InfiniteDiffusion requires an embedded PointSource, not surface illumination."
            )

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

        return ScalarField(
            mesh=mesh,
            values=values * ureg.watt / ureg.meter**2,
            name="Fluence rate",
            metadata={
                "transport": "infinite-medium steady-state P1 diffusion",
                "wavelength": source.wavelength,
                "source_position": source.position,
            },
        )


class Diffusion:
    """Steady-state P1 diffusion with linear tetrahedral finite elements.

    All six boundary faces must be configured explicitly. Material properties
    are constant per element, sampled at its centroid. Interfaces not aligned
    to element faces are approximated and require mesh refinement. Continuous
    fluence assumes equal refractive indices across all internal regions;
    internal index jumps are rejected rather than silently ignored.
    """

    def __init__(self, *, boundaries):
        from .boundary import FACES, ReflectingBoundary, RobinBoundary

        if set(boundaries) != set(FACES):
            raise ValueError(f"Configure exactly these boundary faces: {FACES}.")

        if any(
            not isinstance(condition, (RobinBoundary, ReflectingBoundary))
            for condition in boundaries.values()
        ):
            raise TypeError("Each face requires RobinBoundary or ReflectingBoundary.")

        self.boundaries = dict(boundaries)

    def solve(self, *, medium, source, mesh):
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import spsolve

        from .boundary import RobinBoundary
        from .fem import (
            boundary_triangles,
            element_geometry,
            source_irradiance,
            source_load,
            tetrahedralize,
            triangle_areas,
            triangle_quadrature,
        )
        from .fields import BoundaryField, VectorField
        from .medium import StructuredMedium
        from .results import SimulationResult
        from .source import SurfaceSource, UniformSurfaceSource

        if not isinstance(medium, StructuredMedium):
            raise TypeError("Diffusion requires a finite StructuredMedium.")

        if not isinstance(source, UniformSurfaceSource):
            raise TypeError("Diffusion requires SurfaceSource or UniformSurfaceSource.")

        medium = medium.at(wavelength=source.wavelength)

        if medium.background is None:
            raise ValueError("Configure the sample background before solving.")

        materials = [
            medium.background,
            *(shape.material for shape in medium.structures),
        ]

        if any(
            material.refractive_index != medium.background.refractive_index
            for material in materials
        ):
            raise ValueError(
                "Internal refractive-index jumps are not supported by this continuous-fluence FEM solver."
            )

        lower = medium.origin.to("m").magnitude

        upper = (medium.origin + medium.size).to("m").magnitude

        mesh_lower = mesh.origin.to("m").magnitude

        mesh_upper = (
            mesh_lower + (np.asarray(mesh.shape) - 1) * mesh.spacing.to("m").magnitude
        )

        if not np.allclose(
            mesh_lower, lower, rtol=1e-12, atol=1e-15
        ) or not np.allclose(mesh_upper, upper, rtol=1e-12, atol=1e-15):
            raise ValueError(
                "The mesh must span exactly the finite sample; use CartesianMesh.from_sample()."
            )

        condition = self.boundaries[source.face]

        if not isinstance(condition, RobinBoundary):
            raise ValueError("The illuminated face requires a RobinBoundary.")

        if isinstance(source, SurfaceSource):
            axis = "xyz".index(source.face[0])

            position = source.position.to("m").magnitude

            face_coordinate = (
                lower[axis] if source.face.endswith("min") else upper[axis]
            )

            if not np.isclose(position[axis], face_coordinate, rtol=1e-12, atol=1e-15):
                raise ValueError("Source position must lie on its configured face.")

            in_plane = [i for i in range(3) if i != axis]

            radius = source.radius.to("m").magnitude

            if np.any(position[in_plane] - radius < lower[in_plane]) or np.any(
                position[in_plane] + radius > upper[in_plane]
            ):
                raise ValueError("The entire illumination disk must fit on its face.")

        points = mesh.points.to("m").magnitude.reshape(-1, 3)

        tetrahedra = tetrahedralize(mesh=mesh)

        volumes, gradients = element_geometry(points=points, tetrahedra=tetrahedra)

        centers = points[tetrahedra].mean(axis=1)

        labels = medium.material_labels(points=centers * ureg.meter)

        absorption = np.array(
            [material.absorption.to("1/m").magnitude for material in materials]
        )[labels]

        scattering = np.array(
            [material.reduced_scattering.to("1/m").magnitude for material in materials]
        )[labels]

        coefficient = 1 / (3 * (absorption + scattering))

        stiffness = (
            volumes[:, None, None]
            * coefficient[:, None, None]
            * np.einsum("tik,tjk->tij", gradients, gradients)
        )

        mass = (
            volumes[:, None, None]
            * absorption[:, None, None]
            * (np.ones((4, 4)) + np.eye(4))
            / 20
        )

        rows = np.broadcast_to(tetrahedra[:, :, None], stiffness.shape).ravel()

        columns = np.broadcast_to(tetrahedra[:, None, :], stiffness.shape).ravel()

        row_blocks = [rows]

        column_blocks = [columns]

        data_blocks = [(stiffness + mass).ravel()]

        faces = boundary_triangles(tetrahedra=tetrahedra, shape=mesh.shape)

        nodal_leakage = {}

        for face, triangles in faces.items():
            areas = triangle_areas(points=points, triangles=triangles)

            beta = self.boundaries[face].leakage_coefficient

            boundary_mass = (
                beta * areas[:, None, None] * (np.ones((3, 3)) + np.eye(3)) / 12
            )

            row_blocks.append(
                np.broadcast_to(triangles[:, :, None], boundary_mass.shape).ravel()
            )

            column_blocks.append(
                np.broadcast_to(triangles[:, None, :], boundary_mass.shape).ravel()
            )

            data_blocks.append(boundary_mass.ravel())

            weights = np.zeros(len(points))

            np.add.at(weights, triangles.ravel(), np.repeat(beta * areas / 3, 3))

            nodal_leakage[face] = weights

        matrix = coo_matrix(
            (
                np.concatenate(data_blocks),
                (np.concatenate(row_blocks), np.concatenate(column_blocks)),
            ),
            shape=(len(points), len(points)),
        ).tocsr()

        local_load, sampled_area = source_load(
            points=points, triangles=faces[source.face], source=source
        )

        load = np.zeros(len(points))

        np.add.at(
            load,
            faces[source.face].ravel(),
            local_load.ravel() * condition.source_factor,
        )

        values = spsolve(matrix, load)

        if not np.all(np.isfinite(values)):
            raise RuntimeError("The diffusion solve returned nonfinite fluence.")

        relative_residual = np.linalg.norm(matrix @ values - load) / np.linalg.norm(
            load
        )

        if relative_residual > 1e-8:
            raise RuntimeError(
                f"Diffusion linear residual is too large: {relative_residual:g}."
            )

        absorbed = np.sum(absorption * volumes * values[tetrahedra].mean(axis=1))

        power = source.optical_power.to("W").magnitude

        escaping = {
            face: float(weights @ values) for face, weights in nodal_leakage.items()
        }

        # Remove the incident-field contribution to the gross Robin leakage.
        # This distinguishes outgoing transmitted partial current from net current.
        escaping[source.face] -= 2 * condition.leakage_coefficient * power

        boundary_fields = {}

        for face, triangles in faces.items():
            quadrature, barycentric, weights = triangle_quadrature(
                points=points, triangles=triangles
            )

            face_fluence = np.einsum("qa,ta->tq", barycentric, values[triangles])

            incoming = np.zeros_like(face_fluence)

            if face == source.face:
                incoming, _ = source_irradiance(
                    quadrature=quadrature, weights=weights, source=source
                )

            beta = self.boundaries[face].leakage_coefficient

            boundary_fields[face] = BoundaryField(
                face=face,
                positions=quadrature * ureg.meter,
                area_weights=weights * ureg.meter**2,
                escaping_irradiance=beta
                * (face_fluence - 2 * incoming)
                * ureg.watt
                / ureg.meter**2,
                sample_origin=medium.origin,
                sample_size=medium.size,
                wavelength=source.wavelength,
            )

        flux_values = -coefficient[:, None] * np.einsum(
            "ti,tij->tj", values[tetrahedra], gradients
        )

        fluence = ScalarField(
            mesh=mesh,
            values=values.reshape(mesh.shape) * ureg.watt / ureg.meter**2,
            name="Fluence rate",
            metadata={
                "transport": "steady-state tetrahedral P1 FEM diffusion",
                "wavelength": source.wavelength,
                "coefficient_sampling": "element centroid",
                "source_face": source.face,
            },
        )

        flux = VectorField(
            positions=centers * ureg.meter,
            values=flux_values * ureg.watt / ureg.meter**2,
            name="Diffusive flux",
        )

        return SimulationResult(
            fluence=fluence,
            flux=flux,
            energy_balance={
                "incoming_power": source.optical_power,
                "absorbed_power": absorbed * ureg.watt,
                "escaping_power": sum(escaping.values()) * ureg.watt,
                "residual_power": (power - absorbed - sum(escaping.values()))
                * ureg.watt,
            },
            outgoing_power_by_face={
                face: value * ureg.watt for face, value in escaping.items()
            },
            diagnostics={
                "linear_relative_residual": relative_residual,
                "source_quadrature_area": sampled_area * ureg.meter**2,
                "tetrahedron_count": len(tetrahedra),
                "minimum_fluence": values.min() * ureg.watt / ureg.meter**2,
            },
            boundary_fields=boundary_fields,
        )
