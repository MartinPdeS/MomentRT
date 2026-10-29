"""Steady-state SP3 transport with matched-index Marshak boundaries."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg

import MomentRT


class SP3:
    """Solve two coupled simplified P3 moments using tetrahedral elements.

    Parameters
    ----------
    boundaries : mapping of str to RobinBoundary or ReflectingBoundary
        Conditions on all six faces. Robin boundaries must have zero effective
        reflectance, corresponding to a matched-index vacuum boundary.

    Notes
    -----
    Each material supplies its phase function. Henyey–Greenstein materials
    require anisotropy; isotropic, Rayleigh, and tabulated functions provide
    their own first three normalized Legendre moments.
    Illumination is isotropic over the incoming hemisphere. Fluence is the
    zeroth moment, not the first composite unknown. The implementation assumes
    continuous composite moments across internal interfaces and rejects index
    jumps. SP3 is an approximation, and agreement with SP3 slab references does
    not establish accuracy against radiative transfer for a particular sensor.
    """

    def __init__(
        self,
        *,
        boundaries: Mapping[
            MomentRT.boundary.Face,
            MomentRT.boundary.RobinBoundary | MomentRT.boundary.ReflectingBoundary,
        ],
    ) -> None:
        """Configure the six sample boundary conditions.

        Parameters
        ----------
        boundaries : mapping of str to RobinBoundary or ReflectingBoundary
            Six face conditions; Robin reflectance must be zero.
        """
        validated = MomentRT.diffusion.Diffusion(boundaries=boundaries)

        if any(
            isinstance(boundary, MomentRT.boundary.RobinBoundary)
            and boundary.effective_reflectance != 0
            for boundary in validated.boundaries.values()
        ):
            raise ValueError("SP3 supports only zero-reflectance Robin boundaries.")

        self.boundaries = validated.boundaries

    def _prepare(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh,
    ) -> _PreparedSP3:
        validated = SP3(boundaries=self.boundaries)

        return _PreparedSP3(
            medium=medium,
            source=source,
            mesh=mesh,
            boundaries=validated.boundaries,
        )

    def solve(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh,
    ) -> MomentRT.results.SimulationResult:
        """Solve the finite sample for one source wavelength.

        Parameters
        ----------
        medium : StructuredMedium
            Finite sample with a background and equal internal indices.
        source : UniformSurfaceSource or SurfaceSource
            Incoming diffuse power, uniformly distributed over the inward
            hemisphere at each point of the illuminated footprint.
        mesh : CartesianMesh
            Nodes spanning the finite sample.

        Returns
        -------
        SimulationResult
            Zeroth-moment fluence, first-moment flux, boundary fields, energy
            balance, and composite fields in ``moments``.
        """
        prepared = self._prepare(medium=medium, source=source, mesh=mesh)

        return prepared.solve(medium=medium, source=source, mesh=mesh)


class _PreparedSP3(MomentRT.diffusion._PreparedDiffusion):
    """Reuse the finite geometry and source quadrature during a sweep."""

    def __init__(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh,
        boundaries: Mapping[
            MomentRT.boundary.Face,
            MomentRT.boundary.RobinBoundary | MomentRT.boundary.ReflectingBoundary,
        ],
    ) -> None:
        super().__init__(medium=medium, source=source, mesh=mesh, boundaries=boundaries)

    def solve(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh,
    ) -> MomentRT.results.SimulationResult:
        MomentRT.units.validate_units(
            value=source.optical_power,
            unit="W",
            name="optical_power",
            shape=(),
            positive=True,
        )

        materials = MomentRT.diffusion._finite_materials_at(
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
                "SP3 requires anisotropy on every material using Henyey–Greenstein; missing for: "
                + ", ".join(missing)
            )

        phase_moments = np.array([material._phase_moments() for material in materials])

        anisotropy = phase_moments[:, 1]

        factors = np.array([material._transport_factors() for material in materials])

        mu_a = np.array(
            [material.absorption.to("1/m").magnitude for material in materials]
        )

        mu_s1 = np.array(
            [material.reduced_scattering.to("1/m").magnitude for material in materials]
        )

        mu_1 = mu_a + mu_s1

        mu_2 = mu_a + mu_s1 * factors[:, 2]

        mu_3 = mu_a + mu_s1 * factors[:, 3]

        first = self.boundary_matrix.copy()

        second = self.boundary_matrix * (7 / 12)

        coupling = self.boundary_matrix * (-1 / 4)

        for label, (stiffness, mass) in self.region_matrices.items():
            first += stiffness / (3 * mu_1[label]) + mu_a[label] * mass

            second += (
                stiffness / (7 * mu_3[label])
                + (4 * mu_a[label] + 5 * mu_2[label]) / 9 * mass
            )

            coupling += (-2 / 3) * mu_a[label] * mass

        matrix = sparse.bmat([[first, coupling], [coupling, second]], format="csr")

        power = source.optical_power.to("W").magnitude

        drive = self.unit_load * power

        load = np.concatenate((drive, -drive / 4))

        solution = sparse_linalg.spsolve(matrix, load)

        if not np.all(np.isfinite(solution)):
            raise RuntimeError("The SP3 solve returned nonfinite moments.")

        relative_residual = np.linalg.norm(matrix @ solution - load) / np.linalg.norm(
            load
        )

        if relative_residual > 1e-8:
            raise RuntimeError(
                f"SP3 linear residual is too large: {relative_residual:g}."
            )

        u, v = np.split(solution, 2)

        scalar = u - 2 * v / 3

        geometry = self.geometry

        unit = MomentRT.units.ureg.watt / MomentRT.units.ureg.meter**2

        metadata = {
            "transport": "steady-state tetrahedral SP3 FEM",
            "wavelength": source.wavelength.copy(),
            "anisotropy_by_region": tuple(anisotropy),
            "coefficient_sampling": "element centroid",
            "source_face": self.source.face,
            "sample_origin": self.medium.origin.copy(),
            "sample_size": self.medium.size.copy(),
        }

        result_mesh = MomentRT.mesh.CartesianMesh(
            shape=self.mesh.shape,
            spacing=self.mesh.spacing.copy(),
            origin=self.mesh.origin.copy(),
        )

        def field(*, values: np.ndarray, name: str) -> MomentRT.fields.ScalarField:
            return MomentRT.fields.ScalarField(
                mesh=result_mesh,
                values=values.reshape(self.mesh.shape) * unit,
                name=name,
                metadata=metadata,
            )

        absorbed = np.sum(
            mu_a[self.labels]
            * geometry.volumes
            * scalar[geometry.tetrahedra].mean(axis=1)
        )

        boundary_fields: dict[str, MomentRT.fields.BoundaryField] = {}

        escaping = {}

        for face, triangles in geometry.faces.items():
            quadrature, barycentric, weights = geometry.quadrature[face]

            if isinstance(self.boundaries[face], MomentRT.boundary.ReflectingBoundary):
                outgoing = np.zeros(weights.shape)
            else:
                u_face = np.einsum("qa,ta->tq", barycentric, u[triangles])

                v_face = np.einsum("qa,ta->tq", barycentric, v[triangles])

                incoming = (
                    self.unit_irradiance * power if face == self.source.face else 0
                )

                outgoing = u_face / 2 - v_face / 8 - incoming

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

        flux_values = -np.einsum(
            "ti,tij->tj", u[geometry.tetrahedra], geometry.gradients
        ) / (3 * mu_1[self.labels, None])

        reused = self.solve_count > 0

        self.solve_count += 1

        return MomentRT.results.SimulationResult(
            fluence=field(values=scalar, name="Fluence rate"),
            flux=MomentRT.fields.VectorField(
                positions=geometry.centers * MomentRT.units.ureg.meter,
                values=flux_values * unit,
                name="SP3 flux",
            ),
            moments={
                "composite_1": field(values=u, name="SP3 composite moment 1"),
                "composite_2": field(values=v, name="SP3 composite moment 2"),
            },
            energy_balance={
                "incoming_power": source.optical_power.copy(),
                "absorbed_power": absorbed * MomentRT.units.ureg.watt,
                "escaping_power": sum(escaping.values()) * MomentRT.units.ureg.watt,
                "residual_power": (power - absorbed - sum(escaping.values()))
                * MomentRT.units.ureg.watt,
            },
            outgoing_power_by_face={
                face: value * MomentRT.units.ureg.watt
                for face, value in escaping.items()
            },
            diagnostics={
                "linear_relative_residual": relative_residual,
                "source_quadrature_area": self.source_sampled_area
                * MomentRT.units.ureg.meter**2,
                "source_geometric_area": self.source_geometric_area
                * MomentRT.units.ureg.meter**2,
                "tetrahedron_count": len(geometry.tetrahedra),
                "minimum_fluence": scalar.min() * unit,
                "preparation_reused": reused,
                "anisotropy_by_region": tuple(anisotropy),
            },
            boundary_fields=boundary_fields,
        )
