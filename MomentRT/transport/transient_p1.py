"""Time-dependent P1 finite elements with elementwise current relaxation."""

from __future__ import annotations
from collections.abc import Callable, Mapping
import numpy as np
import TypedUnit
from numpy.typing import ArrayLike
from pint import Quantity
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg
import MomentRT


class TransientP1:
    """Time-dependent P1 transport on a finite structured sample.

    Parameters
    ----------
    boundaries : mapping of str to RobinBoundary or ReflectingBoundary
        Explicit conditions on all six faces. The illuminated face must be Robin.
    times : pint.Quantity
        Strictly increasing times in seconds or compatible units, shape (N,),
        N >= 2, beginning at zero.
    power_profile : array_like of float
        Nonnegative dimensionless multipliers for source optical power, shape
        (N,), beginning at zero. Values are linearly interpolated between times.
    energy_speed : pint.Quantity
        Positive uniform transport energy speed with length/time units. It is
        supplied explicitly rather than inferred from refractive index.

    Notes
    -----
    Initial fluence and current are zero. Linear tetrahedra represent nodal
    fluence; element-centered P1 currents obey their own relaxation equation.
    Absorption and reduced scattering are sampled at element centroids, so
    layered or heterogeneous media require interface and mesh convergence.
    Time integration uses implicit midpoint with piecewise-linear input. Matched
    internal refractive indices are required. P1 is not full radiative transfer,
    especially at early times or close to a source. Instantaneous boundary
    partial currents and Lambertian detector acceptance remain approximations.
    """

    model_notes = (
        "Transient P1 uses implicit midpoint and element-centered current relaxation; it is not time-dependent diffusion or full radiative transfer.",
        "Initial fluence and current are zero; source samples define a nonnegative piecewise-linear coupled power profile.",
        "The uniform energy speed is explicit and internal refractive-index jumps are unsupported.",
        "P1 may be inaccurate at early times and near sources; check mesh and time convergence against an independent transport reference.",
    )

    def __init__(
        self,
        *,
        boundaries: Mapping[
            MomentRT.boundary.Face,
            MomentRT.boundary.RobinBoundary | MomentRT.boundary.ReflectingBoundary,
        ],
        times: Quantity,
        power_profile: ArrayLike,
        energy_speed: Quantity,
    ) -> None:
        """Configure boundary conditions, time samples, pulse, and energy speed.

        Parameters
        ----------
        boundaries : mapping of str to RobinBoundary or ReflectingBoundary
            Conditions on all six sample faces.
        times : pint.Quantity
            Increasing time samples, shape (N,), starting at zero.
        power_profile : array_like of float
            Nonnegative, dimensionless source power multipliers, shape (N,),
            starting at zero.
        energy_speed : pint.Quantity
            Positive scalar uniform speed with length/time units.
        """
        validated = MomentRT.transport.p1.P1(boundaries=boundaries)
        TypedUnit.validate_quantity(
            value=times,
            unit="s",
            name="times",
            registry=MomentRT.units.ureg,
            nonnegative=True,
        )
        if times.ndim != 1 or len(times) < 2:
            raise ValueError("times must contain at least two samples.")
        values = times.to("s").magnitude
        if values[0] != 0 or np.any(np.diff(values) <= 0):
            raise ValueError("times must start at zero and strictly increase.")
        if isinstance(power_profile, Quantity):
            raise TypeError("power_profile must contain plain dimensionless numbers.")
        profile = np.asarray(power_profile)
        if (
            profile.dtype.kind not in "iuf"
            or profile.shape != times.shape
            or not np.all(np.isfinite(profile))
            or np.any(profile < 0)
            or profile[0] != 0
        ):
            raise ValueError(
                "power_profile must be finite, nonnegative, match times, and start at zero."
            )
        TypedUnit.validate_quantity(
            value=energy_speed,
            unit="m/s",
            name="energy_speed",
            registry=MomentRT.units.ureg,
            shape=(),
            positive=True,
        )
        self.boundaries = validated.boundaries
        self.times = times.copy()
        self.power_profile = profile.astype(float, copy=True)
        self.energy_speed = energy_speed.copy()

    def solve(
        self,
        *,
        medium: MomentRT.medium.StructuredMedium,
        source: MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> MomentRT.results.TransientResult:
        """Integrate coupled P1 fluence and current on a fixed finite mesh.

        Parameters
        ----------
        medium : StructuredMedium
            Finite sample with explicit materials and matched internal indices.
        source : UniformSurfaceSource or SurfaceSource
            Monochromatic diffuse source with explicit coupled optical power.
        mesh : CartesianMesh or TetrahedralMesh
            Finite sample mesh used for spatial discretization.

        Returns
        -------
        TransientResult
            Per-time nodal fluence, element current, boundary fields, and
            interval energy diagnostics. Detector signals are added by Simulation.
        """
        if not isinstance(medium, MomentRT.medium.StructuredMedium):
            raise TypeError("TransientP1 requires a finite StructuredMedium.")
        if not isinstance(source, MomentRT.source.UniformSurfaceSource):
            raise TypeError("TransientP1 requires a surface source.")
        if source.spectrum is not None:
            raise ValueError("TransientP1 requires a monochromatic source.")
        preparation = MomentRT.transport.p1.P1(boundaries=self.boundaries)._prepare(
            medium=medium,
            source=source,
            mesh=mesh,
        )
        materials = MomentRT.transport.p1._finite_materials_at(
            medium=preparation.medium,
            wavelength=source.wavelength,
        )
        absorption = np.array(
            [material.absorption.to("1/m").magnitude for material in materials]
        )
        transport = np.array(
            [
                (material.absorption + material.reduced_scattering).to("1/m").magnitude
                for material in materials
            ]
        )
        geometry = preparation.geometry
        matrices = preparation.region_matrices
        mass = sum(
            (mass for _, mass in matrices.values()),
            start=sparse.csr_matrix(preparation.boundary_matrix.shape),
        )
        absorption_matrix = sum(
            (absorption[label] * pair[1] for label, pair in matrices.items()),
            start=sparse.csr_matrix(preparation.boundary_matrix.shape),
        )
        velocity = float(self.energy_speed.to("m/s").magnitude)
        times = self.times.to("s").magnitude
        powers = self.power_profile * source.optical_power.to("W").magnitude
        values = np.zeros(len(geometry.points))
        current = np.zeros((len(geometry.tetrahedra), 3))
        frames = []
        stored = []
        absorbed = []
        escaping = []
        ones = np.ones(len(values))
        storage_weights = np.asarray(mass @ ones).ravel() / velocity
        factors: dict[float, Callable[[np.ndarray], np.ndarray]] = {}
        for index, power in enumerate(powers):
            relative_residual = 0.0
            if index:
                step = float(times[index] - times[index - 1])
                coefficient = step * velocity * transport / 2
                diffusion = step * velocity / (6 * (1 + coefficient))
                operator = preparation.boundary_matrix + absorption_matrix
                for label, (stiffness, _) in matrices.items():
                    operator = operator + diffusion[label] * stiffness
                left = mass / (velocity * step) + operator / 2
                right = (mass / (velocity * step) - operator / 2) @ values
                relaxed = current / (1 + coefficient[preparation.labels, None])
                local = geometry.volumes[:, None] * np.einsum(
                    "tij,tj->ti", geometry.gradients, relaxed
                )
                np.add.at(right, geometry.tetrahedra.ravel(), local.ravel())
                right += preparation.unit_load * ((powers[index - 1] + power) / 2)
                if step not in factors:
                    factors[step] = sparse_linalg.factorized(left.tocsc())
                new_values = factors[step](right)
                if not np.all(np.isfinite(new_values)):
                    raise RuntimeError("TransientP1 returned nonfinite fluence.")
                residual_norm = np.linalg.norm(left @ new_values - right)
                right_norm = np.linalg.norm(right)
                relative_residual = (
                    float(residual_norm / right_norm)
                    if right_norm
                    else float(residual_norm)
                )
                if relative_residual > 1e-8:
                    raise RuntimeError(
                        f"TransientP1 linear residual is too large: {relative_residual:g}."
                    )
                gradient = np.einsum(
                    "ti,tij->tj",
                    (new_values + values)[geometry.tetrahedra],
                    geometry.gradients,
                )
                current = (
                    (1 - coefficient[preparation.labels])
                    / (1 + coefficient[preparation.labels])
                )[:, None] * current - diffusion[preparation.labels, None] * gradient
                values = new_values
            frame, energy, absorption_power, escaped = self._frame(
                preparation=preparation,
                source=source,
                values=values,
                current=current,
                power=float(power),
                storage_weights=storage_weights,
                absorption_matrix=absorption_matrix,
                time=float(times[index]),
                relative_residual=relative_residual,
            )
            frames.append(frame)
            stored.append(energy)
            absorbed.append(absorption_power)
            escaping.append(escaped)
        stored_values = np.asarray(stored)
        net = np.asarray(absorbed) + np.asarray(escaping) - powers
        residual = np.diff(stored_values) / np.diff(times) + (net[:-1] + net[1:]) / 2
        unit = MomentRT.units.ureg
        return MomentRT.results.TransientResult(
            times=self.times,
            frames=tuple(frames),
            stored_energy=stored_values * unit.joule,
            absorbed_power=np.asarray(absorbed) * unit.watt,
            incoming_power=powers * unit.watt,
            escaping_power=np.asarray(escaping) * unit.watt,
            interval_energy_residual=residual * unit.watt,
        )

    def _frame(
        self,
        *,
        preparation: MomentRT.transport.p1._PreparedDiffusion,
        source: MomentRT.source.UniformSurfaceSource,
        values: np.ndarray,
        current: np.ndarray,
        power: float,
        storage_weights: np.ndarray,
        absorption_matrix: sparse.csr_matrix,
        time: float,
        relative_residual: float,
    ) -> tuple[MomentRT.results.SimulationResult, float, float, float]:
        """Construct one field snapshot and integrated SI diagnostics.

        Parameters
        ----------
        preparation : _PreparedDiffusion
            Fixed mesh, material labels, and boundary quadrature.
        source : UniformSurfaceSource
            Monochromatic illumination defining the source face and wavelength.
        values : ndarray of float, shape (N,)
            Nodal fluence in W/m².
        current : ndarray of float, shape (T, 3)
            Element-centered P1 current in W/m².
        power : float
            Instantaneous incoming coupled optical power in watts.
        storage_weights : ndarray of float, shape (N,)
            Nodal integration weights divided by energy speed, in m² s.
        absorption_matrix : scipy.sparse.csr_matrix
            Absorption-weighted mass matrix, in m².
        time : float
            Snapshot time in seconds.
        relative_residual : float
            Dimensionless linear-system relative residual for this step.

        Returns
        -------
        tuple of SimulationResult and float
            Snapshot, stored energy in joules, absorbed power in watts, and
            signed total escaping power in watts.
        """
        unit = MomentRT.units.ureg
        geometry = preparation.geometry
        boundary_fields: dict[str, MomentRT.fields.BoundaryField] = {}
        for face, triangles in geometry.faces.items():
            positions, barycentric, weights = geometry.quadrature[face]
            face_fluence = np.einsum("qa,ta->tq", barycentric, values[triangles])
            incoming = preparation.unit_irradiance * power if face == source.face else 0
            beta = preparation.boundaries[face].leakage_coefficient
            boundary_fields[face] = MomentRT.fields.BoundaryField(
                face=face,
                positions=positions * unit.meter,
                area_weights=weights * unit.meter**2,
                escaping_irradiance=beta
                * (face_fluence - 2 * incoming)
                * unit.watt
                / unit.meter**2,
                sample_origin=preparation.medium.origin,
                sample_size=preparation.medium.size,
                wavelength=source.wavelength,
            )
        field = MomentRT.fields.ScalarField(
            mesh=preparation.mesh,
            values=values.reshape(preparation.mesh.shape).copy()
            * unit.watt
            / unit.meter**2,
            name="Fluence rate",
            metadata={
                "transport": "time-dependent P1 tetrahedral FEM",
                "time": time * unit.second,
                "wavelength": source.wavelength.copy(),
                "energy_speed": self.energy_speed.copy(),
                "time_integrator": "implicit midpoint",
                "coefficient_sampling": "element centroid",
                "source_face": source.face,
                "sample_origin": preparation.medium.origin.copy(),
                "sample_size": preparation.medium.size.copy(),
            },
        )
        flux = MomentRT.fields.VectorField(
            positions=geometry.centers * unit.meter,
            values=current.copy() * unit.watt / unit.meter**2,
            name="P1 current",
        )
        escaped = sum(
            field.escaping_power.to("W").magnitude for field in boundary_fields.values()
        )
        frame = MomentRT.results.SimulationResult(
            fluence=field,
            flux=flux,
            boundary_fields=boundary_fields,
            outgoing_power_by_face={
                face: boundary.escaping_power
                for face, boundary in boundary_fields.items()
            },
            diagnostics={
                "time": time * unit.second,
                "linear_relative_residual": relative_residual,
                "source_quadrature_area": preparation.source_sampled_area
                * unit.meter**2,
                "source_geometric_area": preparation.source_geometric_area
                * unit.meter**2,
                "minimum_fluence": values.min() * unit.watt / unit.meter**2,
            },
        )
        return (
            frame,
            float(storage_weights @ values),
            float(np.sum(absorption_matrix @ values)),
            float(escaped),
        )
