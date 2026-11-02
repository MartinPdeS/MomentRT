"""Results retain physical fields and their provenance."""

from __future__ import annotations
import copy
from collections.abc import Mapping, Sequence
from typing import cast
import numpy as np
from pint import Quantity
from scipy import integrate
import MomentRT


def _energy_accounting(
    *,
    incoming_power: Quantity,
    absorbed_watts: float,
    escaping_watts: Mapping[MomentRT.boundary.Face, float],
) -> tuple[dict[str, Quantity], dict[str, Quantity]]:
    """Package finite-solver power accounting without changing signed outputs.

    Parameters
    ----------
    incoming_power : pint.Quantity
        Incoming coupled power, with power units.
    absorbed_watts : float
        Integrated absorption in watts.
    escaping_watts : mapping of str to float
        Signed escaping power by face in watts.

    Returns
    -------
    tuple of dict
        Energy balance and per-face escaping powers, both with power units.
    """
    watt = MomentRT.units.ureg.watt

    escaping_total = sum(escaping_watts.values())

    balance = {
        "incoming_power": incoming_power.copy(),
        "absorbed_power": absorbed_watts * watt,
        "escaping_power": escaping_total * watt,
        "residual_power": (
            incoming_power.to("W").magnitude - absorbed_watts - escaping_total
        )
        * watt,
    }

    by_face: dict[str, Quantity] = {
        face: value * watt for face, value in escaping_watts.items()
    }

    return balance, by_face


class SimulationResult:
    """Monochromatic transport fields, energy accounting, and measurements.

    Parameters
    ----------
    fluence : ScalarField
        Nodal fluence rate, with power-per-area units.
    flux : VectorField or None, optional
        Element-centered diffuse current density, with power-per-area units.
        Default is None.
    moments : mapping of str to ScalarField, or None, optional
        Additional solver-specific nodal moment fields.
    energy_balance : mapping of str to pint.Quantity, or None, optional
        Incoming, absorbed, net outward, and residual power accounting.
        Default is None.
    outgoing_power_by_face : mapping of str to pint.Quantity, or None, optional
        Escaping transmitted power per complete sample face; not a finite
        aperture measurement.
        Default is None.
    diagnostics : mapping of str to object, or None, optional
        Numerical provenance, residuals, geometry statistics, or timing quantities.
        Default is None.
    boundary_fields : mapping of str to BoundaryField, or None, optional
        Per-face snapshots of escaping-irradiance quadrature samples.
        Default is None.
    detectors : mapping of str to DetectorResult, or None, optional
        Named finite-aperture measurements already evaluated on this solution.
        Default is None.

    Attributes
    ----------
    fluence : ScalarField
        Nodal fluence rate with power-per-area units.
    flux : VectorField or None
        Element-centered current density, absent for the infinite Green function.
    boundary_fields : dict of str to BoundaryField
        Quadrature snapshots for finite-sample contact measurements.
    detectors : dict of str to DetectorResult
        Evaluated measurements, keyed by detector name.
    diagnostics : dict of str to object
        Solver provenance, residuals, and resolution statistics.
    """

    def __init__(
        self,
        *,
        fluence: MomentRT.fields.ScalarField,
        flux: MomentRT.fields.VectorField | None = None,
        moments: Mapping[str, MomentRT.fields.ScalarField] | None = None,
        energy_balance: Mapping[str, Quantity] | None = None,
        outgoing_power_by_face: Mapping[str, Quantity] | None = None,
        diagnostics: Mapping[str, object] | None = None,
        boundary_fields: Mapping[str, MomentRT.fields.BoundaryField] | None = None,
        detectors: Mapping[str, MomentRT.results.DetectorResult] | None = None,
    ) -> None:
        """Initialize the SimulationResult configuration.

        Parameters
        ----------
        fluence : ScalarField
            Nodal fluence rate, with power-per-area units.
        flux : VectorField or None, optional
            Element-centered diffuse current density, with power-per-area units.
            Default is None.
        moments : mapping of str to ScalarField, or None, optional
            Additional solver-specific nodal moment fields.
        energy_balance : mapping of str to pint.Quantity, or None, optional
            Incoming, absorbed, net outward, and residual power accounting.
            Default is None.
        outgoing_power_by_face : mapping of str to pint.Quantity, or None, optional
            Escaping transmitted power per complete sample face; not a finite
            aperture measurement.
            Default is None.
        diagnostics : mapping of str to object, or None, optional
            Numerical provenance, residuals, geometry statistics, or timing quantities.
            Default is None.
        boundary_fields : mapping of str to BoundaryField, or None, optional
            Per-face snapshots of escaping-irradiance quadrature samples.
            Default is None.
        detectors : mapping of str to DetectorResult, or None, optional
            Named finite-aperture measurements already evaluated on this solution.
            Default is None.
        """

        self.fluence = fluence

        self.flux = flux

        self.moments = dict(moments or {})

        self.energy_balance = dict(energy_balance or {})

        self.outgoing_power_by_face = dict(outgoing_power_by_face or {})

        self.diagnostics = dict(diagnostics or {})

        self.boundary_fields = dict(boundary_fields or {})

        self.detectors = dict(detectors or {})

    @property
    def source_sampled_area(self) -> Quantity | None:
        """Return the sampled source footprint area, when available.

        Returns
        -------
        pint.Quantity or None
            Sampled area with physical units, or None for an infinite solution.
        """
        return cast(Quantity | None, self.diagnostics.get("source_quadrature_area"))

    @property
    def source_geometric_area(self) -> Quantity | None:
        """Return the exact source footprint area, when available.

        Returns
        -------
        pint.Quantity or None
            Geometric area with physical units, or None for an infinite solution.
        """
        return cast(Quantity | None, self.diagnostics.get("source_geometric_area"))

    @property
    def linear_relative_residual(self) -> float | None:
        """Return the relative linear-system residual, when available.

        Returns
        -------
        float or None
            Relative residual, or None for an analytical solution.
        """
        return cast(float | None, self.diagnostics.get("linear_relative_residual"))

    @property
    def minimum_fluence(self) -> Quantity | None:
        """Return the minimum nodal fluence, when recorded by the solver.

        Returns
        -------
        pint.Quantity or None
            Minimum fluence rate; negative values are not clipped.
        """
        return cast(Quantity | None, self.diagnostics.get("minimum_fluence"))

    def summary(
        self, *, convergence: MomentRT.results.DetectorConvergenceResult | None = None
    ) -> str:
        """Format detector powers, sampling, energy balance, and diagnostics.

        Parameters
        ----------
        convergence : DetectorConvergenceResult or None, optional
            Separately supplied convergence study to summarize. The caller must
            verify it applies to this result; None means none was supplied.
            Default is None.

        Returns
        -------
        str
            Human-readable result summary with the stored physical units.
        """
        lines = ["Transport result"]

        transport = self.fluence.metadata.get("transport")

        if transport is not None:
            lines.append(f"Model: {transport}")

        wavelength = self.fluence.metadata.get("wavelength")

        if isinstance(wavelength, Quantity):
            lines.append(f"Wavelength: {wavelength:~P}")

        lines.append(f"Mesh shape: {self.fluence.mesh.shape}")

        if convergence is None:
            lines.append(
                "Convergence study: not supplied; convergence is not established."
            )
        else:
            lines.append(
                f"Supplied convergence study: {'converged' if convergence.converged else 'not converged'}; "
                f"recommended mesh: {convergence.recommended_shape}. "
                "Verify that the study applies to this result."
            )

        for name, measurement in self.detectors.items():
            lines.append(
                f"Detector {name}: {measurement.collected_power:~P} (sampled area {measurement.sampled_area:~P}; geometric area {measurement.geometric_area:~P})"
            )

            if measurement.minimum_escaping_irradiance.magnitude < 0:
                lines.append(
                    f"Warning: detector {name} samples negative escaping irradiance."
                )

            area_error = abs(
                (measurement.sampled_area / measurement.geometric_area).to("").magnitude
                - 1
            )

            if area_error > 0.02:
                lines.append(
                    f"Warning: detector {name} sampled area differs from geometric area by {area_error:.1%}."
                )

        if (
            self.source_sampled_area is not None
            and self.source_geometric_area is not None
        ):
            lines.append(
                f"Source sampled area: {self.source_sampled_area:~P}; geometric area: {self.source_geometric_area:~P}"
            )

            area_error = abs(
                (self.source_sampled_area / self.source_geometric_area).to("").magnitude
                - 1
            )

            if area_error > 0.02:
                lines.append(
                    f"Warning: source sampled area differs from geometric area by {area_error:.1%}."
                )

        for name, power in self.energy_balance.items():
            lines.append(f"{name.replace('_', ' ').capitalize()}: {power:~P}")

        if self.linear_relative_residual is not None:
            lines.append(
                f"Linear relative residual: {self.linear_relative_residual:.3g}"
            )

        if self.minimum_fluence is not None and self.minimum_fluence.magnitude < 0:
            lines.append(
                f"Warning: minimum fluence is negative: {self.minimum_fluence:~P}"
            )

        return "\n".join(lines)

    @property
    def detector_power(self) -> dict[str, Quantity]:
        """Return collected powers keyed by detector name.

        Returns
        -------
        dict of str to pint.Quantity
            Collected optical powers keyed by detector name.
        """

        return {
            name: measurement.collected_power
            for name, measurement in self.detectors.items()
        }


class DetectorResult:
    """Finite-aperture diffuse collection and quadrature provenance.

    Parameters
    ----------
    collected_power : pint.Quantity
        Integrated accepted diffuse optical power, with power units.
    sampled_area : pint.Quantity
        Aperture area represented by selected quadrature weights.
        It is not normalized to the geometric area.
    geometric_area : pint.Quantity
        Exact circular aperture area, with area units.
    angular_fraction : float
        Fraction of Lambertian escaping irradiance accepted by the cone.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.
    minimum_escaping_irradiance : pint.Quantity
        Smallest sampled escaping irradiance inside the aperture.
        Negative values flag discretization artifacts and remain unchanged.
    detector : SurfaceDetector
        Snapshot of the aperture and angular collection configuration.
    photocurrent : pint.Quantity or None, optional
        Linear electrical current in amperes when responsivity is configured.
        Default is None.
    """

    def __init__(
        self,
        *,
        collected_power: Quantity,
        sampled_area: Quantity,
        geometric_area: Quantity,
        angular_fraction: float,
        wavelength: Quantity,
        minimum_escaping_irradiance: Quantity,
        detector: MomentRT.detector.SurfaceDetector,
        photocurrent: Quantity | None = None,
    ) -> None:
        """Initialize the DetectorResult configuration.

        Parameters
        ----------
        collected_power : pint.Quantity
            Integrated accepted diffuse optical power, with power units.
        sampled_area : pint.Quantity
            Aperture area represented by selected quadrature weights.
            It is not normalized to the geometric area.
        geometric_area : pint.Quantity
            Exact circular aperture area, with area units.
        angular_fraction : float
            Fraction of Lambertian escaping irradiance accepted by the cone.
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.
        minimum_escaping_irradiance : pint.Quantity
            Smallest sampled escaping irradiance inside the aperture.
            Negative values flag discretization artifacts and remain unchanged.
        detector : SurfaceDetector
            Snapshot of the aperture and angular collection configuration.
        photocurrent : pint.Quantity or None, optional
            Linear electrical current in amperes when responsivity is configured.
            Default is None.
        """

        self.photocurrent = photocurrent

        self.collected_power = collected_power

        self.sampled_area = sampled_area

        self.geometric_area = geometric_area

        self.angular_fraction = angular_fraction

        self.wavelength = wavelength.copy()

        self.minimum_escaping_irradiance = minimum_escaping_irradiance

        self.detector = detector

        self.angular_model = "Lambertian escaping radiance"


class DetectorConvergenceResult:
    """Detector signals and footprint sampling on successively refined meshes.

    Parameters
    ----------
    meshes : sequence of CartesianMesh or TetrahedralMesh
        Coarse-to-fine meshes used in the study.
    results : sequence of SimulationResult or None
        Per-level solutions; None records an unresolved source footprint.
    detectors : mapping of str to SurfaceDetector
        Named passive contact detectors.
    source_geometric_area : pint.Quantity
        Exact source footprint area, with area units.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.
    level_errors : sequence of mapping of str to str
        Per-level resolution failures and negative-field diagnostics.
    power_rtol : float
        Nonnegative relative power tolerance against adjacent and finest meshes.
    area_rtol : float
        Nonnegative relative tolerance of sampled source and detector areas
        against their geometric areas.
    power_atol : pint.Quantity
        Nonnegative absolute power tolerance, with power units.
    required_consecutive : int
        Positive number of consecutive passing refinements required at the
        fine end of the mesh sequence.

    Attributes
    ----------
    node_counts : numpy.ndarray of int
        Number of nodes at each refinement level.
    tetrahedron_counts : numpy.ndarray of int
        Number of tetrahedra at each refinement level.
    source_sampled_area : pint.Quantity
        Sampled illumination areas in square meters, one per level.
    source_area_relative_error : numpy.ndarray of float
        Absolute relative errors against the exact source footprint area.
    detector_power : dict of str to pint.Quantity
        Per-detector collected power vectors in watts. Unavailable powers are NaN.
    detector_sampled_area : dict of str to pint.Quantity
        Per-detector sampled aperture areas in square meters. Unsampled apertures
        have zero sampled area; unresolved sources produce NaN.
    detector_geometric_area : dict of str to pint.Quantity
        Exact aperture area for each detector.
    detector_area_relative_error : dict of str to numpy.ndarray
        Absolute relative aperture-area errors at each mesh level.
    power_relative_change : dict of str to numpy.ndarray
        Relative power differences against the preceding level. The first
        difference and unavailable values are NaN.
    power_relative_to_finest : dict of str to numpy.ndarray
        Relative power differences against the finest level; unavailable
        values are NaN.
    meets_tolerances : numpy.ndarray of bool
        Per-level power, area, and resolution checks.
    stable_refinements : numpy.ndarray of int
        Number of consecutive passing refinements ending at each level.
    recommended_index : int or None
        Index of the recommended mesh; None if the final sequence fails.

    Notes
    -----
    Only a stable final sequence can support a recommendation. Sampled
    source and detector areas must approximate their exact geometric areas,
    and powers must meet adjacent-level and finest-level tolerances.
    Unresolved levels and negative fields remain visible; passing these
    checks is not proof of accuracy of the P1 approximation.
    """

    def __init__(
        self,
        *,
        meshes: Sequence[MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh],
        results: Sequence[MomentRT.results.SimulationResult | None],
        detectors: Mapping[str, MomentRT.detector.SurfaceDetector],
        source_geometric_area: Quantity,
        wavelength: Quantity,
        level_errors: Sequence[Mapping[str, str]],
        power_rtol: float,
        area_rtol: float,
        power_atol: Quantity,
        required_consecutive: int,
    ) -> None:
        """Initialize the DetectorConvergenceResult configuration.

        Parameters
        ----------
        meshes : sequence of CartesianMesh or TetrahedralMesh
            Coarse-to-fine meshes used in the study.
        results : sequence of SimulationResult or None
            Per-level solutions; None records an unresolved source footprint.
        detectors : mapping of str to SurfaceDetector
            Named passive contact detectors.
        source_geometric_area : pint.Quantity
            Exact source footprint area, with area units.
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.
        level_errors : sequence of mapping of str to str
            Per-level resolution failures and negative-field diagnostics.
        power_rtol : float
            Nonnegative relative power tolerance against adjacent and finest meshes.
        area_rtol : float
            Nonnegative relative tolerance of sampled source and detector areas
            against their geometric areas.
        power_atol : pint.Quantity
            Nonnegative absolute power tolerance, with power units.
        required_consecutive : int
            Positive number of consecutive passing refinements required at the
            fine end of the mesh sequence.
        """

        self.meshes = tuple(meshes)

        self.mesh_shapes = tuple(mesh.shape for mesh in meshes)

        self.results = tuple(results)

        self.detectors = dict(detectors)

        self.wavelength = wavelength.copy()

        self.level_errors = tuple(dict(errors) for errors in level_errors)

        self.power_rtol = float(power_rtol)

        self.area_rtol = float(area_rtol)

        self.power_atol = power_atol.copy()

        self.required_consecutive = int(required_consecutive)

        self.node_counts = np.array([np.prod(mesh.shape) for mesh in meshes], dtype=int)

        self.tetrahedron_counts = np.array(
            [
                len(mesh.tetrahedra)
                if isinstance(mesh, MomentRT.mesh.TetrahedralMesh)
                else 6 * np.prod(np.asarray(mesh.shape) - 1)
                for mesh in meshes
            ],
            dtype=int,
        )

        self.mesh_resolutions = tuple(mesh.resolution for mesh in meshes)

        self.source_geometric_area = source_geometric_area.copy()

        self.source_sampled_area = (
            np.array(
                [
                    0
                    if result is None
                    else cast(Quantity, result.diagnostics["source_quadrature_area"])
                    .to("m^2")
                    .magnitude
                    for result in results
                ]
            )
            * MomentRT.units.ureg.meter**2
        )

        self.source_area_relative_error = np.abs(
            (self.source_sampled_area / self.source_geometric_area).to("").magnitude - 1
        )

        self.detector_power = {}

        self.detector_sampled_area = {}

        self.detector_geometric_area = {}

        self.detector_area_relative_error = {}

        self.power_relative_change = {}

        self.power_within_tolerance = {}

        self.power_relative_to_finest = {}

        within_finest_tolerance = {}

        for name, detector in detectors.items():
            measurements = [
                None if result is None else result.detectors.get(name)
                for result in results
            ]

            powers = np.array(
                [
                    np.nan if value is None else value.collected_power.to("W").magnitude
                    for value in measurements
                ],
                dtype=float,
            )

            self.detector_power[name] = powers * MomentRT.units.ureg.watt

            self.detector_sampled_area[name] = (
                np.array(
                    [
                        np.nan
                        if result is None
                        else (
                            0
                            if value is None
                            else value.sampled_area.to("m^2").magnitude
                        )
                        for result, value in zip(results, measurements)
                    ]
                )
                * MomentRT.units.ureg.meter**2
            )

            self.detector_geometric_area[name] = np.pi * detector.radius**2

            self.detector_area_relative_error[name] = np.abs(
                (self.detector_sampled_area[name] / self.detector_geometric_area[name])
                .to("")
                .magnitude
                - 1
            )

            difference = np.abs(np.diff(powers))

            denominator = np.abs(powers[1:])

            finite_pair = np.isfinite(powers[1:]) & np.isfinite(powers[:-1])

            relative = np.divide(
                difference,
                denominator,
                out=np.full_like(difference, np.inf),
                where=denominator != 0,
            )

            relative[(denominator == 0) & (difference == 0)] = 0

            relative[~finite_pair] = np.nan

            self.power_relative_change[name] = np.concatenate(([np.nan], relative))

            limit = self.power_atol.to("W").magnitude + self.power_rtol * denominator

            self.power_within_tolerance[name] = np.concatenate(
                ([False], finite_pair & (difference <= limit))
            )

            finest_difference = np.abs(powers - powers[-1])

            reference = abs(powers[-1])

            finest_relative = np.divide(
                finest_difference,
                reference,
                out=np.full_like(powers, np.inf),
                where=reference != 0,
            )

            finest_relative[(reference == 0) & (finest_difference == 0)] = 0

            finest_relative[~np.isfinite(powers) | ~np.isfinite(reference)] = np.nan

            self.power_relative_to_finest[name] = finest_relative

            within_finest_tolerance[name] = (
                np.isfinite(powers)
                & np.isfinite(reference)
                & (
                    finest_difference
                    <= self.power_atol.to("W").magnitude + self.power_rtol * reference
                )
            )

        acceptable = np.isfinite(self.source_area_relative_error) & (
            self.source_area_relative_error <= self.area_rtol
        )

        for name in detectors:
            area_error = self.detector_area_relative_error[name]

            acceptable &= (
                self.power_within_tolerance[name]
                & np.isfinite(area_error)
                & (area_error <= self.area_rtol)
            )

        acceptable &= np.array([not errors for errors in self.level_errors])

        self.meets_tolerances = acceptable

        self.stable_refinements = np.zeros(len(meshes), dtype=int)

        self.recommended_index = None

        for index, passed in enumerate(acceptable):
            previous = self.stable_refinements[index - 1] if index else 0

            self.stable_refinements[index] = previous + 1 if passed else 0

        # Only recommend from the final uninterrupted sequence of passing checks.
        # A later failure invalidates an apparent plateau on coarser meshes.
        final_streak = self.stable_refinements[-1]

        if final_streak >= self.required_consecutive:
            first_candidate = len(meshes) - final_streak + self.required_consecutive - 1

            for index in range(first_candidate, len(meshes)):
                if all(within_finest_tolerance[name][index] for name in detectors):
                    self.recommended_index = index

                    break

    @property
    def converged(self) -> bool:
        """Report whether the final refinement sequence supports a mesh recommendation.

        Returns
        -------
        bool
            Whether the final stable refinement sequence supports a recommendation.
        """

        return self.recommended_index is not None

    @property
    def recommended_mesh(
        self,
    ) -> MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh | None:
        """Return the recommended mesh when the convergence criteria are met.

        Returns
        -------
        CartesianMesh, TetrahedralMesh, or None
            Recommended existing mesh, or None when criteria are not met.
        """

        return (
            None
            if self.recommended_index is None
            else self.meshes[self.recommended_index]
        )

    @property
    def recommended_shape(self) -> tuple[int, int, int] | None:
        """Return the node shape of the recommended mesh, if any.

        Returns
        -------
        tuple of int or None
            Recommended node counts, or None when no mesh is recommended.
        """

        if not isinstance(self.recommended_mesh, MomentRT.mesh.CartesianMesh):
            return None

        return self.recommended_mesh.shape

    @property
    def recommended_resolution(self) -> Quantity | None:
        """Return the target size of the recommended generated tetrahedral mesh.

        Returns
        -------
        pint.Quantity or None
            Scalar target element size, or None for Cartesian or imported meshes.
        """

        if not isinstance(self.recommended_mesh, MomentRT.mesh.TetrahedralMesh):
            return None

        return self.recommended_mesh.resolution

    @property
    def records(self) -> list[dict[str, object]]:
        """Return per-level records with quantities intact and explicit failure flags.

        Returns
        -------
        list of dict
            Per-level quantities, detector diagnostics, and explicit failure flags.
        """
        return [
            {
                "mesh_shape": shape,
                "mesh_resolution": self.mesh_resolutions[index],
                "node_count": int(self.node_counts[index]),
                "source_sampled_area": self.source_sampled_area[index],
                "source_area_relative_error": float(
                    self.source_area_relative_error[index]
                ),
                "detectors": {
                    name: {
                        "collected_power": self.detector_power[name][index],
                        "sampled_area": self.detector_sampled_area[name][index],
                        "area_relative_error": float(
                            self.detector_area_relative_error[name][index]
                        ),
                        "power_relative_change": float(
                            self.power_relative_change[name][index]
                        ),
                        "power_relative_to_finest": float(
                            self.power_relative_to_finest[name][index]
                        ),
                    }
                    for name in self.detectors
                },
                "meets_tolerances": bool(self.meets_tolerances[index]),
                "stable_refinements": int(self.stable_refinements[index]),
                "errors": dict(self.level_errors[index]),
            }
            for index, shape in enumerate(self.mesh_shapes)
        ]

    def plot(
        self, *, power_unit: str = "mW", area_unit: str = "mm^2", show: bool = True
    ) -> MomentRT.plotting.Figure:
        """Compare powers, aperture areas, and source area against mesh node count.

        Parameters
        ----------
        power_unit : str, optional
            Power unit used for plotted detector signals.
            Default is 'mW'.
        area_unit : str, optional
            Area unit used for plotted sampled and geometric areas.
            Default is 'mm^2'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable interactive figure, also displayed when show is True.
        """

        go = MomentRT.plotting.go

        figure = MomentRT.plotting.make_subplots(
            rows=1,
            cols=3,
            subplot_titles=(
                "Collected power",
                "Detector sampled area",
                "Source sampled area",
            ),
        )

        for number, name in enumerate(self.detectors):
            color = MomentRT.plotting.qualitative.Plotly[
                number % len(MomentRT.plotting.qualitative.Plotly)
            ]

            for column, values in (
                (1, self.detector_power[name].to(power_unit).magnitude),
                (2, self.detector_sampled_area[name].to(area_unit).magnitude),
            ):
                figure.add_trace(
                    go.Scatter(
                        x=self.node_counts,
                        y=values,
                        mode="lines+markers",
                        name=name,
                        legendgroup=f"detector:{name}",
                        showlegend=column == 1,
                        line={"color": color},
                    ),
                    row=1,
                    col=column,
                )

            figure.add_trace(
                go.Scatter(
                    x=self.node_counts,
                    y=np.full(
                        len(self.meshes),
                        self.detector_geometric_area[name].to(area_unit).magnitude,
                    ),
                    mode="lines",
                    name=f"{name} geometric area",
                    legendgroup=f"detector:{name}",
                    showlegend=False,
                    line={"color": color, "dash": "dash"},
                ),
                row=1,
                col=2,
            )

        for label, values, dash in (
            (
                "Source sampled area",
                self.source_sampled_area.to(area_unit).magnitude,
                "solid",
            ),
            (
                "Source geometric area",
                np.full(
                    len(self.meshes), self.source_geometric_area.to(area_unit).magnitude
                ),
                "dash",
            ),
        ):
            figure.add_trace(
                go.Scatter(
                    x=self.node_counts,
                    y=values,
                    mode="lines+markers",
                    name=label,
                    legendgroup="source",
                    showlegend=dash == "solid",
                    line={"color": "crimson", "dash": dash},
                ),
                row=1,
                col=3,
            )

        figure.update_xaxes(
            title_text="Mesh nodes",
            type="log",
            tickmode="array",
            tickvals=self.node_counts,
            ticktext=[f"{int(count):,}" for count in self.node_counts],
            tickangle=-35,
        )

        figure.update_yaxes(title_text=f"Power [{power_unit}]", row=1, col=1)

        figure.update_yaxes(title_text=f"Area [{area_unit}]", row=1, col=2)

        figure.update_yaxes(title_text=f"Area [{area_unit}]", row=1, col=3)

        recommendation = (
            "no mesh meets the criteria"
            if not self.converged
            else f"recommended {self.recommended_shape}"
        )

        figure.update_layout(
            template="plotly_white",
            title={
                "text": f"Detector refinement at {self.wavelength:~P} — {recommendation}",
                "y": 0.98,
                "yanchor": "top",
                "x": 0.06,
            },
            legend={
                "orientation": "h",
                "y": 1.08,
                "yanchor": "bottom",
                "x": 0,
            },
            margin={"t": 120},
        )

        if show:
            figure.show()

        return figure


class SpectralConvergenceResult:
    """Detector power and photocurrent stability on nested wavelength grids.

    Parameters
    ----------
    wavelength_grids : sequence of pint.Quantity
        Increasing-resolution wavelength vectors over one fixed source interval.
    results : sequence of BroadbandSimulationResult or None
        Per-level predictions. None retains a failed refinement.
    detectors : mapping of str to SurfaceDetector
        Fixed detector geometry and optional electrical response.
    level_errors : sequence of mapping of str to str
        Per-level solve failures and negative-field diagnostics.
    power_rtol : float
        Nonnegative relative optical-power tolerance.
    photocurrent_rtol : float
        Nonnegative relative electrical-current tolerance.
    power_atol : pint.Quantity
        Nonnegative absolute optical-power tolerance.
    photocurrent_atol : pint.Quantity
        Nonnegative absolute electrical-current tolerance.
    required_consecutive : int
        Consecutive passing refinements required at the fine end.

    Attributes
    ----------
    wavelength_grids : tuple of pint.Quantity
        Actual nested integration grids after required table knots are merged.
    sample_counts : numpy.ndarray of int
        Number of wavelengths at each level.
    detector_power : dict of str to pint.Quantity
        Integrated optical-power vectors in watts, with NaN for failed levels.
    detector_photocurrent : dict of str to pint.Quantity
        Integrated current vectors in amperes for detectors with responsivity.
    power_relative_change : dict of str to numpy.ndarray
        Adjacent power changes divided by the larger absolute adjacent value.
    photocurrent_relative_change : dict of str to numpy.ndarray
        Adjacent current changes using the same normalization as power changes.
    power_relative_to_finest : dict of str to numpy.ndarray
        Absolute power differences divided by the finest absolute power.
    photocurrent_relative_to_finest : dict of str to numpy.ndarray
        Current differences relative to the finest absolute current.
    meets_tolerances : numpy.ndarray of bool
        Per-level stability checks including failures and negative outputs.
    recommended_index : int or None
        Recommended grid index, or None after an unresolved final sequence.

    Notes
    -----
    A recommendation requires adjacent-level stability and agreement with the
    finest result for every measured quantity. Failed or negative outputs are
    retained and prevent a recommendation from a failed final sequence.
    These criteria concern wavelength integration, not spatial or model accuracy.
    """

    def __init__(
        self,
        *,
        wavelength_grids: Sequence[Quantity],
        results: Sequence[MomentRT.results.BroadbandSimulationResult | None],
        detectors: Mapping[str, MomentRT.detector.SurfaceDetector],
        level_errors: Sequence[Mapping[str, str]],
        power_rtol: float,
        photocurrent_rtol: float,
        power_atol: Quantity,
        photocurrent_atol: Quantity,
        required_consecutive: int,
    ) -> None:
        """Collect spectral refinements and evaluate final-sequence stability.

        Parameters
        ----------
        wavelength_grids : sequence of pint.Quantity
            Nested wavelength vectors over one source interval.
        results : sequence of BroadbandSimulationResult or None
            Predictions or failed refinement placeholders.
        detectors : mapping of str to SurfaceDetector
            Fixed apertures and optional responsivity tables.
        level_errors : sequence of mapping of str to str
            Failures and numerical artifacts for each refinement.
        power_rtol : float
            Nonnegative relative power tolerance.
        photocurrent_rtol : float
            Nonnegative relative current tolerance.
        power_atol : pint.Quantity
            Nonnegative absolute optical-power tolerance.
        photocurrent_atol : pint.Quantity
            Nonnegative absolute electrical-current tolerance.
        required_consecutive : int
            Consecutive passing refinements required at the fine end.
        """
        if (
            len(wavelength_grids) < 2
            or len(results) != len(wavelength_grids)
            or len(level_errors) != len(results)
        ):
            raise ValueError(
                "Supply matching grids, results, and errors for at least two levels."
            )
        self.wavelength_grids = tuple(grid.copy() for grid in wavelength_grids)
        self.sample_counts = np.array([len(grid) for grid in wavelength_grids])
        self.results = tuple(results)
        self.detectors = copy.deepcopy(dict(detectors))
        self.level_errors = tuple(dict(errors) for errors in level_errors)
        self.power_rtol = power_rtol
        self.photocurrent_rtol = photocurrent_rtol
        self.power_atol = power_atol.copy()
        self.photocurrent_atol = photocurrent_atol.copy()
        self.required_consecutive = required_consecutive
        self.detector_power: dict[str, Quantity] = {}
        self.detector_photocurrent: dict[str, Quantity] = {}
        self.power_relative_change: dict[str, np.ndarray] = {}
        self.power_relative_to_finest: dict[str, np.ndarray] = {}
        self.photocurrent_relative_change: dict[str, np.ndarray] = {}
        self.photocurrent_relative_to_finest: dict[str, np.ndarray] = {}
        acceptable = np.array([not errors for errors in self.level_errors])
        finest_checks = []
        for name, detector in self.detectors.items():
            power = np.array(
                [
                    np.nan
                    if result is None
                    else result.detectors[name].collected_power.to("W").magnitude
                    for result in results
                ]
            )
            self.detector_power[name] = power * MomentRT.units.ureg.watt
            change, finest, adjacent_ok, finest_ok = self._compare(
                values=power,
                rtol=power_rtol,
                atol=power_atol.to("W").magnitude,
            )
            self.power_relative_change[name] = change
            self.power_relative_to_finest[name] = finest
            acceptable &= adjacent_ok
            finest_checks.append(finest_ok)
            if detector.responsivity is not None:
                current = np.array(
                    [
                        np.nan
                        if result is None or result.detectors[name].photocurrent is None
                        else cast(Quantity, result.detectors[name].photocurrent)
                        .to("A")
                        .magnitude
                        for result in results
                    ]
                )
                self.detector_photocurrent[name] = current * MomentRT.units.ureg.ampere
                change, finest, adjacent_ok, finest_ok = self._compare(
                    values=current,
                    rtol=photocurrent_rtol,
                    atol=photocurrent_atol.to("A").magnitude,
                )
                self.photocurrent_relative_change[name] = change
                self.photocurrent_relative_to_finest[name] = finest
                acceptable &= adjacent_ok
                finest_checks.append(finest_ok)
        self.meets_tolerances = acceptable
        self.stable_refinements = np.zeros(len(results), dtype=int)
        self.recommended_index: int | None = None
        for index, passed in enumerate(acceptable):
            previous = self.stable_refinements[index - 1] if index else 0
            self.stable_refinements[index] = previous + 1 if passed else 0
        final_streak = self.stable_refinements[-1]
        if final_streak >= required_consecutive:
            first_candidate = len(results) - final_streak + required_consecutive - 1
            for index in range(first_candidate, len(results)):
                if all(check[index] for check in finest_checks):
                    self.recommended_index = index
                    break

    @staticmethod
    def _compare(
        *,
        values: np.ndarray,
        rtol: float,
        atol: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compare signed finite outputs against adjacent and finest levels."""
        difference = np.abs(np.diff(values))
        scale = np.maximum(np.abs(values[:-1]), np.abs(values[1:]))
        relative = np.divide(
            difference, scale, out=np.zeros_like(difference), where=scale != 0
        )
        valid = np.isfinite(values) & (values >= 0)
        pair_valid = valid[:-1] & valid[1:]
        relative[~pair_valid] = np.nan
        change = np.concatenate(([np.nan], relative))
        adjacent_ok = np.concatenate(
            ([False], pair_valid & (difference <= atol + rtol * scale))
        )
        finest_difference = np.abs(values - values[-1])
        finest_scale = abs(values[-1])
        finest = np.divide(
            finest_difference,
            finest_scale,
            out=np.full_like(values, np.inf),
            where=finest_scale != 0,
        )
        finest[(finest_scale == 0) & (finest_difference == 0)] = 0
        finest[~valid | ~valid[-1]] = np.nan
        finest_ok = (
            valid & valid[-1] & (finest_difference <= atol + rtol * finest_scale)
        )
        return change, finest, adjacent_ok, finest_ok

    @property
    def converged(self) -> bool:
        """Report whether the final sequence supports a wavelength recommendation.

        Returns
        -------
        bool
            True when all detector quantities meet the final stability criteria.
        """
        return self.recommended_index is not None

    @property
    def recommended_wavelengths(self) -> Quantity | None:
        """Return a copy of the recommended grid when tolerances are met.

        Returns
        -------
        pint.Quantity or None
            Increasing vacuum wavelengths, or None after an unresolved sequence.
        """
        return (
            None
            if self.recommended_index is None
            else self.wavelength_grids[self.recommended_index].copy()
        )

    @property
    def records(self) -> list[dict[str, object]]:
        """Return per-level powers, currents, changes, and explicit failures.

        Returns
        -------
        list of dict
            Records preserving physical units and all failed refinement levels.
        """
        return [
            {
                "wavelengths": self.wavelength_grids[index].copy(),
                "sample_count": int(self.sample_counts[index]),
                "detectors": {
                    name: {
                        "collected_power": power[index],
                        "photocurrent": self.detector_photocurrent[name][index]
                        if name in self.detector_photocurrent
                        else None,
                        "power_relative_change": float(
                            self.power_relative_change[name][index]
                        ),
                        "power_relative_to_finest": float(
                            self.power_relative_to_finest[name][index]
                        ),
                        "photocurrent_relative_change": float(
                            self.photocurrent_relative_change[name][index]
                        )
                        if name in self.detector_photocurrent
                        else None,
                        "photocurrent_relative_to_finest": float(
                            self.photocurrent_relative_to_finest[name][index]
                        )
                        if name in self.detector_photocurrent
                        else None,
                    }
                    for name, power in self.detector_power.items()
                },
                "meets_tolerances": bool(self.meets_tolerances[index]),
                "stable_refinements": int(self.stable_refinements[index]),
                "errors": dict(self.level_errors[index]),
            }
            for index in range(len(self.results))
        ]

    def summary(self) -> str:
        """Describe wavelength refinement outcomes and the final recommendation.

        Returns
        -------
        str
            Labeled sample counts, passing levels, failures, and recommendation.
        """
        lines = []
        for index, count in enumerate(self.sample_counts):
            lines.append(
                f"Level {index}: {count} wavelengths; meets tolerances: {bool(self.meets_tolerances[index])}; errors: {self.level_errors[index]}"
            )
            for name, power in self.detector_power.items():
                line = (
                    f"  {name}: collected power {power[index]:~P}; "
                    f"relative change {self.power_relative_change[name][index]:.3g}"
                )
                if name in self.detector_photocurrent:
                    line += (
                        f"; photocurrent {self.detector_photocurrent[name][index]:~P}; "
                        f"relative change {self.photocurrent_relative_change[name][index]:.3g}"
                    )
                lines.append(line)
        lines.append(
            f"Converged: {self.converged}; recommended index: {self.recommended_index}"
        )
        return "\n".join(lines)

    def plot(
        self,
        *,
        power_unit: str = "uW",
        current_unit: str = "uA",
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Plot integrated measurements against wavelength sample count.

        Parameters
        ----------
        power_unit : str, optional
            Optical-power display unit. Default is 'uW'.
        current_unit : str, optional
            Electrical-current display unit. Default is 'uA'.
        show : bool, optional
            Display immediately when True. Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable power and optional photocurrent curves with failed gaps.
        """
        figure = MomentRT.plotting.make_subplots(
            rows=1, cols=2 if self.detector_photocurrent else 1
        )
        for name, power in self.detector_power.items():
            figure.add_trace(
                MomentRT.plotting.go.Scatter(
                    x=self.sample_counts,
                    y=power.to(power_unit).magnitude,
                    mode="lines+markers",
                    name=f"{name}: optical power",
                ),
                row=1,
                col=1,
            )
        figure.update_yaxes(title_text=f"Collected power ({power_unit})", row=1, col=1)
        for name, current in self.detector_photocurrent.items():
            figure.add_trace(
                MomentRT.plotting.go.Scatter(
                    x=self.sample_counts,
                    y=current.to(current_unit).magnitude,
                    mode="lines+markers",
                    name=f"{name}: photocurrent",
                ),
                row=1,
                col=2,
            )
        if self.detector_photocurrent:
            figure.update_yaxes(
                title_text=f"Photocurrent ({current_unit})", row=1, col=2
            )
        figure.update_xaxes(title_text="Wavelength samples")
        figure.update_layout(title="Wavelength integration convergence")
        if show:
            figure.show()
        return figure


class BroadbandDetectorResult:
    """Collected spectral density and its explicit wavelength integrals.

    Parameters
    ----------
    wavelengths : pint.Quantity
        Strictly increasing positive vacuum wavelengths, shape (n,), n >= 2.
    spectral_power : pint.Quantity
        Collected power density, shape (n,), with W / length units. Negative
        numerical outputs are retained. Integration uses the supplied endpoints.
    responsivity : DetectorResponsivity or None, optional
        Linear electrical response. None leaves photocurrent unavailable.
        Default is None.
    reference_measurements : sequence of DetectorResult or None, optional
        Measurements on the same aperture and mesh, one per wavelength.
        Supplies the sampled and geometric areas. Default is None.

    Attributes
    ----------
    spectral_power : pint.Quantity
        Signed collected power density in its supplied W / length units.
    collected_power : pint.Quantity
        Trapezoidal optical integral in watts.
    photocurrent : pint.Quantity or None
        Trapezoidal electrical integral in amperes when responsivity is supplied.
    responsivity : DetectorResponsivity or None
        Snapshot of the linear response used for electrical conversion.
    sampled_area : pint.Quantity or None
        Shared sampled aperture area. None when no geometry was supplied.
    geometric_area : pint.Quantity or None
        Shared exact aperture area. None when no geometry was supplied.
    """

    def __init__(
        self,
        *,
        wavelengths: Quantity,
        spectral_power: Quantity,
        responsivity: MomentRT.detector.DetectorResponsivity | None = None,
        reference_measurements: Sequence[MomentRT.results.DetectorResult] | None = None,
    ) -> None:
        """Integrate collected optical density and optional electrical response.

        Parameters
        ----------
        wavelengths : pint.Quantity
            Strictly increasing positive vacuum wavelengths, shape (n,), n >= 2.
        spectral_power : pint.Quantity
            Signed collected spectral power density, shape (n,), in W / length.
        responsivity : DetectorResponsivity or None, optional
            Linear spectral electrical response. Default is None.
        reference_measurements : sequence of DetectorResult or None, optional
            Same-aperture measurements on the wavelength grid. Default is None.
        """
        MomentRT.units.validate_units(
            value=wavelengths, unit="m", name="wavelengths", positive=True
        )
        if wavelengths.ndim != 1 or len(wavelengths) < 2:
            raise ValueError("wavelengths must contain at least two samples.")
        grid = wavelengths.to("m").magnitude
        if np.any(np.diff(grid) <= 0):
            raise ValueError("wavelengths must be strictly increasing.")
        MomentRT.units.validate_units(
            value=spectral_power,
            unit="W/m",
            name="spectral_power",
            shape=wavelengths.shape,
        )
        self.wavelengths = wavelengths.copy()
        self.spectral_power = spectral_power.copy()
        self.responsivity = copy.deepcopy(responsivity)
        self.reference_measurements = tuple(copy.deepcopy(reference_measurements or ()))
        self.sampled_area: Quantity | None = None
        self.geometric_area: Quantity | None = None
        self.angular_fraction: float | None = None
        self.detector: MomentRT.detector.SurfaceDetector | None = None
        self.angular_model = "Lambertian escaping radiance"
        if self.reference_measurements:
            if len(self.reference_measurements) != len(wavelengths):
                raise ValueError("Supply one reference measurement per wavelength.")
            first = self.reference_measurements[0]
            for index, measurement in enumerate(self.reference_measurements):
                if not np.isclose(
                    measurement.wavelength.to("m").magnitude,
                    wavelengths[index].to("m").magnitude,
                    rtol=1e-12,
                    atol=0,
                ):
                    raise ValueError(
                        "Reference measurement wavelengths must match the grid."
                    )
                if (
                    not np.isclose(
                        measurement.sampled_area.to("m^2").magnitude,
                        first.sampled_area.to("m^2").magnitude,
                        rtol=1e-12,
                        atol=0,
                    )
                    or not np.isclose(
                        measurement.geometric_area.to("m^2").magnitude,
                        first.geometric_area.to("m^2").magnitude,
                        rtol=1e-12,
                        atol=0,
                    )
                    or measurement.angular_fraction != first.angular_fraction
                ):
                    raise ValueError(
                        "Broadband measurements must use the same aperture quadrature."
                    )
            self.sampled_area = first.sampled_area.copy()
            self.geometric_area = first.geometric_area.copy()
            self.angular_fraction = first.angular_fraction
            self.detector = copy.deepcopy(first.detector)
        self.collected_power = (
            integrate.trapezoid(spectral_power.to("W/m").magnitude, x=grid)
            * MomentRT.units.ureg.watt
        )
        self.spectral_photocurrent: Quantity | None = None
        self.photocurrent: Quantity | None = None
        if self.responsivity is not None:
            self.spectral_photocurrent = spectral_power * self.responsivity.at(
                wavelength=wavelengths
            )
            self.photocurrent = (
                integrate.trapezoid(
                    self.spectral_photocurrent.to("A/m").magnitude,
                    x=grid,
                )
                * MomentRT.units.ureg.ampere
            )


class BroadbandSimulationResult:
    """Broadband detector predictions and the underlying unit-power solves.

    Parameters
    ----------
    spectrum : SourceSpectrum
        Source spectral density and the integration wavelength grid.
    reference_results : SpectralSimulationResult
        Independent monochromatic solves on the source grid, each excited by
        reference_optical_power. Their fields retain those reference units.
    reference_optical_power : pint.Quantity
        Positive scalar optical power used for every reference solve.
    """

    def __init__(
        self,
        *,
        spectrum: MomentRT.source.SourceSpectrum,
        reference_results: MomentRT.results.SpectralSimulationResult,
        reference_optical_power: Quantity,
    ) -> None:
        """Scale unit responses by source density and integrate detector outputs.

        Parameters
        ----------
        spectrum : SourceSpectrum
            Incoming source density with increasing vacuum wavelengths.
        reference_results : SpectralSimulationResult
            Monochromatic reference-power solutions on exactly the same grid.
        reference_optical_power : pint.Quantity
            Positive scalar reference optical power, with power units.
        """
        MomentRT.units.validate_units(
            value=reference_optical_power,
            unit="W",
            name="reference_optical_power",
            shape=(),
            positive=True,
        )
        if (
            spectrum.wavelengths.shape != reference_results.wavelengths.shape
            or not np.allclose(
                spectrum.wavelengths.to("m").magnitude,
                reference_results.wavelengths.to("m").magnitude,
                rtol=8 * np.finfo(float).eps,
                atol=0,
            )
        ):
            raise ValueError("Reference wavelengths must match the source spectrum.")
        self.wavelengths = spectrum.wavelengths.copy()
        self.source_spectrum = copy.deepcopy(spectrum)
        self.reference_results = reference_results
        self.reference_optical_power = reference_optical_power.copy()
        self.detectors = {
            name: BroadbandDetectorResult(
                wavelengths=self.wavelengths,
                spectral_power=(power / reference_optical_power)
                * spectrum.spectral_power,
                responsivity=reference_results.detectors[name][0].detector.responsivity,
                reference_measurements=reference_results.detectors[name],
            )
            for name, power in reference_results.detector_power.items()
        }
        self.diagnostics = {
            **reference_results.diagnostics,
            "integration_rule": "trapezoidal wavelength quadrature",
            "integration_interval": self.wavelengths[[0, -1]].copy(),
            "incoming_optical_power": spectrum.optical_power,
        }

    @property
    def detector_power(self) -> dict[str, Quantity]:
        """Return the wavelength-integrated optical powers for named detectors.

        Returns
        -------
        dict of str to pint.Quantity
            Accepted broadband optical powers in watts.
        """
        return {name: result.collected_power for name, result in self.detectors.items()}


class SpectralSimulationResult:
    """Ordered independent monochromatic solves on one instrument geometry.

    Parameters
    ----------
    wavelengths : pint.Quantity
        Nonempty one-dimensional vector of positive vacuum wavelengths.
    results : sequence of SimulationResult
        One solution for each wavelength, in the same order. All solutions
        must contain the same detector names.
    diagnostics : mapping of str to object, or None, optional
        Numerical provenance, residuals, geometry statistics, or timing quantities.
        Default is None.

    Attributes
    ----------
    wavelengths : pint.Quantity
        Supplied vacuum wavelengths in their original units and order.
    results : tuple of SimulationResult
        One independent solution at each wavelength.
    diagnostics : dict of str to object
        Sweep preparation count and unit-bearing preparation, solve, and total
        timings when created by Simulation.sweep.
    """

    def __init__(
        self,
        *,
        wavelengths: Quantity,
        results: Sequence[MomentRT.results.SimulationResult],
        diagnostics: Mapping[str, object] | None = None,
    ) -> None:
        """Initialize the SpectralSimulationResult configuration.

        Parameters
        ----------
        wavelengths : pint.Quantity
            Nonempty one-dimensional vector of positive vacuum wavelengths.
        results : sequence of SimulationResult
            One solution for each wavelength, in the same order. All solutions
            must contain the same detector names.
        diagnostics : mapping of str to object, or None, optional
            Numerical provenance, residuals, geometry statistics, or timing quantities.
            Default is None.
        """

        MomentRT.units.validate_units(
            value=wavelengths, unit="m", name="wavelengths", positive=True
        )

        if wavelengths.ndim != 1 or not len(wavelengths):
            raise ValueError("wavelengths must be a nonempty vector.")

        results = tuple(results)

        if len(results) != len(wavelengths) or any(
            not isinstance(result, SimulationResult) for result in results
        ):
            raise ValueError("Supply one SimulationResult for each wavelength.")

        names = set(results[0].detectors)

        if any(set(result.detectors) != names for result in results):
            raise ValueError(
                "All spectral results must contain the same detector names."
            )

        self.wavelengths = wavelengths.copy()

        self.results = results

        self.diagnostics = dict(diagnostics or {})

    @property
    def detectors(self) -> dict[str, tuple[MomentRT.results.DetectorResult, ...]]:
        """Return each detector's measurements in supplied wavelength order.

        Returns
        -------
        dict of str to tuple of DetectorResult
            Detector names mapped to monochromatic measurements, including
            collected power and sampled and geometric aperture areas. These
            are independent solves, not a broadband integral.
        """
        return {
            name: tuple(result.detectors[name] for result in self.results)
            for name in self.results[0].detectors
        }

    @property
    def detector_power(self) -> dict[str, Quantity]:
        """Return collected powers keyed by detector name.

        Returns
        -------
        dict of str to pint.Quantity
            Collected optical powers keyed by detector name.
        """

        return {
            name: np.array(
                [
                    result.detector_power[name].to("W").magnitude
                    for result in self.results
                ]
            )
            * MomentRT.units.ureg.watt
            for name in self.results[0].detectors
        }

    def plot(
        self, *, wavelength_unit: str = "nm", power_unit: str = "mW", show: bool = True
    ) -> MomentRT.plotting.Figure:
        """Plot collected power for each detector; no spectral integration is applied.

        Parameters
        ----------
        wavelength_unit : str, optional
            Length unit used for plotted vacuum wavelengths.
            Default is 'nm'.
        power_unit : str, optional
            Power unit used for plotted detector signals.
            Default is 'mW'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable interactive figure, also displayed when show is True.
        """

        powers = self.detector_power

        if not powers:
            raise ValueError("Configure detectors to plot a collected-power spectrum.")

        go = MomentRT.plotting.go

        figure = go.Figure()

        for name, power in powers.items():
            figure.add_trace(
                go.Scatter(
                    x=self.wavelengths.to(wavelength_unit).magnitude,
                    y=power.to(power_unit).magnitude,
                    mode="lines+markers",
                    name=name,
                )
            )

        figure.update_layout(
            template="plotly_white",
            title="Collected power across wavelengths",
            xaxis_title=f"Wavelength [{wavelength_unit}]",
            yaxis_title=f"Collected power [{power_unit}]",
        )

        if show:
            figure.show()

        return figure
