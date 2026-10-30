"""Results retain physical fields and their provenance."""

from __future__ import annotations
from collections.abc import Mapping, Sequence
from typing import cast
import numpy as np
from pint import Quantity
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
        """

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
    meshes : sequence of CartesianMesh
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
        meshes: Sequence[MomentRT.mesh.CartesianMesh],
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
        meshes : sequence of CartesianMesh
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
            [6 * np.prod(np.asarray(mesh.shape) - 1) for mesh in meshes], dtype=int
        )

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
    def recommended_mesh(self) -> MomentRT.mesh.CartesianMesh | None:
        """Return the recommended mesh when the convergence criteria are met.

        Returns
        -------
        CartesianMesh or None
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

        return None if self.recommended_mesh is None else self.recommended_mesh.shape

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

        go = MomentRT.plotting.require_plotly()

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

        go = MomentRT.plotting.require_plotly()

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
