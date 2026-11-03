"""Convergence result containers."""

from __future__ import annotations
import copy
from collections.abc import Mapping, Sequence
from typing import cast
import numpy as np
from pint import Quantity
import MomentRT


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

    def assessment(self) -> str:
        """Explain the final convergence decision and any blocking checks.

        Returns
        -------
        str
            Final passing-refinement count, requested tolerances, recommendation,
            and reasons the relevant refinement failed, when applicable.
        """
        final_streak = int(self.stable_refinements[-1])
        lines = [
            f"Converged: {self.converged}",
            f"Consecutive passing refinements at the fine end: "
            f"{final_streak} of {self.required_consecutive} required.",
            f"Tolerances: power {self.power_rtol:.2%} relative + "
            f"{self.power_atol.to('uW').magnitude:.4g} µW absolute; "
            f"footprint area {self.area_rtol:.2%} relative.",
        ]

        if self.converged:
            if self.recommended_resolution is not None:
                recommendation = f"{self.recommended_resolution:~P} target resolution"
            elif self.recommended_shape is not None:
                recommendation = f"Cartesian shape {self.recommended_shape}"
            else:
                recommendation = (
                    f"mesh with {self.node_counts[self.recommended_index]} nodes"
                )
            lines.append(f"Recommended mesh: {recommendation}.")
            return "\n".join(lines)

        lines.append("No mesh is recommended.")
        failed_index = len(self.meshes) - 1 - final_streak

        if failed_index <= 0:
            lines.append(
                "Not enough successive passing refinements are available; "
                "add finer mesh levels and repeat the study."
            )
            if failed_index < 0:
                return "\n".join(lines)

        if failed_index > 0:
            lines.append(f"Most recent failing refinement: level {failed_index}.")

        source_error = self.source_area_relative_error[failed_index]
        if not np.isfinite(source_error):
            lines.append("  Source footprint area is unavailable.")
        elif source_error > self.area_rtol:
            lines.append(
                f"  Source footprint area error {source_error:.2%} exceeds "
                f"{self.area_rtol:.2%}."
            )

        for name in self.detectors:
            area_error = self.detector_area_relative_error[name][failed_index]
            if not np.isfinite(area_error):
                lines.append(f"  Detector {name}: footprint area is unavailable.")
            elif area_error > self.area_rtol:
                lines.append(
                    f"  Detector {name}: footprint area error {area_error:.2%} "
                    f"exceeds {self.area_rtol:.2%}."
                )

            if failed_index > 0 and not self.power_within_tolerance[name][failed_index]:
                relative = self.power_relative_change[name][failed_index]
                if np.isnan(relative):
                    lines.append(f"  Detector {name}: power comparison is unavailable.")
                else:
                    delta = (
                        abs(
                            self.detector_power[name][failed_index]
                            - self.detector_power[name][failed_index - 1]
                        )
                        .to("uW")
                        .magnitude
                    )
                    current = (
                        abs(self.detector_power[name][failed_index]).to("uW").magnitude
                    )
                    limit = (
                        self.power_atol.to("uW").magnitude + self.power_rtol * current
                    )
                    lines.append(
                        f"  Detector {name}: power changed by {delta:.4g} µW; "
                        f"allowed change is {limit:.4g} µW."
                    )

        for label, message in self.level_errors[failed_index].items():
            lines.append(f"  {label}: {message}")

        if final_streak < self.required_consecutive and failed_index > 0:
            lines.append("Refine the mesh further and repeat the study.")
        elif final_streak >= self.required_consecutive:
            lines.append(
                "The final passing sequence did not satisfy the comparison "
                "against the finest mesh. Refine further and repeat the study."
            )

        return "\n".join(lines)

    def summary(self) -> str:
        """Describe mesh levels, detector powers, failures, and the assessment.

        Returns
        -------
        str
            Labeled per-level quantities, footprint errors, relative power
            changes, and the final mesh recommendation, if available.
        """
        lines = [f"Detector mesh convergence at {self.wavelength:~P}"]

        for index, mesh in enumerate(self.meshes):
            description = (
                f"shape {mesh.shape}"
                if isinstance(mesh, MomentRT.mesh.CartesianMesh)
                else f"target resolution {self.mesh_resolutions[index]:~P}"
                if self.mesh_resolutions[index] is not None
                else "imported tetrahedral mesh"
            )

            lines.append(
                f"Level {index}: {description}; {self.node_counts[index]} nodes; "
                f"source area error {self.source_area_relative_error[index]:.2%}; "
                f"meets tolerances: {bool(self.meets_tolerances[index])}"
            )

            for name, powers in self.detector_power.items():
                power = powers[index].to("uW").magnitude
                relative_change = self.power_relative_change[name][index]
                change_text = (
                    f"{relative_change:.2%}" if np.isfinite(relative_change) else "—"
                )
                lines.append(
                    f"  {name}: power {power:.4g} µW; "
                    f"area error {self.detector_area_relative_error[name][index]:.2%}; "
                    f"relative power change {change_text}"
                )

            if self.level_errors[index]:
                for label, message in self.level_errors[index].items():
                    lines.append(f"  {label}: {message}")

        lines.append(self.assessment())

        return "\n".join(lines)

    def plot(
        self, *, power_unit: str = "mW", area_unit: str = "mm^2", show: bool = True
    ) -> MomentRT.plotting.Figure:
        """Compare powers, aperture areas, and source area against mesh node count.

        Shaded bands mark power tolerance around each finest-mesh detector
        value and area tolerance around each exact footprint area. Power bands
        include the absolute tolerance and are omitted for unresolved finest
        values. Passing these bands alone does not establish convergence: the
        study also checks adjacent refinements and diagnostic failures.

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
                f"Power: {self.power_rtol:.1%} + absolute tolerance",
                f"Detector area: ±{self.area_rtol:.1%}",
                f"Source area: ±{self.area_rtol:.1%}",
            ),
        )

        for number, name in enumerate(self.detectors):
            color = MomentRT.plotting.qualitative.Plotly[
                number % len(MomentRT.plotting.qualitative.Plotly)
            ]

            finest_power = self.detector_power[name][-1].to(power_unit).magnitude

            if np.isfinite(finest_power):
                power_tolerance = (
                    self.power_rtol * abs(finest_power)
                    + self.power_atol.to(power_unit).magnitude
                )

                figure.add_hrect(
                    y0=finest_power - power_tolerance,
                    y1=finest_power + power_tolerance,
                    fillcolor=color,
                    opacity=0.12,
                    line_width=0,
                    layer="below",
                    exclude_empty_subplots=False,
                    row=1,
                    col=1,
                )

            geometric_area = self.detector_geometric_area[name].to(area_unit).magnitude

            figure.add_hrect(
                y0=geometric_area * (1 - self.area_rtol),
                y1=geometric_area * (1 + self.area_rtol),
                fillcolor=color,
                opacity=0.12,
                line_width=0,
                layer="below",
                exclude_empty_subplots=False,
                row=1,
                col=2,
            )

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

        source_area = self.source_geometric_area.to(area_unit).magnitude

        figure.add_hrect(
            y0=source_area * (1 - self.area_rtol),
            y1=source_area * (1 + self.area_rtol),
            fillcolor="crimson",
            opacity=0.12,
            line_width=0,
            layer="below",
            exclude_empty_subplots=False,
            row=1,
            col=3,
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

        # Plotly does not include layout shapes in autoranging. Invisible
        # points keep the full tolerance bands visible on each y axis.
        for band in figure.layout.shapes:
            column = {"y": 1, "y2": 2, "y3": 3}[band.yref]
            figure.add_trace(
                go.Scatter(
                    x=[self.node_counts[0], self.node_counts[-1]],
                    y=[band.y0, band.y1],
                    mode="markers",
                    marker={"opacity": 0},
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=1,
                col=column,
            )

        figure.update_xaxes(
            title_text="Mesh nodes",
            type="log",
            tickmode="array",
            tickvals=self.node_counts,
            ticktext=[f"{int(count):,}" for count in self.node_counts],
            tickangle=-35,
        )

        figure.update_yaxes(
            title_text=f"Power [{power_unit}]; band around finest mesh", row=1, col=1
        )

        figure.update_yaxes(title_text=f"Area [{area_unit}]", row=1, col=2)

        figure.update_yaxes(title_text=f"Area [{area_unit}]", row=1, col=3)

        if not self.converged:
            recommendation = "no mesh meets the criteria"
        elif self.recommended_resolution is not None:
            recommendation = f"recommended {self.recommended_resolution:~P}"
        else:
            recommendation = f"recommended {self.recommended_shape}"

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
