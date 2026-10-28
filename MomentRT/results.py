"""Results retain physical fields and their provenance."""

import numpy as np

import MomentRT


class SimulationResult:
    def __init__(
        self,
        *,
        fluence,
        flux=None,
        energy_balance=None,
        outgoing_power_by_face=None,
        diagnostics=None,
        boundary_fields=None,
        detectors=None,
    ):
        self.fluence = fluence

        self.flux = flux

        self.energy_balance = dict(energy_balance or {})

        self.outgoing_power_by_face = dict(outgoing_power_by_face or {})

        self.diagnostics = dict(diagnostics or {})

        self.boundary_fields = dict(boundary_fields or {})

        self.detectors = dict(detectors or {})

    @property
    def detector_power(self):
        return {
            name: measurement.collected_power
            for name, measurement in self.detectors.items()
        }


class DetectorResult:
    """Finite-aperture diffuse collection and quadrature provenance."""

    def __init__(
        self,
        *,
        collected_power,
        sampled_area,
        geometric_area,
        angular_fraction,
        wavelength,
        minimum_escaping_irradiance,
        detector,
    ):
        self.collected_power = collected_power

        self.sampled_area = sampled_area

        self.geometric_area = geometric_area

        self.angular_fraction = angular_fraction

        self.wavelength = wavelength.copy()

        self.minimum_escaping_irradiance = minimum_escaping_irradiance

        self.detector = detector

        self.angular_model = "Lambertian escaping radiance"


class DetectorConvergenceResult:
    """Detector signals and footprint sampling on successively refined meshes."""

    def __init__(
        self,
        *,
        meshes,
        results,
        detectors,
        source_geometric_area,
        wavelength,
        level_errors,
        power_rtol,
        area_rtol,
        power_atol,
        required_consecutive,
    ):
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
                    else result.diagnostics["source_quadrature_area"]
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
    def converged(self):
        return self.recommended_index is not None

    @property
    def recommended_mesh(self):
        return (
            None
            if self.recommended_index is None
            else self.meshes[self.recommended_index]
        )

    @property
    def recommended_shape(self):
        return None if self.recommended_mesh is None else self.recommended_mesh.shape

    @property
    def records(self):
        """Return per-level records with quantities intact and explicit failure flags."""
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

    def plot(self, *, power_unit="mW", area_unit="mm^2", show=True):
        """Compare powers, aperture areas, and source area against mesh node count."""

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
    """Ordered independent monochromatic solves on one instrument geometry."""

    def __init__(self, *, wavelengths, results, diagnostics=None):
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
    def detector_power(self):
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

    def plot(self, *, wavelength_unit="nm", power_unit="mW", show=True):
        """Plot collected power for each detector; no spectral integration is applied."""

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
