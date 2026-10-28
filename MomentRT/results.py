"""Results retain physical fields and their provenance."""


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


class SpectralSimulationResult:
    """Ordered independent monochromatic solves on one instrument geometry."""

    def __init__(self, *, wavelengths, results):
        from .units import validate_units

        validate_units(value=wavelengths, unit="m", name="wavelengths", positive=True)

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

    @property
    def detector_power(self):
        import numpy as np

        from .units import ureg

        return {
            name: np.array(
                [
                    result.detector_power[name].to("W").magnitude
                    for result in self.results
                ]
            )
            * ureg.watt
            for name in self.results[0].detectors
        }

    def plot(self, *, wavelength_unit="nm", power_unit="mW", show=True):
        """Plot collected power for each detector; no spectral integration is applied."""
        from .plotting import require_plotly

        powers = self.detector_power

        if not powers:
            raise ValueError("Configure detectors to plot a collected-power spectrum.")

        go = require_plotly()

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
