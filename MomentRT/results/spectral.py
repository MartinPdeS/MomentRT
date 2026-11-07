"""Spectral result containers."""

from __future__ import annotations
import copy
from collections.abc import Mapping, Sequence
import numpy as np
import TypedUnit
from pint import Quantity
from scipy import integrate
import MomentRT


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
        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=wavelengths,
            unit="m",
            name="wavelengths",
            positive=True,
        )
        if wavelengths.ndim != 1 or len(wavelengths) < 2:
            raise ValueError("wavelengths must contain at least two samples.")
        grid = wavelengths.to("m").magnitude
        if np.any(np.diff(grid) <= 0):
            raise ValueError("wavelengths must be strictly increasing.")
        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
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
        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
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
            name: MomentRT.results.BroadbandDetectorResult(
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

        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=wavelengths,
            unit="m",
            name="wavelengths",
            positive=True,
        )

        if wavelengths.ndim != 1 or not len(wavelengths):
            raise ValueError("wavelengths must be a nonempty vector.")

        results = tuple(results)

        if len(results) != len(wavelengths) or any(
            not isinstance(result, MomentRT.results.SimulationResult)
            for result in results
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
