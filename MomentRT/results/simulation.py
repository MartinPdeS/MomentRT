"""Simulation result containers."""

from __future__ import annotations
from collections.abc import Mapping
from typing import cast
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
