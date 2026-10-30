"""Non-solving validation report for simulation configurations."""

from __future__ import annotations
from collections.abc import Mapping
from pint import Quantity


class ValidationReport:
    """Configuration errors, numerical cautions, and sampled footprint areas.

    Parameters
    ----------
    errors : tuple of str
        Configuration errors that must be corrected before solving.
    warnings : tuple of str
        Numerical cautions that do not prevent a solve.
    source_sampled_area : pint.Quantity or None
        Boundary quadrature area sampled by the source, when available.
    detector_sampled_areas : mapping of str to pint.Quantity
        Boundary quadrature area sampled by each valid detector.

    Notes
    -----
    Passing validation does not establish mesh convergence or physical model
    accuracy. This report is a snapshot and never caches a prepared operator.
    """

    def __init__(
        self,
        *,
        errors: tuple[str, ...],
        warnings: tuple[str, ...],
        source_sampled_area: Quantity | None,
        detector_sampled_areas: Mapping[str, Quantity],
    ) -> None:
        """Store the independently collected validation findings.

        Parameters
        ----------
        errors : tuple of str
            Configuration errors that prevent solving.
        warnings : tuple of str
            Numerical cautions that do not prevent solving.
        source_sampled_area : pint.Quantity or None
            Boundary quadrature area sampled by the source, when available.
        detector_sampled_areas : mapping of str to pint.Quantity
            Sampled boundary quadrature areas for valid detectors.
        """
        self.errors = errors

        self.warnings = warnings

        self.source_sampled_area = source_sampled_area

        self.detector_sampled_areas = dict(detector_sampled_areas)

    @property
    def valid(self) -> bool:
        """Report whether any configuration errors were found.

        Returns
        -------
        bool
            True when no errors were found; warnings may still be present.
        """
        return not self.errors

    def summary(self) -> str:
        """Format the validation findings for display.

        Returns
        -------
        str
            Human-readable findings and sampled areas with physical units.
        """
        lines = ["Configuration valid" if self.valid else "Configuration has errors"]

        lines.extend(f"Error: {error}" for error in self.errors)

        lines.extend(f"Warning: {warning}" for warning in self.warnings)

        if self.source_sampled_area is not None:
            lines.append(f"Source sampled area: {self.source_sampled_area:~P}")

        lines.extend(
            f"Detector {name} sampled area: {area:~P}"
            for name, area in self.detector_sampled_areas.items()
        )

        return "\n".join(lines)
