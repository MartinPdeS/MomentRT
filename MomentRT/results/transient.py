"""Time-resolved P1 fields and energy diagnostics."""

from __future__ import annotations
import numpy as np
import TypedUnit
from pint import Quantity
import MomentRT


class TransientResult:
    """Time-resolved fields and power accounting for a finite P1 simulation.

    Parameters
    ----------
    times : pint.Quantity
        Simulation times with time units, shape (N,).
    frames : tuple of SimulationResult
        Fields and per-face boundary snapshots at each time, length N.
    stored_energy : pint.Quantity
        Integrated radiation energy, shape (N,), with energy units.
    absorbed_power : pint.Quantity
        Integrated absorption, shape (N,), with power units.
    incoming_power : pint.Quantity
        Incoming coupled power, shape (N,), with power units.
    escaping_power : pint.Quantity
        Signed total escaping power, shape (N,), with power units.
    interval_energy_residual : pint.Quantity
        Midpoint energy-balance residual for each interval, shape (N-1,),
        with power units. This is a numerical conservation diagnostic, not a
        measure of transport accuracy.
    """

    def __init__(
        self,
        *,
        times: Quantity,
        frames: tuple[MomentRT.results.SimulationResult, ...],
        stored_energy: Quantity,
        absorbed_power: Quantity,
        incoming_power: Quantity,
        escaping_power: Quantity,
        interval_energy_residual: Quantity,
    ) -> None:
        """Store the time-resolved field snapshots and energy diagnostics.

        Parameters
        ----------
        times : pint.Quantity
            Simulation times, shape (N,), with time units.
        frames : tuple of SimulationResult
            Per-time fluence, element current, and boundary fields.
        stored_energy : pint.Quantity
            Integrated radiation energy, shape (N,), with energy units.
        absorbed_power : pint.Quantity
            Integrated absorption, shape (N,), with power units.
        incoming_power : pint.Quantity
            Incoming coupled power, shape (N,), with power units.
        escaping_power : pint.Quantity
            Signed escaping power, shape (N,), with power units.
        interval_energy_residual : pint.Quantity
            Midpoint energy-balance residual, shape (N-1,), with power units.
        """
        TypedUnit.validate_quantity(
            value=times,
            unit="s",
            name="times",
            registry=MomentRT.units.ureg,
            nonnegative=True,
        )
        if (
            times.ndim != 1
            or len(times) < 2
            or np.any(np.diff(times.to("s").magnitude) <= 0)
        ):
            raise ValueError(
                "times must contain at least two strictly increasing samples."
            )
        if len(frames) != len(times) or any(
            not isinstance(frame, MomentRT.results.SimulationResult) for frame in frames
        ):
            raise ValueError("frames must contain one SimulationResult per time.")
        for name, value, unit, shape in (
            ("stored_energy", stored_energy, "J", times.shape),
            ("absorbed_power", absorbed_power, "W", times.shape),
            ("incoming_power", incoming_power, "W", times.shape),
            ("escaping_power", escaping_power, "W", times.shape),
            (
                "interval_energy_residual",
                interval_energy_residual,
                "W",
                (len(times) - 1,),
            ),
        ):
            TypedUnit.validate_quantity(
                value=value,
                unit=unit,
                name=name,
                registry=MomentRT.units.ureg,
                shape=shape,
            )
        self.times = times.copy()
        self.frames = frames
        self.stored_energy = stored_energy.copy()
        self.absorbed_power = absorbed_power.copy()
        self.incoming_power = incoming_power.copy()
        self.escaping_power = escaping_power.copy()
        self.interval_energy_residual = interval_energy_residual.copy()

    @property
    def detector_power(self) -> dict[str, Quantity]:
        """Return each detector's time-resolved collected power.

        Returns
        -------
        dict of str to pint.Quantity
            Detector names mapped to signed power arrays of shape (N,).
            Aperture integration uses the existing Lambertian angular model.
        """
        if not self.frames:
            return {}
        return {
            name: np.array(
                [
                    frame.detectors[name].collected_power.to("W").magnitude
                    for frame in self.frames
                ]
            )
            * MomentRT.units.ureg.watt
            for name in self.frames[0].detectors
        }

    @property
    def net_outward_power(self) -> Quantity:
        """Return escaping minus incoming power over the complete boundary.

        Returns
        -------
        pint.Quantity
            Signed net outward power, shape (N,), with power units.
        """
        return self.escaping_power - self.incoming_power
