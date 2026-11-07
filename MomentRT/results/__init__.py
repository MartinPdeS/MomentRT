"""Result containers for transport simulations and detectors."""

from __future__ import annotations
from MomentRT.results.convergence import (
    DetectorConvergenceResult,
    SpectralConvergenceResult,
)
from MomentRT.results.detector import DetectorResult
from MomentRT.results.simulation import SimulationResult, _energy_accounting
from MomentRT.results.spectral import (
    BroadbandDetectorResult,
    BroadbandSimulationResult,
    SpectralSimulationResult,
)
from MomentRT.results.transient import TransientResult

__all__ = [
    "BroadbandDetectorResult",
    "BroadbandSimulationResult",
    "DetectorConvergenceResult",
    "DetectorResult",
    "SimulationResult",
    "SpectralConvergenceResult",
    "SpectralSimulationResult",
    "TransientResult",
    "_energy_accounting",
]
