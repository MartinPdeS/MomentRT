"""Result classes retain their public import paths after organization."""

from __future__ import annotations
import MomentRT
import MomentRT.results
import MomentRT.results.convergence
import MomentRT.results.detector
import MomentRT.results.simulation
import MomentRT.results.spectral


def test_result_import_paths_share_classes() -> None:
    """Public imports resolve to the classes in their respective modules."""
    modules = {
        "SimulationResult": MomentRT.results.simulation,
        "DetectorResult": MomentRT.results.detector,
        "DetectorConvergenceResult": MomentRT.results.convergence,
        "SpectralConvergenceResult": MomentRT.results.convergence,
        "BroadbandDetectorResult": MomentRT.results.spectral,
        "BroadbandSimulationResult": MomentRT.results.spectral,
        "SpectralSimulationResult": MomentRT.results.spectral,
    }

    for name, module in modules.items():
        assert getattr(MomentRT, name) is getattr(MomentRT.results, name)
        assert getattr(MomentRT.results, name) is getattr(module, name)
