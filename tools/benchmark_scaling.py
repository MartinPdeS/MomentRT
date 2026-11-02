"""Measure tetrahedral mesh generation, finite solves, and process peak memory."""

from __future__ import annotations
import argparse
import json
import platform
import subprocess
import sys
import time
from typing import Any
import numpy as np
import scipy
import MomentRT

resource: Any
try:
    import resource as resource_module
except ImportError:
    resource = None
else:
    resource = resource_module


def benchmark_case(*, solver_name: str, resolution_mm: float) -> dict[str, object]:
    """Run one isolated tetrahedral finite-element scaling case."""
    unit = MomentRT.ureg
    material = MomentRT.Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=0.03 / unit.mm,
        reduced_scattering=1.0 / unit.mm,
        refractive_index=1.4,
        anisotropy=0.8,
    )
    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * unit.mm,
        size=np.array([6, 6, 4]) * unit.mm,
    )
    sample.add_background(material=material)

    mesh_started = time.perf_counter()
    mesh = MomentRT.TetrahedralMesh.from_sample(
        sample=sample,
        resolution=resolution_mm * unit.mm,
    )
    mesh_seconds = time.perf_counter() - mesh_started

    source = MomentRT.UniformSurfaceSource(
        face="z_min",
        wavelength=850 * unit.nm,
        optical_power=1 * unit.mW,
    )
    boundaries = MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0)
    )
    solver_types: dict[str, Any] = {
        "p1": MomentRT.Diffusion,
        "sp3": MomentRT.SP3,
        "p3": MomentRT.P3,
    }
    solver_type = solver_types[solver_name]
    simulation = MomentRT.Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=solver_type(boundaries=boundaries),
    )

    solve_started = time.perf_counter()
    simulation.run()
    solve_seconds = time.perf_counter() - solve_started

    peak_rss_bytes = None
    if resource is not None:
        peak_rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != "darwin":
            peak_rss_bytes *= 1024

    return {
        "solver": solver_name,
        "target_resolution_mm": resolution_mm,
        "nodes": len(mesh.points),
        "tetrahedra": len(mesh.tetrahedra),
        "mesh_generation_seconds": mesh_seconds,
        "solve_seconds": solve_seconds,
        "peak_process_rss_bytes": peak_rss_bytes,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }


def _worker(*, solver: str, resolution_mm: float) -> None:
    """Print one case as JSON for the parent process."""
    print(json.dumps(benchmark_case(solver_name=solver, resolution_mm=resolution_mm)))


def main() -> None:
    """Run isolated solver and tetrahedral target-size combinations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solvers", nargs="+", choices=("p1", "sp3", "p3"), default=("p1", "sp3", "p3")
    )
    parser.add_argument(
        "--resolutions-mm", nargs="+", type=float, default=(2.0, 1.5, 1.0)
    )
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--solver", choices=("p1", "sp3", "p3"), help=argparse.SUPPRESS)
    parser.add_argument("--resolution-mm", type=float, help=argparse.SUPPRESS)
    arguments = parser.parse_args()

    if arguments.worker:
        if arguments.solver is None or arguments.resolution_mm is None:
            parser.error("Worker mode requires --solver and --resolution-mm.")
        _worker(solver=arguments.solver, resolution_mm=arguments.resolution_mm)
        return

    if any(not np.isfinite(value) or value <= 0 for value in arguments.resolutions_mm):
        parser.error("Resolutions must be finite positive millimeter values.")
    if any(
        previous <= current
        for previous, current in zip(
            arguments.resolutions_mm[:-1], arguments.resolutions_mm[1:]
        )
    ):
        parser.error("Resolutions must decrease from coarse to fine.")

    cases = []
    for solver in arguments.solvers:
        for resolution_mm in arguments.resolutions_mm:
            completed = subprocess.run(
                [
                    sys.executable,
                    __file__,
                    "--worker",
                    "--solver",
                    solver,
                    "--resolution-mm",
                    str(resolution_mm),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            cases.append(json.loads(completed.stdout))

    print(json.dumps({"sample_size_mm": [6, 6, 4], "cases": cases}, indent=2))


if __name__ == "__main__":
    main()
