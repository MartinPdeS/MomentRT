"""Measure finite-sweep time and traced allocation peaks for reproducible cases."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import tracemalloc

import numpy as np
import scipy

import MomentRT

try:
    import resource
except ImportError:
    resource = None


def main() -> None:
    """Run one configured spectral sweep and print machine-readable diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--shape", type=int, nargs=3, default=(9, 9, 9))

    parser.add_argument("--wavelengths", type=int, default=3)

    parser.add_argument("--solver", choices=("p1", "sp3", "p3"), default="p1")

    args = parser.parse_args()

    if min(args.shape) < 3 or args.wavelengths < 1:
        parser.error("Each mesh dimension must be at least 3 and wavelengths positive.")

    unit = MomentRT.ureg

    material = MomentRT.SpectralMaterial(
        wavelengths=np.array([700, 900]) * unit.nm,
        absorption=np.array([0.03, 0.02]) / unit.mm,
        reduced_scattering=np.array([1.2, 0.9]) / unit.mm,
        refractive_index=1.4,
        anisotropy=np.array([0.75, 0.8]),
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )

    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * unit.mm,
        size=np.full(3, 4) * unit.mm,
    )

    sample.add_background(material=material)

    mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=tuple(args.shape))

    source = MomentRT.UniformSurfaceSource(
        face="z_min", wavelength=800 * unit.nm, optical_power=1 * unit.mW
    )

    boundary = MomentRT.RobinBoundary(effective_reflectance=0)

    boundaries = {face: boundary for face in MomentRT.boundary.FACES}

    solver = {"p1": MomentRT.Diffusion, "sp3": MomentRT.SP3, "p3": MomentRT.P3}[
        args.solver
    ](boundaries=boundaries)

    simulation = MomentRT.Simulation(
        medium=sample, mesh=mesh, source=source, solver=solver
    )

    tracemalloc.start()

    result = simulation.sweep(
        wavelengths=np.linspace(700, 900, args.wavelengths) * unit.nm
    )

    _, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()

    diagnostics = result.diagnostics

    peak_rss = None

    if resource is not None:
        peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

        if sys.platform != "darwin":
            peak_rss *= 1024

    print(
        json.dumps(
            {
                "solver": args.solver,
                "shape": args.shape,
                "wavelength_count": args.wavelengths,
                "preparation_seconds": diagnostics["preparation_time"]
                .to("s")
                .magnitude,
                "solve_seconds": diagnostics["solve_times"].to("s").magnitude.tolist(),
                "total_seconds": diagnostics["total_time"].to("s").magnitude,
                "geometry_preparations": diagnostics["geometry_preparations"],
                "peak_traced_bytes": peak,
                "peak_process_rss_bytes": peak_rss,
                "python": platform.python_version(),
                "platform": platform.platform(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
