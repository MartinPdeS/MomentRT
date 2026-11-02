"""Compare finite-moment sensor predictions with independent photon transport."""

from __future__ import annotations
import argparse
import json
import platform
import time
from pathlib import Path
from typing import Any
import matplotlib.pyplot as plt
import numpy as np
import scipy
from sensor_reference import LayerOptics, simulate_layered_slab_reference
import MomentRT

LAYERS_MM = (0.25, 0.65, 1.10)
DETECTOR_OFFSETS_MM = (0.75, 1.5, 2.5)
MESH_RESOLUTIONS_MM = (1.0, 0.8, 0.6)
CONTACT_RESOLUTIONS_MM = (0.2, 0.15, 0.1)
CONTACT_DISTANCE_MM = 0.35
SAMPLE_WIDTH_MM = 6.0

CASES = {
    "moderate_absorption": (
        LayerOptics(0.25, 0.18, 1.50, 0.85),
        LayerOptics(0.65, 0.05, 1.10, 0.80),
        LayerOptics(1.10, 0.02, 0.65, 0.75),
    ),
    "higher_absorption": (
        LayerOptics(0.25, 0.35, 1.20, 0.85),
        LayerOptics(0.65, 0.14, 0.90, 0.80),
        LayerOptics(1.10, 0.08, 0.60, 0.75),
    ),
}


def build_case(
    *, layers: tuple[LayerOptics, ...]
) -> tuple[
    MomentRT.medium.StructuredMedium,
    MomentRT.source.SurfaceSource,
    dict[str, MomentRT.detector.SurfaceDetector],
]:
    """Construct the matching MomentRT sample, source, and detector geometry."""
    unit = MomentRT.ureg
    materials = [
        MomentRT.Material(
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
            absorption=layer.absorption_per_mm / unit.mm,
            reduced_scattering=layer.reduced_scattering_per_mm / unit.mm,
            refractive_index=1.4,
            anisotropy=layer.anisotropy,
        )
        for layer in layers
    ]
    sample = MomentRT.StructuredMedium(
        origin=(-SAMPLE_WIDTH_MM / 2, -SAMPLE_WIDTH_MM / 2, 0) * unit.mm,
        size=(SAMPLE_WIDTH_MM, SAMPLE_WIDTH_MM, sum(LAYERS_MM)) * unit.mm,
    )
    sample.add_background(material=materials[-1])
    for index, (thickness, material) in enumerate(zip(LAYERS_MM[:-1], materials[:-1])):
        sample.add_layer(
            name=f"Layer {index + 1}",
            thickness=thickness * unit.mm,
            material=material,
        )

    source = MomentRT.SurfaceSource(
        face="z_min",
        position=np.zeros(3) * unit.mm,
        radius=0.2 * unit.mm,
        wavelength=850 * unit.nm,
        optical_power=1 * unit.mW,
        name="Diffuse source",
    )
    detectors = {
        f"separation_{offset:g}_mm": MomentRT.SurfaceDetector.full_hemisphere_on_face(
            sample=sample,
            face="z_min",
            coordinates=np.array([offset, 0]) * unit.mm,
            radius=0.2 * unit.mm,
            name=f"Separation {offset:g} mm",
        )
        for offset in DETECTOR_OFFSETS_MM
    }
    return sample, source, detectors


def run_case(
    *, name: str, layers: tuple[LayerOptics, ...], photon_count: int, seed: int
) -> dict[str, Any]:
    """Run independent Monte Carlo and P1/SP3/P3 mesh sequences for one tissue."""
    unit = MomentRT.ureg
    monte_carlo = simulate_layered_slab_reference(
        layers=layers,
        sample_width_mm=SAMPLE_WIDTH_MM,
        sample_height_mm=sum(LAYERS_MM),
        source_radius_mm=0.2,
        detector_offsets_mm=DETECTOR_OFFSETS_MM,
        detector_radius_mm=0.2,
        photon_count=photon_count,
        seed=seed,
        optical_power_w=1e-3,
    )
    sample, source, detectors = build_case(layers=layers)
    boundaries = MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0)
    )

    meshes = []
    mesh_generation_seconds = []
    for resolution_mm, contact_resolution_mm in zip(
        MESH_RESOLUTIONS_MM, CONTACT_RESOLUTIONS_MM
    ):
        started = time.perf_counter()
        mesh = MomentRT.TetrahedralMesh.from_sample(
            sample=sample,
            resolution=resolution_mm * unit.mm,
        )
        mesh.refine_near(
            objects=[source, *detectors.values()],
            resolution=contact_resolution_mm * unit.mm,
            distance=CONTACT_DISTANCE_MM * unit.mm,
        )
        meshes.append(mesh)
        mesh_generation_seconds.append(time.perf_counter() - started)

    solver_types: dict[str, Any] = {
        "P1": MomentRT.Diffusion,
        "SP3": MomentRT.SP3,
        "P3": MomentRT.P3,
    }
    finite_element: dict[str, list[dict[str, Any]]] = {}
    for solver_name, solver_type in solver_types.items():
        levels: list[dict[str, Any]] = []
        for mesh, mesh_time, target_mm in zip(
            meshes, mesh_generation_seconds, MESH_RESOLUTIONS_MM
        ):
            simulation = MomentRT.Simulation(
                medium=sample,
                source=source,
                mesh=mesh,
                solver=solver_type(boundaries=boundaries),
                detectors=detectors,
            )
            report = simulation.validate()
            if not report.valid:
                raise ValueError(f"{name} {solver_name} preflight: {report.summary()}")
            started = time.perf_counter()
            result = simulation.run()
            solve_seconds = time.perf_counter() - started
            levels.append(
                {
                    "target_resolution_mm": target_mm,
                    "contact_resolution_mm": CONTACT_RESOLUTIONS_MM[len(levels)],
                    "nodes": len(mesh.points),
                    "tetrahedra": len(mesh.tetrahedra),
                    "mesh_generation_seconds": mesh_time,
                    "solve_seconds": solve_seconds,
                    "source_sampled_area_mm2": result.source_sampled_area.to(
                        "mm^2"
                    ).magnitude,
                    "source_geometric_area_mm2": result.source_geometric_area.to(
                        "mm^2"
                    ).magnitude,
                    "detectors": {
                        detector_name: {
                            "power_w": result.detector_power[detector_name]
                            .to("W")
                            .magnitude,
                            "sampled_area_mm2": result.detectors[detector_name]
                            .sampled_area.to("mm^2")
                            .magnitude,
                            "geometric_area_mm2": result.detectors[detector_name]
                            .geometric_area.to("mm^2")
                            .magnitude,
                            "sampled_area_fraction": (
                                result.detectors[detector_name]
                                .sampled_area.to("mm^2")
                                .magnitude
                                / result.detectors[detector_name]
                                .geometric_area.to("mm^2")
                                .magnitude
                            ),
                        }
                        for detector_name in detectors
                    },
                }
            )
        finite_element[solver_name] = levels

    mc_power = monte_carlo["detector_power_w"]
    mc_error = monte_carlo["standard_error_w"]
    for solver_levels in finite_element.values():
        for level in solver_levels:
            for detector_index, detector_name in enumerate(detectors):
                value = level["detectors"][detector_name]["power_w"]
                reference = mc_power[detector_index]
                level["detectors"][detector_name]["relative_error_percent"] = (
                    100 * (value - reference) / reference if reference else None
                )
                level["detectors"][detector_name]["mc_95_percent"] = (
                    196 * mc_error[detector_index] / reference if reference else None
                )

    return {
        "name": name,
        "monte_carlo": monte_carlo,
        "finite_element": finite_element,
    }


def plot_results(*, report: dict[str, Any], destination: Path) -> None:
    """Plot detector power errors and Monte Carlo confidence intervals."""
    case_reports = report["cases"]
    figure, axes = plt.subplots(
        len(case_reports),
        len(DETECTOR_OFFSETS_MM),
        figsize=(13, 7),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    colors = {"P1": "#0072B2", "SP3": "#D55E00", "P3": "#009E73"}
    for row, case in enumerate(case_reports):
        for column, offset in enumerate(DETECTOR_OFFSETS_MM):
            axis = axes[row, column]
            detector_name = f"separation_{offset:g}_mm"
            for solver_name, levels in case["finite_element"].items():
                x = np.array([level["tetrahedra"] for level in levels])
                errors = np.array(
                    [
                        level["detectors"][detector_name]["relative_error_percent"]
                        for level in levels
                    ]
                )
                axis.plot(
                    x,
                    errors,
                    marker="o",
                    color=colors[solver_name],
                    label=solver_name,
                )
            mc_band = case["monte_carlo"]["relative_standard_error"][column] * 196
            axis.axhspan(-mc_band, mc_band, color="0.8", zorder=0)
            axis.axhline(0, color="black", linewidth=0.8)
            axis.set_xscale("log")
            axis.set_title(
                f"{case['name'].replace('_', ' ')}, separation {offset:g} mm"
            )
            axis.grid(True, alpha=0.25)
            if row == len(case_reports) - 1:
                axis.set_xlabel("Tetrahedra")
            if column == 0:
                axis.set_ylabel("Detector power error vs MC (%)")
    axes[0, 0].legend()
    figure.suptitle("Layered contact-sensor transport validation")
    figure.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> None:
    """Run the reproducible sensor validation matrix and save JSON plus plot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--photons", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/source/_static/images/sensor_validation"),
    )
    arguments = parser.parse_args()
    if arguments.photons < 2:
        parser.error("--photons must be at least 2.")

    started = time.perf_counter()
    cases = [
        run_case(
            name=name,
            layers=layers,
            photon_count=arguments.photons,
            seed=arguments.seed + index,
        )
        for index, (name, layers) in enumerate(CASES.items())
    ]
    report = {
        "monte_carlo_method": "independent layered-slab photon transport",
        "photon_count_per_case": arguments.photons,
        "mesh_resolutions_mm": MESH_RESOLUTIONS_MM,
        "contact_resolutions_mm": CONTACT_RESOLUTIONS_MM,
        "contact_distance_mm": CONTACT_DISTANCE_MM,
        "sample_width_mm": SAMPLE_WIDTH_MM,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "momentrt": MomentRT.__version__,
        "elapsed_seconds": time.perf_counter() - started,
        "cases": cases,
    }
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = arguments.output_dir / "sensor_validation.json"
    json_path.write_text(json.dumps(report, indent=2))
    plot_path = arguments.output_dir / "sensor_validation.png"
    plot_results(report=report, destination=plot_path)
    print(f"Validation data: {json_path}")
    print(f"Error plot: {plot_path}")
    print(f"Elapsed: {report['elapsed_seconds']:.1f} s")


if __name__ == "__main__":
    main()
