"""Compare finite-moment sensor predictions with independent photon transport."""

from __future__ import annotations
import argparse
import json
import platform
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any
import matplotlib.pyplot as plt
import numpy as np
import scipy
from sensor_reference import LayerOptics, simulate_layered_slab_reference
import MomentRT

DETECTOR_OFFSETS_MM = (0.75, 1.25, 1.5, 2.0, 2.5)
MESH_RESOLUTIONS_MM = (1.0, 0.8, 0.6)
CONTACT_DISTANCE_MM = 0.35
SAMPLE_WIDTH_MM = 6.0


@dataclass(frozen=True)
class SensorCase:
    """One fully specified optical and contact-sensor benchmark case."""

    layers: tuple[LayerOptics, ...]
    source_radius_mm: float = 0.2
    detector_radius_mm: float = 0.2
    detector_offsets_mm: tuple[float, ...] = DETECTOR_OFFSETS_MM

    @property
    def layer_thicknesses_mm(self) -> tuple[float, ...]:
        """Return the configured thickness of each optical layer."""
        return tuple(layer.thickness_mm for layer in self.layers)


BASE_LAYERS = (
    LayerOptics(0.25, 0.18, 1.50, 0.85),
    LayerOptics(0.65, 0.05, 1.10, 0.80),
    LayerOptics(1.10, 0.02, 0.65, 0.75),
)
HIGH_ABSORPTION_LAYERS = (
    LayerOptics(0.25, 0.35, 1.20, 0.85),
    LayerOptics(0.65, 0.14, 0.90, 0.80),
    LayerOptics(1.10, 0.08, 0.60, 0.75),
)
CASES = {
    "moderate_absorption": SensorCase(BASE_LAYERS),
    "higher_absorption": SensorCase(HIGH_ABSORPTION_LAYERS),
    "thin_top_layer": SensorCase(
        tuple(
            replace(layer, thickness_mm=thickness)
            for layer, thickness in zip(BASE_LAYERS, (0.10, 0.65, 1.25))
        )
    ),
    "thick_top_layer": SensorCase(
        tuple(
            replace(layer, thickness_mm=thickness)
            for layer, thickness in zip(BASE_LAYERS, (0.45, 0.65, 0.90))
        )
    ),
    "small_source": SensorCase(BASE_LAYERS, source_radius_mm=0.10),
    "large_source": SensorCase(BASE_LAYERS, source_radius_mm=0.35),
    "small_detector": SensorCase(BASE_LAYERS, detector_radius_mm=0.10),
    "large_detector": SensorCase(BASE_LAYERS, detector_radius_mm=0.30),
}


def build_case(
    *, case: SensorCase
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
        for layer in case.layers
    ]
    sample = MomentRT.StructuredMedium(
        origin=(-SAMPLE_WIDTH_MM / 2, -SAMPLE_WIDTH_MM / 2, 0) * unit.mm,
        size=(
            SAMPLE_WIDTH_MM,
            SAMPLE_WIDTH_MM,
            sum(case.layer_thicknesses_mm),
        )
        * unit.mm,
    )
    sample.add_background(material=materials[-1])
    for index, (thickness, material) in enumerate(
        zip(case.layer_thicknesses_mm[:-1], materials[:-1])
    ):
        sample.add_layer(
            name=f"Layer {index + 1}",
            thickness=thickness * unit.mm,
            material=material,
        )

    source = MomentRT.SurfaceSource(
        face="z_min",
        position=np.zeros(3) * unit.mm,
        radius=case.source_radius_mm * unit.mm,
        wavelength=850 * unit.nm,
        optical_power=1 * unit.mW,
        name="Diffuse source",
    )
    detectors = {
        f"separation_{offset:g}_mm": MomentRT.SurfaceDetector.full_hemisphere_on_face(
            sample=sample,
            face="z_min",
            coordinates=np.array([offset, 0]) * unit.mm,
            radius=case.detector_radius_mm * unit.mm,
            name=f"Separation {offset:g} mm",
        )
        for offset in case.detector_offsets_mm
    }
    return sample, source, detectors


def run_case(
    *, name: str, case: SensorCase, photon_count: int, seed: int
) -> dict[str, Any]:
    """Run Monte Carlo and P1/SP3/P3 mesh sequences for one sensor case."""
    unit = MomentRT.ureg
    sample_height_mm = sum(case.layer_thicknesses_mm)
    contact_resolutions_mm = tuple(
        min(resolution, min(case.source_radius_mm, case.detector_radius_mm))
        for resolution in (0.2, 0.15, 0.1)
    )
    monte_carlo = simulate_layered_slab_reference(
        layers=case.layers,
        sample_width_mm=SAMPLE_WIDTH_MM,
        sample_height_mm=sample_height_mm,
        source_radius_mm=case.source_radius_mm,
        detector_offsets_mm=case.detector_offsets_mm,
        detector_radius_mm=case.detector_radius_mm,
        photon_count=photon_count,
        seed=seed,
        optical_power_w=1e-3,
    )
    sample, source, detectors = build_case(case=case)
    boundaries = MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0)
    )

    meshes = []
    mesh_generation_seconds = []
    for resolution_mm, contact_resolution_mm in zip(
        MESH_RESOLUTIONS_MM, contact_resolutions_mm
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
                    "contact_resolution_mm": contact_resolutions_mm[len(levels)],
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
        "configuration": {
            "sample_width_mm": SAMPLE_WIDTH_MM,
            "layer_thicknesses_mm": case.layer_thicknesses_mm,
            "layers": [asdict(layer) for layer in case.layers],
            "source_radius_mm": case.source_radius_mm,
            "detector_radius_mm": case.detector_radius_mm,
            "detector_offsets_mm": case.detector_offsets_mm,
            "source_power_w": 1e-3,
            "wavelength_nm": 850,
            "boundary_model": "index-matched transmitting faces; no Fresnel reflection",
            "detector_angular_model": "full outward hemisphere",
        },
        "monte_carlo": monte_carlo,
        "finite_element": finite_element,
    }


def plot_results(*, report: dict[str, Any], destination: Path) -> None:
    """Plot detector power errors and Monte Carlo confidence intervals."""
    all_case_reports = report["cases"]
    base_case_reports = [
        case
        for case in all_case_reports
        if case["name"] in ("moderate_absorption", "higher_absorption")
    ]
    case_reports = (
        base_case_reports
        if len(all_case_reports) > 2 and len(base_case_reports) == 2
        else all_case_reports
    )
    figure, axes = plt.subplots(
        len(case_reports),
        len(DETECTOR_OFFSETS_MM),
        figsize=(18, max(7, 3.5 * len(case_reports))),
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
            if row == 0:
                axis.set_title(f"Separation {offset:g} mm")
            axis.grid(True, alpha=0.25)
            if row == len(case_reports) - 1:
                axis.set_xlabel("Tetrahedra")
            if column == 0:
                axis.set_ylabel(
                    f"{case['name'].replace('_', ' ')}\nDetector power error vs MC (%)"
                )
    axes[0, 0].legend()
    figure.suptitle("Layered contact-sensor transport validation")
    figure.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def plot_case_matrix(*, report: dict[str, Any], destination: Path) -> None:
    """Plot finest-mesh relative errors for every case and separation."""
    case_reports = report["cases"]
    solver_names = ("P1", "SP3", "P3")
    detector_names = [f"separation_{offset:g}_mm" for offset in DETECTOR_OFFSETS_MM]
    errors = np.array(
        [
            [
                [
                    case["finite_element"][solver_name][-1]["detectors"][detector_name][
                        "relative_error_percent"
                    ]
                    for detector_name in detector_names
                ]
                for case in case_reports
            ]
            for solver_name in solver_names
        ],
        dtype=float,
    )
    color_limit = max(100, min(200, np.ceil(np.nanmax(np.abs(errors)) / 25) * 25))
    figure, axes = plt.subplots(
        1,
        len(solver_names),
        figsize=(13, max(6, 0.75 * len(case_reports))),
        sharey=True,
        constrained_layout=True,
    )
    image = None
    for index, (axis, solver_name) in enumerate(zip(axes, solver_names)):
        image = axis.imshow(
            np.clip(errors[index], -color_limit, color_limit),
            cmap="coolwarm",
            vmin=-color_limit,
            vmax=color_limit,
            aspect="auto",
        )
        axis.set_title(solver_name)
        axis.set_xticks(
            np.arange(len(DETECTOR_OFFSETS_MM)),
            labels=[f"{offset:g}" for offset in DETECTOR_OFFSETS_MM],
        )
        axis.set_xlabel("Source-detector separation (mm)")
        axis.set_yticks(
            np.arange(len(case_reports)),
            labels=[case["name"].replace("_", " ") for case in case_reports],
        )
        for row in range(len(case_reports)):
            for column in range(len(DETECTOR_OFFSETS_MM)):
                value = errors[index, row, column]
                axis.text(
                    column,
                    row,
                    f"{value:+.0f}%",
                    ha="center",
                    va="center",
                    color="white" if abs(value) > color_limit * 0.55 else "black",
                    fontsize=8,
                )
    assert image is not None
    axes[0].set_ylabel("Optical and contact geometry case")
    figure.colorbar(
        image,
        ax=axes.tolist(),
        label=f"Relative detector-power error (%, clipped at ±{color_limit}%)",
        shrink=0.8,
    )
    figure.suptitle("Finest-mesh detector-power error across validation cases")
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> None:
    """Run the reproducible sensor validation matrix and save JSON plus plot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--photons", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=tuple(CASES),
        default=tuple(CASES),
        help="Subset of named cases to run (default: all cases).",
    )
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
            case=case,
            photon_count=arguments.photons,
            seed=arguments.seed + index,
        )
        for index, name in enumerate(arguments.cases)
        for case in (CASES[name],)
    ]
    report = {
        "monte_carlo_method": "independent layered-slab photon transport",
        "photon_count_per_case": arguments.photons,
        "mesh_resolutions_mm": MESH_RESOLUTIONS_MM,
        "base_contact_resolution_schedule_mm": (0.2, 0.15, 0.1),
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
    case_matrix_path = arguments.output_dir / "sensor_validation_cases.png"
    plot_case_matrix(report=report, destination=case_matrix_path)
    print(f"Validation data: {json_path}")
    print(f"Error plot: {plot_path}")
    print(f"Case matrix plot: {case_matrix_path}")
    print(f"Elapsed: {report['elapsed_seconds']:.1f} s")


if __name__ == "__main__":
    main()
