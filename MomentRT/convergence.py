"""Detector convergence studies on an explicit sequence of finite meshes."""

from __future__ import annotations

import copy
from collections.abc import Iterable
from typing import cast

import numpy as np
from pint import Quantity

import MomentRT


def _run_study(
    *,
    simulation: MomentRT.simulation.Simulation,
    mesh_shapes: Iterable[tuple[int, int, int]],
    power_rtol: float,
    area_rtol: float,
    power_atol: Quantity | None,
    required_consecutive: int,
) -> MomentRT.results.DetectorConvergenceResult:
    """Run the finite-detector refinement study on independent mesh snapshots.

    Parameters
    ----------
    simulation : Simulation
        Finite Diffusion configuration with at least one named detector.
    mesh_shapes : iterable of tuple of int
        At least two coarse-to-fine node shapes. Each shape must refine at
        least one axis and coarsen none.
    power_rtol : float
        Nonnegative relative power tolerance against adjacent and finest meshes.
    area_rtol : float
        Nonnegative relative tolerance of sampled source and detector areas
        against their geometric areas.
    power_atol : pint.Quantity or None
        Nonnegative absolute power tolerance. None means zero watts.
    required_consecutive : int
        Positive number of consecutive passing refinements required at the
        fine end of the mesh sequence.

    Returns
    -------
    DetectorConvergenceResult
        Completed study including unresolved levels and refinement diagnostics.
    """

    if not isinstance(
        simulation.medium, MomentRT.medium.StructuredMedium
    ) or not isinstance(
        simulation.solver, (MomentRT.diffusion.Diffusion, MomentRT.sp3.SP3)
    ):
        raise TypeError(
            "Detector convergence requires a finite Diffusion or SP3 simulation."
        )

    if not simulation.detectors:
        raise ValueError("Configure at least one detector for a convergence study.")

    if not isinstance(simulation.source, MomentRT.source.UniformSurfaceSource):
        raise TypeError("Detector convergence requires a finite surface source.")

    for name, value in (("power_rtol", power_rtol), ("area_rtol", area_rtol)):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, float, np.integer, np.floating))
            or not np.isfinite(value)
            or value < 0
        ):
            raise ValueError(f"{name} must be a finite nonnegative plain number.")

    if (
        isinstance(required_consecutive, (bool, np.bool_))
        or not isinstance(required_consecutive, (int, np.integer))
        or required_consecutive < 1
    ):
        raise ValueError("required_consecutive must be a positive integer.")

    if power_atol is None:
        power_atol = 0 * MomentRT.units.ureg.watt

    MomentRT.units.validate_units(
        value=power_atol, unit="W", name="power_atol", shape=(), nonnegative=True
    )

    medium = copy.deepcopy(simulation.medium)

    source = copy.deepcopy(simulation.source)

    solver = copy.deepcopy(simulation.solver)

    detectors = copy.deepcopy(simulation.detectors)

    MomentRT.diffusion._finite_materials_at(medium=medium, wavelength=source.wavelength)

    for detector in detectors.values():
        detector._validate_sample(sample_origin=medium.origin, sample_size=medium.size)

    meshes = tuple(
        MomentRT.mesh.CartesianMesh.from_sample(sample=medium, shape=shape)
        for shape in mesh_shapes
    )

    if len(meshes) < 2:
        raise ValueError("Supply at least two mesh_shapes from coarse to fine.")

    for previous, current in zip(meshes[:-1], meshes[1:]):
        if (
            np.any(np.asarray(current.shape) < previous.shape)
            or current.shape == previous.shape
        ):
            raise ValueError(
                "Each mesh must refine at least one axis and coarsen none."
            )

    if isinstance(source, MomentRT.source.SurfaceSource):
        source_geometric_area = np.pi * source.radius**2
    else:
        axis = "xyz".index(source.face[0])

        source_geometric_area = np.prod(
            medium.size[[index for index in range(3) if index != axis]]
        )

    results: list[MomentRT.results.SimulationResult | None] = []

    level_errors = []

    for mesh in meshes:
        candidate = copy.copy(simulation)

        candidate.medium = medium

        candidate.source = source

        candidate.mesh = mesh

        candidate.solver = solver

        candidate.detectors = {}

        errors = {}

        try:
            result = candidate.run()
        except MomentRT.source.SourceResolutionError as error:
            results.append(None)

            level_errors.append({"source": str(error)})

            continue

        for name, detector in detectors.items():
            try:
                measurement = detector.measure(result=result)
            except MomentRT.detector.DetectorResolutionError as error:
                errors[f"detector:{name}"] = str(error)

                continue

            result.detectors[name] = measurement

            if measurement.minimum_escaping_irradiance.magnitude < 0:
                errors[f"detector:{name}"] = (
                    "Negative escaping irradiance in the aperture."
                )

        if cast(Quantity, result.diagnostics["minimum_fluence"]).magnitude < 0:
            errors["fluence"] = "Negative nodal fluence."

        results.append(result)

        level_errors.append(errors)

    return MomentRT.results.DetectorConvergenceResult(
        meshes=meshes,
        results=results,
        detectors=detectors,
        source_geometric_area=source_geometric_area,
        wavelength=source.wavelength,
        level_errors=level_errors,
        power_rtol=power_rtol,
        area_rtol=area_rtol,
        power_atol=power_atol,
        required_consecutive=required_consecutive,
    )
