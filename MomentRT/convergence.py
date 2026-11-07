"""Detector mesh-refinement workflow implementation."""

from __future__ import annotations
import copy
import logging
import time
from collections.abc import Iterable
from typing import cast
import numpy as np
import TypedUnit
from pint import Quantity
import MomentRT

_LOGGER = logging.getLogger(__name__)


class MeshConvergence:
    """Configure a detector-power mesh-refinement study at one wavelength.

    Parameters
    ----------
    simulation : Simulation
        Finite simulation with a surface source and named detectors. Its source
        wavelength is used for every mesh.
    mesh_shapes : iterable of tuple of int or None, optional
        Cartesian node shapes, ordered from coarse to fine. Supply exactly one
        of this and ``mesh_resolutions``. Default is None.
    mesh_resolutions : iterable of pint.Quantity or None, optional
        Strictly decreasing tetrahedral target sizes with length units.
        Default is None.
    contact_resolution : pint.Quantity or None, optional
        Tetrahedral target size near source and detector disks. Supply with
        ``contact_distance``. Default is None.
    contact_distance : pint.Quantity or None, optional
        Distance from each contact defining local refinement. Default is None.
    power_rtol : float, optional
        Relative detector-power tolerance. Default is 0.01.
    area_rtol : float, optional
        Relative footprint-area tolerance. Default is 0.02.
    power_atol : pint.Quantity or None, optional
        Absolute detector-power tolerance, with power units. None uses zero.
        Default is None.
    required_consecutive : int, optional
        Number of passing refinements required at the fine end. Default is 2.
    progress : bool, optional
        Print a live table as mesh levels finish. Default is True.

    Notes
    -----
    A passing mesh study does not establish the accuracy of the transport model.
    The simulation is sampled at run time and is not mutated by the study.
    Set the ``MomentRT.convergence`` logger to INFO for progress and values.
    """

    def __init__(
        self,
        *,
        simulation: MomentRT.simulation.Simulation,
        mesh_shapes: Iterable[tuple[int, int, int]] | None = None,
        mesh_resolutions: Iterable[Quantity] | None = None,
        contact_resolution: Quantity | None = None,
        contact_distance: Quantity | None = None,
        power_rtol: float = 0.01,
        area_rtol: float = 0.02,
        power_atol: Quantity | None = None,
        required_consecutive: int = 2,
        progress: bool = True,
    ) -> None:
        """Store the simulation and the explicit coarse-to-fine study settings.

        Parameters
        ----------
        simulation : Simulation
            Finite detector simulation sampled when the study runs.
        mesh_shapes : iterable of tuple of int or None, optional
            Cartesian node shapes. Supply this or ``mesh_resolutions``.
            Default is None.
        mesh_resolutions : iterable of pint.Quantity or None, optional
            Tetrahedral target sizes with length units. Default is None.
        contact_resolution : pint.Quantity or None, optional
            Tetrahedral target size near contacts. Default is None.
        contact_distance : pint.Quantity or None, optional
            Radius of contact refinement zones. Default is None.
        power_rtol : float, optional
            Relative detector-power tolerance. Default is 0.01.
        area_rtol : float, optional
            Relative footprint-area tolerance. Default is 0.02.
        power_atol : pint.Quantity or None, optional
            Absolute power tolerance, with power units. Default is None.
        required_consecutive : int, optional
            Number of passing refinements required at the fine end. Default is 2.
        progress : bool, optional
            Print a live table as mesh levels finish. Default is True.
        """
        self.simulation = simulation
        self.mesh_shapes = None if mesh_shapes is None else tuple(mesh_shapes)
        self.mesh_resolutions = (
            None if mesh_resolutions is None else tuple(mesh_resolutions)
        )
        self.contact_resolution = contact_resolution
        self.contact_distance = contact_distance
        self.power_rtol = power_rtol
        self.area_rtol = area_rtol
        self.power_atol = power_atol
        self.required_consecutive = required_consecutive
        if not isinstance(progress, bool):
            raise TypeError("progress must be a bool.")
        self.progress = progress

    def run(self) -> MomentRT.results.DetectorConvergenceResult:
        """Run the mesh study with the current simulation configuration.

        Returns
        -------
        DetectorConvergenceResult
            Per-mesh detector powers, footprint diagnostics, and convergence
            assessment. Unresolved levels and negative outputs are retained.
        """
        return MomentRT.convergence._DetectorConvergenceStudy(
            simulation=self.simulation,
        ).run(
            mesh_shapes=self.mesh_shapes,
            mesh_resolutions=self.mesh_resolutions,
            contact_resolution=self.contact_resolution,
            contact_distance=self.contact_distance,
            power_rtol=self.power_rtol,
            area_rtol=self.area_rtol,
            power_atol=self.power_atol,
            required_consecutive=self.required_consecutive,
            progress=self.progress,
        )


class _SpectralConvergenceStudy:
    """Refine wavelength quadrature for a fixed broadband simulation."""

    def __init__(self, *, simulation: MomentRT.simulation.Simulation) -> None:
        self.simulation = simulation

    def run(
        self,
        *,
        wavelength_grids: Iterable[Quantity],
        include_table_knots: bool,
        power_rtol: float,
        photocurrent_rtol: float,
        power_atol: Quantity | None,
        photocurrent_atol: Quantity | None,
        required_consecutive: int,
    ) -> MomentRT.results.SpectralConvergenceResult:
        """Run nested quadrature grids without changing geometry or source data."""
        if self.simulation.source.spectrum is None:
            raise ValueError("Spectral convergence requires a source spectrum.")
        if not self.simulation.detectors:
            raise ValueError(
                "Configure at least one detector for spectral convergence."
            )
        for name, value in (
            ("power_rtol", power_rtol),
            ("photocurrent_rtol", photocurrent_rtol),
        ):
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
        if not isinstance(include_table_knots, bool):
            raise TypeError("include_table_knots must be a bool.")
        power_atol = 0 * MomentRT.units.ureg.watt if power_atol is None else power_atol
        photocurrent_atol = (
            0 * MomentRT.units.ureg.ampere
            if photocurrent_atol is None
            else photocurrent_atol
        )
        for name, absolute_tolerance, unit in (
            ("power_atol", power_atol, "W"),
            ("photocurrent_atol", photocurrent_atol, "A"),
        ):
            TypedUnit.validate_quantity(
                registry=MomentRT.units.ureg,
                value=absolute_tolerance,
                unit=unit,
                name=name,
                shape=(),
                nonnegative=True,
            )

        configuration = copy.copy(self.simulation)
        configuration.medium = copy.deepcopy(self.simulation.medium)
        configuration.source = copy.deepcopy(self.simulation.source)
        configuration.detectors = copy.deepcopy(self.simulation.detectors)
        configuration.solver = copy.deepcopy(self.simulation.solver)
        configuration.mesh = self.simulation.mesh.copy()
        spectrum = configuration.source.spectrum
        assert spectrum is not None
        knots = [spectrum.wavelengths]
        if include_table_knots:
            medium = configuration.medium
            materials = (
                [medium.material]
                if isinstance(medium, MomentRT.medium.HomogeneousMedium)
                else [
                    medium.background,
                    *(shape.material for shape in medium.structures),
                ]
            )
            knots.extend(
                material.wavelengths
                for material in materials
                if isinstance(material, MomentRT.material.SpectralMaterial)
            )
            knots.extend(
                detector.responsivity.wavelengths
                for detector in configuration.detectors.values()
                if detector.responsivity is not None
            )
        grid_unit = spectrum.wavelengths.units
        lower, upper = spectrum.wavelengths.magnitude[[0, -1]]
        tolerance = 8 * np.finfo(float).eps * max(abs(lower), abs(upper))
        fixed_knots = np.concatenate([grid.to(grid_unit).magnitude for grid in knots])
        fixed_knots = fixed_knots[(fixed_knots >= lower) & (fixed_knots <= upper)]
        grids: list[Quantity] = []
        for requested in wavelength_grids:
            TypedUnit.validate_quantity(
                registry=MomentRT.units.ureg,
                value=requested,
                unit="m",
                name="wavelength grid",
                positive=True,
            )
            if requested.ndim != 1 or len(requested) < 2:
                raise ValueError(
                    "Each wavelength grid must contain at least two samples."
                )
            values = requested.to(grid_unit).magnitude
            if np.any(np.diff(values) <= 0):
                raise ValueError("Each wavelength grid must be strictly increasing.")
            if (
                abs(values[0] - lower) > tolerance
                or abs(values[-1] - upper) > tolerance
            ):
                raise ValueError("Every grid must span the original source interval.")
            values = np.sort(np.concatenate((values, fixed_knots)))
            values = values[np.concatenate(([True], np.diff(values) > tolerance))]
            values[[0, -1]] = [lower, upper]
            if grids:
                previous = grids[-1].magnitude
                if len(values) <= len(previous) or not np.all(
                    np.any(
                        np.abs(previous[:, None] - values[None, :]) <= tolerance, axis=1
                    )
                ):
                    raise ValueError(
                        "Grids must be nested and add wavelengths after table knots are included."
                    )
            grids.append(values * grid_unit)
        if len(grids) < 2:
            raise ValueError("Supply at least two coarse-to-fine wavelength grids.")

        results: list[MomentRT.results.BroadbandSimulationResult | None] = []
        level_errors = []
        for grid in grids:
            candidate = copy.copy(configuration)
            candidate.source = copy.copy(configuration.source)
            candidate.source.spectrum = spectrum.resample(wavelengths=grid)
            errors = {}
            try:
                result = candidate.run()
                assert isinstance(result, MomentRT.results.BroadbandSimulationResult)
                for name, measurement in result.detectors.items():
                    if measurement.collected_power.magnitude < 0 or (
                        measurement.photocurrent is not None
                        and measurement.photocurrent.magnitude < 0
                    ):
                        errors[f"detector:{name}"] = (
                            "Negative integrated detector output."
                        )
                for reference in result.reference_results.results:
                    minimum = reference.minimum_fluence
                    if minimum is not None and minimum.magnitude < 0:
                        errors["fluence"] = (
                            "Negative nodal fluence in a reference solve."
                        )
                    for name, reference_measurement in reference.detectors.items():
                        if (
                            reference_measurement.minimum_escaping_irradiance.magnitude
                            < 0
                        ):
                            errors[f"detector:{name}"] = (
                                "Negative escaping irradiance in the aperture."
                            )
                results.append(result)
            except (ValueError, TypeError, RuntimeError) as error:
                results.append(None)
                errors["simulation"] = str(error)
            level_errors.append(errors)
        return MomentRT.results.SpectralConvergenceResult(
            wavelength_grids=grids,
            results=results,
            detectors=configuration.detectors,
            level_errors=level_errors,
            power_rtol=power_rtol,
            photocurrent_rtol=photocurrent_rtol,
            power_atol=power_atol,
            photocurrent_atol=photocurrent_atol,
            required_consecutive=required_consecutive,
        )


class _DetectorConvergenceStudy:
    """Run a detector convergence study for a simulation snapshot.

    Parameters
    ----------
    simulation : Simulation
        Configuration to copy for independent refinement levels.
    """

    def __init__(self, *, simulation: MomentRT.simulation.Simulation) -> None:
        self.simulation = simulation

    def run(
        self,
        *,
        mesh_shapes: Iterable[tuple[int, int, int]] | None,
        mesh_resolutions: Iterable[Quantity] | None,
        contact_resolution: Quantity | None,
        contact_distance: Quantity | None,
        power_rtol: float,
        area_rtol: float,
        power_atol: Quantity | None,
        required_consecutive: int,
        progress: bool = False,
    ) -> MomentRT.results.DetectorConvergenceResult:
        """Run the finite-detector refinement study on independent mesh snapshots.

        Parameters
        ----------
        mesh_shapes : iterable of tuple of int or None
            At least two coarse-to-fine Cartesian node shapes. Each shape must
            refine at least one axis and coarsen none. Supply this or
            ``mesh_resolutions``, but not both.
        mesh_resolutions : iterable of pint.Quantity or None
            At least two strictly decreasing global tetrahedral target sizes.
            Supply this or ``mesh_shapes``, but not both.
        contact_resolution : pint.Quantity or None
            Optional local target element size near sources and detectors.
            Required together with ``contact_distance`` and only with
            ``mesh_resolutions``.
        contact_distance : pint.Quantity or None
            Positive distance defining each local contact refinement zone.
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

        progress : bool, optional
            Print a live table as levels finish. Default is False.

        Returns
        -------
        DetectorConvergenceResult
            Completed study including unresolved levels and refinement diagnostics.
        """

        study_started = time.perf_counter()

        if not isinstance(
            self.simulation.medium, MomentRT.medium.StructuredMedium
        ) or not isinstance(
            self.simulation.solver,
            (
                MomentRT.transport.diffusion.Diffusion,
                MomentRT.transport.sp3.SP3,
                MomentRT.transport.p3.P3,
            ),
        ):
            raise TypeError(
                "Detector convergence requires a finite Diffusion, SP3, or P3 simulation."
            )

        if not self.simulation.detectors:
            raise ValueError("Configure at least one detector for a convergence study.")

        if not isinstance(
            self.simulation.source,
            (
                MomentRT.source.UniformSurfaceSource,
                MomentRT.source.SurfaceSourceGroup,
            ),
        ):
            raise TypeError("Detector convergence requires a finite surface source.")

        if (mesh_shapes is None) == (mesh_resolutions is None):
            raise ValueError("Supply exactly one of mesh_shapes or mesh_resolutions.")

        if (contact_resolution is None) != (contact_distance is None):
            raise ValueError(
                "contact_resolution and contact_distance must be supplied together."
            )

        if mesh_shapes is not None and contact_resolution is not None:
            raise ValueError("Contact refinement applies only to tetrahedral meshes.")

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

        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=power_atol,
            unit="W",
            name="power_atol",
            shape=(),
            nonnegative=True,
        )

        if self.simulation.source.spectrum is not None:
            raise ValueError(
                "Detector convergence studies require a monochromatic source; "
                "check representative wavelengths before broadband integration."
            )

        medium = copy.deepcopy(self.simulation.medium)

        source = copy.deepcopy(self.simulation.source)

        solver = copy.deepcopy(self.simulation.solver)

        detectors = copy.deepcopy(self.simulation.detectors)

        MomentRT.transport.diffusion._finite_materials_at(
            medium=medium, wavelength=source.wavelength
        )

        for detector in detectors.values():
            detector._validate_sample(
                sample_origin=medium.origin, sample_size=medium.size
            )

        meshes: list[MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh] = []
        contacts: list[object] = []

        _LOGGER.info(
            "Starting detector mesh convergence at %s: power rtol %.2g, "
            "power atol %s, area rtol %.2g, required consecutive refinements %d",
            f"{source.wavelength:~P}",
            power_rtol,
            f"{power_atol:~P}",
            area_rtol,
            required_consecutive,
        )

        if mesh_shapes is not None:
            shapes = tuple(mesh_shapes)

            if len(shapes) < 2:
                raise ValueError("Supply at least two mesh_shapes from coarse to fine.")

            for previous, current in zip(shapes[:-1], shapes[1:]):
                if np.any(np.asarray(current) < previous) or current == previous:
                    raise ValueError(
                        "Each mesh must refine at least one axis and coarsen none."
                    )

            requests: tuple[tuple[int, int, int] | Quantity, ...] = shapes
        else:
            resolutions = tuple(mesh_resolutions or ())

            if len(resolutions) < 2:
                raise ValueError(
                    "Supply at least two mesh_resolutions from coarse to fine."
                )

            for index, resolution in enumerate(resolutions):
                MomentRT.mesh.TetrahedralMesh._validate_resolution(
                    resolution=resolution
                )

                if index and not resolutions[index - 1] > resolution:
                    raise ValueError(
                        "mesh_resolutions must decrease strictly from coarse to fine."
                    )

            if contact_resolution is not None:
                assert contact_distance is not None
                MomentRT.mesh.TetrahedralMesh._validate_resolution(
                    resolution=contact_resolution
                )

                TypedUnit.validate_quantity(
                    registry=MomentRT.units.ureg,
                    value=contact_distance,
                    unit="m",
                    name="contact_distance",
                    shape=(),
                    positive=True,
                )

                if any(contact_resolution > resolution for resolution in resolutions):
                    raise ValueError(
                        "contact_resolution must not exceed any global mesh resolution."
                    )

            if isinstance(self.simulation.source, MomentRT.source.SurfaceSource):
                contacts.append(source)

            contacts.extend(detectors.values())

            requests = resolutions

        if isinstance(source, MomentRT.source.SurfaceSourceGroup):
            source_geometric_area = sum(
                np.pi * member.radius**2 for member in source.sources
            )
        elif isinstance(source, MomentRT.source.SurfaceSource):
            source_geometric_area = np.pi * source.radius**2
        else:
            axis = "xyz".index(source.face[0])

            source_geometric_area = np.prod(
                medium.size[[index for index in range(3) if index != axis]]
            )

        results: list[MomentRT.results.SimulationResult | None] = []

        level_errors = []

        if progress:
            columns = [
                f"{'Mesh':>13}",
                f"{'Nodes':>9}",
                f"{'Time (s)':>9}",
                f"{'Source err':>10}",
            ]
            for name in detectors:
                columns.extend(
                    [f"{name + ' µW':>13}", f"{'Δ power':>9}", f"{'Area err':>9}"]
                )
            header = "  ".join(columns)
            print(f"Detector mesh convergence at {source.wavelength:~P}", flush=True)
            print(header, flush=True)
            print("-" * len(header), flush=True)

        for index, request in enumerate(requests):
            level_started = time.perf_counter()
            started = level_started
            if isinstance(request, tuple):
                _LOGGER.info(
                    "Preparing mesh %d/%d: shape %s", index + 1, len(requests), request
                )
                mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh = (
                    MomentRT.mesh.CartesianMesh.from_sample(
                        sample=medium, shape=request
                    )
                )
            else:
                _LOGGER.info(
                    "Preparing mesh %d/%d: target resolution %s",
                    index + 1,
                    len(requests),
                    f"{request:~P}",
                )
                mesh = MomentRT.mesh.TetrahedralMesh.from_sample(
                    sample=medium, resolution=request
                )
                if contact_resolution is not None:
                    assert contact_distance is not None
                    mesh.refine_near(
                        objects=contacts,
                        resolution=contact_resolution,
                        distance=contact_distance,
                    )
                _LOGGER.info(
                    "Prepared mesh %d/%d: %d nodes, %d tetrahedra in %.2f s",
                    index + 1,
                    len(requests),
                    len(mesh.points),
                    len(mesh.tetrahedra),
                    time.perf_counter() - started,
                )

            meshes.append(mesh)
            started = time.perf_counter()
            node_count = (
                int(np.prod(mesh.shape))
                if isinstance(mesh, MomentRT.mesh.CartesianMesh)
                else len(mesh.points)
            )
            _LOGGER.info(
                "Solving mesh %d/%d: %d nodes", index + 1, len(requests), node_count
            )
            candidate = copy.copy(self.simulation)

            candidate.medium = medium

            candidate.source = source

            candidate.mesh = mesh

            candidate.solver = solver

            candidate.detectors = {}

            errors = {}

            try:
                result = cast(MomentRT.results.SimulationResult, candidate.run())
            except MomentRT.source.SourceResolutionError as error:
                results.append(None)

                level_errors.append({"source": str(error)})
                if progress:
                    label = (
                        str(request) if isinstance(request, tuple) else f"{request:~P}"
                    )
                    print(
                        f"{label:>13}  {node_count:>9,}  {time.perf_counter() - level_started:>9.2f}  "
                        f"{'—':>10}  unresolved source: {error}",
                        flush=True,
                    )
                _LOGGER.info(
                    "Mesh %d/%d unresolved source: %s", index + 1, len(requests), error
                )

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

            _LOGGER.info(
                "Mesh %d/%d complete in %.2f s: source sampled area %s, geometric area %s",
                index + 1,
                len(requests),
                time.perf_counter() - started,
                f"{result.source_sampled_area:~P}",
                f"{source_geometric_area:~P}",
            )

            for name, measurement in result.detectors.items():
                _LOGGER.info(
                    "Mesh %d/%d detector %s: collected power %s, sampled area %s, geometric area %s",
                    index + 1,
                    len(requests),
                    name,
                    f"{measurement.collected_power:~P}",
                    f"{measurement.sampled_area:~P}",
                    f"{measurement.geometric_area:~P}",
                )

            if errors:
                _LOGGER.info(
                    "Mesh %d/%d diagnostics: %s", index + 1, len(requests), errors
                )

            if progress:
                label = str(request) if isinstance(request, tuple) else f"{request:~P}"
                source_error = abs(
                    (cast(Quantity, result.source_sampled_area) / source_geometric_area)
                    .to("")
                    .magnitude
                    - 1
                )
                cells = [
                    f"{label:>13}",
                    f"{node_count:>9,}",
                    f"{time.perf_counter() - level_started:>9.2f}",
                    f"{source_error:>9.2%}",
                ]
                previous_result = results[-2] if len(results) > 1 else None
                for name in detectors:
                    current_measurement = result.detectors.get(name)
                    if current_measurement is None:
                        cells.extend([f"{'—':>13}", f"{'—':>9}", f"{'—':>9}"])
                        continue
                    power = current_measurement.collected_power.to("uW").magnitude
                    area_error = abs(
                        (
                            current_measurement.sampled_area
                            / current_measurement.geometric_area
                        )
                        .to("")
                        .magnitude
                        - 1
                    )
                    previous_measurement = (
                        None
                        if previous_result is None
                        else previous_result.detectors.get(name)
                    )
                    if previous_measurement is None:
                        change = "—"
                    else:
                        previous_power = previous_measurement.collected_power.to(
                            "uW"
                        ).magnitude
                        change = (
                            f"{abs(power - previous_power) / abs(power):.2%}"
                            if power != 0
                            else "0.00%"
                            if previous_power == 0
                            else "∞"
                        )
                    cells.extend(
                        [f"{power:>13.4g}", f"{change:>9}", f"{area_error:>8.2%}"]
                    )
                print("  ".join(cells), flush=True)
                for label, message in errors.items():
                    print(f"  {label}: {message}", flush=True)

        report = MomentRT.results.DetectorConvergenceResult(
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

        _LOGGER.info("Mesh convergence complete. Converged: %s", report.converged)

        if progress:
            print(report.assessment(), flush=True)
            print(
                f"Total time: {time.perf_counter() - study_started:.2f} s", flush=True
            )

        return report
