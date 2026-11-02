"""Detector mesh-refinement workflow implementation."""

from __future__ import annotations
import copy
from collections.abc import Iterable
from typing import cast
import numpy as np
from pint import Quantity
import MomentRT


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
            MomentRT.units.validate_units(
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
            MomentRT.units.validate_units(
                value=requested, unit="m", name="wavelength grid", positive=True
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

        Returns
        -------
        DetectorConvergenceResult
            Completed study including unresolved levels and refinement diagnostics.
        """

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

        MomentRT.units.validate_units(
            value=power_atol, unit="W", name="power_atol", shape=(), nonnegative=True
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

        meshes: tuple[MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh, ...]

        if mesh_shapes is not None:
            shapes = tuple(mesh_shapes)

            if len(shapes) < 2:
                raise ValueError("Supply at least two mesh_shapes from coarse to fine.")

            for previous, current in zip(shapes[:-1], shapes[1:]):
                if np.any(np.asarray(current) < previous) or current == previous:
                    raise ValueError(
                        "Each mesh must refine at least one axis and coarsen none."
                    )

            meshes = tuple(
                MomentRT.mesh.CartesianMesh.from_sample(sample=medium, shape=shape)
                for shape in shapes
            )
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

                MomentRT.units.validate_units(
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

            contacts: list[object] = []

            if isinstance(self.simulation.source, MomentRT.source.SurfaceSource):
                contacts.append(source)

            contacts.extend(detectors.values())

            tetrahedral_meshes: list[MomentRT.mesh.TetrahedralMesh] = []

            for resolution in resolutions:
                tetrahedral_mesh = MomentRT.mesh.TetrahedralMesh.from_sample(
                    sample=medium,
                    resolution=resolution,
                )

                if contact_resolution is not None:
                    assert contact_distance is not None
                    tetrahedral_mesh.refine_near(
                        objects=contacts,
                        resolution=contact_resolution,
                        distance=contact_distance,
                    )

                tetrahedral_meshes.append(tetrahedral_mesh)

            meshes = tuple(tetrahedral_meshes)

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

        for mesh in meshes:
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
