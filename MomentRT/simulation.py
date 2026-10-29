"""Composition of medium, illumination, and transport."""

from __future__ import annotations

import copy
import time
from collections.abc import Iterable, Mapping
from typing import Any, Protocol

import numpy as np
from pint import Quantity

import MomentRT


class TransportSolver(Protocol):
    """Structural interface for solvers accepted by :class:`Simulation`.

    Implementations choose their supported medium and source types and validate
    those configurations at runtime. No inheritance from this protocol is needed.
    """

    def solve(
        self,
        *,
        medium: Any,
        source: Any,
        mesh: MomentRT.mesh.CartesianMesh,
    ) -> MomentRT.results.SimulationResult | MomentRT.fields.ScalarField:
        """Solve the transport problem on the supplied Cartesian mesh.

        Parameters
        ----------
        medium : HomogeneousMedium or StructuredMedium
            Physical configuration supported by the implementation.
        source : PointSource or UniformSurfaceSource
            Illumination supported by the implementation.
        mesh : CartesianMesh
            Sampling nodes or finite-element grid.

        Returns
        -------
        SimulationResult or ScalarField
            Transport solution, or a fluence field wrapped by Simulation.run.
        """
        ...


class Simulation:
    """Compose physical geometry, illumination, a transport solver, and detectors.

    Parameters
    ----------
    medium : HomogeneousMedium or StructuredMedium
        Physical material configuration, compatible with the chosen solver.
    source : PointSource or UniformSurfaceSource
        Illumination configuration; SurfaceSource is a finite disk subclass.
    mesh : CartesianMesh
        Cartesian nodes used for sampling or finite-element discretization.
    solver : TransportSolver
        Transport model implementing the keyword-only solve protocol.
    detectors : mapping of str to SurfaceDetector, or None, optional
        Nonempty names mapped to passive contact detectors. None disables
        detector measurements.
        Default is None.
    """

    def __init__(
        self,
        *,
        medium: MomentRT.medium.HomogeneousMedium | MomentRT.medium.StructuredMedium,
        source: MomentRT.source.PointSource | MomentRT.source.UniformSurfaceSource,
        mesh: MomentRT.mesh.CartesianMesh,
        solver: MomentRT.simulation.TransportSolver,
        detectors: Mapping[str, MomentRT.detector.SurfaceDetector] | None = None,
    ) -> None:
        """Initialize the Simulation configuration.

        Parameters
        ----------
        medium : HomogeneousMedium or StructuredMedium
            Physical material configuration, compatible with the chosen solver.
        source : PointSource or UniformSurfaceSource
            Illumination configuration; SurfaceSource is a finite disk subclass.
        mesh : CartesianMesh
            Cartesian nodes used for sampling or finite-element discretization.
        solver : TransportSolver
            Transport model implementing the keyword-only solve protocol.
        detectors : mapping of str to SurfaceDetector, or None, optional
            Nonempty names mapped to passive contact detectors. None disables
            detector measurements.
            Default is None.
        """

        if detectors is not None and not isinstance(detectors, Mapping):
            raise TypeError("detectors must map names to SurfaceDetector objects.")

        detectors = dict(detectors or {})

        if any(
            not isinstance(name, str)
            or not name.strip()
            or not isinstance(detector, MomentRT.detector.SurfaceDetector)
            for name, detector in detectors.items()
        ):
            raise ValueError(
                "detectors must map nonempty names to SurfaceDetector objects."
            )

        self.medium = medium

        self.source = source

        self.mesh = mesh

        self.solver = solver

        self.detectors = detectors

    def validate(
        self, *, wavelengths: Quantity | None = None
    ) -> MomentRT.validation.ValidationReport:
        """Inspect configuration, material coverage, and mesh footprint sampling.

        Parameters
        ----------
        wavelengths : pint.Quantity or None, optional
            One-dimensional sweep wavelengths with length units. None checks
            the configured source wavelength.

        Returns
        -------
        ValidationReport
            Independent errors, cautions, and sampled source and detector areas.
            Validation does not solve transport or establish mesh convergence.
        """
        return MomentRT.validation._validate(simulation=self, wavelengths=wavelengths)

    def run(self) -> MomentRT.results.SimulationResult:
        """Solve the configured transport problem and evaluate all named detectors.

        Returns
        -------
        SimulationResult
            Solved fields, numerical diagnostics, and configured detector measurements.
        """

        solution = self.solver.solve(
            medium=self.medium, source=self.source, mesh=self.mesh
        )

        if not isinstance(solution, MomentRT.results.SimulationResult):
            solution = MomentRT.results.SimulationResult(fluence=solution)

        measurements = {
            name: detector.measure(result=solution)
            for name, detector in self.detectors.items()
        }

        solution.detectors.update(measurements)

        return solution

    def sweep(
        self, *, wavelengths: Quantity, optical_power: Quantity | None = None
    ) -> MomentRT.results.SpectralSimulationResult:
        """Solve wavelengths in supplied order without mutating this configuration.

        optical_power may be a scalar or a same-length vector of powers. When
        omitted, every solve uses the configured source's optical_power.
        This is a set of monochromatic solves, not broadband spectral integration.

        Parameters
        ----------
        wavelengths : pint.Quantity
            Nonempty one-dimensional vector of positive vacuum wavelengths.
        optical_power : pint.Quantity or None, optional
            Positive scalar power or vector matching wavelengths. None retains
            the configured source power for every monochromatic solve.
            Default is None.

        Returns
        -------
        SpectralSimulationResult
            Ordered independent monochromatic solutions. Finite sweeps reuse
            geometry but rebuild the wavelength-dependent operator and factorization.

        Notes
        -----
        These are independent monochromatic powers, not a broadband spectral
        integral. Finite sweeps reuse geometry, material labels, source footprint,
        and boundary quadrature once per call. Each wavelength gets a new bulk
        operator and sparse factorization.
        """

        started = time.perf_counter()

        MomentRT.units.validate_units(
            value=wavelengths, unit="m", name="wavelengths", positive=True
        )

        if wavelengths.ndim != 1 or not len(wavelengths):
            raise ValueError("wavelengths must be a nonempty vector.")

        if optical_power is not None:
            MomentRT.units.validate_units(
                value=optical_power, unit="W", name="optical_power", positive=True
            )

            if optical_power.shape not in ((), wavelengths.shape):
                raise ValueError(
                    "optical_power must be scalar or shaped like wavelengths."
                )

        # Resolve every table first, so an out-of-range sweep fails before solving.
        wavelength: Quantity

        for wavelength in wavelengths:
            self.medium.at(wavelength=wavelength)

        prepare = getattr(self.solver, "_prepare", None)

        preparation_started = time.perf_counter()

        solver = (
            self.solver
            if prepare is None
            else prepare(medium=self.medium, source=self.source, mesh=self.mesh)
        )

        preparation_time = time.perf_counter() - preparation_started

        results = []

        solve_times = []

        for index, wavelength in enumerate(wavelengths):
            source = copy.copy(self.source)

            source.wavelength = wavelength.copy()

            if optical_power is not None:
                source.optical_power = (
                    optical_power if optical_power.shape == () else optical_power[index]
                ).copy()

            simulation = Simulation(
                medium=self.medium,
                source=source,
                mesh=self.mesh,
                solver=solver,
                detectors=self.detectors,
            )

            solve_started = time.perf_counter()

            results.append(simulation.run())

            solve_times.append(time.perf_counter() - solve_started)

        return MomentRT.results.SpectralSimulationResult(
            wavelengths=wavelengths,
            results=results,
            diagnostics={
                "geometry_preparations": int(prepare is not None),
                "preparation_time": preparation_time * MomentRT.units.ureg.second,
                "solve_times": np.array(solve_times) * MomentRT.units.ureg.second,
                "total_time": (time.perf_counter() - started)
                * MomentRT.units.ureg.second,
            },
        )

    def convergence_study(
        self,
        *,
        mesh_shapes: Iterable[tuple[int, int, int]],
        power_rtol: float = 0.01,
        area_rtol: float = 0.02,
        power_atol: Quantity | None = None,
        required_consecutive: int = 2,
    ) -> MomentRT.results.DetectorConvergenceResult:
        """Compare detector signals and footprint sampling across explicit meshes.

        Parameters
        ----------
        mesh_shapes : iterable of tuple of int
            At least two coarse-to-fine node shapes. Each shape must refine at
            least one axis and coarsen none.
        power_rtol : float, optional
            Nonnegative relative power tolerance against adjacent and finest meshes.
            Default is 0.01.
        area_rtol : float, optional
            Nonnegative relative tolerance of sampled source and detector areas
            against their geometric areas.
            Default is 0.02.
        power_atol : pint.Quantity or None, optional
            Nonnegative absolute power tolerance. None means zero watts.
            Default is None.
        required_consecutive : int, optional
            Positive number of consecutive passing refinements required at the
            fine end of the mesh sequence.
            Default is 2.

        Returns
        -------
        DetectorConvergenceResult
            Per-level fields, powers, areas, failures, and a mesh recommendation
            only when the final refinement sequence meets all requested criteria.
        """

        return MomentRT.convergence._run_study(
            simulation=self,
            mesh_shapes=mesh_shapes,
            power_rtol=power_rtol,
            area_rtol=area_rtol,
            power_atol=power_atol,
            required_consecutive=required_consecutive,
        )

    def plot_system(
        self, *, length_unit: str = "mm", show: bool = True
    ) -> MomentRT.plotting.Figure:
        """Display the physical sample, source, and detector geometry in 3D.

        Parameters
        ----------
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Sample construction geometry with source and detector overlays.
        """

        medium = self.medium.at(wavelength=self.source.wavelength)

        figure = medium.plot(mesh=self.mesh, length_unit=length_unit, show=False)

        context = {
            "length_unit": length_unit,
            "sample_origin": getattr(medium, "origin", None),
            "sample_size": getattr(medium, "size", None),
        }

        MomentRT.plotting._add_overlays(
            figure=figure,
            overlays=[self.source],
            context=context,
        )

        MomentRT.plotting._add_overlays(
            figure=figure,
            overlays={
                f"Detector: {name}": detector
                for name, detector in self.detectors.items()
            },
            context=context,
        )

        if show:
            figure.show()

        return figure
