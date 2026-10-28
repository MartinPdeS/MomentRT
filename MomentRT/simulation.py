"""Composition of medium, illumination, and transport."""

import copy
import time
from collections.abc import Mapping

import numpy as np

import MomentRT


class Simulation:
    def __init__(self, *, medium, source, mesh, solver, detectors=None):
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

    def run(self):
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

    def sweep(self, *, wavelengths, optical_power=None):
        """Solve wavelengths in supplied order without mutating this configuration.

        optical_power may be a scalar or a same-length vector of powers. When
        omitted, every solve uses the configured source's optical_power.
        This is a set of monochromatic solves, not broadband spectral integration.
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
        mesh_shapes,
        power_rtol=0.01,
        area_rtol=0.02,
        power_atol=None,
        required_consecutive=2,
    ):
        """Compare detector signals and footprint sampling across explicit meshes."""

        return MomentRT.convergence._run_study(
            simulation=self,
            mesh_shapes=mesh_shapes,
            power_rtol=power_rtol,
            area_rtol=area_rtol,
            power_atol=power_atol,
            required_consecutive=required_consecutive,
        )

    def plot_system(self, *, length_unit="mm", show=True):
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
