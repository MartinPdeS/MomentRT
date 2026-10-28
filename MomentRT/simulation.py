"""Composition of medium, illumination, and transport."""

import numpy as np

from .results import SimulationResult, SpectralSimulationResult


class Simulation:
    def __init__(self, *, medium, source, mesh, solver, detectors=None):
        from collections.abc import Mapping

        from .detector import SurfaceDetector

        if detectors is not None and not isinstance(detectors, Mapping):
            raise TypeError("detectors must map names to SurfaceDetector objects.")

        detectors = dict(detectors or {})

        if any(
            not isinstance(name, str)
            or not name.strip()
            or not isinstance(detector, SurfaceDetector)
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

        if not isinstance(solution, SimulationResult):
            solution = SimulationResult(fluence=solution)

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
        from copy import copy

        from .units import validate_units

        validate_units(value=wavelengths, unit="m", name="wavelengths", positive=True)

        if wavelengths.ndim != 1 or not len(wavelengths):
            raise ValueError("wavelengths must be a nonempty vector.")

        if optical_power is not None:
            validate_units(
                value=optical_power, unit="W", name="optical_power", positive=True
            )

            if optical_power.shape not in ((), wavelengths.shape):
                raise ValueError(
                    "optical_power must be scalar or shaped like wavelengths."
                )

        # Resolve every table first, so an out-of-range sweep fails before solving.
        for wavelength in wavelengths:
            self.medium.at(wavelength=wavelength)

        results = []

        for index, wavelength in enumerate(wavelengths):
            source = copy(self.source)

            source.wavelength = wavelength.copy()

            if optical_power is not None:
                source.optical_power = (
                    optical_power if optical_power.shape == () else optical_power[index]
                ).copy()

            simulation = Simulation(
                medium=self.medium,
                source=source,
                mesh=self.mesh,
                solver=self.solver,
                detectors=self.detectors,
            )

            results.append(simulation.run())

        return SpectralSimulationResult(wavelengths=wavelengths, results=results)

    def plot_system(self, *, length_unit="mm", show=True):
        medium = self.medium.at(wavelength=self.source.wavelength)

        figure = medium.plot(mesh=self.mesh, length_unit=length_unit, show=False)

        from .plotting import require_plotly

        go = require_plotly()

        if hasattr(self.source, "position"):
            position = self.source.position.to(length_unit).magnitude
        else:
            position = (
                (self.medium.origin + self.medium.size / 2).to(length_unit).magnitude
            )

            axis = "xyz".index(self.source.face[0])

            position[axis] = self.medium.origin.to(length_unit).magnitude[axis]

            if self.source.face.endswith("max"):
                position[axis] += self.medium.size.to(length_unit).magnitude[axis]

        if hasattr(self.source, "radius"):
            angle = np.linspace(0, 2 * np.pi, 65)

            ring = np.repeat(position[None, :], len(angle), axis=0).astype(float)

            axes = [
                axis for axis in range(3) if axis != "xyz".index(self.source.face[0])
            ]

            ring[:, axes[0]] += self.source.radius.to(length_unit).magnitude * np.cos(
                angle
            )

            ring[:, axes[1]] += self.source.radius.to(length_unit).magnitude * np.sin(
                angle
            )

            figure.add_trace(
                go.Scatter3d(
                    x=ring[:, 0],
                    y=ring[:, 1],
                    z=ring[:, 2],
                    mode="lines",
                    line={"color": "crimson", "width": 5},
                    name="Illumination spot",
                )
            )

        if hasattr(self.source, "face") and not hasattr(self.source, "radius"):
            axes = [
                axis for axis in range(3) if axis != "xyz".index(self.source.face[0])
            ]

            lower = self.medium.origin.to(length_unit).magnitude

            upper = (self.medium.origin + self.medium.size).to(length_unit).magnitude

            horizontal, vertical = np.meshgrid(
                [lower[axes[0]], upper[axes[0]]],
                [lower[axes[1]], upper[axes[1]]],
            )

            surface = np.broadcast_to(position, (2, 2, 3)).copy()

            surface[..., axes[0]] = horizontal

            surface[..., axes[1]] = vertical

            figure.add_trace(
                go.Surface(
                    x=surface[..., 0],
                    y=surface[..., 1],
                    z=surface[..., 2],
                    surfacecolor=np.zeros((2, 2)),
                    colorscale=[[0, "crimson"], [1, "crimson"]],
                    opacity=0.4,
                    showscale=False,
                    showlegend=True,
                    name="Full-face illumination",
                )
            )

        figure.add_trace(
            go.Scatter3d(
                x=[position[0]],
                y=[position[1]],
                z=[position[2]],
                mode="markers",
                marker={"color": "crimson", "size": 7},
                name="Point source"
                if not hasattr(self.source, "face")
                else "Surface source",
            )
        )

        for name, detector in self.detectors.items():
            position = detector.position.to(length_unit).magnitude

            angle = np.linspace(0, 2 * np.pi, 65)

            ring = np.repeat(position[None, :], len(angle), axis=0).astype(float)

            axes = [axis for axis in range(3) if axis != "xyz".index(detector.face[0])]

            radius = detector.radius.to(length_unit).magnitude

            ring[:, axes[0]] += radius * np.cos(angle)

            ring[:, axes[1]] += radius * np.sin(angle)

            figure.add_trace(
                go.Scatter3d(
                    x=ring[:, 0],
                    y=ring[:, 1],
                    z=ring[:, 2],
                    mode="lines",
                    line={"color": "royalblue", "width": 5},
                    name=f"Detector: {name}",
                )
            )

            tip = position + 2 * radius * detector.collection_axis

            figure.add_trace(
                go.Scatter3d(
                    x=[position[0], tip[0]],
                    y=[position[1], tip[1]],
                    z=[position[2], tip[2]],
                    mode="lines+markers",
                    line={"color": "royalblue", "width": 3},
                    marker={"size": 3},
                    showlegend=False,
                    hovertemplate=f"Collection axis: {name}<extra></extra>",
                )
            )

        if show:
            figure.show()

        return figure
