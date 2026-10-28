"""Composition of medium, illumination, and transport."""

import numpy as np

from .results import SimulationResult


class Simulation:
    def __init__(self, *, medium, source, mesh, solver):
        self.medium = medium

        self.source = source

        self.mesh = mesh

        self.solver = solver

    def run(self):
        solution = self.solver.solve(
            medium=self.medium, source=self.source, mesh=self.mesh
        )

        if isinstance(solution, SimulationResult):
            return solution

        return SimulationResult(fluence=solution)

    def plot_system(self, *, length_unit="mm", show=True):
        figure = self.medium.plot(mesh=self.mesh, length_unit=length_unit, show=False)

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

        if show:
            figure.show()

        return figure
