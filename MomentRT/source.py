"""Ideal isotropic point illumination inside an infinite medium."""

import numpy as np

import MomentRT


class SourceResolutionError(ValueError):
    """No boundary quadrature point samples the illumination footprint."""


class PointSource:
    def __init__(self, *, position, wavelength, optical_power):
        MomentRT.units.validate_units(
            value=position, unit="m", name="position", shape=(3,)
        )

        MomentRT.units.validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        MomentRT.units.validate_units(
            value=optical_power, unit="W", name="optical_power", shape=(), positive=True
        )

        self.position = position

        self.wavelength = wavelength

        self.optical_power = optical_power

    def _add_to_figure(self, *, figure, context):
        go = MomentRT.plotting.require_plotly()

        position = self.position.to(context["length_unit"]).magnitude

        figure.add_trace(
            go.Scatter3d(
                x=[position[0]],
                y=[position[1]],
                z=[position[2]],
                mode="markers",
                marker={"color": "crimson", "size": 7},
                name=context.get("name") or "Point source",
                legendgroup=context.get("legendgroup"),
            )
        )


class UniformSurfaceSource:
    """Uniform incoming diffuse partial current across one complete sample face.

    optical_power is power already transmitted into the sample, before diffuse
    internal reflection. This does not model collimated LED entry or specular
    Fresnel loss. The boundary condition's source factor is applied by the solver.
    """

    def __init__(self, *, face, wavelength, optical_power):
        if face not in MomentRT.boundary.FACES:
            raise ValueError(f"face must be one of {MomentRT.boundary.FACES}.")

        MomentRT.units.validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        MomentRT.units.validate_units(
            value=optical_power, unit="W", name="optical_power", shape=(), positive=True
        )

        self.face = face

        self.wavelength = wavelength

        self.optical_power = optical_power

    def _add_to_figure(self, *, figure, context):
        origin = context.get("sample_origin")

        size = context.get("sample_size")

        if origin is None or size is None:
            raise ValueError(
                "UniformSurfaceSource overlays require finite-sample bounds; "
                "use a finite Diffusion result or simulation.plot_system()."
            )

        lower = origin.to(context["length_unit"]).magnitude

        upper = (origin + size).to(context["length_unit"]).magnitude

        axis = "xyz".index(self.face[0])

        axes = [coordinate for coordinate in range(3) if coordinate != axis]

        position = (lower + upper) / 2

        position[axis] = lower[axis] if self.face.endswith("min") else upper[axis]

        horizontal, vertical = np.meshgrid(
            [lower[axes[0]], upper[axes[0]]],
            [lower[axes[1]], upper[axes[1]]],
        )

        surface = np.broadcast_to(position, (2, 2, 3)).copy()

        surface[..., axes[0]] = horizontal

        surface[..., axes[1]] = vertical

        go = MomentRT.plotting.require_plotly()

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
                name=context.get("name") or "Full-face illumination",
                legendgroup=context.get("legendgroup"),
            )
        )

        figure.add_trace(
            go.Scatter3d(
                x=[position[0]],
                y=[position[1]],
                z=[position[2]],
                mode="markers",
                marker={"color": "crimson", "size": 7},
                name=context.get("name") or "Surface source",
                legendgroup=context.get("legendgroup"),
                showlegend=False,
            )
        )


class SurfaceSource(UniformSurfaceSource):
    """Circular top-hat incoming diffuse illumination on an axis-aligned face.

    position is an absolute 3D surface coordinate; radius is the spot radius.
    The entire disk must fit on the selected face. Surface quadrature samples
    its footprint and is normalized to optical_power; refine the mesh to
    resolve the disk edge and check the reported quadrature spot area.
    """

    def __init__(self, *, face, position, radius, wavelength, optical_power):
        super().__init__(face=face, wavelength=wavelength, optical_power=optical_power)

        MomentRT.units.validate_units(
            value=position, unit="m", name="position", shape=(3,)
        )

        MomentRT.units.validate_units(
            value=radius, unit="m", name="radius", shape=(), positive=True
        )

        self.position = position.copy()

        self.radius = radius.copy()

    def _add_to_figure(self, *, figure, context):
        go = MomentRT.plotting.require_plotly()

        position = self.position.to(context["length_unit"]).magnitude

        ring = MomentRT.plotting._circle_points(
            position=position,
            radius=self.radius.to(context["length_unit"]).magnitude,
            face=self.face,
        )

        figure.add_trace(
            go.Scatter3d(
                x=ring[:, 0],
                y=ring[:, 1],
                z=ring[:, 2],
                mode="lines",
                line={"color": "crimson", "width": 5},
                name=context.get("name") or "Illumination spot",
                legendgroup=context.get("legendgroup"),
            )
        )

        figure.add_trace(
            go.Scatter3d(
                x=[position[0]],
                y=[position[1]],
                z=[position[2]],
                mode="markers",
                marker={"color": "crimson", "size": 7},
                name=context.get("name") or "Surface source",
                legendgroup=context.get("legendgroup"),
                showlegend=False,
            )
        )
