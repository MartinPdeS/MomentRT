"""Ideal isotropic point illumination inside an infinite medium."""

from __future__ import annotations
from collections.abc import Mapping
from typing import Any
import numpy as np
from pint import Quantity
import MomentRT


class SourceResolutionError(ValueError):
    """No boundary quadrature point samples the illumination footprint."""


class PointSource:
    """Isotropic point illumination embedded in an infinite medium.

    Parameters
    ----------
    position : pint.Quantity
        Absolute Cartesian position, shape (3,), with length units.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.
    optical_power : pint.Quantity
        Positive scalar optical power with power units.

    Notes
    -----
    Use with InfiniteDiffusion. The point-source singularity is not
    regularized, and the Cartesian mesh is an observation window.
    """

    def __init__(
        self, *, position: Quantity, wavelength: Quantity, optical_power: Quantity
    ) -> None:
        """Initialize the PointSource configuration.

        Parameters
        ----------
        position : pint.Quantity
            Absolute Cartesian position, shape (3,), with length units.
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.
        optical_power : pint.Quantity
            Positive scalar optical power with power units.
        """

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

    def _add_to_figure(
        self, *, figure: MomentRT.plotting.Figure, context: Mapping[str, Any]
    ) -> None:
        """Append this instrument geometry to an existing figure in place.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Figure to modify in place.
        context : mapping of str to object
            Coordinate units, legend labels, and physical sample bounds.
        """

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

    Parameters
    ----------
    face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
        Axis-aligned sample face; min and max refer to outward boundaries.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.
    optical_power : pint.Quantity
        Positive scalar optical power with power units.

    Notes
    -----
    Power is incoming diffuse power already coupled into the sample.
    Collimated propagation and specular Fresnel entry losses are not modeled.
    The solver applies the Robin source factor separately.
    """

    def __init__(
        self,
        *,
        face: MomentRT.boundary.Face,
        wavelength: Quantity,
        optical_power: Quantity,
    ) -> None:
        """Initialize the UniformSurfaceSource configuration.

        Parameters
        ----------
        face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
            Axis-aligned sample face; min and max refer to outward boundaries.
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.
        optical_power : pint.Quantity
            Positive scalar optical power with power units.
        """

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

    def _add_to_figure(
        self, *, figure: MomentRT.plotting.Figure, context: Mapping[str, Any]
    ) -> None:
        """Append this instrument geometry to an existing figure in place.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Figure to modify in place.
        context : mapping of str to object
            Coordinate units, legend labels, and physical sample bounds.
        """

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

        direction = np.zeros(3)

        direction[axis] = 1 if self.face.endswith("min") else -1

        length = 0.25 * min((upper - lower)[axes])

        MomentRT.plotting._add_direction_arrow(
            figure=figure,
            start=position - length * direction,
            direction=direction,
            length=length,
            color="crimson",
            label="Inward normal (diffuse illumination)",
            legendgroup=context.get("legendgroup"),
        )


class SurfaceSource(UniformSurfaceSource):
    """Circular top-hat incoming diffuse illumination on an axis-aligned face.

    position is an absolute 3D surface coordinate; radius is the spot radius.
    The entire disk must fit on the selected face. Surface quadrature samples
    its footprint and is normalized to optical_power; refine the mesh to
    resolve the disk edge and check the reported quadrature spot area.

    Parameters
    ----------
    face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
        Axis-aligned sample face; min and max refer to outward boundaries.
    position : pint.Quantity
        Absolute Cartesian position, shape (3,), with length units.
    radius : pint.Quantity
        Positive scalar radius with length units.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.
    optical_power : pint.Quantity
        Positive scalar optical power with power units.

    Notes
    -----
    The entire disk must fit on the selected Robin face. Boundary quadrature
    is normalized to the specified incoming power; energy conservation alone
    does not establish that the disk edge is spatially resolved.
    """

    def __init__(
        self,
        *,
        face: MomentRT.boundary.Face,
        position: Quantity,
        radius: Quantity,
        wavelength: Quantity,
        optical_power: Quantity,
    ) -> None:
        """Initialize the SurfaceSource configuration.

        Parameters
        ----------
        face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
            Axis-aligned sample face; min and max refer to outward boundaries.
        position : pint.Quantity
            Absolute Cartesian position, shape (3,), with length units.
        radius : pint.Quantity
            Positive scalar radius with length units.
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.
        optical_power : pint.Quantity
            Positive scalar optical power with power units.
        """

        super().__init__(face=face, wavelength=wavelength, optical_power=optical_power)

        MomentRT.units.validate_units(
            value=position, unit="m", name="position", shape=(3,)
        )

        MomentRT.units.validate_units(
            value=radius, unit="m", name="radius", shape=(), positive=True
        )

        self.position = position.copy()

        self.radius = radius.copy()

    @classmethod
    def on_face(
        cls,
        *,
        sample: MomentRT.medium.StructuredMedium,
        face: MomentRT.boundary.Face,
        coordinates: Quantity,
        radius: Quantity,
        wavelength: Quantity,
        optical_power: Quantity,
    ) -> SurfaceSource:
        """Place a circular source using absolute coordinates in a sample face.

        Parameters
        ----------
        sample : StructuredMedium
            Physical sample used to determine the face coordinate.
        face : str
            Cartesian face name, such as ``z_min``.
        coordinates : pint.Quantity
            Two absolute in-plane coordinates, ordered by Cartesian axis
            (y, z for x faces; x, z for y faces; x, y for z faces).
        radius : pint.Quantity
            Positive illumination radius with length units.
        wavelength : pint.Quantity
            Positive vacuum wavelength with length units.
        optical_power : pint.Quantity
            Incoming diffuse power already coupled into the sample.

        Returns
        -------
        SurfaceSource
            Source with its absolute three-dimensional position recorded.
        """
        position = MomentRT.boundary._face_position(
            sample=sample, face=face, coordinates=coordinates
        )

        source = cls(
            face=face,
            position=position,
            radius=radius,
            wavelength=wavelength,
            optical_power=optical_power,
        )

        MomentRT.boundary._validate_face_disk(
            sample=sample, face=face, position=source.position, radius=source.radius
        )

        return source

    def _add_to_figure(
        self, *, figure: MomentRT.plotting.Figure, context: Mapping[str, Any]
    ) -> None:
        """Append this instrument geometry to an existing figure in place.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Figure to modify in place.
        context : mapping of str to object
            Coordinate units, legend labels, and physical sample bounds.
        """

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

        MomentRT.plotting._add_footprint_fill(
            figure=figure,
            position=position,
            ring=ring,
            color="crimson",
            legendgroup=context.get("legendgroup"),
        )

        direction = np.zeros(3)

        direction["xyz".index(self.face[0])] = 1 if self.face.endswith("min") else -1

        length = 2 * self.radius.to(context["length_unit"]).magnitude

        MomentRT.plotting._add_direction_arrow(
            figure=figure,
            start=position - length * direction,
            direction=direction,
            length=length,
            color="crimson",
            label="Inward normal (diffuse illumination)",
            legendgroup=context.get("legendgroup"),
        )
