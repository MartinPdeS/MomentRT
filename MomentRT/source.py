"""Monochromatic and broadband point and diffuse surface illumination."""

from __future__ import annotations
from collections.abc import Iterable, Mapping
from html import escape
from typing import Any
import numpy as np
import TypedUnit
from pint import Quantity
from scipy import integrate
import MomentRT


class SourceResolutionError(ValueError):
    """No boundary quadrature point samples the illumination footprint."""


class SourceSpectrum:
    """Tabulated incoming spectral power density for a broadband source.

    Parameters
    ----------
    wavelengths : pint.Quantity
        Strictly increasing positive vacuum wavelengths, shape (n,), n >= 2.
        Their endpoints define the integration interval.
    spectral_power : pint.Quantity
        Nonnegative spectral power density, shape (n,), such as W / nm.
        The trapezoidal integral must be positive. Surface-source power is
        diffuse power already coupled into the sample.
    """

    def __init__(self, *, wavelengths: Quantity, spectral_power: Quantity) -> None:
        """Validate and snapshot a spectral power-density table.

        Parameters
        ----------
        wavelengths : pint.Quantity
            Increasing positive vacuum wavelengths, shape (n,), n >= 2.
        spectral_power : pint.Quantity
            Nonnegative spectral power density, shape (n,), with W / length units.
        """
        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=wavelengths,
            unit="m",
            name="wavelengths",
            positive=True,
        )
        if wavelengths.ndim != 1 or len(wavelengths) < 2:
            raise ValueError("wavelengths must be a vector of at least two samples.")
        if np.any(np.diff(wavelengths.to("m").magnitude) <= 0):
            raise ValueError("wavelengths must be strictly increasing.")
        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=spectral_power,
            unit="W/m",
            name="spectral_power",
            shape=wavelengths.shape,
            nonnegative=True,
        )
        if not np.any(spectral_power.magnitude > 0):
            raise ValueError("spectral_power must have a positive integral.")
        self.wavelengths = wavelengths.copy()
        self.spectral_power = spectral_power.copy()

    @property
    def optical_power(self) -> Quantity:
        """Return the trapezoidal spectral integral over the supplied interval.

        Returns
        -------
        pint.Quantity
            Total incoming optical power in watts.
        """
        return (
            integrate.trapezoid(
                self.spectral_power.to("W/m").magnitude,
                x=self.wavelengths.to("m").magnitude,
            )
            * MomentRT.units.ureg.watt
        )

    def resample(self, *, wavelengths: Quantity) -> SourceSpectrum:
        """Interpolate density onto an explicit grid without extrapolation.

        Parameters
        ----------
        wavelengths : pint.Quantity
            Strictly increasing positive vacuum wavelengths, shape (n,), n >= 2,
            within the original table. Include existing knots when preserving
            the piecewise linear source and its total optical power.

        Returns
        -------
        SourceSpectrum
            Independent interpolated table, retaining spectral density units.

        Raises
        ------
        ValueError
            If the grid is invalid, outside coverage, or loses all source power.
        """
        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=wavelengths,
            unit="m",
            name="wavelengths",
            positive=True,
        )
        if wavelengths.ndim != 1 or len(wavelengths) < 2:
            raise ValueError("wavelengths must contain at least two samples.")
        requested = wavelengths.to(self.wavelengths.units).magnitude
        grid = self.wavelengths.magnitude
        if np.any(np.diff(requested) <= 0):
            raise ValueError("wavelengths must be strictly increasing.")
        tolerance = 8 * np.finfo(float).eps * max(abs(grid[0]), abs(grid[-1]))
        if np.any(requested < grid[0] - tolerance) or np.any(
            requested > grid[-1] + tolerance
        ):
            raise ValueError("Requested wavelengths are outside the source spectrum.")
        values = np.interp(
            np.clip(requested, grid[0], grid[-1]),
            grid,
            self.spectral_power.magnitude,
        )
        return SourceSpectrum(
            wavelengths=wavelengths,
            spectral_power=values * self.spectral_power.units,
        )

    def plot(
        self,
        *,
        wavelength_unit: str = "nm",
        power_unit: str = "uW/nm",
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Plot the supplied spectral power-density samples.

        Parameters
        ----------
        wavelength_unit : str, optional
            Display length unit. Default is 'nm'.
        power_unit : str, optional
            Display spectral power-density unit. Default is 'uW/nm'.
        show : bool, optional
            Display the composable figure immediately. Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Source density curve with physical wavelength and power labels.
        """
        figure = MomentRT.plotting.go.Figure(
            data=[
                MomentRT.plotting.go.Scatter(
                    x=self.wavelengths.to(wavelength_unit).magnitude,
                    y=self.spectral_power.to(power_unit).magnitude,
                    mode="lines+markers",
                    name="Source spectral power",
                )
            ],
        )
        figure.update_layout(
            title="Source emission spectrum",
            xaxis_title=f"Vacuum wavelength ({wavelength_unit})",
            yaxis_title=f"Spectral power ({power_unit})",
        )
        if show:
            figure.show()
        return figure


class _SourceEmission:
    """Shared monochromatic or broadband emission configuration."""

    def _configure_emission(
        self,
        *,
        wavelength: Quantity | None,
        optical_power: Quantity | None,
        spectrum: SourceSpectrum | None,
    ) -> None:
        """Require either a monochromatic pair or a spectral-density table."""
        if spectrum is not None:
            if not isinstance(spectrum, SourceSpectrum):
                raise TypeError("spectrum must be a SourceSpectrum object.")
            if wavelength is not None or optical_power is not None:
                raise ValueError(
                    "Supply spectrum or wavelength and optical_power exclusively."
                )
        else:
            if wavelength is None or optical_power is None:
                raise ValueError(
                    "Supply both wavelength and optical_power, or spectrum."
                )
            TypedUnit.validate_quantity(
                registry=MomentRT.units.ureg,
                value=wavelength,
                unit="m",
                name="wavelength",
                shape=(),
                positive=True,
            )
            TypedUnit.validate_quantity(
                registry=MomentRT.units.ureg,
                value=optical_power,
                unit="W",
                name="optical_power",
                shape=(),
                positive=True,
            )
        self.spectrum = spectrum
        self._wavelength = wavelength
        self._optical_power = optical_power

    @property
    def wavelength(self) -> Quantity:
        """Return the monochromatic vacuum wavelength.

        Returns
        -------
        pint.Quantity
            Positive scalar wavelength with length units.

        Raises
        ------
        ValueError
            If the source is broadband and has no single wavelength.
        """
        if self.spectrum is not None or self._wavelength is None:
            raise ValueError(
                "A broadband source has no single wavelength; use Simulation.run()."
            )
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value: Quantity) -> None:
        if self.spectrum is not None:
            raise ValueError("A broadband source's wavelengths belong to its spectrum.")
        self._wavelength = value

    @property
    def optical_power(self) -> Quantity:
        """Return monochromatic power or the total broadband spectral integral.

        Returns
        -------
        pint.Quantity
            Scalar optical power with power units. Broadband power integrates
            the source's spectral density over its supplied wavelength interval.
        """
        if self.spectrum is not None:
            return self.spectrum.optical_power
        assert self._optical_power is not None
        return self._optical_power

    @optical_power.setter
    def optical_power(self, value: Quantity) -> None:
        if self.spectrum is not None:
            raise ValueError(
                "Set the source spectrum to change broadband optical power."
            )
        self._optical_power = value


class PointSource(_SourceEmission):
    """Isotropic point illumination embedded in an infinite medium.

    Parameters
    ----------
    position : pint.Quantity
        Absolute Cartesian position, shape (3,), with length units.
    wavelength : pint.Quantity or None, optional
        Positive scalar vacuum wavelength, with length units.
    optical_power : pint.Quantity or None, optional
        Positive scalar optical power with power units.
    spectrum : SourceSpectrum or None, optional
        Broadband spectral density. Supply this or both wavelength and
        optical_power. Default is None.
    name : str, optional
        Legend label. Default is ``"Point source"``.
    color : str, optional
        Plotly-compatible color. Default is ``"crimson"``.

    Notes
    -----
    Use with InfiniteDiffusion. The point-source singularity is not
    regularized, and the Cartesian mesh is an observation window.
    """

    def __init__(
        self,
        *,
        position: Quantity,
        wavelength: Quantity | None = None,
        optical_power: Quantity | None = None,
        spectrum: SourceSpectrum | None = None,
        name: str = "Point source",
        color: str = "crimson",
    ) -> None:
        """Initialize the PointSource configuration.

        Parameters
        ----------
        position : pint.Quantity
            Absolute Cartesian position, shape (3,), with length units.
        wavelength : pint.Quantity or None, optional
            Positive scalar vacuum wavelength, with length units.
        optical_power : pint.Quantity or None, optional
            Positive scalar optical power with power units.
        spectrum : SourceSpectrum or None, optional
            Broadband spectral density. Supply this or both wavelength and
            optical_power. Default is None.
        name : str, optional
            Nonempty legend label.
        color : str, optional
            Nonempty Plotly-compatible color string.
        """

        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=position,
            unit="m",
            name="position",
            shape=(3,),
        )

        self._configure_emission(
            wavelength=wavelength,
            optical_power=optical_power,
            spectrum=spectrum,
        )

        if not isinstance(name, str) or not name.strip():
            raise ValueError("name must be a nonempty string.")

        if not isinstance(color, str) or not color.strip():
            raise ValueError("color must be a nonempty Plotly color string.")

        self.position = position

        self.name = name

        self.color = color

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

        go = MomentRT.plotting.go

        position = self.position.to(context["length_unit"]).magnitude

        figure.add_trace(
            go.Scatter3d(
                x=[position[0]],
                y=[position[1]],
                z=[position[2]],
                mode="markers",
                marker={"color": self.color, "size": 7},
                name=context["name"],
                legendgroup=context.get("legendgroup"),
            )
        )


class UniformSurfaceSource(_SourceEmission):
    """Uniform incoming diffuse partial current across one complete sample face.

    optical_power is power already transmitted into the sample, before diffuse
    internal reflection. This does not model collimated LED entry or specular
    Fresnel loss. The boundary condition's source factor is applied by the solver.

    Parameters
    ----------
    face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
        Axis-aligned sample face; min and max refer to outward boundaries.
    wavelength : pint.Quantity or None, optional
        Positive scalar vacuum wavelength, with length units.
    optical_power : pint.Quantity or None, optional
        Positive scalar optical power with power units.
    spectrum : SourceSpectrum or None, optional
        Broadband spectral density. Supply this or both wavelength and
        optical_power. Default is None.
    name : str, optional
        Legend label. Default is ``"Full-face illumination"``.
    color : str, optional
        Plotly-compatible color. Default is ``"crimson"``.

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
        wavelength: Quantity | None = None,
        optical_power: Quantity | None = None,
        spectrum: SourceSpectrum | None = None,
        name: str = "Full-face illumination",
        color: str = "crimson",
    ) -> None:
        """Initialize the UniformSurfaceSource configuration.

        Parameters
        ----------
        face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
            Axis-aligned sample face; min and max refer to outward boundaries.
        wavelength : pint.Quantity or None, optional
            Positive scalar vacuum wavelength, with length units.
        optical_power : pint.Quantity or None, optional
            Positive scalar optical power with power units.
        spectrum : SourceSpectrum or None, optional
            Broadband spectral density. Supply this or both wavelength and
            optical_power. Default is None.
        name : str, optional
            Nonempty legend label.
        color : str, optional
            Nonempty Plotly-compatible color string.
        """

        if face not in MomentRT.boundary.FACES:
            raise ValueError(f"face must be one of {MomentRT.boundary.FACES}.")

        self._configure_emission(
            wavelength=wavelength,
            optical_power=optical_power,
            spectrum=spectrum,
        )

        if not isinstance(name, str) or not name.strip():
            raise ValueError("name must be a nonempty string.")

        if not isinstance(color, str) or not color.strip():
            raise ValueError("color must be a nonempty Plotly color string.")

        self.face = face

        self.name = name

        self.color = color

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

        go = MomentRT.plotting.go

        figure.add_trace(
            go.Surface(
                x=surface[..., 0],
                y=surface[..., 1],
                z=surface[..., 2],
                surfacecolor=np.zeros((2, 2)),
                colorscale=[[0, self.color], [1, self.color]],
                opacity=0.4,
                showscale=False,
                showlegend=True,
                name=context["name"],
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
            color=self.color,
            label=f"Inward normal: {self.name}",
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
    wavelength : pint.Quantity or None, optional
        Positive scalar vacuum wavelength, with length units.
    optical_power : pint.Quantity or None, optional
        Positive scalar optical power with power units.
    spectrum : SourceSpectrum or None, optional
        Broadband spectral density. Supply this or both wavelength and
        optical_power. Default is None.
    name : str, optional
        Legend label. Default is ``"Illumination spot"``.
    color : str, optional
        Plotly-compatible color. Default is ``"crimson"``.

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
        wavelength: Quantity | None = None,
        optical_power: Quantity | None = None,
        spectrum: SourceSpectrum | None = None,
        name: str = "Illumination spot",
        color: str = "crimson",
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
        wavelength : pint.Quantity or None, optional
            Positive scalar vacuum wavelength, with length units.
        optical_power : pint.Quantity or None, optional
            Positive scalar optical power with power units.
        spectrum : SourceSpectrum or None, optional
            Broadband spectral density. Supply this or both wavelength and
            optical_power. Default is None.
        name : str, optional
            Nonempty legend label.
        color : str, optional
            Nonempty Plotly-compatible color string.
        """

        super().__init__(
            face=face,
            wavelength=wavelength,
            optical_power=optical_power,
            spectrum=spectrum,
            name=name,
            color=color,
        )

        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=position,
            unit="m",
            name="position",
            shape=(3,),
        )

        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=radius,
            unit="m",
            name="radius",
            shape=(),
            positive=True,
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
        wavelength: Quantity | None = None,
        optical_power: Quantity | None = None,
        spectrum: SourceSpectrum | None = None,
        name: str = "Illumination spot",
        color: str = "crimson",
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
        wavelength : pint.Quantity or None, optional
            Positive vacuum wavelength with length units.
        optical_power : pint.Quantity or None, optional
            Incoming diffuse power already coupled into the sample.
        spectrum : SourceSpectrum or None, optional
            Broadband spectral density. Supply this or both wavelength and
            optical_power. Default is None.
        name : str, optional
            Legend label. Default is ``"Illumination spot"``.
        color : str, optional
            Plotly-compatible color. Default is ``"crimson"``.

        Returns
        -------
        SurfaceSource
            Source with its absolute three-dimensional position recorded.
        """
        position = MomentRT._face_geometry.face_position(
            sample=sample, face=face, coordinates=coordinates
        )

        source = cls(
            face=face,
            position=position,
            radius=radius,
            wavelength=wavelength,
            optical_power=optical_power,
            spectrum=spectrum,
            name=name,
            color=color,
        )

        MomentRT._face_geometry.validate_face_disk(
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

        go = MomentRT.plotting.go

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
                line={"color": self.color, "width": 5},
                name=context["name"],
                legendgroup=context.get("legendgroup"),
            )
        )

        MomentRT.plotting._add_footprint_fill(
            figure=figure,
            position=position,
            ring=ring,
            color=self.color,
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
            color=self.color,
            label=f"Inward normal: {self.name}",
            legendgroup=context.get("legendgroup"),
        )


class SurfaceSourceGroup(UniformSurfaceSource):
    """Simultaneous illumination from multiple circular surface footprints.

    Each member disk is normalized independently to its configured optical
    power. All disks must share a face and wavelength. The solver assembles
    their summed incoming profile into one right-hand side.

    Parameters
    ----------
    sources : iterable of SurfaceSource
        Nonempty set of disk sources sharing one face and wavelength.
    name : str, optional
        Group name. Default is ``"Surface source group"``.
    color : str, optional
        Fallback Plotly color for the group. Member colors are used for plots.
    """

    def __init__(
        self,
        *,
        sources: Iterable[SurfaceSource],
        name: str = "Surface source group",
        color: str = "crimson",
    ) -> None:
        """Initialize a set of simultaneously active surface sources.

        Parameters
        ----------
        sources : iterable of SurfaceSource
            Nonempty set of disk sources sharing one face and wavelength.
        name : str, optional
            Nonempty group name used for simulation provenance.
        color : str, optional
            Nonempty Plotly-compatible fallback color.
        """

        sources = tuple(sources)

        if not sources or any(
            not isinstance(source, SurfaceSource) for source in sources
        ):
            raise ValueError("sources must contain at least one SurfaceSource.")

        first = sources[0]

        if any(source.spectrum is not None for source in sources):
            raise ValueError(
                "SurfaceSourceGroup currently requires monochromatic members."
            )

        if any(
            source.face != first.face
            or not np.isclose(
                source.wavelength.to("m").magnitude,
                first.wavelength.to("m").magnitude,
                rtol=1e-12,
                atol=0,
            )
            for source in sources[1:]
        ):
            raise ValueError("All grouped sources must share a face and wavelength.")

        optical_power = sum(
            (source.optical_power for source in sources),
            0 * MomentRT.units.ureg.watt,
        )

        super().__init__(
            face=first.face,
            wavelength=first.wavelength,
            optical_power=optical_power,
            name=name,
            color=color,
        )

        self.sources = sources

    def _add_to_figure(
        self, *, figure: MomentRT.plotting.Figure, context: Mapping[str, Any]
    ) -> None:
        """Draw every member footprint with its own label and color.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Figure to modify in place.
        context : mapping of str to object
            Coordinate units, sample bounds, and legend grouping.
        """

        for source in self.sources:
            source._add_to_figure(
                figure=figure,
                context={
                    **context,
                    "name": escape(source.name),
                    "legendgroup": context.get("legendgroup"),
                },
            )
