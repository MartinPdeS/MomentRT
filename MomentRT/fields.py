"""Unit-bearing nodal scalar fields with mesh-aware visualization."""

from __future__ import annotations
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Literal
import numpy as np
from pint import Quantity
import MomentRT


class ScalarField:
    """Unit-bearing scalar samples at Cartesian mesh nodes.

    Parameters
    ----------
    mesh : CartesianMesh
        Cartesian nodes used for sampling or finite-element discretization.
    values : pint.Quantity
        Finite nodal scalar values with shape equal to mesh.shape.
        The supplied physical units are retained.
    name : str
        Nonempty name identifying this object.
    metadata : mapping of str to object, or None, optional
        Field provenance and optional physical sample bounds.
        Default is None.
    """

    def __init__(
        self,
        *,
        mesh: MomentRT.mesh.CartesianMesh,
        values: Quantity,
        name: str,
        metadata: Mapping[str, object] | None = None,
    ) -> None:
        """Initialize the ScalarField configuration.

        Parameters
        ----------
        mesh : CartesianMesh
            Cartesian nodes used for sampling or finite-element discretization.
        values : pint.Quantity
            Finite nodal scalar values with shape equal to mesh.shape.
            The supplied physical units are retained.
        name : str
            Nonempty name identifying this object.
        metadata : mapping of str to object, or None, optional
            Field provenance and optional physical sample bounds.
            Default is None.
        """

        if not isinstance(values, MomentRT.units.ureg.Quantity):
            raise TypeError("values must carry units.")

        if values.shape != mesh.shape or not np.all(np.isfinite(values.magnitude)):
            raise ValueError("Field values must be finite and match the mesh shape.")

        self.mesh = mesh

        self.values = values

        self.name = name

        self.metadata = dict(metadata or {})

    def to_pyvista(self, *, length_unit: str = "mm") -> MomentRT.plotting.ImageData:
        """Convert Cartesian data to a PyVista image grid.

        Parameters
        ----------
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.

        Returns
        -------
        pyvista.ImageData
            Cartesian VTK grid with coordinates converted to length_unit.
        """

        grid = self.mesh.to_pyvista(length_unit=length_unit)

        grid.point_data[self.name] = self.values.magnitude.ravel(order="F")

        grid.field_data["value_units"] = [str(self.values.units)]

        grid.field_data["coordinate_units"] = [length_unit]

        grid.field_data["association"] = ["point"]

        for key, value in self.metadata.items():
            grid.field_data[key] = [str(value)]

        return grid

    def export(self, *, path: str | Path, length_unit: str = "mm") -> None:
        """Write the scalar field and provenance to a VTK image file.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination filename ending in .vti. Existing files may be overwritten.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        """

        if Path(path).suffix != ".vti":
            raise ValueError("Cartesian fields must be exported with a .vti extension.")

        grid = self.to_pyvista(length_unit=length_unit)

        grid.save(path)

    def plot(
        self,
        *,
        view: Literal["volume", "isosurfaces", "slices"] = "volume",
        scale: Literal["linear", "log"] = "linear",
        length_unit: str = "mm",
        opacity: float | None = None,
        isosurface_count: int | None = None,
        overlays: Iterable[object] | None = None,
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Plot the scalar field with optional source and detector geometry.

        Sources, detectors, and layers can be overlaid. Each object's ``name``
        sets its legend label. Full-face sources and layers require physical
        sample bounds recorded by the finite solver.

        Parameters
        ----------
        view : {'volume', 'isosurfaces', 'slices'}, optional
            Three-dimensional scalar-field rendering mode.
            Default is 'volume'.
        scale : {'linear', 'log'}, optional
            Display scaling. Logarithmic scaling requires strictly positive values.
            Default is 'linear'.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        opacity : float or None, optional
            Initial surface opacity in [0, 1]. None selects the rendering default.
            The volume view exposes an interactive opacity slider.
            Default is None.
        isosurface_count : int or None, optional
            Initial positive number of iso levels for volume or isosurface
            rendering. Both views expose a density slider. None uses 8 for
            volume and 6 for isosurfaces; invalid for slices.
            Default is None.
        overlays : iterable of objects or None, optional
            Sources, detectors, or layers implementing the plotting hook.
            Their ``name`` sets the legend label. Full-face sources and layers
            require finite-sample bounds.
            Default is None.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable interactive figure, also displayed when show is True.

        Notes
        -----
        Logarithmic rendering rejects nonpositive data rather than clipping it.
        Full-face source and layer overlays require physical sample bounds
        saved by the finite solver; an infinite observation window is not a
        physical sample.
        """

        return MomentRT.plotting.plot_field(
            field=self,
            view=view,
            scale=scale,
            length_unit=length_unit,
            opacity=opacity,
            isosurface_count=isosurface_count,
            overlays=overlays,
            show=show,
        )


class BoundaryField:
    """Snapshot of escaping irradiance at boundary triangle quadrature points.

    positions, area_weights, and escaping_irradiance share a leading shape.
    Escaping irradiance includes the incident-field correction on a source face.
    It is not an angle-resolved radiance.

    Parameters
    ----------
    face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
        Axis-aligned sample face; min and max refer to outward boundaries.
    positions : pint.Quantity
        Cartesian sampling coordinates with length units and trailing axis (3,).
    area_weights : pint.Quantity
        Positive quadrature area weights with area units; shape matches
        the leading axes of positions.
    escaping_irradiance : pint.Quantity
        Escaping diffuse irradiance with power-per-area units; shape matches
        area_weights. Values are not clipped to zero.
    sample_origin : pint.Quantity
        Physical sample lower corner, shape (3,), with length units.
    sample_size : pint.Quantity
        Positive physical sample extents, shape (3,), with length units.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.
    """

    def __init__(
        self,
        *,
        face: MomentRT.boundary.Face,
        positions: Quantity,
        area_weights: Quantity,
        escaping_irradiance: Quantity,
        sample_origin: Quantity,
        sample_size: Quantity,
        wavelength: Quantity,
    ) -> None:
        """Initialize the BoundaryField configuration.

        Parameters
        ----------
        face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
            Axis-aligned sample face; min and max refer to outward boundaries.
        positions : pint.Quantity
            Cartesian sampling coordinates with length units and trailing axis (3,).
        area_weights : pint.Quantity
            Positive quadrature area weights with area units; shape matches
            the leading axes of positions.
        escaping_irradiance : pint.Quantity
            Escaping diffuse irradiance with power-per-area units; shape matches
            area_weights. Values are not clipped to zero.
        sample_origin : pint.Quantity
            Physical sample lower corner, shape (3,), with length units.
        sample_size : pint.Quantity
            Positive physical sample extents, shape (3,), with length units.
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.
        """

        if face not in MomentRT.boundary.FACES:
            raise ValueError(f"face must be one of {MomentRT.boundary.FACES}.")

        MomentRT.units.validate_units(value=positions, unit="m", name="positions")

        if positions.ndim < 2 or positions.shape[-1] != 3:
            raise ValueError("positions must have a trailing axis of length three.")

        MomentRT.units.validate_units(
            value=area_weights,
            unit="m^2",
            name="area_weights",
            shape=positions.shape[:-1],
            positive=True,
        )

        MomentRT.units.validate_units(
            value=escaping_irradiance,
            unit="W/m^2",
            name="escaping_irradiance",
            shape=positions.shape[:-1],
        )

        MomentRT.units.validate_units(
            value=sample_origin, unit="m", name="sample_origin", shape=(3,)
        )

        MomentRT.units.validate_units(
            value=sample_size, unit="m", name="sample_size", shape=(3,), positive=True
        )

        MomentRT.units.validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        self.face = face

        self.positions = positions.copy()

        self.area_weights = area_weights.copy()

        self.escaping_irradiance = escaping_irradiance.copy()

        self.sample_origin = sample_origin.copy()

        self.sample_size = sample_size.copy()

        self.wavelength = wavelength.copy()

    @property
    def escaping_power(self) -> Quantity:
        """Integrate the escaping irradiance over the complete boundary face.

        Returns
        -------
        pint.Quantity
            Integral of escaping irradiance over this entire face, with power units.
        """

        return (self.escaping_irradiance * self.area_weights).sum()


class VectorField:
    """Element-centered vectors, without smoothing across material interfaces.

    Parameters
    ----------
    positions : pint.Quantity
        Element-center coordinates, shape (N, 3), with length units; N > 0.
    values : pint.Quantity
        Finite vector values, shape (N, 3), with physical units.
        Values are not smoothed across material interfaces.
    name : str
        Nonempty name identifying this object.
    """

    def __init__(self, *, positions: Quantity, values: Quantity, name: str) -> None:
        """Initialize the VectorField configuration.

        Parameters
        ----------
        positions : pint.Quantity
            Element-center coordinates, shape (N, 3), with length units; N > 0.
        values : pint.Quantity
            Finite vector values, shape (N, 3), with physical units.
            Values are not smoothed across material interfaces.
        name : str
            Nonempty name identifying this object.
        """

        MomentRT.units.validate_units(value=positions, unit="m", name="positions")

        if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) == 0:
            raise ValueError("positions must have shape (N, 3).")

        if (
            not isinstance(values, MomentRT.units.ureg.Quantity)
            or values.shape != positions.shape
            or not np.all(np.isfinite(values.magnitude))
        ):
            raise ValueError(
                "values must be finite unit-bearing vectors matching positions."
            )

        self.positions = positions

        self.values = values

        self.name = name

    def plot(
        self, *, max_vectors: int = 500, length_unit: str = "mm", show: bool = True
    ) -> MomentRT.plotting.Figure:
        """Create an interactive visualization of this object.

        Parameters
        ----------
        max_vectors : int, optional
            Positive maximum number of evenly selected element-centered arrows.
            Default is 500.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable interactive figure, also displayed when show is True.
        """

        if (
            isinstance(max_vectors, bool)
            or not isinstance(max_vectors, (int, np.integer))
            or max_vectors < 1
        ):
            raise ValueError("max_vectors must be a positive integer.")

        go = MomentRT.plotting.require_plotly()

        selection = np.linspace(
            0, len(self.values) - 1, min(max_vectors, len(self.values)), dtype=int
        )

        points = self.positions.to(length_unit).magnitude[selection]

        vectors = self.values.magnitude[selection]

        figure = MomentRT.plotting.make_figure(length_unit=length_unit)

        figure.add_trace(
            go.Cone(
                x=points[:, 0],
                y=points[:, 1],
                z=points[:, 2],
                u=vectors[:, 0],
                v=vectors[:, 1],
                w=vectors[:, 2],
                colorscale="Viridis",
                sizemode="scaled",
                sizeref=0.25,
                colorbar={"title": {"text": f"{self.name} [{self.values.units:~}]"}},
            )
        )

        figure.update_layout(
            title=f"{self.name} (element centers; relative arrow lengths)"
        )

        if show:
            figure.show()

        return figure
