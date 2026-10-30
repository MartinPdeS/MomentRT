"""Resolved planar material maps independent of the transport mesh."""

from __future__ import annotations
from typing import Literal, cast
import numpy as np
from pint import Quantity
import MomentRT


class SampleSection:
    """Snapshot of cell-center samples on an axis-aligned plane.

    Arrays use (vertical, horizontal) order. Physical inputs retain their units;
    optical-property arrays are stored in SI. No interpolation is applied.

    Parameters
    ----------
    sample : StructuredMedium
        Finite rectangular sample with an explicit background.
    normal : {'x', 'y', 'z'}
        Axis normal to the section plane.
    position : pint.Quantity
        Absolute scalar coordinate along the plane normal, with length units.
    resolution : tuple of int
        Horizontal and vertical cell counts, each at least two.
        The stored section arrays have shape (vertical, horizontal).
    wavelength : pint.Quantity or None, optional
        Vacuum wavelength used to resolve spectral materials. Required when
        any material is spectral; None retains fixed coefficients.
        Default is None.

    Attributes
    ----------
    labels : numpy.ndarray of int
        Cell-center material labels, shape (vertical, horizontal). Zero denotes
        background; positive labels follow structure insertion order.
    fields : dict of str to pint.Quantity
        Absorption and reduced scattering in inverse meters and dimensionless
        refractive index, with the same shape as labels.
    names : tuple of str
        Background and structure names indexed by labels.
    horizontal : pint.Quantity
        One-dimensional horizontal cell-center coordinates in meters.
    vertical : pint.Quantity
        One-dimensional vertical cell-center coordinates in meters.
    extent : pint.Quantity
        Horizontal and vertical lower/upper bounds, shape (2, 2), in meters.
    axes : tuple of str
        Names of the horizontal and vertical Cartesian axes.
    """

    def __init__(
        self,
        *,
        sample: MomentRT.medium.StructuredMedium,
        normal: Literal["x", "y", "z"],
        position: Quantity,
        resolution: tuple[int, int],
        wavelength: Quantity | None = None,
    ) -> None:
        """Initialize the SampleSection configuration.

        Parameters
        ----------
        sample : StructuredMedium
            Finite rectangular sample with an explicit background.
        normal : {'x', 'y', 'z'}
            Axis normal to the section plane.
        position : pint.Quantity
            Absolute scalar coordinate along the plane normal, with length units.
        resolution : tuple of int
            Horizontal and vertical cell counts, each at least two.
            The stored section arrays have shape (vertical, horizontal).
        wavelength : pint.Quantity or None, optional
            Vacuum wavelength used to resolve spectral materials. Required when
            any material is spectral; None retains fixed coefficients.
            Default is None.
        """

        sample = sample._at_for_display(wavelength=wavelength)

        if normal not in ("x", "y", "z"):
            raise ValueError("normal must be x, y, or z.")

        MomentRT.units.validate_units(
            value=position, unit="m", name="position", shape=()
        )

        if len(resolution) != 2 or any(
            isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2
            for n in resolution
        ):
            raise ValueError("resolution must contain two integers of at least two.")

        normal_axis = "xyz".index(normal)

        plane_axes = [axis for axis in range(3) if axis != normal_axis]

        lower = sample.origin.to("m").magnitude

        size = sample.size.to("m").magnitude

        coordinate = position.to("m").magnitude

        if (
            coordinate < lower[normal_axis]
            or coordinate > lower[normal_axis] + size[normal_axis]
        ):
            raise ValueError("position must lie inside the sample.")

        coordinates = [
            lower[axis] + (np.arange(count) + 0.5) * size[axis] / count
            for axis, count in zip(plane_axes, resolution)
        ]

        horizontal, vertical = np.meshgrid(*coordinates, indexing="xy")

        points = np.empty((*horizontal.shape, 3))

        points[..., normal_axis] = coordinate

        points[..., plane_axes[0]] = horizontal

        points[..., plane_axes[1]] = vertical

        self.labels = sample.material_labels(points=points * MomentRT.units.ureg.meter)

        self.normal = normal

        self.position = position.copy()

        self.wavelength = None if wavelength is None else wavelength.copy()

        self.axes = tuple("xyz"[axis] for axis in plane_axes)

        self.horizontal = coordinates[0] * MomentRT.units.ureg.meter

        self.vertical = coordinates[1] * MomentRT.units.ureg.meter

        self.extent = (
            np.array([[lower[axis], lower[axis] + size[axis]] for axis in plane_axes])
            * MomentRT.units.ureg.meter
        )

        self.names = ("Background", *(shape.name for shape in sample.structures))

        materials = [
            cast(MomentRT.material.Material, sample.background),
            *(shape.material for shape in sample.structures),
        ]

        self.fields = {}

        for name, unit in (
            ("absorption", "1/m"),
            ("reduced_scattering", "1/m"),
            ("refractive_index", "dimensionless"),
        ):
            values = np.array(
                [
                    getattr(material, name).to(unit).magnitude
                    if name != "refractive_index"
                    else material.refractive_index
                    for material in materials
                ]
            )

            self.fields[name] = values[self.labels] * MomentRT.units.ureg(unit)

    def plot(
        self,
        *,
        quantity: Literal[
            "material", "absorption", "reduced_scattering", "refractive_index"
        ] = "material",
        length_unit: str = "mm",
        value_unit: str | None = None,
        scale: Literal["linear", "log"] = "linear",
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Create an interactive visualization of this object.

        Parameters
        ----------
        quantity : {'material', 'absorption', 'reduced_scattering', 'refractive_index'}, optional
            Material labels or optical property to display.
            Default is 'material'.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        value_unit : str or None, optional
            Unit for displayed optical-property values; None uses stored SI units.
            Default is None.
        scale : {'linear', 'log'}, optional
            Display scaling. Logarithmic scaling requires strictly positive values.
            Default is 'linear'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable interactive figure, also displayed when show is True.
        """

        return MomentRT.plotting.plot_section(
            section=self,
            quantity=quantity,
            length_unit=length_unit,
            value_unit=value_unit,
            scale=scale,
            show=show,
        )
