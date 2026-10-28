"""Resolved planar material maps independent of the transport mesh."""

import numpy as np

import MomentRT


class SampleSection:
    """Snapshot of cell-center samples on an axis-aligned plane.

    Arrays use (vertical, horizontal) order. Physical inputs retain their units;
    optical-property arrays are stored in SI. No interpolation is applied.
    """

    def __init__(self, *, sample, normal, position, resolution, wavelength=None):
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
            sample.background,
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
        quantity="material",
        length_unit="mm",
        value_unit=None,
        scale="linear",
        show=True,
    ):
        return MomentRT.plotting.plot_section(
            section=self,
            quantity=quantity,
            length_unit=length_unit,
            value_unit=value_unit,
            scale=scale,
            show=show,
        )
