"""Material-bearing geometry for finite samples; z increases into the sample."""

import numpy as np

import MomentRT


def _validate_region(*, material, name):
    if not isinstance(
        material, (MomentRT.material.Material, MomentRT.material.SpectralMaterial)
    ):
        raise TypeError("material must be a Material or SpectralMaterial.")

    if not isinstance(name, str) or not name.strip():
        raise ValueError("name must be a nonempty string.")


class Box:
    """Axis-aligned box, specified by its lower corner and full size."""

    def __init__(self, *, origin, size, material, name):
        MomentRT.units.validate_units(value=origin, unit="m", name="origin", shape=(3,))

        MomentRT.units.validate_units(
            value=size, unit="m", name="size", shape=(3,), positive=True
        )

        _validate_region(material=material, name=name)

        self.origin = origin.copy()

        self.size = size.copy()

        self.material = material

        self.name = name

    def contains(self, *, points):
        lower = self.origin.to("m").magnitude

        upper = lower + self.size.to("m").magnitude

        coordinates = points.to("m").magnitude

        return np.all((coordinates >= lower) & (coordinates <= upper), axis=-1)


class Sphere:
    """Sphere specified by center and radius."""

    def __init__(self, *, center, radius, material, name):
        MomentRT.units.validate_units(value=center, unit="m", name="center", shape=(3,))

        MomentRT.units.validate_units(
            value=radius, unit="m", name="radius", shape=(), positive=True
        )

        _validate_region(material=material, name=name)

        self.center = center.copy()

        self.radius = radius.copy()

        self.material = material

        self.name = name

    def contains(self, *, points):
        distance = np.linalg.norm((points - self.center).to("m").magnitude, axis=-1)

        return distance <= self.radius.to("m").magnitude


class Layer:
    """Horizontal layer spanning the sample's full x/y extent.

    depth is the absolute z coordinate of the upper interface. The layer
    occupies depth <= z <= depth + thickness; later regions win at interfaces.
    """

    def __init__(self, *, depth, thickness, material, name):
        MomentRT.units.validate_units(value=depth, unit="m", name="depth", shape=())

        MomentRT.units.validate_units(
            value=thickness, unit="m", name="thickness", shape=(), positive=True
        )

        _validate_region(material=material, name=name)

        self.depth = depth.copy()

        self.thickness = thickness.copy()

        self.material = material

        self.name = name

    def contains(self, *, points):
        z = points.to("m").magnitude[..., 2]

        lower = self.depth.to("m").magnitude

        upper = (self.depth + self.thickness).to("m").magnitude

        return (z >= lower) & (z <= upper)
