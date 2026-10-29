"""Material-bearing geometry for finite samples; z increases into the sample."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from pint import Quantity

import MomentRT


def _validate_region(
    *,
    material: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
    name: str,
) -> None:
    """Validate the material type and nonempty region name.

    Parameters
    ----------
    material : Material or SpectralMaterial
        Optical coefficients assigned to this region.
    name : str
        Nonempty name identifying this object.
    """

    if not isinstance(
        material, (MomentRT.material.Material, MomentRT.material.SpectralMaterial)
    ):
        raise TypeError("material must be a Material or SpectralMaterial.")

    if not isinstance(name, str) or not name.strip():
        raise ValueError("name must be a nonempty string.")


class Box:
    """Axis-aligned box, specified by its lower corner and full size.

    Parameters
    ----------
    origin : pint.Quantity
        Lower Cartesian corner, shape (3,), with length units.
    size : pint.Quantity
        Positive full Cartesian extents, shape (3,), with length units.
    material : Material or SpectralMaterial
        Optical coefficients assigned to this region.
    name : str
        Nonempty name identifying this object.
    """

    def __init__(
        self,
        *,
        origin: Quantity,
        size: Quantity,
        material: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
        name: str,
    ) -> None:
        """Initialize the Box configuration.

        Parameters
        ----------
        origin : pint.Quantity
            Lower Cartesian corner, shape (3,), with length units.
        size : pint.Quantity
            Positive full Cartesian extents, shape (3,), with length units.
        material : Material or SpectralMaterial
            Optical coefficients assigned to this region.
        name : str
            Nonempty name identifying this object.
        """

        MomentRT.units.validate_units(value=origin, unit="m", name="origin", shape=(3,))

        MomentRT.units.validate_units(
            value=size, unit="m", name="size", shape=(3,), positive=True
        )

        _validate_region(material=material, name=name)

        self.origin = origin.copy()

        self.size = size.copy()

        self.material = material

        self.name = name

    def contains(self, *, points: Quantity) -> NDArray[np.bool_]:
        """Test point membership, including the region interfaces.

        Parameters
        ----------
        points : pint.Quantity
            Query coordinates with length units and shape (..., 3).

        Returns
        -------
        ndarray of bool
            Membership mask with shape points.shape[:-1]; interfaces are included.
        """

        lower = self.origin.to("m").magnitude

        upper = lower + self.size.to("m").magnitude

        coordinates = points.to("m").magnitude

        return np.all((coordinates >= lower) & (coordinates <= upper), axis=-1)


class Sphere:
    """Sphere specified by center and radius.

    Parameters
    ----------
    center : pint.Quantity
        Absolute Cartesian center, shape (3,), with length units.
    radius : pint.Quantity
        Positive scalar radius with length units.
    material : Material or SpectralMaterial
        Optical coefficients assigned to this region.
    name : str
        Nonempty name identifying this object.
    """

    def __init__(
        self,
        *,
        center: Quantity,
        radius: Quantity,
        material: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
        name: str,
    ) -> None:
        """Initialize the Sphere configuration.

        Parameters
        ----------
        center : pint.Quantity
            Absolute Cartesian center, shape (3,), with length units.
        radius : pint.Quantity
            Positive scalar radius with length units.
        material : Material or SpectralMaterial
            Optical coefficients assigned to this region.
        name : str
            Nonempty name identifying this object.
        """

        MomentRT.units.validate_units(value=center, unit="m", name="center", shape=(3,))

        MomentRT.units.validate_units(
            value=radius, unit="m", name="radius", shape=(), positive=True
        )

        _validate_region(material=material, name=name)

        self.center = center.copy()

        self.radius = radius.copy()

        self.material = material

        self.name = name

    def contains(self, *, points: Quantity) -> NDArray[np.bool_]:
        """Test point membership, including the region interfaces.

        Parameters
        ----------
        points : pint.Quantity
            Query coordinates with length units and shape (..., 3).

        Returns
        -------
        ndarray of bool
            Membership mask with shape points.shape[:-1]; interfaces are included.
        """

        distance = np.linalg.norm((points - self.center).to("m").magnitude, axis=-1)

        return distance <= self.radius.to("m").magnitude


class Layer:
    """Horizontal layer spanning the sample's full x/y extent.

    depth is the absolute z coordinate of the upper interface. The layer
    occupies depth <= z <= depth + thickness; later regions win at interfaces.

    Parameters
    ----------
    depth : pint.Quantity
        Absolute z coordinate of the upper interface, with length units.
    thickness : pint.Quantity
        Positive scalar layer thickness with length units.
    material : Material or SpectralMaterial
        Optical coefficients assigned to this region.
    name : str
        Nonempty name identifying this object.
    """

    def __init__(
        self,
        *,
        depth: Quantity,
        thickness: Quantity,
        material: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
        name: str,
    ) -> None:
        """Initialize the Layer configuration.

        Parameters
        ----------
        depth : pint.Quantity
            Absolute z coordinate of the upper interface, with length units.
        thickness : pint.Quantity
            Positive scalar layer thickness with length units.
        material : Material or SpectralMaterial
            Optical coefficients assigned to this region.
        name : str
            Nonempty name identifying this object.
        """

        MomentRT.units.validate_units(value=depth, unit="m", name="depth", shape=())

        MomentRT.units.validate_units(
            value=thickness, unit="m", name="thickness", shape=(), positive=True
        )

        _validate_region(material=material, name=name)

        self.depth = depth.copy()

        self.thickness = thickness.copy()

        self.material = material

        self.name = name

    def contains(self, *, points: Quantity) -> NDArray[np.bool_]:
        """Test point membership, including the region interfaces.

        Parameters
        ----------
        points : pint.Quantity
            Query coordinates with length units and shape (..., 3).

        Returns
        -------
        ndarray of bool
            Membership mask with shape points.shape[:-1]; interfaces are included.
        """

        z = points.to("m").magnitude[..., 2]

        lower = self.depth.to("m").magnitude

        upper = (self.depth + self.thickness).to("m").magnitude

        return (z >= lower) & (z <= upper)
