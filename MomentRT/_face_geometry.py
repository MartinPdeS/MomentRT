"""Shared internal geometry for positions and footprints on sample faces."""

from __future__ import annotations
import numpy as np
import TypedUnit
from pint import Quantity
import MomentRT


def face_position(
    *,
    sample: MomentRT.medium.StructuredMedium,
    face: MomentRT.boundary.Face,
    coordinates: Quantity,
) -> Quantity:
    """Convert absolute in-plane coordinates to a sample-face position.

    Parameters
    ----------
    sample : StructuredMedium
        Finite rectangular sample providing face bounds.
    face : Face
        Selected Cartesian sample face.
    coordinates : pint.Quantity
        Absolute in-plane coordinates, shape (2,), with length units, ordered
        by the remaining Cartesian axes.

    Returns
    -------
    pint.Quantity
        Cartesian face position, shape (3,), in the coordinate units.
    """
    if not isinstance(sample, MomentRT.medium.StructuredMedium):
        raise TypeError("sample must be a StructuredMedium.")

    if face not in MomentRT.boundary.FACES:
        raise ValueError(f"face must be one of {MomentRT.boundary.FACES}.")

    TypedUnit.validate_quantity(
        registry=MomentRT.units.ureg,
        value=coordinates,
        unit="m",
        name="coordinates",
        shape=(2,),
    )

    axis = "xyz".index(face[0])

    position = np.asarray(
        sample.origin.to(coordinates.units).magnitude, dtype=float
    ).copy()

    position[[i for i in range(3) if i != axis]] = coordinates.magnitude

    if face.endswith("max"):
        position[axis] += sample.size.to(coordinates.units).magnitude[axis]

    return position * coordinates.units


def validate_face_disk(
    *,
    sample: MomentRT.medium.StructuredMedium,
    face: MomentRT.boundary.Face,
    position: Quantity,
    radius: Quantity,
) -> None:
    """Require a circular face footprint to fit within the sample.

    Parameters
    ----------
    sample : StructuredMedium
        Finite rectangular sample providing the face bounds.
    face : Face
        Selected Cartesian sample face.
    position : pint.Quantity
        Cartesian disk center, shape (3,), with length units.
    radius : pint.Quantity
        Scalar disk radius with length units.

    Raises
    ------
    ValueError
        If the disk extends outside the sample's in-plane face bounds.

    Notes
    -----
    This check concerns in-plane bounds only; callers validate the face,
    center, and radius and check placement on the face separately.
    """
    axis = "xyz".index(face[0])

    in_plane = [index for index in range(3) if index != axis]

    lower = sample.origin.to("m").magnitude

    upper = (sample.origin + sample.size).to("m").magnitude

    center = position.to("m").magnitude

    distance = radius.to("m").magnitude

    tolerance = 16 * np.finfo(float).eps * np.maximum(np.abs(lower), np.abs(upper))

    if np.any(
        center[in_plane] - distance < lower[in_plane] - tolerance[in_plane]
    ) or np.any(center[in_plane] + distance > upper[in_plane] + tolerance[in_plane]):
        raise ValueError(
            disk_bounds_message(
                face=face,
                position=position,
                radius=radius,
                lower=sample.origin,
                upper=sample.origin + sample.size,
                label="circular footprint",
            )
        )


def disk_bounds_message(
    *,
    face: MomentRT.boundary.Face,
    position: Quantity,
    radius: Quantity,
    lower: Quantity,
    upper: Quantity,
    label: str,
) -> str:
    """Describe a footprint and the in-plane bounds it must fit within.

    Parameters
    ----------
    face : Face
        Selected Cartesian sample face.
    position : pint.Quantity
        Cartesian disk center, shape (3,), with length units.
    radius : pint.Quantity
        Scalar disk radius with length units.
    lower, upper : pint.Quantity
        Cartesian bounds, each shape (3,), with length units.
    label : str
        Description of the footprint used in the message.

    Returns
    -------
    str
        Error message containing face, center, radius, and in-plane bounds.
    """
    axes = [index for index in range(3) if index != "xyz".index(face[0])]

    unit = position.units

    bounds = ", ".join(
        f"{'xyz'[axis]}=[{lower[axis].to(unit):~P}, {upper[axis].to(unit):~P}]"
        for axis in axes
    )

    return (
        f"The entire {label} must fit on its sample face {face}: "
        f"center={position:~P}, radius={radius.to(unit):~P}; "
        f"in-plane bounds: {bounds}."
    )
