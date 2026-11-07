"""Reusable layered phantom configurations."""

from __future__ import annotations
from collections.abc import Sequence
from pint import Quantity
import MomentRT


def LayeredPhantom(
    *,
    origin: Quantity,
    size: Quantity,
    background: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
    layers: Sequence[
        tuple[
            str,
            Quantity,
            MomentRT.material.Material | MomentRT.material.SpectralMaterial,
        ]
    ],
) -> MomentRT.medium.StructuredMedium:
    """Build a finite sample with an ordered stack of material layers.

    Parameters
    ----------
    origin : pint.Quantity
        Lower Cartesian corner, shape (3,), with length units.
    size : pint.Quantity
        Positive sample extents, shape (3,), with length units.
    background : Material or SpectralMaterial
        Material filling all space not covered by the supplied layers.
    layers : sequence of (str, pint.Quantity, Material or SpectralMaterial)
        Layer name, positive thickness with length units, and material, in
        increasing z order. Layers start at origin.z and must fit the sample.

    Returns
    -------
    StructuredMedium
        Finite sample with the supplied background and layers.

    Notes
    -----
    Positive z points into the sample. Interfaces are not automatically
    aligned to a mesh; numerical convergence must be checked separately.
    """
    sample = MomentRT.medium.StructuredMedium(
        origin=origin,
        size=size,
    )
    sample.add_background(material=background)
    for name, thickness, material in layers:
        sample.add_layer(
            name=name,
            thickness=thickness,
            material=material,
        )
    return sample
