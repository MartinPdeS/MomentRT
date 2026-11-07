"""Reusable skin sample configurations."""

from __future__ import annotations
import TypedUnit
from pint import Quantity
import MomentRT


def TwoLayerSkin(
    *,
    origin: Quantity,
    size: Quantity,
    epidermis: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
    dermis: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
    epidermis_thickness: Quantity,
) -> MomentRT.medium.StructuredMedium:
    """Build an epidermis over a dermis filling the remaining sample depth.

    Parameters
    ----------
    origin : pint.Quantity
        Lower Cartesian corner, shape (3,), with length units.
    size : pint.Quantity
        Positive sample extents, shape (3,), with length units.
    epidermis : Material or SpectralMaterial
        Optical properties of the upper layer.
    dermis : Material or SpectralMaterial
        Optical properties of the lower layer and background.
    epidermis_thickness : pint.Quantity
        Positive thickness with length units, strictly less than sample depth.

    Returns
    -------
    StructuredMedium
        Finite sample with epidermis and dermis in increasing z order.

    Notes
    -----
    This geometry does not prescribe tissue properties or establish their
    accuracy. Wavelength coverage, interface compatibility and mesh
    convergence depend on the selected materials and solver.
    """
    TypedUnit.validate_quantity(
        registry=MomentRT.units.ureg,
        value=size,
        unit="m",
        name="size",
        shape=(3,),
        positive=True,
    )
    TypedUnit.validate_quantity(
        registry=MomentRT.units.ureg,
        value=epidermis_thickness,
        unit="m",
        name="epidermis_thickness",
        shape=(),
        positive=True,
    )
    if epidermis_thickness >= size[2]:
        raise ValueError("epidermis_thickness must be less than the sample depth.")
    return MomentRT.models.phantom.LayeredPhantom(
        origin=origin,
        size=size,
        background=dermis,
        layers=(
            ("Epidermis", epidermis_thickness, epidermis),
            ("Dermis", size[2] - epidermis_thickness, dermis),
        ),
    )


def ShimojoThreeLayer850nm(
    *,
    origin: Quantity,
    size: Quantity,
) -> MomentRT.medium.StructuredMedium:
    """Build the illustrative 850 nm three-layer skin sample.

    Parameters
    ----------
    origin : pint.Quantity
        Lower Cartesian corner, shape (3,), with length units.
    size : pint.Quantity
        Positive sample extents, shape (3,), with length units.

    Returns
    -------
    StructuredMedium
        Finite sample with epidermis, dermis and hypodermis in increasing
        depth order. The hypodermis is also the background material.

    Notes
    -----
    The fixed absorption and reduced scattering coefficients in inverse mm
    are respectively (0.329, 3.75), (0.055, 1.85), and (0.04, 1.33).
    They are approximate values at 850 nm interpolated from one ex vivo skin
    dataset; they are not wavelength-dependent or universal skin properties.
    All layers use refractive index 1.4 and Henyey–Greenstein anisotropy 0.9.
    The requested layer thicknesses are 0.20, 1.5 and 2.35 mm. Layers are
    clipped to the sample depth; beyond their total depth, the hypodermis
    background continues. No wavelength compatibility check is performed.

    References
    ----------
    Shimojo, Y. et al. (2020), Journal of Biomedical Optics, 25, 045002.
    doi:10.1117/1.JBO.25.4.045002.
    """
    unit = MomentRT.units.ureg
    epidermis = MomentRT.material.Material(
        absorption=0.329 / unit.millimeter,
        reduced_scattering=3.75 / unit.millimeter,
        refractive_index=1.4,
        anisotropy=0.9,
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )
    dermis = MomentRT.material.Material(
        absorption=0.055 / unit.millimeter,
        reduced_scattering=1.85 / unit.millimeter,
        refractive_index=1.4,
        anisotropy=0.9,
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )
    hypodermis = MomentRT.material.Material(
        absorption=0.04 / unit.millimeter,
        reduced_scattering=1.33 / unit.millimeter,
        refractive_index=1.4,
        anisotropy=0.9,
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )
    sample = MomentRT.medium.StructuredMedium(origin=origin, size=size)
    sample.add_background(material=hypodermis)
    for name, thickness, material in (
        ("Epidermis", 0.20, epidermis),
        ("Dermis", 1.5, dermis),
        ("Hypodermis (subcutaneous fat)", 2.35, hypodermis),
    ):
        sample.add_layer(
            name=name,
            thickness=thickness * unit.millimeter,
            material=material,
            clip_to_sample=True,
        )
    return sample
