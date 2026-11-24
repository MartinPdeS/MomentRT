"""Tests for named sample configurations."""

import numpy as np
import pytest
import MomentRT
from MomentRT.models import phantom, skin


def test_layered_phantom_stacks_layers_and_retains_background() -> None:
    """Named materials stack from the absolute origin with uncovered background."""
    unit = MomentRT.ureg
    material = MomentRT.Material(
        absorption=0.1 / unit.mm,
        reduced_scattering=1 / unit.mm,
        refractive_index=1.4,
        phase_function=MomentRT.Isotropic(),
    )
    background = MomentRT.Material(
        absorption=0.2 / unit.mm,
        reduced_scattering=2 / unit.mm,
        refractive_index=1.4,
        phase_function=MomentRT.Isotropic(),
    )
    sample = phantom.LayeredPhantom(
        origin=(0, 0, 2) * unit.mm,
        size=(2, 2, 4) * unit.mm,
        background=background,
        layers=(("Top", 0.5 * unit.mm, material), ("Middle", 1 * unit.mm, background)),
    )
    assert sample.background is background
    assert [layer.name for layer in sample.structures] == ["Top", "Middle"]
    np.testing.assert_allclose(
        [layer.depth.to("mm").magnitude for layer in sample.structures],
        [2, 2.5],
    )
    assert sample.structures[0].material is material
    with pytest.raises(ValueError):
        phantom.LayeredPhantom(
            origin=(0, 0, 0) * unit.mm,
            size=(1, 1, 1) * unit.mm,
            background=background,
            layers=(("Too deep", 2 * unit.mm, material),),
        )


def test_two_layer_skin_uses_explicit_materials_and_depth() -> None:
    """Two-layer geometry fills the sample and rejects a missing dermis."""
    unit = MomentRT.ureg
    epidermis = MomentRT.Epidermis(
        melanin_fraction=0.1,
        wavelength=850 * unit.nm,
        refractive_index=1.4,
    )
    dermis = MomentRT.Dermis(
        hydration=0.7,
        wavelength=850 * unit.nm,
        refractive_index=1.4,
    )
    arguments = dict(
        origin=(0, 0, 2) * unit.mm,
        size=(2, 2, 4) * unit.mm,
        epidermis=epidermis,
        dermis=dermis,
    )
    sample = skin.TwoLayerSkin(
        **arguments,
        epidermis_thickness=200 * unit.micrometer,
    )
    assert [layer.name for layer in sample.structures] == ["Epidermis", "Dermis"]
    assert sample.structures[0].material is epidermis
    assert sample.structures[1].material is dermis
    assert sample.background is dermis
    assert sample.structures[1].thickness.to("mm").magnitude == pytest.approx(3.8)
    with pytest.raises(ValueError, match="less than"):
        skin.TwoLayerSkin(**arguments, epidermis_thickness=4 * unit.mm)


def test_three_layer_skin_constructs_the_documented_sample() -> None:
    """The named model retains the example's geometry and optical values."""
    unit = MomentRT.ureg
    sample = skin.ShimojoThreeLayer850nm(
        origin=(-0.6, -0.6, 0) * unit.mm,
        size=(1.2, 1.2, 4) * unit.mm,
    )
    assert isinstance(sample, MomentRT.StructuredMedium)
    assert [layer.name for layer in sample.structures] == [
        "Epidermis",
        "Dermis",
        "Hypodermis (subcutaneous fat)",
    ]
    np.testing.assert_allclose(
        [layer.thickness.to("mm").magnitude for layer in sample.structures],
        [0.2, 1.5, 2.3],
    )
    np.testing.assert_allclose(
        [layer.material.absorption.to("1/mm").magnitude for layer in sample.structures],
        [0.329, 0.055, 0.04],
    )
    assert sample.background is sample.structures[-1].material


def test_three_layer_skin_clips_and_builds_independent_samples() -> None:
    """Thin samples clip layers and each call creates fresh materials."""
    unit = MomentRT.ureg
    first = skin.ShimojoThreeLayer850nm(
        origin=(0, 0, 3) * unit.mm,
        size=(1, 1, 0.1) * unit.mm,
    )
    second = skin.ShimojoThreeLayer850nm(
        origin=(0, 0, 3) * unit.mm,
        size=(1, 1, 0.1) * unit.mm,
    )
    assert len(first.structures) == 1
    assert first.structures[0].depth.to("mm").magnitude == pytest.approx(3)
    assert first.structures[0].thickness.to("mm").magnitude == pytest.approx(0.1)
    assert first.background is not second.background


def test_three_layer_skin_requires_unit_bearing_geometry() -> None:
    """The model retains the sample's quantity validation."""
    with pytest.raises((TypeError, ValueError)):
        skin.ShimojoThreeLayer850nm(origin=(0, 0, 0), size=(1, 1, 1))
