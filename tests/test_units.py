"""Unit validation at MomentRT's public configuration boundary."""

import numpy as np
import pint
import pytest
import TypedUnit
import MomentRT


def test_validation_accepts_compatible_units_without_modifying_input() -> None:
    """Validation accepts conversion and preserves the caller's quantity."""
    value = np.array([1.0, 2.0, 3.0]) * MomentRT.ureg.mm
    MomentRT.units.validate_units(
        value=value,
        unit="m",
        name="position",
        shape=(3,),
        positive=True,
    )
    np.testing.assert_array_equal(value.magnitude, [1, 2, 3])
    assert value.units == MomentRT.ureg.mm


def test_validation_rejects_foreign_registry_and_wrong_dimensions() -> None:
    """Quantities must use MomentRT's registry and compatible dimensions."""
    foreign = pint.UnitRegistry()
    with pytest.raises(ValueError, match="wavelength"):
        MomentRT.units.validate_units(
            value=850 * foreign.nm,
            unit="m",
            name="wavelength",
            shape=(),
            positive=True,
        )
    with pytest.raises(ValueError, match="compatible"):
        MomentRT.units.validate_units(
            value=1 * MomentRT.ureg.second,
            unit="m",
            name="wavelength",
        )


def test_validation_enforces_shape_finiteness_and_sign() -> None:
    """Shape, finite values, and sign retain their independent checks."""
    unit = MomentRT.ureg.mm
    with pytest.raises(ValueError, match="shape"):
        MomentRT.units.validate_units(
            value=np.array([1]) * unit,
            unit="m",
            name="length",
            shape=(),
        )
    with pytest.raises(ValueError, match="finite"):
        MomentRT.units.validate_units(
            value=np.nan * unit,
            unit="m",
            name="length",
        )
    with pytest.raises(ValueError, match="positive"):
        MomentRT.units.validate_units(
            value=0 * unit,
            unit="m",
            name="length",
            positive=True,
        )
    MomentRT.units.validate_units(
        value=0 * unit,
        unit="m",
        name="length",
        nonnegative=True,
    )


def test_water_uses_typedunit_with_momentrt_registry() -> None:
    """Direct validation retains the water model's registry and scalar checks."""
    wavelength = 600 * MomentRT.ureg.nm
    assert MomentRT.Water.absorption_at(wavelength=wavelength).magnitude > 0
    assert wavelength == 600 * MomentRT.ureg.nm

    with pytest.raises(ValueError, match="registry"):
        MomentRT.Water.absorption_at(wavelength=600 * pint.UnitRegistry().nm)
    with pytest.raises(ValueError, match="shape"):
        MomentRT.Water.absorption_at(wavelength=np.array([600]) * MomentRT.ureg.nm)
    with pytest.raises(ValueError, match="positive"):
        MomentRT.Water.absorption_at(wavelength=0 * MomentRT.ureg.nm)


def test_momentrt_and_typedunit_share_registry() -> None:
    """TypedUnit quantities can enter MomentRT without registry conversion."""
    assert MomentRT.ureg is TypedUnit.ureg
    absorption = 0.02 / TypedUnit.ureg.mm
    material = MomentRT.Material(
        absorption=absorption,
        reduced_scattering=1 / TypedUnit.ureg.mm,
        refractive_index=1.4,
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    )
    assert material.absorption is absorption
    assert material.absorption.units == TypedUnit.ureg.mm**-1


@pytest.mark.parametrize(
    "value",
    [
        True,
        np.bool_(True),
        np.nan,
        np.inf,
        1.4 + 0.1j,
        [1.4],
        1.4 * MomentRT.ureg.dimensionless,
        1.4 * pint.UnitRegistry().dimensionless,
    ],
)
def test_material_unitless_scalar_contract(value: object) -> None:
    """Optical indices reject quantities from every registry before conversion."""
    with pytest.raises(ValueError, match="refractive_index"):
        MomentRT.Material(
            absorption=0.02 / MomentRT.ureg.mm,
            reduced_scattering=1 / MomentRT.ureg.mm,
            refractive_index=value,
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        )


def test_foreign_dimensionless_inputs_cannot_bypass_plain_number_validation() -> None:
    """Boundary and phase-function inputs must be bare even in other registries."""
    foreign = pint.UnitRegistry()
    with pytest.raises(ValueError, match="effective_reflectance"):
        MomentRT.RobinBoundary(effective_reflectance=0.2 * foreign.dimensionless)
    phase = MomentRT.phase_function.HenyeyGreenstein()
    with pytest.raises(ValueError, match="cosine"):
        phase.evaluate(cosine=np.array([0.0]) * foreign.dimensionless, anisotropy=0.5)
    with pytest.raises(ValueError, match="anisotropy"):
        phase.moments(anisotropy=0.5 * foreign.dimensionless)
