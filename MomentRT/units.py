"""Physical units shared by all MomentRT objects."""

from __future__ import annotations

import numpy as np
from pint import Quantity, UnitRegistry

ureg = UnitRegistry()


def validate_units(
    *,
    value: Quantity,
    unit: str,
    name: str,
    shape: tuple[int, ...] | None = None,
    positive: bool = False,
    nonnegative: bool = False,
) -> None:
    """Validate registry, dimensional compatibility, finiteness, shape, and sign.

    Parameters
    ----------
    value : pint.Quantity
        Quantity to validate without changing or replacing it.
    unit : str
        Required compatible physical unit.
    name : str
        Parameter label used in validation error messages.
    shape : tuple of int or None, optional
        Required magnitude shape. None leaves shape unrestricted; ()
        requires a scalar.
        Default is None.
    positive : bool, optional
        Require every magnitude to be strictly positive.
        Default is False.
    nonnegative : bool, optional
        Require every magnitude to be nonnegative.
        Default is False.

    Raises
    ------
    TypeError
        If value is not a quantity from MomentRT.ureg.
    pint.DimensionalityError
        If the supplied units are incompatible with unit.
    ValueError
        If magnitudes are nonfinite or violate the requested shape or sign.
    """

    if not isinstance(value, ureg.Quantity):
        raise TypeError(f"{name} must carry units from MomentRT.ureg.")

    magnitude = np.asarray(value.to(unit).magnitude, dtype=float)

    if shape is not None and magnitude.shape != shape:
        raise ValueError(f"{name} must have shape {shape}.")

    if not np.all(np.isfinite(magnitude)):
        raise ValueError(f"{name} must be finite.")

    if positive and np.any(magnitude <= 0):
        raise ValueError(f"{name} must be positive.")

    if nonnegative and np.any(magnitude < 0):
        raise ValueError(f"{name} must be nonnegative.")
