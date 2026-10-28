"""Physical units shared by all MomentRT objects."""

import numpy as np
from pint import UnitRegistry

ureg = UnitRegistry()


def validate_units(*, value, unit, name, shape=None, positive=False, nonnegative=False):
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
