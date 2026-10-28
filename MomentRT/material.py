"""Homogeneous optical properties at the configured source wavelength."""

import numpy as np

from .units import ureg, validate_units


class Material:
    def __init__(self, *, absorption, reduced_scattering, refractive_index):
        validate_units(
            value=absorption, unit="1/m", name="absorption", shape=(), nonnegative=True
        )

        validate_units(
            value=reduced_scattering,
            unit="1/m",
            name="reduced_scattering",
            shape=(),
            positive=True,
        )

        if (
            isinstance(refractive_index, ureg.Quantity)
            or not np.isscalar(refractive_index)
            or not np.isfinite(refractive_index)
            or refractive_index <= 0
        ):
            raise ValueError("refractive_index must be a positive, unitless number.")

        self.absorption = absorption

        self.reduced_scattering = reduced_scattering

        self.refractive_index = refractive_index
