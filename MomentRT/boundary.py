"""Explicit partial-current boundary conditions for steady-state diffusion."""

import numpy as np

from .units import ureg

FACES = ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")


class RobinBoundary:
    """Diffuse internal reflection with an explicitly supplied effective reflectance.

    This is a partial-current coefficient, not normal-incidence Fresnel R.
    No refractive-index-to-reflectance approximation is selected implicitly.
    With R=effective_reflectance, beta=(1-R)/(2(1+R)) and outward net
    current J.n = beta Phi - 2 I_in/(1+R). I_in is incoming diffuse
    partial current already coupled into the sample (W/m^2).
    """

    def __init__(self, *, effective_reflectance):
        value = effective_reflectance

        if (
            isinstance(value, (ureg.Quantity, bool))
            or not np.isscalar(value)
            or not np.isfinite(value)
            or not 0 <= value < 1
        ):
            raise ValueError("effective_reflectance must be a plain number in [0, 1).")

        self.effective_reflectance = value

    @property
    def leakage_coefficient(self):
        return (1 - self.effective_reflectance) / (2 * (1 + self.effective_reflectance))

    @property
    def source_factor(self):
        return 2 / (1 + self.effective_reflectance)


class ReflectingBoundary:
    """Zero normal current; illumination on this face is unsupported."""

    leakage_coefficient = 0.0
