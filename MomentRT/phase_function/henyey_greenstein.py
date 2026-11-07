"""Henyey–Greenstein scattering phase function."""

from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike, NDArray
import MomentRT
from MomentRT.phase_function.base import PhaseFunction, _cosines


class HenyeyGreenstein(PhaseFunction):
    """Henyey–Greenstein phase function parameterized by material anisotropy.

    Supply anisotropy on the material, including a wavelength-dependent table
    for SpectralMaterial. This phase-function object contains no extra value.
    """

    def evaluate(
        self, *, cosine: ArrayLike, anisotropy: float | None = None
    ) -> NDArray[np.float64]:
        """Evaluate the Henyey–Greenstein density.

        Parameters
        ----------
        cosine : array_like of float
            Plain scattering-angle cosines in [-1, 1].
        anisotropy : float or None, optional
            Required mean scattering cosine g in (-1, 1).

        Returns
        -------
        numpy.ndarray
            Phase density per steradian, normalized over 4 pi.
        """
        if anisotropy is None:
            raise ValueError("HenyeyGreenstein requires material anisotropy.")

        g = float(
            MomentRT.material.Material._anisotropy_values(value=anisotropy, shape=())
        )

        mu = _cosines(cosine=cosine)

        return (1 - g**2) / (4 * np.pi * ((1 - g) ** 2 + 2 * g * (1 - mu)) ** 1.5)

    def moments(self, *, anisotropy: float | None = None) -> NDArray[np.float64]:
        """Return the Henyey–Greenstein Legendre moments.

        Parameters
        ----------
        anisotropy : float or None, optional
            Required mean scattering cosine g in (-1, 1).

        Returns
        -------
        numpy.ndarray
            Moments [1, g, g², g³].
        """
        if anisotropy is None:
            raise ValueError("HenyeyGreenstein requires material anisotropy.")

        g = float(
            MomentRT.material.Material._anisotropy_values(value=anisotropy, shape=())
        )

        return np.array([1, g, g**2, g**3])
