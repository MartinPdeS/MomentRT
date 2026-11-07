"""Unpolarized Rayleigh scattering phase function."""

from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike, NDArray
import MomentRT
from MomentRT.phase_function.base import PhaseFunction, _cosines


class Rayleigh(PhaseFunction):
    """Unpolarized scalar Rayleigh scattering; polarization is not modeled."""

    def evaluate(
        self, *, cosine: ArrayLike, anisotropy: float | None = None
    ) -> NDArray[np.float64]:
        """Evaluate the unpolarized Rayleigh phase density.

        Parameters
        ----------
        cosine : array_like of float
            Plain scattering-angle cosines in [-1, 1].
        anisotropy : float or None, optional
            Must be zero or None for Rayleigh scattering.

        Returns
        -------
        numpy.ndarray
            Phase density 3(1 + cosine²)/(16 pi) per steradian.
        """
        self.moments(anisotropy=anisotropy)

        mu = _cosines(cosine=cosine)

        return 3 * (1 + mu**2) / (16 * np.pi)

    def moments(self, *, anisotropy: float | None = None) -> NDArray[np.float64]:
        """Return the unpolarized Rayleigh Legendre moments.

        Parameters
        ----------
        anisotropy : float or None, optional
            Must be zero or None for Rayleigh scattering.

        Returns
        -------
        numpy.ndarray
            Moments [1, 0, 0.1, 0].
        """
        if (
            anisotropy is not None
            and float(
                MomentRT.material.Material._anisotropy_values(
                    value=anisotropy,
                    shape=(),
                )
            )
            != 0
        ):
            raise ValueError("Rayleigh scattering requires anisotropy=0.")

        return np.array([1.0, 0.0, 0.1, 0.0])
