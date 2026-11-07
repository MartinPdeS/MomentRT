"""Isotropic scattering phase function."""

from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike, NDArray
import MomentRT
from MomentRT.phase_function.base import PhaseFunction, _cosines


class Isotropic(PhaseFunction):
    """Equal scalar scattering probability per unit solid angle."""

    def evaluate(
        self, *, cosine: ArrayLike, anisotropy: float | None = None
    ) -> NDArray[np.float64]:
        """Evaluate the isotropic phase density.

        Parameters
        ----------
        cosine : array_like of float
            Plain scattering-angle cosines in [-1, 1].
        anisotropy : float or None, optional
            Must be zero or None for isotropic scattering.

        Returns
        -------
        numpy.ndarray
            Constant phase density 1/(4 pi) per steradian.
        """
        self.moments(anisotropy=anisotropy)

        return np.ones_like(_cosines(cosine=cosine)) / (4 * np.pi)

    def moments(self, *, anisotropy: float | None = None) -> NDArray[np.float64]:
        """Return the isotropic Legendre moments.

        Parameters
        ----------
        anisotropy : float or None, optional
            Must be zero or None for isotropic scattering.

        Returns
        -------
        numpy.ndarray
            Moments [1, 0, 0, 0].
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
            raise ValueError("Isotropic scattering requires anisotropy=0.")

        return np.array([1.0, 0.0, 0.0, 0.0])
