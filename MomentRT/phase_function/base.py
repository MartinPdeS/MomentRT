"""Common phase-function interface and angle validation."""

from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
from numpy.typing import ArrayLike, NDArray
import MomentRT


def _cosines(*, cosine: ArrayLike) -> NDArray[np.float64]:
    """Validate the dimensionless cosines of scattering angles."""
    if isinstance(cosine, MomentRT.units.ureg.Quantity):
        raise TypeError("cosine must contain plain dimensionless numbers.")

    values = np.asarray(cosine)

    if (
        values.dtype.kind not in "iuf"
        or not np.all(np.isfinite(values))
        or np.any(np.abs(values) > 1)
    ):
        raise ValueError("cosine must contain finite plain numbers in [-1, 1].")

    return values.astype(float)


class PhaseFunction(ABC):
    """Axisymmetric scalar phase function normalized over the unit sphere.

    The normalized moment of degree l is the spherical integral of the
    phase function times the Legendre polynomial P_l. Polarization and
    azimuth-dependent scattering are outside this scalar model.
    """

    @abstractmethod
    def evaluate(
        self, *, cosine: ArrayLike, anisotropy: float | None = None
    ) -> NDArray[np.float64]:
        """Evaluate probability density per unit solid angle.

        Parameters
        ----------
        cosine : array_like of float
            Plain cosines of scattering angles, in [-1, 1].
        anisotropy : float or None, optional
            Material anisotropy when required by the phase function.

        Returns
        -------
        numpy.ndarray
            Scalar phase density per steradian, normalized over 4 pi.
        """
        raise NotImplementedError

    @abstractmethod
    def moments(self, *, anisotropy: float | None = None) -> NDArray[np.float64]:
        """Return the normalized Legendre moments from degree zero to three.

        Parameters
        ----------
        anisotropy : float or None, optional
            Material anisotropy when required by the phase function.

        Returns
        -------
        numpy.ndarray
            Four dimensionless values, with the first equal to one.
        """
        raise NotImplementedError
