"""Scalar scattering phase functions and their normalized Legendre moments."""

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

        g = float(MomentRT.material._anisotropy_values(value=anisotropy, shape=()))

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

        g = float(MomentRT.material._anisotropy_values(value=anisotropy, shape=()))

        return np.array([1, g, g**2, g**3])


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
            and float(MomentRT.material._anisotropy_values(value=anisotropy, shape=()))
            != 0
        ):
            raise ValueError("Isotropic scattering requires anisotropy=0.")

        return np.array([1.0, 0.0, 0.0, 0.0])


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
            and float(MomentRT.material._anisotropy_values(value=anisotropy, shape=()))
            != 0
        ):
            raise ValueError("Rayleigh scattering requires anisotropy=0.")

        return np.array([1.0, 0.0, 0.1, 0.0])


class Tabulated(PhaseFunction):
    """Piecewise-linear scalar phase function over scattering-angle cosines.

    Parameters
    ----------
    cosine : array_like of float
        Strictly increasing plain cosines covering exactly [-1, 1].
    values : array_like of float
        Nonnegative relative phase densities on the same grid. The supplied
        shape is normalized over the sphere; an overall scale is immaterial.

    Notes
    -----
    Values are interpolated linearly in cosine, not in angle. The phase
    function is fixed across wavelengths. The first three moments are
    integrated exactly for each piecewise-linear interval.
    """

    def __init__(self, *, cosine: ArrayLike, values: ArrayLike) -> None:
        """Construct and normalize a tabulated phase function.

        Parameters
        ----------
        cosine : array_like of float
            Strictly increasing plain cosines spanning [-1, 1].
        values : array_like of float
            Nonnegative relative density values on the same grid.
        """
        grid = _cosines(cosine=cosine)

        if isinstance(values, MomentRT.units.ureg.Quantity):
            raise TypeError("values must be plain relative densities.")

        density = np.asarray(values)

        if (
            grid.ndim != 1
            or len(grid) < 2
            or grid[0] != -1
            or grid[-1] != 1
            or np.any(np.diff(grid) <= 0)
        ):
            raise ValueError("cosine must strictly increase from -1 to 1.")

        if (
            density.dtype.kind not in "iuf"
            or density.shape != grid.shape
            or not np.all(np.isfinite(density))
            or np.any(density < 0)
        ):
            raise ValueError(
                "values must be finite nonnegative numbers matching cosine."
            )

        integral = float(np.sum(np.diff(grid) * (density[:-1] + density[1:]) / 2))

        if integral <= 0:
            raise ValueError("values must have a positive integral.")

        self.cosine = grid.copy()

        self.values = density.astype(float) / (2 * np.pi * integral)

        self.cosine.flags.writeable = False

        self.values.flags.writeable = False

    def evaluate(
        self, *, cosine: ArrayLike, anisotropy: float | None = None
    ) -> NDArray[np.float64]:
        """Interpolate the normalized phase density in cosine.

        Parameters
        ----------
        cosine : array_like of float
            Plain scattering-angle cosines in [-1, 1].
        anisotropy : float or None, optional
            If supplied, must match the table's first moment.

        Returns
        -------
        numpy.ndarray
            Interpolated phase density per steradian.
        """
        if anisotropy is not None:
            self.moments(anisotropy=anisotropy)

        return np.asarray(np.interp(_cosines(cosine=cosine), self.cosine, self.values))

    def moments(self, *, anisotropy: float | None = None) -> NDArray[np.float64]:
        """Integrate the tabulated phase function against Legendre polynomials.

        Parameters
        ----------
        anisotropy : float or None, optional
            If supplied, must match the table's first moment.

        Returns
        -------
        numpy.ndarray
            Normalized moments of degrees zero through three.
        """
        nodes, weights = np.polynomial.legendre.leggauss(3)

        lower, upper = self.cosine[:-1], self.cosine[1:]

        mu = (lower[:, None] + upper[:, None]) / 2 + (upper - lower)[
            :, None
        ] * nodes / 2

        density = np.interp(mu, self.cosine, self.values)

        weighted = 2 * np.pi * density * (upper - lower)[:, None] * weights / 2

        basis = np.stack(
            (np.ones_like(mu), mu, (3 * mu**2 - 1) / 2, (5 * mu**3 - 3 * mu) / 2)
        )

        moments = np.einsum("lij,ij->l", basis, weighted)

        if anisotropy is not None:
            supplied = float(
                MomentRT.material._anisotropy_values(value=anisotropy, shape=())
            )

            if not np.isclose(supplied, moments[1], atol=1e-10, rtol=1e-10):
                raise ValueError(
                    "anisotropy must match the tabulated phase function's first moment."
                )

        return moments
