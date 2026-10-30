"""Homogeneous optical properties at the configured source wavelength."""

from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike
from pint import Quantity
import MomentRT


class Material:
    """Wavelength-independent optical coefficients for moment models.

    Parameters
    ----------
    absorption : pint.Quantity
        Nonnegative absorption coefficient, with inverse-length units.
    reduced_scattering : pint.Quantity
        Positive reduced scattering coefficient, with inverse-length units.
    refractive_index : float
        Positive dimensionless refractive index, supplied without units.
    anisotropy : float or None, optional
        Mean scattering cosine in (-1, 1), supplied without units. Required
        by Henyey–Greenstein SP3/P3 calculations. Phase functions with their
        own first moment can omit it; P1 does not use it.
    phase_function : PhaseFunction
        Explicit scalar scattering phase function. SP3 and P3 use its moments;
        P1 only needs reduced scattering.

    Notes
    -----
    The diffusion length coefficient is D = 1 / (3 * (absorption +
    reduced_scattering)). No wavelength dependence is inferred.

    Examples
    --------
    >>> import MomentRT
    >>> ureg = MomentRT.ureg
    >>> tissue = MomentRT.Material(
    ...     absorption=0.1 / ureg.millimeter,
    ...     reduced_scattering=1.0 / ureg.millimeter,
    ...     refractive_index=1.4,
    ...     phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    ... )
    """

    @staticmethod
    def _anisotropy_values(
        *,
        value: float | ArrayLike,
        shape: tuple[int, ...],
    ) -> np.ndarray:
        """Validate scalar or tabulated dimensionless scattering anisotropy.

        Parameters
        ----------
        value : float or array_like
            Plain scalar or values to validate.
        shape : tuple of int
            Required tabulated shape. A scalar is accepted for any shape.

        Returns
        -------
        numpy.ndarray
            Validated anisotropy values.

        Raises
        ------
        ValueError
            If values have units, a boolean or nonnumeric dtype, the wrong
            shape, nonfinite entries, or values outside (-1, 1).
        """
        if isinstance(value, (MomentRT.units.ureg.Quantity, bool, np.bool_)):
            raise ValueError("anisotropy must contain plain numbers in (-1, 1).")

        values = np.asarray(value)

        if (
            values.dtype.kind not in "iuf"
            or values.shape not in ((), shape)
            or not np.all(np.isfinite(values))
            or np.any(values <= -1)
            or np.any(values >= 1)
        ):
            raise ValueError(
                "anisotropy must be scalar or match the wavelength grid, with values in (-1, 1)."
            )

        return values

    def __init__(
        self,
        *,
        absorption: Quantity,
        reduced_scattering: Quantity,
        refractive_index: float,
        anisotropy: float | None = None,
        phase_function: MomentRT.phase_function.PhaseFunction,
    ) -> None:
        """Initialize the Material configuration.

        Parameters
        ----------
        absorption : pint.Quantity
            Nonnegative absorption coefficient, with inverse-length units.
        reduced_scattering : pint.Quantity
            Positive reduced scattering coefficient, with inverse-length units.
        refractive_index : float
            Positive dimensionless refractive index, supplied without units.
        anisotropy : float or None, optional
            Mean scattering cosine in (-1, 1). Required for Henyey–Greenstein
            SP3/P3 calculations; unused by P1.
        phase_function : PhaseFunction
            Explicit scalar scattering phase function for SP3 and P3.
        """

        MomentRT.units.validate_units(
            value=absorption, unit="1/m", name="absorption", shape=(), nonnegative=True
        )

        MomentRT.units.validate_units(
            value=reduced_scattering,
            unit="1/m",
            name="reduced_scattering",
            shape=(),
            positive=True,
        )

        if (
            isinstance(refractive_index, MomentRT.units.ureg.Quantity)
            or not np.isscalar(refractive_index)
            or not np.isfinite(refractive_index)
            or refractive_index <= 0
        ):
            raise ValueError("refractive_index must be a positive, unitless number.")

        self.absorption = absorption

        self.reduced_scattering = reduced_scattering

        self.refractive_index = refractive_index

        if not isinstance(phase_function, MomentRT.phase_function.PhaseFunction):
            raise ValueError("phase_function must be a PhaseFunction object.")

        self.phase_function = phase_function

        self.anisotropy = (
            None
            if anisotropy is None
            else float(self._anisotropy_values(value=anisotropy, shape=()))
        )

        if isinstance(
            phase_function, MomentRT.phase_function.PhaseFunction
        ) and not isinstance(phase_function, MomentRT.phase_function.HenyeyGreenstein):
            phase_function.moments(anisotropy=self.anisotropy)

    def _phase_moments(self) -> np.ndarray:
        """Resolve and validate normalized phase-function moments for transport."""
        phase = self.phase_function

        if not isinstance(phase, MomentRT.phase_function.PhaseFunction):
            raise ValueError("phase_function must be a PhaseFunction object.")

        values = np.asarray(phase.moments(anisotropy=self.anisotropy), dtype=float)

        if (
            values.shape != (4,)
            or not np.all(np.isfinite(values))
            or not np.isclose(values[0], 1, rtol=1e-10, atol=1e-10)
            or np.any(np.abs(values[1:]) > 1 + 1e-10)
            or not -1 < values[1] < 1
        ):
            raise ValueError(
                "phase-function moments must be finite, normalized, and physically admissible with first moment in (-1, 1)."
            )

        first = values[1]

        second = (2 * values[2] + 1) / 3

        third = (2 * values[3] + 3 * first) / 5

        for sign in (-1, 1):
            moment_matrix = np.array(
                [
                    [1 + sign * first, first + sign * second],
                    [first + sign * second, second + sign * third],
                ]
            )

            if np.linalg.eigvalsh(moment_matrix).min() < -1e-10:
                raise ValueError(
                    "phase-function moments cannot represent a nonnegative scattering distribution on [-1, 1]."
                )

        return values

    def _transport_factors(self) -> np.ndarray:
        """Return scattering attenuation factors relative to reduced scattering."""
        moments = self._phase_moments()

        if isinstance(self.phase_function, MomentRT.phase_function.HenyeyGreenstein):
            g = moments[1]

            return np.array([0, 1, 1 + g, 1 + g + g**2])

        return np.array(
            [
                0,
                1,
                (1 - moments[2]) / (1 - moments[1]),
                (1 - moments[3]) / (1 - moments[1]),
            ]
        )

    def at(self, *, wavelength: Quantity) -> MomentRT.material.Material:
        """Return this wavelength-independent material at a positive wavelength.

        Parameters
        ----------
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.

        Returns
        -------
        Material
            Monochromatic material at the requested wavelength.
        """
        MomentRT.units.validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        return self


class SpectralMaterial:
    """Tabulated optical coefficients with linear interpolation and no extrapolation.

    wavelengths must increase strictly. Coefficient arrays have the same shape.
    refractive_index is a positive plain scalar or an array on the same grid.

    Parameters
    ----------
    wavelengths : pint.Quantity
        Strictly increasing one-dimensional wavelength grid with at least two
        positive entries and length units.
    absorption : pint.Quantity
        Nonnegative absorption values with inverse-length units, matching
        the wavelength grid shape.
    reduced_scattering : pint.Quantity
        Positive reduced scattering values with inverse-length units,
        matching the wavelength grid shape.
    refractive_index : float or array_like of float
        Positive plain scalar or vector matching the wavelength grid.
        Dimensionless quantities are not accepted.
    anisotropy : float or array_like of float or None, optional
        Plain scalar or wavelength-grid vector of mean scattering cosines in
        (-1, 1). Interpolated linearly; required for Henyey–Greenstein
        SP3/P3 calculations, unused by P1.
    phase_function : PhaseFunction
        Explicit scattering phase function, fixed across the wavelength table.
    """

    def __init__(
        self,
        *,
        wavelengths: Quantity,
        absorption: Quantity,
        reduced_scattering: Quantity,
        refractive_index: float | ArrayLike,
        anisotropy: float | ArrayLike | None = None,
        phase_function: MomentRT.phase_function.PhaseFunction,
    ) -> None:
        """Initialize the SpectralMaterial configuration.

        Parameters
        ----------
        wavelengths : pint.Quantity
            Strictly increasing one-dimensional wavelength grid with at least two
            positive entries and length units.
        absorption : pint.Quantity
            Nonnegative absorption values with inverse-length units, matching
            the wavelength grid shape.
        reduced_scattering : pint.Quantity
            Positive reduced scattering values with inverse-length units,
            matching the wavelength grid shape.
        refractive_index : float or array_like of float
            Positive plain scalar or vector matching the wavelength grid.
            Dimensionless quantities are not accepted.
        anisotropy : float or array_like of float or None, optional
            Plain scalar or wavelength-grid vector of scattering cosines in
            (-1, 1). Required for Henyey–Greenstein SP3/P3; unused by P1.
        phase_function : PhaseFunction
            Explicit scattering phase function, fixed across wavelengths.
        """

        MomentRT.units.validate_units(
            value=wavelengths, unit="m", name="wavelengths", positive=True
        )

        if wavelengths.ndim != 1 or len(wavelengths) < 2:
            raise ValueError("wavelengths must be a vector with at least two entries.")

        if np.any(np.diff(wavelengths.to("m").magnitude) <= 0):
            raise ValueError("wavelengths must increase strictly.")

        MomentRT.units.validate_units(
            value=absorption,
            unit="1/m",
            name="absorption",
            shape=wavelengths.shape,
            nonnegative=True,
        )

        MomentRT.units.validate_units(
            value=reduced_scattering,
            unit="1/m",
            name="reduced_scattering",
            shape=wavelengths.shape,
            positive=True,
        )

        if isinstance(refractive_index, (MomentRT.units.ureg.Quantity, bool)):
            raise ValueError("refractive_index must contain positive plain numbers.")

        indices = np.asarray(refractive_index)

        if (
            indices.dtype.kind not in "iuf"
            or indices.shape not in ((), wavelengths.shape)
            or not np.all(np.isfinite(indices))
            or np.any(indices <= 0)
        ):
            raise ValueError(
                "refractive_index must be positive, scalar or shaped like wavelengths."
            )

        self.wavelengths = wavelengths.copy()

        self.absorption = absorption.copy()

        self.reduced_scattering = reduced_scattering.copy()

        self.refractive_index = indices.copy()

        if not isinstance(phase_function, MomentRT.phase_function.PhaseFunction):
            raise ValueError("phase_function must be a PhaseFunction object.")

        self.phase_function = phase_function

        self.anisotropy = (
            None
            if anisotropy is None
            else Material._anisotropy_values(
                value=anisotropy,
                shape=wavelengths.shape,
            ).copy()
        )

        if isinstance(
            phase_function, MomentRT.phase_function.PhaseFunction
        ) and not isinstance(phase_function, MomentRT.phase_function.HenyeyGreenstein):
            if self.anisotropy is None:
                phase_function.moments()
            else:
                for value in np.atleast_1d(self.anisotropy):
                    phase_function.moments(anisotropy=float(value))

    def at(self, *, wavelength: Quantity) -> MomentRT.material.Material:
        """Evaluate a monochromatic Material without modifying the table.

        Parameters
        ----------
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.

        Returns
        -------
        Material
            Monochromatic material at the requested wavelength.

        Raises
        ------
        ValueError
            If wavelength is nonpositive, nonscalar, or outside the tabulated range.
        """
        MomentRT.units.validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        grid = self.wavelengths.magnitude

        coordinate = float(wavelength.to(self.wavelengths.units).magnitude)

        tolerance = 8 * np.finfo(float).eps * max(abs(grid[0]), abs(grid[-1]))

        if coordinate < grid[0] - tolerance or coordinate > grid[-1] + tolerance:
            raise ValueError("wavelength lies outside the spectral material table.")

        # Endpoint conversion roundoff is allowed; physical extrapolation is not.
        coordinate = min(max(coordinate, grid[0]), grid[-1])

        absorption = (
            np.interp(coordinate, grid, self.absorption.magnitude)
            * self.absorption.units
        )

        reduced_scattering = (
            np.interp(coordinate, grid, self.reduced_scattering.magnitude)
            * self.reduced_scattering.units
        )

        refractive_index = (
            float(self.refractive_index)
            if self.refractive_index.ndim == 0
            else float(np.interp(coordinate, grid, self.refractive_index))
        )

        anisotropy = (
            None
            if self.anisotropy is None
            else float(self.anisotropy)
            if self.anisotropy.ndim == 0
            else float(np.interp(coordinate, grid, self.anisotropy))
        )

        return Material(
            absorption=absorption,
            reduced_scattering=reduced_scattering,
            refractive_index=refractive_index,
            anisotropy=anisotropy,
            phase_function=self.phase_function,
        )
