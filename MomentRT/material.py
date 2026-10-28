"""Homogeneous optical properties at the configured source wavelength."""

import numpy as np

import MomentRT


class Material:
    def __init__(self, *, absorption, reduced_scattering, refractive_index):
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

    def at(self, *, wavelength):
        """Return this wavelength-independent material at a positive wavelength."""
        MomentRT.units.validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        return self


class SpectralMaterial:
    """Tabulated optical coefficients with linear interpolation and no extrapolation.

    wavelengths must increase strictly. Coefficient arrays have the same shape.
    refractive_index is a positive plain scalar or an array on the same grid.
    """

    def __init__(
        self, *, wavelengths, absorption, reduced_scattering, refractive_index
    ):
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

    def at(self, *, wavelength):
        """Evaluate a monochromatic Material without modifying the table."""
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

        return Material(
            absorption=absorption,
            reduced_scattering=reduced_scattering,
            refractive_index=refractive_index,
        )
