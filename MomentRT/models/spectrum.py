"""Analytical wavelength profiles sampled as source spectra."""

from __future__ import annotations
from typing import Literal
import numpy as np
import TypedUnit
from pint import Quantity
from scipy import integrate
import MomentRT


def _spectrum(
    *,
    wavelengths: Quantity,
    center: Quantity,
    fwhm: Quantity,
    optical_power: Quantity,
    profile: Literal["gaussian", "lorentzian"],
) -> MomentRT.source.SourceSpectrum:
    """Sample and normalize a wavelength profile on the supplied interval."""
    TypedUnit.validate_quantity(
        registry=MomentRT.units.ureg,
        value=wavelengths,
        unit="m",
        name="wavelengths",
        positive=True,
    )
    if wavelengths.ndim != 1 or len(wavelengths) < 2:
        raise ValueError("wavelengths must be a vector of at least two samples.")
    TypedUnit.validate_quantity(
        registry=MomentRT.units.ureg,
        value=center,
        unit="m",
        name="center",
        shape=(),
        positive=True,
    )
    TypedUnit.validate_quantity(
        registry=MomentRT.units.ureg,
        value=fwhm,
        unit="m",
        name="fwhm",
        shape=(),
        positive=True,
    )
    TypedUnit.validate_quantity(
        registry=MomentRT.units.ureg,
        value=optical_power,
        unit="W",
        name="optical_power",
        shape=(),
        positive=True,
    )
    grid = wavelengths.to("m").magnitude
    if np.any(np.diff(grid) <= 0):
        raise ValueError("wavelengths must be strictly increasing.")
    offset = (grid - center.to("m").magnitude) / fwhm.to("m").magnitude
    with np.errstate(over="ignore", under="ignore"):
        if profile == "gaussian":
            density = np.exp(-4 * np.log(2) * offset**2)
        else:
            density = 1 / (1 + 4 * offset**2)
    integral = integrate.trapezoid(density, x=grid)
    if not np.isfinite(integral) or integral <= 0:
        raise ValueError("The wavelength grid does not resolve the source profile.")
    spectral_power = (
        density
        / integral
        * optical_power.to("W").magnitude
        * MomentRT.ureg.W
        / MomentRT.ureg.m
    )
    return MomentRT.source.SourceSpectrum(
        wavelengths=wavelengths,
        spectral_power=spectral_power,
    )


def Gaussian(
    *,
    wavelengths: Quantity,
    center: Quantity,
    fwhm: Quantity,
    optical_power: Quantity,
) -> MomentRT.source.SourceSpectrum:
    """Sample a Gaussian wavelength spectrum on an explicit grid.

    Parameters
    ----------
    wavelengths : pint.Quantity
        Strictly increasing positive wavelengths, shape (n,), n >= 2.
    center : pint.Quantity
        Positive center wavelength, with length units.
    fwhm : pint.Quantity
        Positive full width at half maximum in wavelength units.
    optical_power : pint.Quantity
        Positive integrated power over the supplied interval, with power units.

    Returns
    -------
    SourceSpectrum
        Tabulated power per wavelength, normalized by trapezoidal integration.

    Notes
    -----
    The profile is truncated to the grid endpoints and normalized over that
    finite interval. Grid refinement and interval selection affect the sampled
    shape; normalization does not establish spectral convergence. For a surface
    source, the power is diffuse power already coupled into the sample.
    """
    return _spectrum(
        wavelengths=wavelengths,
        center=center,
        fwhm=fwhm,
        optical_power=optical_power,
        profile="gaussian",
    )


def Lorentzian(
    *,
    wavelengths: Quantity,
    center: Quantity,
    fwhm: Quantity,
    optical_power: Quantity,
) -> MomentRT.source.SourceSpectrum:
    """Sample a Lorentzian wavelength spectrum on an explicit grid.

    Parameters
    ----------
    wavelengths : pint.Quantity
        Strictly increasing positive wavelengths, shape (n,), n >= 2.
    center : pint.Quantity
        Positive center wavelength, with length units.
    fwhm : pint.Quantity
        Positive full width at half maximum in wavelength units.
    optical_power : pint.Quantity
        Positive integrated power over the supplied interval, with power units.

    Returns
    -------
    SourceSpectrum
        Tabulated power per wavelength, normalized by trapezoidal integration.

    Notes
    -----
    The profile is truncated to the grid endpoints and normalized over that
    finite interval. Lorentzian tails can require a wide interval. Grid and
    interval convergence must be assessed separately. Surface-source power is
    diffuse power already coupled into the sample.
    """
    return _spectrum(
        wavelengths=wavelengths,
        center=center,
        fwhm=fwhm,
        optical_power=optical_power,
        profile="lorentzian",
    )
