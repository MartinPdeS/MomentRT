"""Parameterized optical models for epidermis and dermis."""

from __future__ import annotations
import numpy as np
from pint import Quantity
import MomentRT


class Epidermis(MomentRT.material.Material):
    """Parameterized epidermal material with melanin-dependent absorption.

    Parameters
    ----------
    melanin_fraction : float
        Unitless epidermal melanin fraction from 0 to 1. This is a model
        fraction, not a direct clinical pigmentation measurement.
    wavelength : pint.Quantity
        Vacuum wavelength from 400 to 1000 nm.
    refractive_index : float
        Positive dimensionless refractive index.

    Notes
    -----
    Absorption combines the melanosome power law and a baseline epidermal
    absorber fit. Reduced scattering uses a two-component Rayleigh/Mie power
    law. Both are simplified literature parameterizations; this is not an
    individualized or validated patient model. A Henyey–Greenstein phase
    function with g=0.8 is supplied for moment solvers that require higher
    phase moments.

    References
    ----------
    Jolivot, R. et al. (2013), International Journal of Biomedical Imaging,
    2013, 978289. doi:10.1155/2013/978289.
    Jacques, S. L. (2013), Physics in Medicine & Biology, 58, R37–R61.
    doi:10.1088/0031-9155/58/11/R37.
    """

    def __init__(
        self,
        *,
        melanin_fraction: float,
        wavelength: Quantity,
        refractive_index: float,
    ) -> None:
        """Create an epidermal material from pigmentation and wavelength.

        Parameters
        ----------
        melanin_fraction : float
            Unitless epidermal melanin fraction from 0 to 1.
        wavelength : pint.Quantity
            Vacuum wavelength from 400 to 1000 nm.
        refractive_index : float
            Positive dimensionless refractive index.
        """
        if (
            isinstance(melanin_fraction, (MomentRT.ureg.Quantity, bool))
            or not np.isscalar(melanin_fraction)
            or not np.isfinite(melanin_fraction)
            or not 0 <= melanin_fraction <= 1
        ):
            raise ValueError("melanin_fraction must be a unitless fraction in [0, 1].")
        MomentRT.units.validate_units(
            value=wavelength,
            unit="m",
            name="wavelength",
            shape=(),
            positive=True,
        )
        if not 400 * MomentRT.ureg.nm <= wavelength.to("nm") <= 1000 * MomentRT.ureg.nm:
            raise ValueError("The skin models support wavelengths from 400 to 1000 nm.")

        if (
            isinstance(refractive_index, (MomentRT.ureg.Quantity, bool))
            or not np.isscalar(refractive_index)
            or not np.isfinite(refractive_index)
            or refractive_index <= 0
        ):
            raise ValueError("refractive_index must be a positive, unitless number.")

        wavelength_nm_value = float(wavelength.to("nm").magnitude)

        melanin_per_cm = 6.6e11 * wavelength_nm_value**-3.33
        baseline_per_cm = 0.244 + 85.3 * np.exp(-(wavelength_nm_value - 164) / 66.2)
        absorption: Quantity = (
            (
                melanin_fraction * melanin_per_cm
                + (1 - melanin_fraction) * baseline_per_cm
            )
            / 10
            / MomentRT.ureg.millimeter
        )
        wavelength_ratio = wavelength_nm_value / 500
        reduced_scattering_per_cm = 66.7 * (
            0.29 * wavelength_ratio**-4 + (1 - 0.29) * wavelength_ratio**-0.689
        )
        reduced_scattering = reduced_scattering_per_cm / 10 / MomentRT.ureg.millimeter

        super().__init__(
            absorption=absorption,
            reduced_scattering=reduced_scattering,
            refractive_index=refractive_index,
            anisotropy=0.8,
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        )

        self.melanin_fraction = melanin_fraction
        self.wavelength = wavelength

    def __repr__(self) -> str:
        """Return model parameters and cited optical-property fits."""
        return (
            f"Epidermis(melanin_fraction={self.melanin_fraction:g}, "
            f"wavelength={self.wavelength!r}, "
            f"refractive_index={self.refractive_index!r}, "
            "model='Jolivot et al. 2013 absorption; Jacques 2013 scattering')"
        )


class Dermis(MomentRT.material.Material):
    """Parameterized dermal material with hydration-dependent absorption.

    Parameters
    ----------
    hydration : float
        Unitless water fraction from 0 to 1 in the non-blood dermal component.
    wavelength : pint.Quantity
        Vacuum wavelength from 400 to 1000 nm.
    refractive_index : float
        Positive dimensionless refractive index.

    Notes
    -----
    Absorption mixes pure-water absorption with a baseline dermal absorber
    using the hydration fraction. Reduced scattering uses a two-component
    Rayleigh/Mie power law. Melanin is modeled in the epidermis, rather than
    the dermis, in this simple layered model. A Henyey–Greenstein phase
    function with g=0.9 is supplied for higher moment calculations. These
    literature fits are approximate and omit blood-vessel structure,
    oxygenation, and patient-specific variation.

    References
    ----------
    Jacques, S. L. (2013), Physics in Medicine & Biology, 58, R37–R61.
    doi:10.1088/0031-9155/58/11/R37.
    Hale, G. M. and Querry, M. R. (1973), Applied Optics, 12(3), 555–563.
    doi:10.1364/AO.12.000555.
    """

    def __init__(
        self,
        *,
        hydration: float,
        wavelength: Quantity,
        refractive_index: float,
    ) -> None:
        """Create a dermal material from hydration and wavelength.

        Parameters
        ----------
        hydration : float
            Unitless water fraction from 0 to 1 in non-blood dermal tissue.
        wavelength : pint.Quantity
            Vacuum wavelength from 400 to 1000 nm.
        refractive_index : float
            Positive dimensionless refractive index.
        """
        if (
            isinstance(hydration, (MomentRT.ureg.Quantity, bool))
            or not np.isscalar(hydration)
            or not np.isfinite(hydration)
            or not 0 <= hydration <= 1
        ):
            raise ValueError("hydration must be a unitless fraction in [0, 1].")
        MomentRT.units.validate_units(
            value=wavelength,
            unit="m",
            name="wavelength",
            shape=(),
            positive=True,
        )
        if not 400 * MomentRT.ureg.nm <= wavelength.to("nm") <= 1000 * MomentRT.ureg.nm:
            raise ValueError("The skin models support wavelengths from 400 to 1000 nm.")

        if (
            isinstance(refractive_index, (MomentRT.ureg.Quantity, bool))
            or not np.isscalar(refractive_index)
            or not np.isfinite(refractive_index)
            or refractive_index <= 0
        ):
            raise ValueError("refractive_index must be a positive, unitless number.")

        water_absorption = MomentRT.materials.liquids.Water.absorption_at(
            wavelength=wavelength,
        )
        baseline_absorption = 0.1 / MomentRT.ureg.millimeter
        absorption: Quantity = (
            hydration * water_absorption + (1 - hydration) * baseline_absorption
        )
        wavelength_nm_value = float(wavelength.to("nm").magnitude)
        wavelength_ratio = wavelength_nm_value / 500
        reduced_scattering_per_cm = 43.6 * (
            0.41 * wavelength_ratio**-4 + (1 - 0.41) * wavelength_ratio**-0.562
        )
        reduced_scattering = reduced_scattering_per_cm / 10 / MomentRT.ureg.millimeter

        super().__init__(
            absorption=absorption,
            reduced_scattering=reduced_scattering,
            refractive_index=refractive_index,
            anisotropy=0.9,
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        )

        self.hydration = hydration
        self.wavelength = wavelength

    def __repr__(self) -> str:
        """Return model parameters and cited optical-property fits."""
        return (
            f"Dermis(hydration={self.hydration:g}, "
            f"wavelength={self.wavelength!r}, "
            f"refractive_index={self.refractive_index!r}, "
            "model='Hale–Querry 1973 water absorption; Jacques 2013 scattering')"
        )
