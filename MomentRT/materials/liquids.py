"""Derived material classes for pure water and Intralipid liquid phantoms."""

from __future__ import annotations
import importlib.resources
import numpy as np
from pint import Quantity
import MomentRT


class Water(MomentRT.material.Material):
    """Pure-water material at one visible wavelength.

    Parameters
    ----------
    wavelength : pint.Quantity
        Vacuum wavelength from 400 to 750 nm.
    refractive_index : float
        Positive dimensionless refractive index at the configured wavelength
        and temperature.

    Notes
    -----
    Absorption is linearly interpolated from Hale and Querry's tabulated data.
    Molecular scattering uses the Buiteveld fit over 340–750 nm and a Rayleigh
    phase function. Pure water is weakly scattering, so MomentRT's diffusion
    approximation is generally unsuitable for this material.

    References
    ----------
    Hale, G. M. and Querry, M. R. (1973), Applied Optics, 12(3), 555–563.
    doi:10.1364/AO.12.000555.
    Buiteveld, H. et al. (1994), Proc. SPIE 2258, 174–183.
    doi:10.1117/12.190060.
    """

    @classmethod
    def absorption_at(cls, *, wavelength: Quantity) -> Quantity:
        """Interpolate Hale–Querry water absorption at a wavelength.

        Parameters
        ----------
        wavelength : pint.Quantity
            Vacuum wavelength from 400 to 1700 nm.

        Returns
        -------
        pint.Quantity
            Absorption coefficient in mm^-1.
        """
        MomentRT.units.validate_units(
            value=wavelength,
            unit="m",
            name="wavelength",
            shape=(),
            positive=True,
        )
        wavelength_nm_value = float(wavelength.to("nm").magnitude)

        data_resource = importlib.resources.files("MomentRT.materials").joinpath(
            "water_absorption.csv",
        )
        with data_resource.open("r", encoding="utf-8") as data_file:
            water_data = np.loadtxt(data_file, delimiter=",", skiprows=1)

        if not (water_data[0, 0] <= wavelength_nm_value <= water_data[-1, 0]):
            raise ValueError(
                "Water absorption data support wavelengths from 400 to 1700 nm."
            )

        absorption_per_cm = np.interp(
            wavelength_nm_value,
            water_data[:, 0],
            water_data[:, 1],
        )

        return (absorption_per_cm / 10) / MomentRT.ureg.millimeter

    def __init__(self, *, wavelength: Quantity, refractive_index: float) -> None:
        """Initialize pure water at the selected wavelength.

        Parameters
        ----------
        wavelength : pint.Quantity
            Vacuum wavelength from 400 to 750 nm.
        refractive_index : float
            Positive dimensionless refractive index.
        """
        MomentRT.units.validate_units(
            value=wavelength,
            unit="m",
            name="wavelength",
            shape=(),
            positive=True,
        )

        if not 400 * MomentRT.ureg.nm <= wavelength.to("nm") <= 750 * MomentRT.ureg.nm:
            raise ValueError("Water supports wavelengths from 400 to 750 nm.")

        wavelength_nm_value = float(wavelength.to("nm").magnitude)

        if (
            isinstance(refractive_index, (MomentRT.ureg.Quantity, bool))
            or not np.isscalar(refractive_index)
            or not np.isfinite(refractive_index)
            or refractive_index <= 0
        ):
            raise ValueError("refractive_index must be a positive, unitless number.")

        reduced_scattering_per_meter = 3.4123664e8 * wavelength_nm_value ** (-4.1537676)

        super().__init__(
            absorption=self.absorption_at(wavelength=wavelength),
            reduced_scattering=reduced_scattering_per_meter
            / 1000
            / MomentRT.ureg.millimeter,
            refractive_index=float(refractive_index),
            phase_function=MomentRT.phase_function.Rayleigh(),
        )

        self.wavelength = wavelength

    def __repr__(self) -> str:
        """Return the configured wavelength and optical model provenance."""
        return (
            f"Water(wavelength={self.wavelength!r}, "
            f"refractive_index={self.refractive_index!r}, "
            "model='Hale–Querry 1973 absorption; Buiteveld et al. 1994 scattering')"
        )


class Intralipid(MomentRT.material.Material):
    """Intralipid-20% phantom using a concentration-dependent scattering model.

    Parameters
    ----------
    concentration : float
        Lipid-particle volume fraction in percent, from 0.227 to 22.7. A value
        of 22.7 represents the undiluted Intralipid-20% stock; 1.0 represents
        a 1% particle volume fraction.
    wavelength : pint.Quantity
        Vacuum wavelength from 400 to 1700 nm and above the model's
        concentration-dependent validity limit.
    refractive_index : float
        Positive dimensionless refractive index of the prepared solution.

    Notes
    -----
    The dependent-scattering equations estimate the scattering coefficient
    and anisotropy from wavelength and particle volume fraction. Absorption is
    approximated by the water fraction of the mixture using Hale–Querry data.
    The model's empirical equations were validated only when
    wavelength_nm >= 560 + 2200 * particle_volume_fraction. This class uses a
    Henyey–Greenstein phase function matched to the modeled anisotropy; it does
    not reproduce the model's full phase function, so SP3/P3 results are an
    approximation. Measure the actual mixture for high-accuracy work.

    References
    ----------
    Attendu, X. et al. (2026), Journal of Biomedical Optics, 31(8), 085003.
    doi:10.1117/1.JBO.31.8.085003.
    """

    def __init__(
        self,
        *,
        concentration: float,
        wavelength: Quantity,
        refractive_index: float,
    ) -> None:
        """Initialize an Intralipid material from concentration and wavelength.

        Parameters
        ----------
        concentration : float
            Particle volume fraction percentage from 0.227 to 22.7.
        wavelength : pint.Quantity
            Vacuum wavelength within the model's concentration-dependent range.
        refractive_index : float
            Positive dimensionless refractive index of the mixture.
        """
        if (
            isinstance(concentration, (MomentRT.ureg.Quantity, bool))
            or not np.isscalar(concentration)
            or not np.isfinite(concentration)
            or not 0.227 <= concentration <= 22.7
        ):
            raise ValueError("concentration must be a percentage in [0.227, 22.7].")

        if (
            isinstance(refractive_index, (MomentRT.ureg.Quantity, bool))
            or not np.isscalar(refractive_index)
            or not np.isfinite(refractive_index)
            or refractive_index <= 0
        ):
            raise ValueError("refractive_index must be a positive, unitless number.")

        MomentRT.units.validate_units(
            value=wavelength,
            unit="m",
            name="wavelength",
            shape=(),
            positive=True,
        )
        wavelength_nm = wavelength.to("nm")
        wavelength_nm_value = float(wavelength_nm.magnitude)
        particle_volume_fraction = concentration / 100
        minimum_wavelength_nm = 560 + 2200 * particle_volume_fraction

        if not (
            minimum_wavelength_nm * MomentRT.ureg.nm
            <= wavelength_nm
            <= 1700 * MomentRT.ureg.nm
        ):
            raise ValueError(
                "Intralipid empirical model supports wavelengths from "
                f"{minimum_wavelength_nm:.1f} to 1700 nm at this concentration."
            )

        packaging = 1.31 + 5.481e-4 * wavelength_nm_value
        independent_scattering = 1.868e9 * wavelength_nm_value ** (-2.59)
        dependent_factor = (1 - particle_volume_fraction) ** (packaging + 1) / (
            1 + particle_volume_fraction * (packaging - 1)
        ) ** (packaging - 1)
        scattering = (
            independent_scattering * particle_volume_fraction / 0.227 * dependent_factor
        )

        a = 1.094
        b = -5.653e-4
        c = 5.3e-3
        f = 0.3516
        h = 0.1933
        d = a * (f - 1) / (b * (h - 1))
        logistic = 1 + np.exp(c * (wavelength_nm_value + d))
        independent_anisotropy = (
            a * ((1 - f) / logistic + f)
            + b * ((1 - h) / logistic + h) * wavelength_nm_value
        )
        concentration_slope = (
            -10
            + 0.05047 * wavelength_nm_value
            - 9.693e-5 * wavelength_nm_value**2
            + 8.449e-8 * wavelength_nm_value**3
            - 3.412e-11 * wavelength_nm_value**4
            + 5.204e-15 * wavelength_nm_value**5
        )
        anisotropy = (
            independent_anisotropy + concentration_slope * particle_volume_fraction
        )

        if not 0 <= anisotropy < 1:
            raise ValueError("The Intralipid model returned inadmissible anisotropy.")

        absorption: Quantity = (1 - particle_volume_fraction) * (
            MomentRT.materials.liquids.Water.absorption_at(wavelength=wavelength)
        )

        super().__init__(
            absorption=absorption,
            reduced_scattering=(scattering * (1 - anisotropy))
            / MomentRT.ureg.millimeter,
            refractive_index=float(refractive_index),
            anisotropy=float(anisotropy),
            phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        )

        self.concentration = float(concentration)
        self.wavelength = wavelength

    def __repr__(self) -> str:
        """Return the configured mixture and scattering model provenance."""
        return (
            f"Intralipid(concentration={self.concentration:g}%, "
            f"wavelength={self.wavelength!r}, "
            f"refractive_index={self.refractive_index!r}, "
            "model='Attendu et al. 2026 dependent-scattering equations; "
            "Henyey–Greenstein phase approximation')"
        )
