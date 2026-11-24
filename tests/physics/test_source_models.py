"""Tests of sampled analytical source spectra."""

import numpy as np
import pytest
import MomentRT
from MomentRT.models import spectrum


@pytest.mark.parametrize("model", [spectrum.Gaussian, spectrum.Lorentzian])
def test_spectrum_has_requested_power_and_fwhm(model) -> None:
    """Both profiles have the specified half-maximum and finite-grid power."""
    unit = MomentRT.ureg
    wavelengths = np.arange(800, 901, 5) * unit.nm
    spectrum = model(
        wavelengths=wavelengths,
        center=850 * unit.nm,
        fwhm=20 * unit.nm,
        optical_power=2 * unit.mW,
    )
    assert isinstance(spectrum, MomentRT.SourceSpectrum)
    assert spectrum.optical_power.to("mW").magnitude == pytest.approx(2)
    density = spectrum.spectral_power.magnitude
    assert density[8] == pytest.approx(density[10] / 2)
    assert density[12] == pytest.approx(density[10] / 2)
    illumination = MomentRT.SurfaceSource(
        face="z_min",
        position=np.zeros(3) * unit.mm,
        radius=1 * unit.mm,
        spectrum=spectrum,
    )
    assert illumination.optical_power.to("mW").magnitude == pytest.approx(2)


@pytest.mark.parametrize("model", [spectrum.Gaussian, spectrum.Lorentzian])
def test_spectrum_rejects_invalid_inputs(model) -> None:
    """Profiles reject unsorted grids, nonpositive widths and unitless power."""
    unit = MomentRT.ureg
    arguments = dict(
        wavelengths=np.array([800, 850, 900]) * unit.nm,
        center=850 * unit.nm,
        fwhm=20 * unit.nm,
        optical_power=1 * unit.mW,
    )
    with pytest.raises(ValueError, match="increasing"):
        model(**{**arguments, "wavelengths": np.array([850, 800]) * unit.nm})
    with pytest.raises(ValueError):
        model(**{**arguments, "fwhm": 0 * unit.nm})
    with pytest.raises((TypeError, ValueError)):
        model(**{**arguments, "optical_power": 1})


def test_gaussian_rejects_unresolved_underflow() -> None:
    """An unsampled extremely narrow Gaussian cannot be normalized."""
    unit = MomentRT.ureg
    with pytest.raises(ValueError, match="does not resolve"):
        spectrum.Gaussian(
            wavelengths=np.array([800, 900]) * unit.nm,
            center=850 * unit.nm,
            fwhm=0.01 * unit.nm,
            optical_power=1 * unit.mW,
        )
