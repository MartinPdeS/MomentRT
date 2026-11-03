"""Detector result containers."""

from __future__ import annotations
from pint import Quantity
import MomentRT


class DetectorResult:
    """Finite-aperture diffuse collection and quadrature provenance.

    Parameters
    ----------
    collected_power : pint.Quantity
        Integrated accepted diffuse optical power, with power units.
    sampled_area : pint.Quantity
        Aperture area represented by selected quadrature weights.
        It is not normalized to the geometric area.
    geometric_area : pint.Quantity
        Exact circular aperture area, with area units.
    angular_fraction : float
        Fraction of Lambertian escaping irradiance accepted by the cone.
    wavelength : pint.Quantity
        Positive scalar vacuum wavelength, with length units.
    minimum_escaping_irradiance : pint.Quantity
        Smallest sampled escaping irradiance inside the aperture.
        Negative values flag discretization artifacts and remain unchanged.
    detector : SurfaceDetector
        Snapshot of the aperture and angular collection configuration.
    photocurrent : pint.Quantity or None, optional
        Linear electrical current in amperes when responsivity is configured.
        Default is None.
    """

    def __init__(
        self,
        *,
        collected_power: Quantity,
        sampled_area: Quantity,
        geometric_area: Quantity,
        angular_fraction: float,
        wavelength: Quantity,
        minimum_escaping_irradiance: Quantity,
        detector: MomentRT.detector.SurfaceDetector,
        photocurrent: Quantity | None = None,
    ) -> None:
        """Initialize the DetectorResult configuration.

        Parameters
        ----------
        collected_power : pint.Quantity
            Integrated accepted diffuse optical power, with power units.
        sampled_area : pint.Quantity
            Aperture area represented by selected quadrature weights.
            It is not normalized to the geometric area.
        geometric_area : pint.Quantity
            Exact circular aperture area, with area units.
        angular_fraction : float
            Fraction of Lambertian escaping irradiance accepted by the cone.
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.
        minimum_escaping_irradiance : pint.Quantity
            Smallest sampled escaping irradiance inside the aperture.
            Negative values flag discretization artifacts and remain unchanged.
        detector : SurfaceDetector
            Snapshot of the aperture and angular collection configuration.
        photocurrent : pint.Quantity or None, optional
            Linear electrical current in amperes when responsivity is configured.
            Default is None.
        """

        self.photocurrent = photocurrent

        self.collected_power = collected_power

        self.sampled_area = sampled_area

        self.geometric_area = geometric_area

        self.angular_fraction = angular_fraction

        self.wavelength = wavelength.copy()

        self.minimum_escaping_irradiance = minimum_escaping_irradiance

        self.detector = detector

        self.angular_model = "Lambertian escaping radiance"
