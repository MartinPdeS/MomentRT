"""Explicit partial-current boundary conditions for steady-state diffusion."""

from __future__ import annotations
import copy
from collections.abc import Iterator, Mapping
from typing import Literal
import numpy as np
import MomentRT

Face = Literal["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"]

FACES: tuple[Face, ...] = ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")


class RobinBoundary:
    """Diffuse internal reflection with an explicitly supplied effective reflectance.

    This is a partial-current coefficient, not normal-incidence Fresnel R.
    No refractive-index-to-reflectance approximation is selected implicitly.
    With R=effective_reflectance, beta=(1-R)/(2(1+R)) and outward net
    current J.n = beta Phi - 2 I_in/(1+R). I_in is incoming diffuse
    partial current already coupled into the sample (W/m^2).

    Parameters
    ----------
    effective_reflectance : float
        Effective diffuse partial-current reflectance in [0, 1).
        This is not the normal-incidence Fresnel reflectance.
    """

    def __init__(self, *, effective_reflectance: float) -> None:
        """Initialize the RobinBoundary configuration.

        Parameters
        ----------
        effective_reflectance : float
            Effective diffuse partial-current reflectance in [0, 1).
            This is not the normal-incidence Fresnel reflectance.
        """

        value = effective_reflectance

        if (
            isinstance(value, (MomentRT.units.ureg.Quantity, bool))
            or not np.isscalar(value)
            or not np.isfinite(value)
            or not 0 <= value < 1
        ):
            raise ValueError("effective_reflectance must be a plain number in [0, 1).")

        self.effective_reflectance = value

    @property
    def leakage_coefficient(self) -> float:
        """Return the dimensionless Robin leakage coefficient.

        Returns
        -------
        float
            Dimensionless beta = (1 - R) / (2 * (1 + R)).
        """

        return (1 - self.effective_reflectance) / (2 * (1 + self.effective_reflectance))

    @property
    def source_factor(self) -> float:
        """Return the dimensionless incoming-current multiplier.

        Returns
        -------
        float
            Dimensionless incoming-current factor 2 / (1 + R).
        """

        return 2 / (1 + self.effective_reflectance)


class ReflectingBoundary:
    """Zero normal current; illumination on this face is unsupported."""

    leakage_coefficient = 0.0


class Boundaries(Mapping[Face, RobinBoundary | ReflectingBoundary]):
    """A complete set of explicitly chosen conditions on the six sample faces.

    Parameters
    ----------
    conditions : mapping of str to RobinBoundary or ReflectingBoundary
        Exactly one condition for each of the six Cartesian sample faces.
    """

    def __init__(
        self, *, conditions: Mapping[Face, RobinBoundary | ReflectingBoundary]
    ) -> None:
        """Validate and snapshot the supplied face conditions.

        Parameters
        ----------
        conditions : mapping of str to RobinBoundary or ReflectingBoundary
            Exactly one condition per sample face.
        """
        if set(conditions) != set(FACES):
            raise ValueError(f"Configure exactly these boundary faces: {FACES}.")

        if any(
            not isinstance(value, (RobinBoundary, ReflectingBoundary))
            for value in conditions.values()
        ):
            raise TypeError("Each face requires RobinBoundary or ReflectingBoundary.")

        self._conditions = {
            face: copy.deepcopy(condition) for face, condition in conditions.items()
        }

    def __getitem__(self, face: Face) -> RobinBoundary | ReflectingBoundary:
        return self._conditions[face]

    def __iter__(self) -> Iterator[Face]:
        return iter(self._conditions)

    def __len__(self) -> int:
        return len(self._conditions)

    @classmethod
    def uniform(cls, *, condition: RobinBoundary | ReflectingBoundary) -> Boundaries:
        """Apply the supplied explicit condition to every face.

        Parameters
        ----------
        condition : RobinBoundary or ReflectingBoundary
            Condition copied independently to each of the six faces.

        Returns
        -------
        Boundaries
            Complete set of face conditions.
        """
        return cls(conditions={face: condition for face in FACES})

    def with_overrides(
        self, *, conditions: Mapping[Face, RobinBoundary | ReflectingBoundary]
    ) -> Boundaries:
        """Return a copy with explicitly selected face conditions replaced.

        Parameters
        ----------
        conditions : mapping of str to RobinBoundary or ReflectingBoundary
            Replacement conditions keyed by face name.

        Returns
        -------
        Boundaries
            New complete face mapping; this set is not modified.
        """
        unknown = set(conditions) - set(FACES)

        if unknown:
            raise ValueError(f"Unknown boundary faces: {sorted(unknown)}.")

        return Boundaries(conditions={**self._conditions, **conditions})
