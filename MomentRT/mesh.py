"""Cartesian observation grids; these do not impose physical boundaries."""

from __future__ import annotations

import numpy as np
from pint import Quantity

import MomentRT


class CartesianMesh:
    """Cartesian node grid for observations or finite-element discretization.

    Parameters
    ----------
    shape : tuple of int
        Node counts (nx, ny, nz), each at least two.
    spacing : pint.Quantity
        Positive node separations (dx, dy, dz), with length units.
    origin : pint.Quantity
        Lower Cartesian corner, shape (3,), with length units.

    Examples
    --------
    >>> import numpy as np
    >>> import MomentRT
    >>> mesh = MomentRT.CartesianMesh(
    ...     shape=(11, 11, 11),
    ...     spacing=np.ones(3) * MomentRT.ureg.millimeter,
    ...     origin=np.zeros(3) * MomentRT.ureg.millimeter,
    ... )
    >>> mesh.points.shape
    (11, 11, 11, 3)
    """

    def __init__(
        self, *, shape: tuple[int, int, int], spacing: Quantity, origin: Quantity
    ) -> None:
        """Initialize the CartesianMesh configuration.

        Parameters
        ----------
        shape : tuple of int
            Node counts (nx, ny, nz), each at least two.
        spacing : pint.Quantity
            Positive node separations (dx, dy, dz), with length units.
        origin : pint.Quantity
            Lower Cartesian corner, shape (3,), with length units.
        """

        if len(shape) != 3 or any(
            isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2
            for n in shape
        ):
            raise ValueError("shape must contain three integers of at least two.")

        MomentRT.units.validate_units(
            value=spacing, unit="m", name="spacing", shape=(3,), positive=True
        )

        MomentRT.units.validate_units(value=origin, unit="m", name="origin", shape=(3,))

        self.shape: tuple[int, int, int] = (shape[0], shape[1], shape[2])

        self.spacing = spacing

        self.origin = origin

    @property
    def points(self) -> Quantity:
        """Generate Cartesian node coordinates using ij array ordering.

        Returns
        -------
        pint.Quantity
            Cartesian coordinates in meters, shape ``(*shape, 3)``, using ij ordering.
        """

        axes = [
            self.origin.to("m").magnitude[i]
            + np.arange(n) * self.spacing.to("m").magnitude[i]
            for i, n in enumerate(self.shape)
        ]

        return (
            np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
            * MomentRT.units.ureg.meter
        )

    def to_pyvista(self, *, length_unit: str = "mm") -> MomentRT.plotting.ImageData:
        """Convert Cartesian data to a PyVista image grid.

        Parameters
        ----------
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.

        Returns
        -------
        pyvista.ImageData
            Cartesian VTK grid with coordinates converted to length_unit.
        """

        pv = MomentRT.plotting.require_pyvista()

        return pv.ImageData(
            dimensions=self.shape,
            spacing=self.spacing.to(length_unit).magnitude,
            origin=self.origin.to(length_unit).magnitude,
        )

    @classmethod
    def from_sample(
        cls,
        *,
        sample: MomentRT.medium.StructuredMedium,
        shape: tuple[int, int, int],
    ) -> MomentRT.mesh.CartesianMesh:
        """Create nodes spanning every face of the finite rectangular sample.

        Parameters
        ----------
        sample : StructuredMedium
            Finite rectangular sample with an explicit background.
        shape : tuple of int
            Node counts (nx, ny, nz), each at least two.

        Returns
        -------
        CartesianMesh
            Nodes spanning all sample faces, with spacing size / (shape - 1).
        """
        # Constructor validates counts before spacing is calculated.
        candidate = cls(shape=shape, spacing=sample.size, origin=sample.origin)

        return cls(
            shape=candidate.shape,
            spacing=sample.size / (np.asarray(candidate.shape) - 1),
            origin=sample.origin.copy(),
        )
