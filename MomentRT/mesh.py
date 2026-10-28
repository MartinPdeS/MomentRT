"""Cartesian observation grids; these do not impose physical boundaries."""

import numpy as np

from .units import ureg, validate_units


class CartesianMesh:
    def __init__(self, *, shape, spacing, origin):
        if len(shape) != 3 or any(
            isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2
            for n in shape
        ):
            raise ValueError("shape must contain three integers of at least two.")

        validate_units(
            value=spacing, unit="m", name="spacing", shape=(3,), positive=True
        )

        validate_units(value=origin, unit="m", name="origin", shape=(3,))

        self.shape = tuple(shape)

        self.spacing = spacing

        self.origin = origin

    @property
    def points(self):
        axes = [
            self.origin.to("m").magnitude[i]
            + np.arange(n) * self.spacing.to("m").magnitude[i]
            for i, n in enumerate(self.shape)
        ]

        return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1) * ureg.meter

    def to_pyvista(self, *, length_unit="mm"):
        from .plotting import require_pyvista

        pv = require_pyvista()

        return pv.ImageData(
            dimensions=self.shape,
            spacing=self.spacing.to(length_unit).magnitude,
            origin=self.origin.to(length_unit).magnitude,
        )

    @classmethod
    def from_sample(cls, *, sample, shape):
        """Create nodes spanning every face of the finite rectangular sample."""
        # Constructor validates counts before spacing is calculated.
        candidate = cls(shape=shape, spacing=sample.size, origin=sample.origin)

        return cls(
            shape=candidate.shape,
            spacing=sample.size / (np.asarray(candidate.shape) - 1),
            origin=sample.origin.copy(),
        )
