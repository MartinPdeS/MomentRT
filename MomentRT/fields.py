"""Unit-bearing nodal scalar fields with mesh-aware visualization."""

import numpy as np

from .units import ureg


class ScalarField:
    def __init__(self, *, mesh, values, name, metadata=None):
        if not isinstance(values, ureg.Quantity):
            raise TypeError("values must carry units.")

        if values.shape != mesh.shape or not np.all(np.isfinite(values.magnitude)):
            raise ValueError("Field values must be finite and match the mesh shape.")

        self.mesh = mesh

        self.values = values

        self.name = name

        self.metadata = dict(metadata or {})

    def to_pyvista(self, *, length_unit="mm"):
        grid = self.mesh.to_pyvista(length_unit=length_unit)

        grid.point_data[self.name] = self.values.magnitude.ravel(order="F")

        grid.field_data["value_units"] = [str(self.values.units)]

        grid.field_data["coordinate_units"] = [length_unit]

        grid.field_data["association"] = ["point"]

        for key, value in self.metadata.items():
            grid.field_data[key] = [str(value)]

        return grid

    def export(self, *, path, length_unit="mm"):
        from pathlib import Path

        if Path(path).suffix != ".vti":
            raise ValueError("Cartesian fields must be exported with a .vti extension.")

        grid = self.to_pyvista(length_unit=length_unit)

        grid.save(path)

    def plot(
        self,
        *,
        view="volume",
        scale="linear",
        length_unit="mm",
        opacity=None,
        isosurface_count=None,
        show=True,
    ):
        from .plotting import plot_field

        return plot_field(
            field=self,
            view=view,
            scale=scale,
            length_unit=length_unit,
            opacity=opacity,
            isosurface_count=isosurface_count,
            show=show,
        )


class BoundaryField:
    """Snapshot of escaping irradiance at boundary triangle quadrature points.

    positions, area_weights, and escaping_irradiance share a leading shape.
    Escaping irradiance includes the incident-field correction on a source face.
    It is not an angle-resolved radiance.
    """

    def __init__(
        self,
        *,
        face,
        positions,
        area_weights,
        escaping_irradiance,
        sample_origin,
        sample_size,
        wavelength,
    ):
        from .boundary import FACES
        from .units import validate_units

        if face not in FACES:
            raise ValueError(f"face must be one of {FACES}.")

        validate_units(value=positions, unit="m", name="positions")

        if positions.ndim < 2 or positions.shape[-1] != 3:
            raise ValueError("positions must have a trailing axis of length three.")

        validate_units(
            value=area_weights,
            unit="m^2",
            name="area_weights",
            shape=positions.shape[:-1],
            positive=True,
        )

        validate_units(
            value=escaping_irradiance,
            unit="W/m^2",
            name="escaping_irradiance",
            shape=positions.shape[:-1],
        )

        validate_units(value=sample_origin, unit="m", name="sample_origin", shape=(3,))

        validate_units(
            value=sample_size, unit="m", name="sample_size", shape=(3,), positive=True
        )

        validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        self.face = face

        self.positions = positions.copy()

        self.area_weights = area_weights.copy()

        self.escaping_irradiance = escaping_irradiance.copy()

        self.sample_origin = sample_origin.copy()

        self.sample_size = sample_size.copy()

        self.wavelength = wavelength.copy()

    @property
    def escaping_power(self):
        return (self.escaping_irradiance * self.area_weights).sum()


class VectorField:
    """Element-centered vectors, without smoothing across material interfaces."""

    def __init__(self, *, positions, values, name):
        from .units import validate_units

        validate_units(value=positions, unit="m", name="positions")

        if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) == 0:
            raise ValueError("positions must have shape (N, 3).")

        if (
            not isinstance(values, ureg.Quantity)
            or values.shape != positions.shape
            or not np.all(np.isfinite(values.magnitude))
        ):
            raise ValueError(
                "values must be finite unit-bearing vectors matching positions."
            )

        self.positions = positions

        self.values = values

        self.name = name

    def plot(self, *, max_vectors=500, length_unit="mm", show=True):
        from .plotting import make_figure, require_plotly

        if (
            isinstance(max_vectors, bool)
            or not isinstance(max_vectors, (int, np.integer))
            or max_vectors < 1
        ):
            raise ValueError("max_vectors must be a positive integer.")

        go = require_plotly()

        selection = np.linspace(
            0, len(self.values) - 1, min(max_vectors, len(self.values)), dtype=int
        )

        points = self.positions.to(length_unit).magnitude[selection]

        vectors = self.values.magnitude[selection]

        figure = make_figure(length_unit=length_unit)

        figure.add_trace(
            go.Cone(
                x=points[:, 0],
                y=points[:, 1],
                z=points[:, 2],
                u=vectors[:, 0],
                v=vectors[:, 1],
                w=vectors[:, 2],
                colorscale="Viridis",
                sizemode="scaled",
                sizeref=0.25,
                colorbar={"title": {"text": f"{self.name} [{self.values.units:~}]"}},
            )
        )

        figure.update_layout(
            title=f"{self.name} (element centers; relative arrow lengths)"
        )

        if show:
            figure.show()

        return figure
