"""Homogeneous material viewed through a finite observation window."""


class HomogeneousMedium:
    def __init__(self, *, material):
        self.material = material

    def plot(self, *, mesh, length_unit="mm", show=True):
        from .plotting import plot_medium

        return plot_medium(
            medium=self,
            mesh=mesh,
            length_unit=length_unit,
            show=show,
        )


class StructuredMedium:
    """Finite sample with an explicit background and ordered material regions.

    origin is the lower x/y/z corner; size is the full sample extent.
    Builders mutate in place and return None. Later structures replace earlier
    materials, including at shared interfaces. Shapes must fit inside the sample.
    This geometry is independent of the numerical discretization.
    """

    def __init__(self, *, origin, size):
        from .units import validate_units

        validate_units(value=origin, unit="m", name="origin", shape=(3,))

        validate_units(value=size, unit="m", name="size", shape=(3,), positive=True)

        self.origin = origin.copy()

        self.size = size.copy()

        self.background = None

        self._structures = ()

    @property
    def structures(self):
        return self._structures

    def add_background(self, *, material):
        """Replace the background while retaining existing structures."""
        from .material import Material

        if not isinstance(material, Material):
            raise TypeError("material must be a Material.")

        self.background = material

    def add_structures(self, *structures):
        """Append named shapes in argument order; validate the whole batch first."""
        import numpy as np

        from .geometry import Box, Layer, Sphere

        lower = self.origin.to("m").magnitude

        upper = (self.origin + self.size).to("m").magnitude

        names = {shape.name for shape in self.structures}

        for shape in structures:
            if not isinstance(shape, (Box, Sphere, Layer)):
                raise TypeError("structures must be Box, Sphere, or Layer objects.")

            if shape.name in names:
                raise ValueError(f"Duplicate structure name: {shape.name}")

            names.add(shape.name)

            if isinstance(shape, Box):
                shape_lower = shape.origin.to("m").magnitude

                shape_upper = (shape.origin + shape.size).to("m").magnitude
            elif isinstance(shape, Sphere):
                shape_lower = (shape.center - shape.radius).to("m").magnitude

                shape_upper = (shape.center + shape.radius).to("m").magnitude
            else:
                shape_lower = lower.copy()

                shape_upper = upper.copy()

                shape_lower[2] = shape.depth.to("m").magnitude

                shape_upper[2] = (shape.depth + shape.thickness).to("m").magnitude

            # Small conversion roundoff at exactly shared faces is harmless.
            tolerance = (
                np.finfo(float).eps * 16 * np.maximum(np.abs(lower), np.abs(upper))
            )

            if np.any(shape_lower < lower - tolerance) or np.any(
                shape_upper > upper + tolerance
            ):
                raise ValueError(
                    f"Structure {shape.name!r} extends outside the sample."
                )

        self._structures = (*self.structures, *structures)

    def material_labels(self, *, points):
        """Return 0 for background and 1..N for ordered structures.

        Points must lie inside the finite sample. The trailing axis is x/y/z.
        """
        import numpy as np

        from .units import validate_units

        validate_units(value=points, unit="m", name="points")

        if points.ndim < 1 or points.shape[-1] != 3:
            raise ValueError("points must have a trailing axis of length three.")

        if self.background is None:
            raise ValueError("Configure an explicit background with add_background().")

        coordinates = points.to("m").magnitude

        lower = self.origin.to("m").magnitude

        upper = (self.origin + self.size).to("m").magnitude

        if np.any(coordinates < lower) or np.any(coordinates > upper):
            raise ValueError("All query points must lie inside the sample.")

        labels = np.zeros(points.shape[:-1], dtype=int)

        for region_number, shape in enumerate(self.structures, start=1):
            labels = np.where(shape.contains(points=points), region_number, labels)

        return labels

    def optical_fields(self, *, mesh):
        """Sample the material map onto mesh nodes, returning unit-bearing fields."""
        import numpy as np

        from .fields import ScalarField
        from .units import ureg

        labels = self.material_labels(points=mesh.points)

        materials = [self.background, *(shape.material for shape in self.structures)]

        fields = {}

        for name, unit in (
            ("absorption", "1/m"),
            ("reduced_scattering", "1/m"),
            ("refractive_index", "dimensionless"),
        ):
            values = np.array(
                [
                    getattr(material, name).to(unit).magnitude
                    if name != "refractive_index"
                    else material.refractive_index
                    for material in materials
                ]
            )

            fields[name] = ScalarField(
                mesh=mesh,
                values=values[labels] * ureg(unit),
                name=name,
                metadata={"association": "node", "overlap": "last structure wins"},
            )

        return fields

    def plot(self, *, mesh=None, length_unit="mm", opacity=0.45, show=True):
        """Plot named region surfaces; mesh is optional and does not define geometry.

        Transparent primitive surfaces show the construction geometry; overlap
        precedence follows insertion order, not optical mixing or CSG subtraction.
        """
        from .plotting import plot_sample

        if self.background is None:
            raise ValueError("Configure an explicit background with add_background().")

        return plot_sample(
            sample=self, length_unit=length_unit, opacity=opacity, show=show
        )

    def add_layer(self, *, name, thickness, material):
        """Append a layer below the deepest existing layer; return None.

        The first layer starts at origin.z. Boxes and spheres do not change
        stacking depth. Explicit Layer objects added through add_structures
        also participate. New layers follow the usual insertion precedence.
        """
        from .geometry import Layer

        depth = self.origin[2]

        for structure in self.structures:
            if isinstance(structure, Layer):
                depth = max(depth, structure.depth + structure.thickness)

        layer = Layer(
            name=name,
            depth=depth,
            thickness=thickness,
            material=material,
        )

        self.add_structures(layer)

    def section(self, *, normal, position, resolution=(200, 200)):
        """Sample a plane at an absolute coordinate, resolving region overlaps.

        resolution gives (horizontal, vertical) cell counts. Cell-center
        samples represent a discrete map, not an exact geometric intersection.
        """
        from .sections import SampleSection

        return SampleSection(
            sample=self,
            normal=normal,
            position=position,
            resolution=resolution,
        )

    def plot_section(
        self,
        *,
        normal,
        position,
        quantity="material",
        resolution=(200, 200),
        length_unit="mm",
        value_unit=None,
        scale="linear",
        show=True,
    ):
        """Plot a resolved cross-section, with depth increasing downwards."""
        section = self.section(
            normal=normal,
            position=position,
            resolution=resolution,
        )

        return section.plot(
            quantity=quantity,
            length_unit=length_unit,
            value_unit=value_unit,
            scale=scale,
            show=show,
        )
