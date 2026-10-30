"""Homogeneous material viewed through a finite observation window."""

from __future__ import annotations
import copy
from collections.abc import Iterable
from typing import Literal, cast
import numpy as np
from numpy.typing import NDArray
from pint import Quantity
import MomentRT


class HomogeneousMedium:
    """Uniform optical material in an unbounded physical medium.

    Parameters
    ----------
    material : Material or SpectralMaterial
        Optical coefficients assigned to this region.
    """

    def __init__(
        self,
        *,
        material: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
    ) -> None:
        """Initialize the HomogeneousMedium configuration.

        Parameters
        ----------
        material : Material or SpectralMaterial
            Optical coefficients assigned to this region.
        """

        self.material = material

    def at(self, *, wavelength: Quantity) -> MomentRT.medium.HomogeneousMedium:
        """Resolve optical coefficients at a positive vacuum wavelength.

        Parameters
        ----------
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.

        Returns
        -------
        HomogeneousMedium
            Geometry-preserving monochromatic medium; the original is unchanged.
        """

        return HomogeneousMedium(material=self.material.at(wavelength=wavelength))

    def plot(
        self,
        *,
        mesh: MomentRT.mesh.CartesianMesh,
        length_unit: str = "mm",
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Create an interactive visualization of this object.

        Parameters
        ----------
        mesh : CartesianMesh
            Cartesian nodes used for sampling or finite-element discretization.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable interactive figure, also displayed when show is True.
        """

        return MomentRT.plotting.plot_medium(
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

    Parameters
    ----------
    origin : pint.Quantity
        Lower Cartesian corner, shape (3,), with length units.
    size : pint.Quantity
        Positive full Cartesian extents, shape (3,), with length units.
    """

    def __init__(self, *, origin: Quantity, size: Quantity) -> None:
        """Initialize the StructuredMedium configuration.

        Parameters
        ----------
        origin : pint.Quantity
            Lower Cartesian corner, shape (3,), with length units.
        size : pint.Quantity
            Positive full Cartesian extents, shape (3,), with length units.
        """

        MomentRT.units.validate_units(value=origin, unit="m", name="origin", shape=(3,))

        MomentRT.units.validate_units(
            value=size, unit="m", name="size", shape=(3,), positive=True
        )

        self.origin = origin.copy()

        self.size = size.copy()

        self.background: (
            MomentRT.material.Material | MomentRT.material.SpectralMaterial | None
        ) = None

        self._structures: tuple[
            MomentRT.geometry.Box | MomentRT.geometry.Sphere | MomentRT.geometry.Layer,
            ...,
        ] = ()

    @property
    def structures(
        self,
    ) -> tuple[
        MomentRT.geometry.Box | MomentRT.geometry.Sphere | MomentRT.geometry.Layer, ...
    ]:
        """Return the immutable ordered collection of material-bearing regions.

        Returns
        -------
        tuple of Box, Sphere, or Layer
            Regions in insertion order; later regions replace earlier materials.
        """

        return self._structures

    def add_background(
        self,
        *,
        material: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
    ) -> None:
        """Replace the background while retaining existing structures.

        Parameters
        ----------
        material : Material or SpectralMaterial
            Optical coefficients assigned to this region.
        """

        if not isinstance(
            material, (MomentRT.material.Material, MomentRT.material.SpectralMaterial)
        ):
            raise TypeError("material must be a Material or SpectralMaterial.")

        self.background = material

    def at(self, *, wavelength: Quantity) -> MomentRT.medium.StructuredMedium:
        """Resolve all materials on a geometry-preserving monochromatic copy.

        Parameters
        ----------
        wavelength : pint.Quantity
            Positive scalar vacuum wavelength, with length units.

        Returns
        -------
        StructuredMedium
            Geometry-preserving monochromatic medium; the original is unchanged.
        """

        if self.background is None:
            raise ValueError("Configure an explicit background with add_background().")

        sample = StructuredMedium(origin=self.origin, size=self.size)

        sample.add_background(material=self.background.at(wavelength=wavelength))

        structures = []

        for structure in self.structures:
            resolved = copy.copy(structure)

            for name, value in vars(structure).items():
                if isinstance(value, MomentRT.units.ureg.Quantity):
                    setattr(resolved, name, value.copy())

            resolved.material = structure.material.at(wavelength=wavelength)

            structures.append(resolved)

        sample.add_structures(*structures)

        return sample

    def _at_for_display(
        self, *, wavelength: Quantity | None
    ) -> MomentRT.medium.StructuredMedium:
        """Resolve spectral materials when needed for property sampling or display.

        Parameters
        ----------
        wavelength : pint.Quantity or None
            Vacuum wavelength used to resolve spectral materials. Required when
            any material is spectral; None retains fixed coefficients.

        Returns
        -------
        StructuredMedium
            Sample with fixed monochromatic coefficients for display.
        """

        if wavelength is not None:
            return self.at(wavelength=wavelength)

        materials = [self.background, *(shape.material for shape in self.structures)]

        if any(
            isinstance(material, MomentRT.material.SpectralMaterial)
            for material in materials
        ):
            raise ValueError("Supply wavelength when displaying spectral properties.")

        return self

    def add_structures(
        self,
        *structures: MomentRT.geometry.Box
        | MomentRT.geometry.Sphere
        | MomentRT.geometry.Layer,
    ) -> None:
        """Append named shapes in argument order; validate the whole batch first.

        Parameters
        ----------
        *structures : Box, Sphere, or Layer
            Named regions appended in argument order. Later regions replace
            earlier ones, including on shared interfaces; the batch is validated
            before any mutation.
        """

        lower = self.origin.to("m").magnitude

        upper = (self.origin + self.size).to("m").magnitude

        names = {shape.name for shape in self.structures}

        for shape in structures:
            if not isinstance(
                shape,
                (
                    MomentRT.geometry.Box,
                    MomentRT.geometry.Sphere,
                    MomentRT.geometry.Layer,
                ),
            ):
                raise TypeError("structures must be Box, Sphere, or Layer objects.")

            if shape.name in names:
                raise ValueError(f"Duplicate structure name: {shape.name}")

            names.add(shape.name)

            if isinstance(shape, MomentRT.geometry.Box):
                shape_lower = shape.origin.to("m").magnitude

                shape_upper = (shape.origin + shape.size).to("m").magnitude
            elif isinstance(shape, MomentRT.geometry.Sphere):
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

    def material_labels(self, *, points: Quantity) -> NDArray[np.int64]:
        """Return 0 for background and 1..N for ordered structures.

        Points must lie inside the finite sample. The trailing axis is x/y/z.

        Parameters
        ----------
        points : pint.Quantity
            Query coordinates with length units and shape (..., 3).

        Returns
        -------
        ndarray of int
            Labels with shape points.shape[:-1]: zero for background and 1..N
            for ordered structures. Later structures win at overlaps.

        Raises
        ------
        ValueError
            If no background is configured, points have invalid shape, or any
            query lies outside the finite sample.
        """

        MomentRT.units.validate_units(value=points, unit="m", name="points")

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

    def optical_fields(
        self,
        *,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
        wavelength: Quantity | None = None,
    ) -> dict[str, MomentRT.fields.ScalarField]:
        """Sample the material map onto mesh nodes, returning unit-bearing fields.

        Parameters
        ----------
        mesh : CartesianMesh or TetrahedralMesh
            Cartesian nodes used for sampling or finite-element discretization.
        wavelength : pint.Quantity or None, optional
            Vacuum wavelength used to resolve spectral materials. Required when
            any material is spectral; None retains fixed coefficients.
            Default is None.

        Returns
        -------
        dict of str to ScalarField
            Nodal absorption and reduced scattering in inverse meters, plus
            dimensionless refractive index. This is sampling, not interface fitting.
        """

        sample = self._at_for_display(wavelength=wavelength)

        labels = sample.material_labels(points=mesh.points)

        materials = [
            cast(MomentRT.material.Material, sample.background),
            *(shape.material for shape in sample.structures),
        ]

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

            fields[name] = MomentRT.fields.ScalarField(
                mesh=mesh,
                values=values[labels] * MomentRT.units.ureg(unit),
                name=name,
                metadata={"association": "node", "overlap": "last structure wins"},
            )

        return fields

    def plot(
        self,
        *,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh | None = None,
        wavelength: Quantity | None = None,
        length_unit: str = "mm",
        opacity: float = 0.45,
        overlays: Iterable[object] | None = None,
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Plot named region surfaces; mesh is optional and does not define geometry.

        Transparent primitive surfaces show the construction geometry; overlap
        precedence follows insertion order, not optical mixing or CSG subtraction.

        Parameters
        ----------
        mesh : CartesianMesh, TetrahedralMesh or None, optional
            Optional compatibility argument; construction surfaces use the
            physical sample bounds and do not depend on a numerical mesh.
            Default is None.
        wavelength : pint.Quantity or None, optional
            Vacuum wavelength used to resolve spectral materials. Required when
            any material is spectral; None retains fixed coefficients.
            Default is None.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        opacity : float, optional
            Construction-surface opacity in [0, 1].
            Default is 0.45.
        overlays : iterable of objects or None, optional
            Sources, detectors, and layers to draw on top of the sample
            geometry. Objects must implement the
            internal plotting hook used by the field visualizations.
            Default is None.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable interactive figure, also displayed when show is True.
        """

        if self.background is None:
            raise ValueError("Configure an explicit background with add_background().")

        return MomentRT.plotting.plot_sample(
            sample=self._at_for_display(wavelength=wavelength),
            length_unit=length_unit,
            opacity=opacity,
            overlays=overlays,
            show=show,
        )

    def add_layer(
        self,
        *,
        name: str,
        thickness: Quantity,
        material: MomentRT.material.Material | MomentRT.material.SpectralMaterial,
        clip_to_sample: bool = False,
    ) -> None:
        """Append a layer below the deepest existing layer; return None.

        The first layer starts at origin.z. Boxes and spheres do not change
        stacking depth. Explicit Layer objects added through add_structures
        also participate. New layers follow the usual insertion precedence.
        When clip_to_sample is True, truncate a layer at the sample's maximum z
        boundary when its requested thickness extends beyond the sample. A
        layer with no overlap with the sample is omitted.

        Parameters
        ----------
        name : str
            Nonempty name identifying this object.
        thickness : pint.Quantity
            Positive scalar layer thickness with length units.
        material : Material or SpectralMaterial
            Optical coefficients assigned to this region.
        clip_to_sample : bool, optional
            Whether to truncate the layer to the remaining depth up to the
            sample's maximum z coordinate.
            A fully out-of-bounds layer is omitted when True. Default is False,
            which raises an error if the layer extends outside.
        """

        if not isinstance(clip_to_sample, (bool, np.bool_)):
            raise TypeError("clip_to_sample must be a boolean.")

        MomentRT.units.validate_units(
            value=thickness,
            unit="m",
            name="thickness",
            shape=(),
            positive=True,
        )

        depth = self.origin[2]

        for structure in self.structures:
            if isinstance(structure, MomentRT.geometry.Layer):
                depth = max(depth, structure.depth + structure.thickness)

        if clip_to_sample:
            remaining_depth = (self.origin[2] + self.size[2] - depth).to(
                thickness.units
            )

            if remaining_depth <= 0 * remaining_depth.units:
                return

            if thickness > remaining_depth:
                thickness = remaining_depth

        layer = MomentRT.geometry.Layer(
            name=name,
            depth=depth,
            thickness=thickness,
            material=material,
        )

        self.add_structures(layer)

    def section(
        self,
        *,
        normal: Literal["x", "y", "z"],
        position: Quantity,
        resolution: tuple[int, int] = (200, 200),
        wavelength: Quantity | None = None,
    ) -> MomentRT.sections.SampleSection:
        """Sample a plane at an absolute coordinate, resolving region overlaps.

        resolution gives (horizontal, vertical) cell counts. Cell-center
        samples represent a discrete map, not an exact geometric intersection.

        Parameters
        ----------
        normal : {'x', 'y', 'z'}
            Axis normal to the section plane.
        position : pint.Quantity
            Absolute scalar coordinate along the plane normal, with length units.
        resolution : tuple of int, optional
            Horizontal and vertical cell counts, each at least two.
            The stored section arrays have shape (vertical, horizontal).
            Default is (200, 200).
        wavelength : pint.Quantity or None, optional
            Vacuum wavelength used to resolve spectral materials. Required when
            any material is spectral; None retains fixed coefficients.
            Default is None.

        Returns
        -------
        SampleSection
            Snapshot of cell-center labels and SI optical-property arrays.
        """

        return MomentRT.sections.SampleSection(
            sample=self,
            normal=normal,
            position=position,
            resolution=resolution,
            wavelength=wavelength,
        )

    def plot_section(
        self,
        *,
        normal: Literal["x", "y", "z"],
        position: Quantity,
        quantity: Literal[
            "material", "absorption", "reduced_scattering", "refractive_index"
        ] = "material",
        resolution: tuple[int, int] = (200, 200),
        wavelength: Quantity | None = None,
        length_unit: str = "mm",
        value_unit: str | None = None,
        scale: Literal["linear", "log"] = "linear",
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Plot a resolved cross-section, with depth increasing downwards.

        Parameters
        ----------
        normal : {'x', 'y', 'z'}
            Axis normal to the section plane.
        position : pint.Quantity
            Absolute scalar coordinate along the plane normal, with length units.
        quantity : {'material', 'absorption', 'reduced_scattering', 'refractive_index'}, optional
            Material labels or optical property to display.
            Default is 'material'.
        resolution : tuple of int, optional
            Horizontal and vertical cell counts, each at least two.
            The stored section arrays have shape (vertical, horizontal).
            Default is (200, 200).
        wavelength : pint.Quantity or None, optional
            Vacuum wavelength used to resolve spectral materials. Required when
            any material is spectral; None retains fixed coefficients.
            Default is None.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        value_unit : str or None, optional
            Unit for displayed optical-property values; None uses stored SI units.
            Default is None.
        scale : {'linear', 'log'}, optional
            Display scaling. Logarithmic scaling requires strictly positive values.
            Default is 'linear'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Resolved cross-section figure with depth increasing downwards.
        """
        section = self.section(
            normal=normal,
            position=position,
            resolution=resolution,
            wavelength=wavelength,
        )

        return section.plot(
            quantity=quantity,
            length_unit=length_unit,
            value_unit=value_unit,
            scale=scale,
            show=show,
        )
