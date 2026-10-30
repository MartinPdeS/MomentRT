"""Cartesian observation grids and conforming finite-sample tetrahedral meshes."""

from __future__ import annotations
import copy
from collections.abc import Iterable, Mapping
from threading import RLock
from typing import Any, cast
from uuid import uuid4
import numpy as np
from numpy.typing import ArrayLike, NDArray
from pint import Quantity
import MomentRT

_gmsh_import_error: OSError | None = None

try:
    import gmsh
except OSError as error:
    gmsh = None
    _gmsh_import_error = error


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

        self._resolution: Quantity | None = None

    @property
    def bounds(self) -> Quantity:
        """Return the lower and upper node coordinates.

        Returns
        -------
        pint.Quantity
            Cartesian bounds with shape (2, 3) and length units.
        """

        return cast(
            Quantity,
            np.stack(
                [
                    self.origin,
                    self.origin + (np.asarray(self.shape) - 1) * self.spacing,
                ],
            ),
        )

    @property
    def name(self) -> str:
        """Return the default overlay name for this mesh.

        Returns
        -------
        str
            Display name used in plot legends.
        """
        return "Cartesian mesh"

    @property
    def color(self) -> str:
        """Return the default color used when this mesh is an overlay.

        Returns
        -------
        str
            Slate color used for mesh edges.
        """
        return "#64748b"

    @property
    def resolution(self) -> Quantity | None:
        """Return the requested maximum spacing, when supplied to from_sample.

        Returns
        -------
        pint.Quantity or None
            Independent copy of the requested scalar or (3,) length quantity.
            None when the mesh was configured using explicit node counts.
            Actual node separations are available through spacing.
        """

        return None if self._resolution is None else self._resolution.copy()

    def copy(self) -> MomentRT.mesh.CartesianMesh:
        """Copy the mesh coordinates and requested resolution independently.

        Returns
        -------
        CartesianMesh
            Independent mesh snapshot with identical nodes and resolution.
        """

        mesh = type(self)(
            shape=self.shape,
            spacing=self.spacing.copy(),
            origin=self.origin.copy(),
        )

        mesh._resolution = self.resolution

        return mesh

    def plot(
        self, *, length_unit: str = "mm", max_edges: int = 20000, show: bool = True
    ) -> MomentRT.plotting.Figure:
        """Plot the Cartesian grid and its translucent bounding surface.

        Parameters
        ----------
        length_unit : str, optional
            Coordinate unit. Default is 'mm'.
        max_edges : int, optional
            Maximum number of evenly selected grid edges. Default is 20000.
        show : bool, optional
            Display immediately when True. Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable grid plot with physical coordinates and depth direction.
        """
        if (
            isinstance(max_edges, bool)
            or not isinstance(max_edges, (int, np.integer))
            or max_edges < 1
        ):
            raise ValueError("max_edges must be a positive integer.")

        edge_shapes = [
            tuple(count - (index == axis) for index, count in enumerate(self.shape))
            for axis in range(3)
        ]
        edge_counts = np.asarray([np.prod(shape) for shape in edge_shapes])
        cumulative_counts = np.cumsum(edge_counts)
        total_edges = int(cumulative_counts[-1])
        selected = np.linspace(
            0, total_edges - 1, min(max_edges, total_edges), dtype=np.intp
        )
        edge_axes = np.searchsorted(cumulative_counts, selected, side="right")
        coordinates = np.empty((len(selected), 3, 3))
        origin = self.origin.to(length_unit).magnitude
        spacing = self.spacing.to(length_unit).magnitude

        for axis in range(3):
            mask = edge_axes == axis

            if not np.any(mask):
                continue

            offset = 0 if axis == 0 else cumulative_counts[axis - 1]
            starts = np.column_stack(
                np.unravel_index(selected[mask] - offset, edge_shapes[axis])
            )
            ends = starts.copy()
            ends[:, axis] += 1
            coordinates[mask, 0] = origin + starts * spacing
            coordinates[mask, 1] = origin + ends * spacing

        lines = np.full((len(selected), 3, 3), np.nan)
        lines[:, :2] = coordinates[:, :2]
        line_points = lines.reshape(-1, 3)
        bounds = self.bounds.to(length_unit).magnitude
        vertices = np.asarray(
            [
                [x, y, z]
                for x in bounds[:, 0]
                for y in bounds[:, 1]
                for z in bounds[:, 2]
            ]
        )
        triangles = np.asarray(
            [
                [0, 1, 3],
                [0, 3, 2],
                [4, 6, 7],
                [4, 7, 5],
                [0, 4, 5],
                [0, 5, 1],
                [2, 3, 7],
                [2, 7, 6],
                [0, 2, 6],
                [0, 6, 4],
                [1, 5, 7],
                [1, 7, 3],
            ]
        )
        go = MomentRT.plotting.require_plotly()
        figure = MomentRT.plotting.make_figure(length_unit=length_unit)
        figure.add_trace(
            go.Mesh3d(
                x=vertices[:, 0],
                y=vertices[:, 1],
                z=vertices[:, 2],
                i=triangles[:, 0],
                j=triangles[:, 1],
                k=triangles[:, 2],
                color="lightsteelblue",
                opacity=0.16,
                name="Mesh bounds",
                hoverinfo="skip",
            )
        )
        figure.add_trace(
            go.Scatter3d(
                x=line_points[:, 0],
                y=line_points[:, 1],
                z=line_points[:, 2],
                mode="lines",
                line={"color": "slategray", "width": 1},
                name="Cartesian grid edges",
            )
        )
        figure.update_layout(
            title=f"Cartesian mesh {self.shape}; {len(selected)}/{total_edges} edges shown",
            scene_zaxis={"autorange": "reversed"},
        )

        if show:
            figure.show()

        return figure

    def _add_to_figure(
        self,
        *,
        figure: MomentRT.plotting.Figure,
        context: Mapping[str, Any],
    ) -> None:
        """Add Cartesian grid edges as an overlay to a Plotly figure.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Figure modified in place.
        context : mapping of str to object
            Shared display units, name, color, and legend group.
        """
        max_edges = context.get("max_edges", 2000)
        edge_shapes = [
            tuple(count - (index == axis) for index, count in enumerate(self.shape))
            for axis in range(3)
        ]
        edge_counts = np.asarray([np.prod(shape) for shape in edge_shapes])
        cumulative_counts = np.cumsum(edge_counts)
        total_edges = int(cumulative_counts[-1])
        selected = np.linspace(
            0, total_edges - 1, min(max_edges, total_edges), dtype=np.intp
        )
        edge_axes = np.searchsorted(cumulative_counts, selected, side="right")
        coordinates = np.empty((len(selected), 3, 3))
        origin = self.origin.to(context["length_unit"]).magnitude
        spacing = self.spacing.to(context["length_unit"]).magnitude

        for axis in range(3):
            mask = edge_axes == axis

            if not np.any(mask):
                continue

            offset = 0 if axis == 0 else cumulative_counts[axis - 1]
            starts = np.column_stack(
                np.unravel_index(selected[mask] - offset, edge_shapes[axis])
            )
            ends = starts.copy()
            ends[:, axis] += 1
            coordinates[mask, 0] = origin + starts * spacing
            coordinates[mask, 1] = origin + ends * spacing

        lines = np.full((len(selected), 3, 3), np.nan)
        lines[:, :2] = coordinates[:, :2]
        line_points = lines.reshape(-1, 3)
        figure.add_trace(
            MomentRT.plotting.require_plotly().Scatter3d(
                x=line_points[:, 0],
                y=line_points[:, 1],
                z=line_points[:, 2],
                mode="lines",
                line={"color": context["overlay_color"], "width": 1},
                opacity=0.55,
                name=context["name"],
                legendgroup=context["legendgroup"],
                hoverinfo="skip",
            )
        )

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
        resolution: Quantity | None = None,
        shape: tuple[int, int, int] | None = None,
    ) -> MomentRT.mesh.CartesianMesh:
        """Create nodes spanning every face of the finite rectangular sample.

        Parameters
        ----------
        sample : StructuredMedium
            Finite rectangular sample with an explicit background.
        resolution : pint.Quantity or None, optional
            Positive maximum node spacing with length units. A scalar applies
            to all axes; a quantity of shape (3,) specifies (rx, ry, rz).
            Actual uniform spacing on each axis is adjusted downward to fit
            the sample exactly, within floating-point precision. Default is
            None. Supply exactly one of resolution and shape.
        shape : tuple of int or None, optional
            Node counts (nx, ny, nz), each at least two. Default is None.
            Supply exactly one of resolution and shape.

        Returns
        -------
        CartesianMesh
            Nodes spanning all sample faces, with spacing size / (shape - 1).
        """

        if (resolution is None) == (shape is None):
            raise ValueError("Supply exactly one of resolution or shape.")

        if resolution is not None:
            MomentRT.units.validate_units(
                value=resolution,
                unit="m",
                name="resolution",
                positive=True,
            )

            if np.shape(resolution.magnitude) not in ((), (3,)):
                raise ValueError("resolution must be scalar or have shape (3,).")

            # Dimensionless ratios become integer interval counts here.
            ratios = (sample.size / resolution).to("").magnitude

            nearest = np.rint(ratios)

            ratios = np.where(
                np.isclose(ratios, nearest, rtol=16 * np.finfo(float).eps, atol=0),
                nearest,
                ratios,
            )

            intervals = np.maximum(1, np.ceil(ratios))

            if np.any(~np.isfinite(intervals)) or np.any(
                intervals >= np.iinfo(np.intp).max
            ):
                raise ValueError("resolution requires too many mesh intervals.")

            shape = (
                int(intervals[0]) + 1,
                int(intervals[1]) + 1,
                int(intervals[2]) + 1,
            )

        assert shape is not None

        # Constructor validates counts before spacing is calculated.
        candidate = cls(shape=shape, spacing=sample.size, origin=sample.origin)

        mesh = cls(
            shape=candidate.shape,
            spacing=sample.size / (np.asarray(candidate.shape) - 1),
            origin=sample.origin.copy(),
        )

        mesh._resolution = None if resolution is None else resolution.copy()

        return mesh


class TetrahedralMesh:
    """Conforming linear tetrahedra filling a finite rectangular sample.

    Generation and local refinement use the required Gmsh dependency.
    Imported connectivity and numerical solves do not require Gmsh. Resolution
    is a target element size, not a guaranteed upper bound on every edge.

    Parameters
    ----------
    points : pint.Quantity
        Node positions with shape (N, 3) and length units.
    tetrahedra : array_like of int
        Tetrahedral node indices with shape (T, 4). Cells must fill the bounding
        rectangle without cavities, duplicate cells, or nonmanifold faces.
    """

    _lock = RLock()

    def __init__(self, *, points: Quantity, tetrahedra: ArrayLike) -> None:
        """Validate and snapshot a tetrahedral discretization.

        Parameters
        ----------
        points : pint.Quantity
            Finite node positions with shape (N, 3) and length units.
        tetrahedra : array_like of int
            Valid node indices with shape (T, 4); N >= 4 and T >= 1.
        """
        MomentRT.units.validate_units(value=points, unit="m", name="points")

        coordinates = np.asarray(points.to("m").magnitude, dtype=float)
        cells = np.asarray(tetrahedra)

        if coordinates.ndim != 2 or coordinates.shape[1] != 3 or len(coordinates) < 4:
            raise ValueError("points must have shape (N, 3), with N >= 4.")

        if (
            cells.ndim != 2
            or cells.shape[1] != 4
            or len(cells) == 0
            or cells.dtype.kind not in "iu"
            or np.any(cells < 0)
            or np.any(cells >= len(coordinates))
        ):
            raise ValueError(
                "tetrahedra must have shape (T, 4) with valid integer indices."
            )

        if len(np.unique(np.sort(cells, axis=1), axis=0)) != len(cells):
            raise ValueError("Duplicate tetrahedra are not allowed.")

        if len(np.unique(cells)) != len(coordinates):
            raise ValueError("Every node must belong to a tetrahedron.")

        vertices = coordinates[cells]
        edges = vertices[:, 1:] - vertices[:, :1]
        determinants = np.abs(np.linalg.det(edges))
        scales = np.max(np.linalg.norm(edges, axis=2), axis=1)

        if np.any(determinants <= 64 * np.finfo(float).eps * scales**3):
            raise ValueError(
                "Tetrahedra must have nonzero, numerically resolved volumes."
            )

        self._points = points.copy()
        self._tetrahedra = cells.astype(np.int64, copy=True)
        self._resolution: Quantity | None = None
        self._sample: MomentRT.medium.StructuredMedium | None = None
        self._region_refinements: dict[str, Quantity] = {}
        self._contact_refinements: list[
            tuple[
                MomentRT.source.SurfaceSource | MomentRT.detector.SurfaceDetector,
                Quantity,
                Quantity,
            ]
        ] = []
        self._faces = self._boundary_faces(coordinates=coordinates)

        extent = coordinates.max(axis=0) - coordinates.min(axis=0)

        if np.any(extent <= 0) or not np.isclose(
            determinants.sum() / 6,
            np.prod(extent),
            rtol=1e-8,
            atol=0,
        ):
            raise ValueError("Tetrahedra must fill the finite rectangular bounds.")

    @property
    def points(self) -> Quantity:
        """Return an independent copy of the unstructured node positions.

        Returns
        -------
        pint.Quantity
            Positions with shape (N, 3) and length units.
        """
        return self._points.copy()

    @property
    def name(self) -> str:
        """Return the default overlay name for this mesh.

        Returns
        -------
        str
            Display name used in plot legends.
        """
        return "Tetrahedral mesh"

    @property
    def color(self) -> str:
        """Return the default color used when this mesh is an overlay.

        Returns
        -------
        str
            Slate color used for mesh edges.
        """
        return "#64748b"

    @property
    def tetrahedra(self) -> NDArray[np.int64]:
        """Return an independent copy of element connectivity.

        Returns
        -------
        ndarray of int
            Tetrahedral node indices with shape (T, 4).
        """
        return self._tetrahedra.copy()

    @property
    def faces(self) -> dict[MomentRT.boundary.Face, NDArray[np.int64]]:
        """Return exterior triangles grouped by the six physical sample faces.

        Returns
        -------
        dict of str to ndarray of int
            Independent triangle connectivity arrays with shape (F, 3).
        """
        return {name: triangles.copy() for name, triangles in self._faces.items()}

    @property
    def shape(self) -> tuple[int]:
        """Return the storage shape of an unstructured nodal scalar field.

        Returns
        -------
        tuple of int
            (N,), where N is the number of nodes; these are not grid axes.
        """
        return (len(self._points),)

    @property
    def bounds(self) -> Quantity:
        """Return the physical lower and upper node coordinates.

        Returns
        -------
        pint.Quantity
            Bounds with shape (2, 3) and length units.
        """
        return cast(
            Quantity, np.stack([self._points.min(axis=0), self._points.max(axis=0)])
        )

    @property
    def resolution(self) -> Quantity | None:
        """Return the requested global target element size, when generated.

        Returns
        -------
        pint.Quantity or None
            Scalar length quantity, or None for imported connectivity. Actual
            edges may exceed the target; local rules can request smaller sizes.
        """
        return None if self._resolution is None else self._resolution.copy()

    def copy(self) -> MomentRT.mesh.TetrahedralMesh:
        """Copy geometry and refinement rules without invoking the mesher.

        Returns
        -------
        TetrahedralMesh
            Independent snapshot unaffected by subsequent refinement.
        """
        return copy.deepcopy(self)

    def _boundary_faces(
        self, *, coordinates: NDArray[np.float64]
    ) -> dict[MomentRT.boundary.Face, NDArray[np.int64]]:
        """Validate face incidence and label exterior triangles by coordinates."""
        candidates = np.concatenate(
            [
                self._tetrahedra[:, face]
                for face in ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])
            ],
        )
        unique, counts = np.unique(
            np.sort(candidates, axis=1), axis=0, return_counts=True
        )

        if np.any(counts > 2):
            raise ValueError("Nonmanifold tetrahedral faces are not allowed.")

        boundary = unique[counts == 1]
        lower, upper = coordinates.min(axis=0), coordinates.max(axis=0)
        tolerance = np.max(upper - lower) * 1e-10
        labels = np.zeros(len(boundary), dtype=int)
        faces = {}

        for face in MomentRT.boundary.FACES:
            axis = "xyz".index(face[0])
            value = lower[axis] if face.endswith("min") else upper[axis]
            selected = np.all(
                np.abs(coordinates[boundary, axis] - value) <= tolerance, axis=1
            )
            faces[face] = boundary[selected]
            labels += selected

            if not np.any(selected):
                raise ValueError(f"Mesh has no exterior triangles on {face}.")

            vertices = coordinates[faces[face]]
            area = (
                np.linalg.norm(
                    np.cross(
                        vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]
                    ),
                    axis=1,
                ).sum()
                / 2
            )
            expected = np.prod(
                (upper - lower)[[index for index in range(3) if index != axis]]
            )

            if not np.isclose(area, expected, rtol=1e-8, atol=0):
                raise ValueError(f"Exterior triangles must cover {face} exactly.")

        if np.any(labels != 1):
            raise ValueError(
                "Exterior triangles must lie on the rectangular sample faces."
            )

        return faces

    @classmethod
    def from_sample(
        cls, *, sample: MomentRT.medium.StructuredMedium, resolution: Quantity
    ) -> MomentRT.mesh.TetrahedralMesh:
        """Generate a mesh conforming to the sample and all region interfaces.

        Box, layer, and sphere surfaces are fragmented together by Gmsh. Final
        material assignment retains the sample's insertion precedence.

        Parameters
        ----------
        sample : StructuredMedium
            Finite rectangular geometry. A snapshot is retained for refinement.
        resolution : pint.Quantity
            Positive scalar length setting the global target element size.

        Returns
        -------
        TetrahedralMesh
            Generated linear tetrahedra with shared nodes at region interfaces.
        """
        cls._validate_resolution(resolution=resolution)

        mesh = cls.__new__(cls)
        mesh._sample = copy.deepcopy(sample)
        checked = MomentRT.medium.StructuredMedium(
            origin=sample.origin, size=sample.size
        )
        checked.add_structures(*mesh._sample.structures)
        mesh._resolution = resolution.copy()
        mesh._region_refinements = {}
        mesh._contact_refinements = []
        mesh._generate()

        return mesh

    @staticmethod
    def _validate_resolution(*, resolution: Quantity) -> None:
        """Validate one physical scalar target size."""
        MomentRT.units.validate_units(
            value=resolution,
            unit="m",
            name="resolution",
            shape=(),
            positive=True,
        )

    def refine_region(self, *, name: str, resolution: Quantity) -> None:
        """Regenerate with a smaller target size in a named region; return None.

        The rule covers the region's geometric support, including parts later
        overridden by other regions. Generation succeeds before mutation.

        Parameters
        ----------
        name : str
            Name of a Box, Sphere, or Layer in the original sample snapshot.
        resolution : pint.Quantity
            Positive scalar target length no larger than the global resolution.
        """
        self._validate_refinement(resolution=resolution)
        assert self._sample is not None

        if name not in {region.name for region in self._sample.structures}:
            raise ValueError(f"Unknown mesh region: {name!r}.")

        candidate = self.copy()
        candidate._region_refinements[name] = resolution.copy()
        candidate._generate()
        self.__dict__.update(candidate.__dict__)

    def refine_near(
        self, *, objects: Iterable[object], resolution: Quantity, distance: Quantity
    ) -> None:
        """Regenerate with refinement near physical contact footprints.

        The fine target applies within distance of each disk, including its
        interior and radius. It transitions to the global size over another
        distance. This changes the mesh, not source or detector quadrature.

        Parameters
        ----------
        objects : iterable of objects
            SurfaceSource, SurfaceSourceGroup, or SurfaceDetector objects.
            Geometry is snapshotted; contacts must fit on their sample face.
        resolution : pint.Quantity
            Positive scalar target length no larger than the global resolution.
        distance : pint.Quantity
            Positive scalar length of the fine neighborhood and transition band.
        """
        self._validate_refinement(resolution=resolution)
        MomentRT.units.validate_units(
            value=distance,
            unit="m",
            name="distance",
            shape=(),
            positive=True,
        )
        assert self._sample is not None
        contacts = []

        for item in objects:
            members = (
                item.sources
                if isinstance(item, MomentRT.source.SurfaceSourceGroup)
                else (item,)
            )

            for member in members:
                if not isinstance(
                    member,
                    (MomentRT.source.SurfaceSource, MomentRT.detector.SurfaceDetector),
                ):
                    raise TypeError(
                        "refine_near requires circular surface sources or detectors."
                    )

                MomentRT.boundary._validate_face_disk(
                    sample=self._sample,
                    face=member.face,
                    position=member.position,
                    radius=member.radius,
                )
                axis = "xyz".index(member.face[0])
                bounds = (
                    self._sample.origin
                    if member.face.endswith("min")
                    else self._sample.origin + self._sample.size
                )

                if not np.isclose(
                    member.position[axis].to("m").magnitude,
                    bounds[axis].to("m").magnitude,
                    rtol=1e-12,
                    atol=1e-15,
                ):
                    raise ValueError(
                        "Refinement contacts must lie on their configured sample face."
                    )

                contacts.append(
                    (copy.deepcopy(member), resolution.copy(), distance.copy())
                )

        if not contacts:
            raise ValueError("objects must contain at least one circular contact.")

        candidate = self.copy()
        candidate._contact_refinements.extend(contacts)
        candidate._generate()
        self.__dict__.update(candidate.__dict__)

    def _validate_refinement(self, *, resolution: Quantity) -> None:
        """Check that refinement can regenerate a geometry snapshot."""
        self._validate_resolution(resolution=resolution)

        if self._sample is None or self._resolution is None:
            raise ValueError("Refinement requires TetrahedralMesh.from_sample().")

        if resolution > self._resolution:
            raise ValueError(
                "Refinement resolution must not exceed the global resolution."
            )

    def _validate_sample(self, *, sample: MomentRT.medium.StructuredMedium) -> None:
        """Reject stale fitted interfaces after the generation geometry changes."""
        if self._sample is None:
            return

        matches = len(sample.structures) == len(self._sample.structures)

        for original, current in zip(self._sample.structures, sample.structures):
            if type(original) is not type(current) or original.name != current.name:
                matches = False
                break

            for attribute in (
                "origin",
                "size",
                "center",
                "radius",
                "depth",
                "thickness",
            ):
                if hasattr(original, attribute) and not np.allclose(
                    getattr(original, attribute).to("m").magnitude,
                    getattr(current, attribute).to("m").magnitude,
                    rtol=1e-12,
                    atol=1e-15,
                ):
                    matches = False

        if not matches:
            raise ValueError(
                "Sample geometry changed after mesh generation; regenerate the fitted mesh."
            )

    def _generate(self) -> None:
        """Build a private Gmsh model and adopt its validated connectivity."""
        if gmsh is None:
            raise RuntimeError(
                "Gmsh is installed, but its native runtime could not be loaded. "
                "Install the system OpenGL utility library (libGLU; "
                "libglu1-mesa on Ubuntu)."
            ) from _gmsh_import_error

        assert self._sample is not None and self._resolution is not None
        lower = self._sample.origin.to("m").magnitude
        size = self._sample.size.to("m").magnitude
        scale = float(size.max())
        base = float(self._resolution.to("m").magnitude / scale)

        with self._lock:
            owned = not gmsh.isInitialized()

            if owned:
                gmsh.initialize(readConfigFiles=False, interruptible=False)

            previous = gmsh.model.getCurrent()
            options = {
                "General.Terminal": 0,
                "Mesh.MeshSizeMax": base,
                "Mesh.MeshSizeMin": 0,
                "Mesh.MeshSizeFromPoints": 0,
                "Mesh.MeshSizeFromCurvature": 0,
                "Mesh.MeshSizeExtendFromBoundary": 0,
                "Mesh.MeshSizeFactor": 1,
                "Mesh.ElementOrder": 1,
                "Mesh.Algorithm3D": 1,
            }
            saved = {key: gmsh.option.getNumber(key) for key in options}
            gmsh.model.add(f"MomentRT-{uuid4().hex}")

            try:
                for key, value in options.items():
                    gmsh.option.setNumber(key, value)

                domain = gmsh.model.occ.addBox(0, 0, 0, *(size / scale))
                tools = []

                for region in self._sample.structures:
                    if isinstance(region, MomentRT.geometry.Sphere):
                        center = (region.center.to("m").magnitude - lower) / scale
                        tag = gmsh.model.occ.addSphere(
                            *center, float(region.radius.to("m").magnitude / scale)
                        )
                    else:
                        region_lower, region_size = self._region_bounds(region=region)
                        tag = gmsh.model.occ.addBox(
                            *((region_lower - lower) / scale), *(region_size / scale)
                        )

                    tools.append((3, tag))

                if tools:
                    gmsh.model.occ.fragment([(3, domain)], tools)

                disks = []
                contact_geometries = set()

                for contact, _, _ in self._contact_refinements:
                    center = (contact.position.to("m").magnitude - lower) / scale
                    radius = float(contact.radius.to("m").magnitude / scale)
                    contact_key = (contact.face, tuple(center), radius)

                    if contact_key in contact_geometries:
                        continue

                    contact_geometries.add(contact_key)
                    disk = gmsh.model.occ.addDisk(*center, radius, radius)

                    if contact.face[0] == "x":
                        gmsh.model.occ.rotate([(2, disk)], *center, 0, 1, 0, np.pi / 2)
                    elif contact.face[0] == "y":
                        gmsh.model.occ.rotate([(2, disk)], *center, 1, 0, 0, np.pi / 2)

                    disks.append((2, disk))

                if disks:
                    gmsh.model.occ.fragment(gmsh.model.occ.getEntities(3), disks)

                gmsh.model.occ.synchronize()
                fields = []

                for region in self._sample.structures:
                    if region.name not in self._region_refinements:
                        continue

                    fine = float(
                        self._region_refinements[region.name].to("m").magnitude / scale
                    )

                    if isinstance(region, MomentRT.geometry.Sphere):
                        field = gmsh.model.mesh.field.add("Ball")
                        center = (region.center.to("m").magnitude - lower) / scale

                        for axis, coordinate in zip("XYZ", center):
                            gmsh.model.mesh.field.setNumber(
                                field, f"{axis}Center", coordinate
                            )

                        gmsh.model.mesh.field.setNumber(
                            field,
                            "Radius",
                            float(region.radius.to("m").magnitude / scale),
                        )
                    else:
                        field = gmsh.model.mesh.field.add("Box")
                        region_lower, region_size = self._region_bounds(region=region)

                        for axis, start, extent in zip(
                            "XYZ", (region_lower - lower) / scale, region_size / scale
                        ):
                            gmsh.model.mesh.field.setNumber(field, f"{axis}Min", start)
                            gmsh.model.mesh.field.setNumber(
                                field, f"{axis}Max", start + extent
                            )

                    gmsh.model.mesh.field.setNumber(field, "VIn", fine)
                    gmsh.model.mesh.field.setNumber(field, "VOut", base)
                    gmsh.model.mesh.field.setNumber(field, "Thickness", base)
                    fields.append(field)

                for contact, resolution, distance in self._contact_refinements:
                    center = (contact.position.to("m").magnitude - lower) / scale
                    radius = float(contact.radius.to("m").magnitude / scale)
                    normal = "xyz".index(contact.face[0])
                    radial = "+".join(
                        f"({axis}-{center[index]:.17g})^2"
                        for index, axis in enumerate("xyz")
                        if index != normal
                    )
                    gap = f"Sqrt(Max(0,Sqrt({radial})-{radius:.17g})^2+({'xyz'[normal]}-{center[normal]:.17g})^2)"
                    fine = float(resolution.to("m").magnitude / scale)
                    band = float(distance.to("m").magnitude / scale)
                    field = gmsh.model.mesh.field.add("MathEval")
                    gmsh.model.mesh.field.setString(
                        field,
                        "F",
                        f"{fine:.17g}+({base - fine:.17g})*Min(1,Max(0,({gap}-{band:.17g})/{band:.17g}))",
                    )
                    fields.append(field)

                if fields:
                    field = gmsh.model.mesh.field.add("Min")
                    gmsh.model.mesh.field.setNumbers(field, "FieldsList", fields)
                    gmsh.model.mesh.field.setAsBackgroundMesh(field)

                gmsh.model.mesh.generate(3)
                tags, coordinates, _ = gmsh.model.mesh.getNodes()
                _, element_tags = gmsh.model.mesh.getElementsByType(4)
                used = np.unique(element_tags)
                order = np.argsort(tags)
                points = np.asarray(coordinates).reshape(-1, 3)[order][
                    np.searchsorted(tags[order], used)
                ]
                cells = np.searchsorted(used, element_tags).reshape(-1, 4)
                points = lower + points * scale

                if np.any(points < lower - scale * 1e-10) or np.any(
                    points > lower + size + scale * 1e-10
                ):
                    raise RuntimeError("Gmsh generated nodes outside the sample.")

                # Remove only coordinate roundoff at the exact physical faces.
                points = np.clip(points, lower, lower + size)
                validated = type(self)(
                    points=points * MomentRT.ureg.m, tetrahedra=cells
                )
                self._points = validated._points
                self._tetrahedra = validated._tetrahedra
                self._faces = validated._faces
            except Exception as error:
                raise RuntimeError(f"Gmsh mesh generation failed: {error}") from error
            finally:
                gmsh.model.remove()

                for key, value in saved.items():
                    gmsh.option.setNumber(key, value)

                if previous in gmsh.model.list():
                    gmsh.model.setCurrent(previous)

                if owned:
                    gmsh.finalize()

    def _region_bounds(
        self, *, region: MomentRT.geometry.Box | MomentRT.geometry.Layer
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return SI box bounds for a layer or box in the sample."""
        assert self._sample is not None

        if isinstance(region, MomentRT.geometry.Box):
            return region.origin.to("m").magnitude, region.size.to("m").magnitude

        lower = self._sample.origin.to("m").magnitude.copy()
        size = self._sample.size.to("m").magnitude.copy()
        lower[2] = region.depth.to("m").magnitude
        size[2] = region.thickness.to("m").magnitude

        return lower, size

    def plot(
        self, *, length_unit: str = "mm", max_edges: int = 20000, show: bool = True
    ) -> MomentRT.plotting.Figure:
        """Display tetrahedral edges to inspect local mesh resolution.

        Parameters
        ----------
        length_unit : str, optional
            Coordinate unit. Default is 'mm'.
        max_edges : int, optional
            Maximum number of evenly selected edges. Default is 20000.
            The title reports when only a subset is displayed.
        show : bool, optional
            Display immediately when True. Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Composable edge plot with physical coordinates and depth direction.
        """
        if (
            isinstance(max_edges, bool)
            or not isinstance(max_edges, (int, np.integer))
            or max_edges < 1
        ):
            raise ValueError("max_edges must be a positive integer.")

        pairs = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
        edges = np.unique(
            np.sort(self._tetrahedra[:, pairs].reshape(-1, 2), axis=1), axis=0
        )
        selection = np.linspace(
            0, len(edges) - 1, min(max_edges, len(edges)), dtype=int
        )
        coordinates = self._points.to(length_unit).magnitude
        lines = np.full((len(selection), 3, 3), np.nan)
        lines[:, :2] = coordinates[edges[selection]]
        line_points = lines.reshape(-1, 3)
        boundary_triangles = np.concatenate(list(self._faces.values()))
        go = MomentRT.plotting.require_plotly()
        figure = MomentRT.plotting.make_figure(length_unit=length_unit)
        figure.add_trace(
            go.Mesh3d(
                x=coordinates[:, 0],
                y=coordinates[:, 1],
                z=coordinates[:, 2],
                i=boundary_triangles[:, 0],
                j=boundary_triangles[:, 1],
                k=boundary_triangles[:, 2],
                color="lightsteelblue",
                opacity=0.16,
                name="Sample boundary",
                hoverinfo="skip",
            )
        )
        figure.add_trace(
            go.Scatter3d(
                x=line_points[:, 0],
                y=line_points[:, 1],
                z=line_points[:, 2],
                mode="lines",
                line={"color": "slategray", "width": 1},
                name="Tetrahedral edges",
            )
        )
        figure.update_layout(
            title=f"{len(coordinates)} nodes, {len(self._tetrahedra)} tetrahedra; {len(selection)}/{len(edges)} edges shown",
            scene_zaxis={"autorange": "reversed"},
        )

        if show:
            figure.show()

        return figure

    def _add_to_figure(
        self,
        *,
        figure: MomentRT.plotting.Figure,
        context: Mapping[str, Any],
    ) -> None:
        """Add tetrahedral edges as an overlay to a Plotly figure.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Figure modified in place.
        context : mapping of str to object
            Shared display units, name, color, and legend group.
        """
        max_edges = context.get("max_edges", 2000)
        pairs = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
        edges = np.unique(
            np.sort(self._tetrahedra[:, pairs].reshape(-1, 2), axis=1), axis=0
        )
        selection = np.linspace(
            0, len(edges) - 1, min(max_edges, len(edges)), dtype=int
        )
        coordinates = self._points.to(context["length_unit"]).magnitude
        lines = np.full((len(selection), 3, 3), np.nan)
        lines[:, :2] = coordinates[edges[selection]]
        line_points = lines.reshape(-1, 3)
        figure.add_trace(
            MomentRT.plotting.require_plotly().Scatter3d(
                x=line_points[:, 0],
                y=line_points[:, 1],
                z=line_points[:, 2],
                mode="lines",
                line={"color": context["overlay_color"], "width": 1},
                opacity=0.55,
                name=context["name"],
                legendgroup=context["legendgroup"],
                hoverinfo="skip",
            )
        )

    def to_pyvista(
        self, *, length_unit: str = "mm"
    ) -> MomentRT.plotting.UnstructuredGrid:
        """Convert the actual tetrahedra to an optional PyVista grid.

        Parameters
        ----------
        length_unit : str, optional
            Coordinate unit. Default is 'mm'.

        Returns
        -------
        pyvista.UnstructuredGrid
            Linear tetrahedra with original node ordering and connectivity.
        """
        pv = MomentRT.plotting.require_pyvista()
        cells = np.column_stack([np.full(len(self._tetrahedra), 4), self._tetrahedra])

        return pv.UnstructuredGrid(
            cells.ravel(),
            np.full(len(cells), 10, dtype=np.uint8),
            self._points.to(length_unit).magnitude,
        )

    def _interpolate(self, *, points: Quantity, values: np.ndarray) -> np.ndarray:
        """Interpolate P1 nodal values using the original tetrahedral connectivity."""
        MomentRT.units.validate_units(value=points, unit="m", name="points")
        queries = np.asarray(points.to("m").magnitude)

        if queries.ndim < 1 or queries.shape[-1] != 3:
            raise ValueError(
                "points must have a trailing coordinate axis of length three."
            )

        output = np.full(queries.shape[:-1], np.nan)
        flat = queries.reshape(-1, 3)
        coordinates = self._points.to("m").magnitude
        vertices = coordinates[self._tetrahedra]
        inverse = np.linalg.inv(vertices[:, 1:] - vertices[:, :1])
        lower, upper = vertices.min(axis=1), vertices.max(axis=1)
        centers = vertices.mean(axis=1)
        tolerance = np.max(coordinates.max(axis=0) - coordinates.min(axis=0)) * 1e-10
        stack = [(np.arange(len(vertices)), np.arange(len(flat)))]

        while stack:
            elements, indices = stack.pop()
            selected = np.all(
                (flat[indices] >= lower[elements].min(axis=0) - tolerance)
                & (flat[indices] <= upper[elements].max(axis=0) + tolerance),
                axis=1,
            )
            indices = indices[selected]

            if len(indices) == 0:
                continue

            if len(elements) > 16:
                axis = np.argmax(np.ptp(centers[elements], axis=0))
                elements = elements[np.argsort(centers[elements, axis])]
                middle = len(elements) // 2
                stack.extend(
                    [(elements[:middle], indices), (elements[middle:], indices)]
                )
                continue

            for element in elements:
                local = (flat[indices] - vertices[element, 0]) @ inverse[element]
                inside = np.all(local >= -1e-10, axis=1) & (
                    local.sum(axis=1) <= 1 + 1e-10
                )
                weights = np.column_stack(
                    [1 - local[inside].sum(axis=1), local[inside]]
                )
                output.ravel()[indices[inside]] = (
                    weights @ values[self._tetrahedra[element]]
                )

        if np.any(~np.isfinite(output)):
            raise ValueError(
                "Interpolation points must lie inside the tetrahedral mesh."
            )

        return output
