"""Contact-aperture detection with an explicit Lambertian escape assumption."""

from __future__ import annotations
from collections.abc import Mapping
from typing import Any
import numpy as np
from numpy.typing import ArrayLike
from pint import Quantity
from scipy import integrate
import MomentRT


class DetectorResolutionError(ValueError):
    """No boundary quadrature point samples the detector aperture."""


class SurfaceDetector:
    """Circular aperture in a sample face with a cone of accepted outgoing rays.

    collection_axis points along collected ray propagation, out of the sample.
    The aperture lies in the configured face even when the collection axis is
    tilted. acceptance_angle is the cone half-angle, with explicit angular units.
    This is a passive contact measurement with Lambertian escaping radiance;
    it does not alter the transport boundary or model refraction and lenses.

    Parameters
    ----------
    face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
        Axis-aligned sample face; min and max refer to outward boundaries.
    position : pint.Quantity
        Absolute Cartesian position, shape (3,), with length units.
    radius : pint.Quantity
        Positive scalar radius with length units.
    collection_axis : array_like of float, shape (3,)
        Nonzero dimensionless direction along accepted outgoing rays.
        It is normalized internally and must point into the outward hemisphere.
    acceptance_angle : pint.Quantity
        Acceptance-cone half-angle, from 0 to 90 degrees, with angular units.
    name : str, optional
        Legend label. Default is ``"Surface detector"``.
    color : str, optional
        Plotly-compatible color. Default is ``"royalblue"``.

    Notes
    -----
    This is a passive contact measurement assuming Lambertian escaping
    radiance. It does not change the transport boundary or model lenses,
    refraction, or angle-resolved Fresnel transmission. For a cone wholly
    inside the outward hemisphere, the accepted fraction is
    cos(tilt) * sin(acceptance_angle)**2.
    """

    def __init__(
        self,
        *,
        face: MomentRT.boundary.Face,
        position: Quantity,
        radius: Quantity,
        collection_axis: ArrayLike,
        acceptance_angle: Quantity,
        name: str = "Surface detector",
        color: str = "royalblue",
    ) -> None:
        """Initialize the SurfaceDetector configuration.

        Parameters
        ----------
        face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
            Axis-aligned sample face; min and max refer to outward boundaries.
        position : pint.Quantity
            Absolute Cartesian position, shape (3,), with length units.
        radius : pint.Quantity
            Positive scalar radius with length units.
        collection_axis : array_like of float, shape (3,)
            Nonzero dimensionless direction along accepted outgoing rays.
            It is normalized internally and must point into the outward hemisphere.
        acceptance_angle : pint.Quantity
            Acceptance-cone half-angle, from 0 to 90 degrees, with angular units.
        name : str, optional
            Nonempty legend label.
        color : str, optional
            Nonempty Plotly-compatible color string.
        """

        if face not in MomentRT.boundary.FACES:
            raise ValueError(f"face must be one of {MomentRT.boundary.FACES}.")

        MomentRT.units.validate_units(
            value=position, unit="m", name="position", shape=(3,)
        )

        MomentRT.units.validate_units(
            value=radius, unit="m", name="radius", shape=(), positive=True
        )

        MomentRT.units.validate_units(
            value=acceptance_angle,
            unit="radian",
            name="acceptance_angle",
            shape=(),
            nonnegative=True,
        )

        if acceptance_angle.to("radian").magnitude > np.pi / 2:
            raise ValueError("acceptance_angle must be at most 90 degrees.")

        if isinstance(collection_axis, MomentRT.units.ureg.Quantity):
            raise TypeError("collection_axis must be a plain dimensionless vector.")

        direction = np.asarray(collection_axis)

        if (
            direction.shape != (3,)
            or direction.dtype.kind not in "iuf"
            or not np.all(np.isfinite(direction))
        ):
            raise ValueError("collection_axis must contain three finite plain numbers.")

        magnitude = np.linalg.norm(direction)

        if not np.isfinite(magnitude) or magnitude == 0:
            raise ValueError("collection_axis must be nonzero and finite.")

        direction = direction / magnitude

        normal = np.zeros(3)

        normal["xyz".index(face[0])] = -1 if face.endswith("min") else 1

        if np.dot(direction, normal) <= 0:
            raise ValueError("collection_axis must point into the outward hemisphere.")

        if not isinstance(name, str) or not name.strip():
            raise ValueError("name must be a nonempty string.")

        if not isinstance(color, str) or not color.strip():
            raise ValueError("color must be a nonempty Plotly color string.")

        self.face = face

        self.position = position.copy()

        self.radius = radius.copy()

        self.collection_axis = direction.copy()

        self.acceptance_angle = acceptance_angle.copy()

        self.name = name

        self.color = color

    @classmethod
    def on_face(
        cls,
        *,
        sample: MomentRT.medium.StructuredMedium,
        face: MomentRT.boundary.Face,
        coordinates: Quantity,
        radius: Quantity,
        acceptance_angle: Quantity,
        orientation: str,
        name: str = "Surface detector",
        color: str = "royalblue",
    ) -> SurfaceDetector:
        """Place a detector on a face with an explicitly selected outward normal.

        Parameters
        ----------
        sample : StructuredMedium
            Physical sample used to determine the face position.
        face : str
            Cartesian face name, such as ``z_min``.
        coordinates : pint.Quantity
            Two absolute in-plane coordinates, ordered by Cartesian axis
            (y, z for x faces; x, z for y faces; x, y for z faces).
        radius : pint.Quantity
            Positive aperture radius with length units.
        acceptance_angle : pint.Quantity
            Acceptance cone half-angle with angular units, from zero to 90 degrees.
        orientation : {'outward_normal'}
            Explicitly select the outward face normal as the collection axis.
            Use the main constructor for a tilted collection axis.
        name : str, optional
            Legend label. Default is ``"Surface detector"``.
        color : str, optional
            Plotly-compatible color. Default is ``"royalblue"``.

        Returns
        -------
        SurfaceDetector
            Detector with the derived 3D position and collection axis recorded.
        """
        if orientation != "outward_normal":
            raise ValueError("orientation must be 'outward_normal'.")

        position = MomentRT.boundary._face_position(
            sample=sample, face=face, coordinates=coordinates
        )

        axis = np.zeros(3)

        axis["xyz".index(face[0])] = -1 if face.endswith("min") else 1

        detector = cls(
            face=face,
            position=position,
            radius=radius,
            collection_axis=axis,
            acceptance_angle=acceptance_angle,
            name=name,
            color=color,
        )

        detector._validate_sample(sample_origin=sample.origin, sample_size=sample.size)

        return detector

    @classmethod
    def full_hemisphere_on_face(
        cls,
        *,
        sample: MomentRT.medium.StructuredMedium,
        face: MomentRT.boundary.Face,
        coordinates: Quantity,
        radius: Quantity,
        name: str = "Surface detector",
        color: str = "royalblue",
    ) -> SurfaceDetector:
        """Collect all escaping irradiance across a circular contact aperture.

        Parameters
        ----------
        sample : StructuredMedium
            Physical sample used to determine and validate the face position.
        face : {'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'}
            Axis-aligned sample face.
        coordinates : pint.Quantity
            Absolute in-plane coordinates, shape (2,), with length units, ordered
            by Cartesian axis.
        radius : pint.Quantity
            Positive aperture radius with length units.
        name : str, optional
            Legend label. Default is ``"Surface detector"``.
        color : str, optional
            Plotly-compatible color. Default is ``"royalblue"``.

        Returns
        -------
        SurfaceDetector
            A detector with outward-normal axis and 90-degree acceptance half-angle.

        Notes
        -----
        This ideal contact measurement accepts the entire outward hemisphere.
        It does not model external optics or change the transport boundary.
        """
        return cls.on_face(
            sample=sample,
            face=face,
            coordinates=coordinates,
            radius=radius,
            acceptance_angle=90 * MomentRT.units.ureg.degree,
            orientation="outward_normal",
            name=name,
            color=color,
        )

    @property
    def angular_fraction(self) -> float:
        """Fraction of Lambertian escaping irradiance inside the acceptance cone.

        Returns
        -------
        float
            Accepted fraction in [0, 1] under Lambertian escaping radiance.
        """

        axis = "xyz".index(self.face[0])

        cosine_tilt = self.collection_axis[axis] * (
            -1 if self.face.endswith("min") else 1
        )

        half_angle = float(self.acceptance_angle.to("radian").magnitude)

        tilt = np.arccos(np.clip(cosine_tilt, 0, 1))

        if half_angle == 0:
            return 0.0

        if tilt + half_angle <= np.pi / 2:
            return float(cosine_tilt * np.sin(half_angle) ** 2)

        sine_tilt = np.sqrt(max(0.0, 1 - cosine_tilt**2))

        if sine_tilt < 1e-14:
            return float(np.sin(half_angle) ** 2)

        def integrand(cosine: float) -> float:
            denominator = sine_tilt * np.sqrt(max(0.0, 1 - cosine**2))

            threshold = (np.cos(half_angle) - cosine_tilt * cosine) / denominator

            azimuth = 2 * np.arccos(np.clip(threshold, -1, 1))

            return cosine * azimuth / np.pi

        transitions = [
            np.cos(tilt + half_angle),
            np.cos(tilt - half_angle),
        ]

        transitions = sorted({value for value in transitions if 0 < value < 1})

        fraction, _ = integrate.quad(
            integrand, 0, 1, points=transitions, epsabs=1e-10, epsrel=1e-10
        )

        return float(fraction)

    def _add_to_figure(
        self, *, figure: MomentRT.plotting.Figure, context: Mapping[str, Any]
    ) -> None:
        """Append this instrument geometry to an existing figure in place.

        Parameters
        ----------
        figure : plotly.graph_objects.Figure
            Figure to modify in place.
        context : mapping of str to object
            Coordinate units, legend labels, and physical sample bounds.
        """

        go = MomentRT.plotting.require_plotly()

        position = self.position.to(context["length_unit"]).magnitude

        radius = self.radius.to(context["length_unit"]).magnitude

        ring = MomentRT.plotting._circle_points(
            position=position, radius=radius, face=self.face
        )

        name = context["name"]

        figure.add_trace(
            go.Scatter3d(
                x=ring[:, 0],
                y=ring[:, 1],
                z=ring[:, 2],
                mode="lines",
                line={"color": self.color, "width": 5},
                name=name,
                legendgroup=context.get("legendgroup"),
            )
        )

        MomentRT.plotting._add_footprint_fill(
            figure=figure,
            position=position,
            ring=ring,
            color=self.color,
            legendgroup=context.get("legendgroup"),
        )

        MomentRT.plotting._add_direction_arrow(
            figure=figure,
            start=position,
            direction=self.collection_axis,
            length=2 * radius,
            color=self.color,
            label=f"Collection axis: {name}",
            legendgroup=context.get("legendgroup"),
        )

    def _validate_sample(
        self, *, sample_origin: Quantity, sample_size: Quantity
    ) -> None:
        """Require the detector disk to lie completely on its configured face.

        Parameters
        ----------
        sample_origin : pint.Quantity
            Physical sample lower corner, shape (3,), with length units.
        sample_size : pint.Quantity
            Positive physical sample extents, shape (3,), with length units.
        """

        axis = "xyz".index(self.face[0])

        in_plane = [coordinate for coordinate in range(3) if coordinate != axis]

        lower = sample_origin.to("m").magnitude

        upper = lower + sample_size.to("m").magnitude

        position = self.position.to("m").magnitude

        radius = self.radius.to("m").magnitude

        face_coordinate = lower[axis] if self.face.endswith("min") else upper[axis]

        if not np.isclose(position[axis], face_coordinate, rtol=1e-12, atol=1e-15):
            raise ValueError(
                "Detector position must lie on its configured sample face."
            )

        tolerance = 16 * np.finfo(float).eps * np.maximum(np.abs(lower), np.abs(upper))

        if np.any(
            position[in_plane] - radius < lower[in_plane] - tolerance[in_plane]
        ) or np.any(
            position[in_plane] + radius > upper[in_plane] + tolerance[in_plane]
        ):
            raise ValueError(
                MomentRT.boundary._disk_bounds_message(
                    face=self.face,
                    position=self.position,
                    radius=self.radius,
                    lower=sample_origin,
                    upper=sample_origin + sample_size,
                    label="detector aperture",
                )
            )

    def measure(
        self, *, result: MomentRT.results.SimulationResult
    ) -> MomentRT.results.DetectorResult:
        """Integrate the solved escaping field over the sampled disk aperture.

        Parameters
        ----------
        result : SimulationResult
            Finite-sample solution with escaping-irradiance boundary fields.

        Returns
        -------
        DetectorResult
            Collected power, sampled and exact areas, angular fraction, and
            minimum sampled irradiance. Quadrature area is not renormalized.

        Raises
        ------
        ValueError
            If boundary data are missing or the aperture does not fit on its face.
        DetectorResolutionError
            If no quadrature point samples the aperture; refine the mesh.
        """

        if self.face not in result.boundary_fields:
            raise ValueError(
                "SurfaceDetector requires finite-sample boundary fields from Diffusion, SP3, or P3."
            )

        field = result.boundary_fields[self.face]

        self._validate_sample(
            sample_origin=field.sample_origin,
            sample_size=field.sample_size,
        )

        axis = "xyz".index(self.face[0])

        in_plane = [coordinate for coordinate in range(3) if coordinate != axis]

        position = self.position.to("m").magnitude

        radius = self.radius.to("m").magnitude

        displacement = (
            field.positions.to("m").magnitude[..., in_plane] - position[in_plane]
        )

        selected = np.linalg.norm(displacement, axis=-1) <= radius

        sampled_area = field.area_weights[selected].sum()

        if sampled_area.to("m^2").magnitude == 0:
            raise DetectorResolutionError(
                "The mesh does not resolve the detector aperture; refine it."
            )

        fraction = self.angular_fraction

        power = (
            field.escaping_irradiance[selected] * field.area_weights[selected]
        ).sum() * fraction

        return MomentRT.results.DetectorResult(
            collected_power=power,
            sampled_area=sampled_area,
            geometric_area=np.pi * self.radius**2,
            angular_fraction=fraction,
            wavelength=field.wavelength,
            minimum_escaping_irradiance=field.escaping_irradiance[selected].min(),
            detector=SurfaceDetector(
                face=self.face,
                position=self.position,
                radius=self.radius,
                collection_axis=self.collection_axis,
                acceptance_angle=self.acceptance_angle,
            ),
        )
