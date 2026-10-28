"""Contact-aperture detection with an explicit Lambertian escape assumption."""

import numpy as np
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
    """

    def __init__(self, *, face, position, radius, collection_axis, acceptance_angle):
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

        self.face = face

        self.position = position.copy()

        self.radius = radius.copy()

        self.collection_axis = direction.copy()

        self.acceptance_angle = acceptance_angle.copy()

    @property
    def angular_fraction(self):
        """Fraction of Lambertian escaping irradiance inside the acceptance cone."""

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

        def integrand(cosine):
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

    def _add_to_figure(self, *, figure, context):
        go = MomentRT.plotting.require_plotly()

        position = self.position.to(context["length_unit"]).magnitude

        radius = self.radius.to(context["length_unit"]).magnitude

        ring = MomentRT.plotting._circle_points(
            position=position, radius=radius, face=self.face
        )

        name = context.get("name") or "Surface detector"

        figure.add_trace(
            go.Scatter3d(
                x=ring[:, 0],
                y=ring[:, 1],
                z=ring[:, 2],
                mode="lines",
                line={"color": "royalblue", "width": 5},
                name=name,
                legendgroup=context.get("legendgroup"),
            )
        )

        tip = position + 2 * radius * self.collection_axis

        figure.add_trace(
            go.Scatter3d(
                x=[position[0], tip[0]],
                y=[position[1], tip[1]],
                z=[position[2], tip[2]],
                mode="lines+markers",
                line={"color": "royalblue", "width": 3},
                marker={"size": 3},
                showlegend=False,
                hovertemplate=f"Collection axis: {name}<extra></extra>",
                legendgroup=context.get("legendgroup"),
            )
        )

    def _validate_sample(self, *, sample_origin, sample_size):
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
                "The entire detector aperture must fit on its sample face."
            )

    def measure(self, *, result):
        """Integrate the solved escaping field over the sampled disk aperture."""

        if self.face not in result.boundary_fields:
            raise ValueError(
                "SurfaceDetector requires finite-sample boundary fields from Diffusion."
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
