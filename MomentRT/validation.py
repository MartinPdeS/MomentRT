"""Non-solving validation of simulation configuration and footprint sampling."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from pint import Quantity

import MomentRT


class ValidationReport:
    """Configuration errors, numerical cautions, and sampled footprint areas.

    Parameters
    ----------
    errors : tuple of str
        Configuration errors that must be corrected before solving.
    warnings : tuple of str
        Numerical cautions that do not prevent a solve.
    source_sampled_area : pint.Quantity or None
        Boundary quadrature area sampled by the source, when available.
    detector_sampled_areas : mapping of str to pint.Quantity
        Boundary quadrature area sampled by each valid detector.

    Notes
    -----
    Passing validation does not establish mesh convergence or physical model
    accuracy. This report is a snapshot and never caches a prepared operator.
    """

    def __init__(
        self,
        *,
        errors: tuple[str, ...],
        warnings: tuple[str, ...],
        source_sampled_area: Quantity | None,
        detector_sampled_areas: Mapping[str, Quantity],
    ) -> None:
        """Store the independently collected validation findings.

        Parameters
        ----------
        errors : tuple of str
            Configuration errors that prevent solving.
        warnings : tuple of str
            Numerical cautions that do not prevent solving.
        source_sampled_area : pint.Quantity or None
            Boundary quadrature area sampled by the source, when available.
        detector_sampled_areas : mapping of str to pint.Quantity
            Sampled boundary quadrature areas for valid detectors.
        """
        self.errors = errors

        self.warnings = warnings

        self.source_sampled_area = source_sampled_area

        self.detector_sampled_areas = dict(detector_sampled_areas)

    @property
    def valid(self) -> bool:
        """Report whether any configuration errors were found.

        Returns
        -------
        bool
            True when no errors were found; warnings may still be present.
        """
        return not self.errors

    def summary(self) -> str:
        """Format the validation findings for display.

        Returns
        -------
        str
            Human-readable findings and sampled areas with physical units.
        """
        lines = ["Configuration valid" if self.valid else "Configuration has errors"]

        lines.extend(f"Error: {error}" for error in self.errors)

        lines.extend(f"Warning: {warning}" for warning in self.warnings)

        if self.source_sampled_area is not None:
            lines.append(f"Source sampled area: {self.source_sampled_area:~P}")

        lines.extend(
            f"Detector {name} sampled area: {area:~P}"
            for name, area in self.detector_sampled_areas.items()
        )

        return "\n".join(lines)


def _validate(
    *, simulation: MomentRT.simulation.Simulation, wavelengths: Quantity | None
) -> ValidationReport:
    """Inspect supported simulation configurations without solving a PDE."""
    errors: list[str] = []

    warnings: list[str] = []

    detector_areas: dict[str, Quantity] = {}

    source_area = None

    medium = simulation.medium

    source = simulation.source

    solver = simulation.solver

    mesh = simulation.mesh

    finite = isinstance(
        solver, (MomentRT.diffusion.Diffusion, MomentRT.sp3.SP3, MomentRT.p3.P3)
    )

    infinite = isinstance(solver, MomentRT.diffusion.InfiniteDiffusion)

    if not finite and not infinite:
        errors.append(
            "Preflight validation supports Diffusion, SP3, P3, and InfiniteDiffusion solvers."
        )

    if finite and not isinstance(medium, MomentRT.medium.StructuredMedium):
        errors.append("Finite solvers require a StructuredMedium.")

    if finite and not isinstance(source, MomentRT.source.UniformSurfaceSource):
        errors.append("Finite solvers require a surface source.")

    if infinite and not isinstance(medium, MomentRT.medium.HomogeneousMedium):
        errors.append("InfiniteDiffusion requires a HomogeneousMedium.")

    if infinite and not isinstance(source, MomentRT.source.PointSource):
        errors.append("InfiniteDiffusion requires a PointSource.")

    if infinite and simulation.detectors:
        errors.append("Contact detectors require finite-sample boundary fields.")

    if not isinstance(mesh, MomentRT.mesh.CartesianMesh):
        errors.append("mesh must be a CartesianMesh.")

    if isinstance(
        source, (MomentRT.source.PointSource, MomentRT.source.UniformSurfaceSource)
    ):
        for name, unit in (("wavelength", "m"), ("optical_power", "W")):
            try:
                MomentRT.units.validate_units(
                    value=getattr(source, name),
                    unit=unit,
                    name=name,
                    shape=(),
                    positive=True,
                )
            except (ValueError, TypeError) as error:
                errors.append(f"Source {name}: {error}")

    selected_wavelengths = []

    if wavelengths is None:
        if hasattr(source, "wavelength"):
            selected_wavelengths = [source.wavelength]
    else:
        try:
            MomentRT.units.validate_units(
                value=wavelengths, unit="m", name="wavelengths", positive=True
            )

            if wavelengths.ndim != 1 or not len(wavelengths):
                raise ValueError("wavelengths must be a nonempty vector.")

            selected_wavelengths = list(wavelengths)
        except (ValueError, TypeError) as error:
            errors.append(f"Wavelengths: {error}")

    if isinstance(medium, MomentRT.medium.StructuredMedium):
        materials = [
            ("background", medium.background),
            *((shape.name, shape.material) for shape in medium.structures),
        ]

        if medium.background is None:
            errors.append("Sample background is not configured.")

        if isinstance(solver, (MomentRT.sp3.SP3, MomentRT.p3.P3)):
            for name, material in materials:
                if material is not None and material.anisotropy is None:
                    errors.append(
                        f"{type(solver).__name__} requires anisotropy for material {name!r}."
                    )

                if (
                    material is not None
                    and material.phase_function != "henyey_greenstein"
                ):
                    errors.append(
                        f"{type(solver).__name__} does not support the phase function for material {name!r}."
                    )

        for wavelength in selected_wavelengths:
            indices = []

            for name, material in materials:
                if material is None:
                    continue

                try:
                    resolved = material.at(wavelength=wavelength)

                    indices.append(resolved.refractive_index)

                    if (
                        isinstance(solver, (MomentRT.sp3.SP3, MomentRT.p3.P3))
                        and resolved.anisotropy is not None
                    ):
                        MomentRT.material._anisotropy_values(
                            value=resolved.anisotropy, shape=()
                        )
                except (ValueError, TypeError) as error:
                    errors.append(f"Material {name!r} at {wavelength:~P}: {error}")

            if indices and any(index != indices[0] for index in indices):
                errors.append(
                    f"Internal refractive-index jumps are unsupported at {wavelength:~P}."
                )
    elif isinstance(medium, MomentRT.medium.HomogeneousMedium):
        for wavelength in selected_wavelengths:
            try:
                medium.at(wavelength=wavelength)
            except (ValueError, TypeError) as error:
                errors.append(f"Material at {wavelength:~P}: {error}")

    if isinstance(
        solver, (MomentRT.diffusion.Diffusion, MomentRT.sp3.SP3, MomentRT.p3.P3)
    ):
        try:
            MomentRT.diffusion.Diffusion(boundaries=solver.boundaries)

            if isinstance(solver, MomentRT.sp3.SP3):
                MomentRT.sp3.SP3(boundaries=solver.boundaries)

            if isinstance(solver, MomentRT.p3.P3):
                MomentRT.p3.P3(boundaries=solver.boundaries)
        except (ValueError, TypeError, AttributeError) as error:
            errors.append(f"Boundaries: {error}")

        if isinstance(source, MomentRT.source.UniformSurfaceSource) and not isinstance(
            solver.boundaries.get(source.face), MomentRT.boundary.RobinBoundary
        ):
            errors.append("The illuminated face requires a RobinBoundary.")

    if (
        finite
        and isinstance(medium, MomentRT.medium.StructuredMedium)
        and isinstance(mesh, MomentRT.mesh.CartesianMesh)
    ):
        mesh_lower = mesh.origin.to("m").magnitude

        mesh_upper = (
            mesh_lower + (np.asarray(mesh.shape) - 1) * mesh.spacing.to("m").magnitude
        )

        lower = medium.origin.to("m").magnitude

        upper = (medium.origin + medium.size).to("m").magnitude

        mesh_valid = np.allclose(
            mesh_lower, lower, rtol=1e-12, atol=1e-15
        ) and np.allclose(mesh_upper, upper, rtol=1e-12, atol=1e-15)

        if not mesh_valid:
            errors.append("The mesh must span the finite sample.")

        valid_detectors = {}

        for name, detector in simulation.detectors.items():
            try:
                detector._validate_sample(
                    sample_origin=medium.origin, sample_size=medium.size
                )

                valid_detectors[name] = detector
            except (ValueError, TypeError) as error:
                errors.append(f"Detector {name!r}: {error}")

        source_valid = isinstance(source, MomentRT.source.UniformSurfaceSource)

        if isinstance(source, MomentRT.source.SurfaceSource):
            try:
                axis = "xyz".index(source.face[0])

                face_coordinate = (
                    lower[axis] if source.face.endswith("min") else upper[axis]
                )

                if not np.isclose(
                    source.position.to("m").magnitude[axis],
                    face_coordinate,
                    rtol=1e-12,
                    atol=1e-15,
                ):
                    raise ValueError("Source position must lie on its configured face.")

                MomentRT.boundary._validate_face_disk(
                    sample=medium,
                    face=source.face,
                    position=source.position,
                    radius=source.radius,
                )
            except (ValueError, TypeError) as error:
                source_valid = False

                errors.append(f"Source: {error}")

        if mesh_valid:
            geometry = MomentRT.fem._MeshGeometry(mesh=mesh)

            if source_valid and isinstance(
                source, MomentRT.source.UniformSurfaceSource
            ):
                try:
                    quadrature, _, weights = geometry.quadrature[source.face]

                    _, sampled = MomentRT.fem.source_irradiance(
                        quadrature=quadrature, weights=weights, source=source
                    )

                    source_area = sampled * MomentRT.units.ureg.meter**2

                    axis = "xyz".index(source.face[0])

                    geometric = (
                        np.pi * source.radius**2
                        if isinstance(source, MomentRT.source.SurfaceSource)
                        else np.prod(
                            medium.size[[index for index in range(3) if index != axis]]
                        )
                    )

                    if abs((source_area / geometric).to("").magnitude - 1) > 0.02:
                        warnings.append(
                            "Source sampled area differs from geometric area by more than 2%; refine the mesh."
                        )
                except (ValueError, TypeError, KeyError) as error:
                    errors.append(f"Source sampling: {error}")

            for name, detector in valid_detectors.items():
                quadrature, _, weights = geometry.quadrature[detector.face]

                axis = "xyz".index(detector.face[0])

                in_plane = [index for index in range(3) if index != axis]

                distance = np.linalg.norm(
                    quadrature[..., in_plane]
                    - detector.position.to("m").magnitude[in_plane],
                    axis=-1,
                )

                sampled = float(
                    np.sum(weights[distance <= detector.radius.to("m").magnitude])
                )

                detector_areas[name] = sampled * MomentRT.units.ureg.meter**2

                if sampled == 0:
                    errors.append(f"Detector {name!r} is not sampled; refine the mesh.")
                elif (
                    abs(sampled / (np.pi * detector.radius.to("m").magnitude ** 2) - 1)
                    > 0.02
                ):
                    warnings.append(
                        f"Detector {name!r} sampled area differs from geometric area by more than 2%; refine the mesh."
                    )

            warnings.append(
                "Footprint sampling alone does not establish mesh convergence or transport accuracy."
            )

    if (
        infinite
        and isinstance(source, MomentRT.source.PointSource)
        and isinstance(mesh, MomentRT.mesh.CartesianMesh)
    ):
        if np.any(
            np.linalg.norm(
                mesh.points.to("m").magnitude - source.position.to("m").magnitude,
                axis=-1,
            )
            == 0
        ):
            errors.append("A mesh node coincides with the singular point source.")

    return ValidationReport(
        errors=tuple(errors),
        warnings=tuple(warnings),
        source_sampled_area=source_area,
        detector_sampled_areas=detector_areas,
    )
