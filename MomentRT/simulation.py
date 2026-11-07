"""Composition of medium, illumination, and transport."""

from __future__ import annotations
import copy
import time
from collections.abc import Iterable, Mapping
from typing import Any, Protocol, cast
import numpy as np
import TypedUnit
from pint import Quantity
import MomentRT


class TransportSolver(Protocol):
    """Structural interface for solvers accepted by :class:`Simulation`.

    Implementations choose their supported medium and source types and validate
    those configurations at runtime. No inheritance from this protocol is needed.
    """

    def solve(
        self,
        *,
        medium: Any,
        source: Any,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> (
        MomentRT.results.SimulationResult
        | MomentRT.results.TransientResult
        | MomentRT.fields.ScalarField
    ):
        """Solve the transport problem on a mesh supported by the implementation.

        Parameters
        ----------
        medium : HomogeneousMedium or StructuredMedium
            Physical configuration supported by the implementation.
        source : PointSource or UniformSurfaceSource
            Illumination supported by the implementation, including groups of
            simultaneous surface sources.
        mesh : CartesianMesh or TetrahedralMesh
            Sampling nodes or finite-element grid.

        Returns
        -------
        SimulationResult, TransientResult, or ScalarField
            Transport solution, or a fluence field wrapped by Simulation.run.
        """
        ...


class Simulation:
    """Compose physical geometry, illumination, a transport solver, and detectors.

    Parameters
    ----------
    medium : HomogeneousMedium or StructuredMedium
        Physical material configuration, compatible with the chosen solver.
    source : PointSource, UniformSurfaceSource, or iterable of SurfaceSource
        Illumination configuration. An iterable of disk sources is converted
        automatically to a SurfaceSourceGroup for one simultaneous solve.
    mesh : CartesianMesh or TetrahedralMesh
        Nodal coordinates for sampling or finite-element discretization.
    solver : TransportSolver
        Transport model implementing the keyword-only solve protocol.
    detectors : mapping of str to SurfaceDetector, or None, optional
        Nonempty names mapped to passive contact detectors. None disables
        detector measurements.
        Default is None.
    """

    def __init__(
        self,
        *,
        medium: MomentRT.medium.HomogeneousMedium | MomentRT.medium.StructuredMedium,
        source: (
            MomentRT.source.PointSource
            | MomentRT.source.UniformSurfaceSource
            | Iterable[MomentRT.source.SurfaceSource]
        ),
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
        solver: MomentRT.simulation.TransportSolver,
        detectors: Mapping[str, MomentRT.detector.SurfaceDetector] | None = None,
    ) -> None:
        """Initialize the Simulation configuration.

        Parameters
        ----------
        medium : HomogeneousMedium or StructuredMedium
            Physical material configuration, compatible with the chosen solver.
        source : PointSource, UniformSurfaceSource, or iterable of SurfaceSource
            Illumination configuration. An iterable of disk sources is
            converted automatically to a SurfaceSourceGroup.
        mesh : CartesianMesh or TetrahedralMesh
            Nodal coordinates for sampling or finite-element discretization.
        solver : TransportSolver
            Transport model implementing the keyword-only solve protocol.
        detectors : mapping of str to SurfaceDetector, or None, optional
            Nonempty names mapped to passive contact detectors. None disables
            detector measurements.
            Default is None.
        """

        if detectors is not None and not isinstance(detectors, Mapping):
            raise TypeError("detectors must map names to SurfaceDetector objects.")

        if not isinstance(
            source, (MomentRT.source.PointSource, MomentRT.source.UniformSurfaceSource)
        ):
            if isinstance(source, (str, bytes, Mapping)):
                raise TypeError(
                    "source must be a source object or an iterable of SurfaceSource objects."
                )

            try:
                source = MomentRT.source.SurfaceSourceGroup(sources=source)
            except TypeError as error:
                raise TypeError(
                    "source must be a source object or an iterable of SurfaceSource objects."
                ) from error

        detectors = dict(detectors or {})

        if any(
            not isinstance(name, str)
            or not name.strip()
            or not isinstance(detector, MomentRT.detector.SurfaceDetector)
            for name, detector in detectors.items()
        ):
            raise ValueError(
                "detectors must map nonempty names to SurfaceDetector objects."
            )

        self.medium = medium

        self.source = source

        self.mesh = mesh

        self.solver = solver

        self.detectors = detectors

    def validate(
        self, *, wavelengths: Quantity | None = None
    ) -> MomentRT.validation.ValidationReport:
        """Inspect configuration, material coverage, and mesh footprint sampling.

        Parameters
        ----------
        wavelengths : pint.Quantity or None, optional
            One-dimensional sweep wavelengths with length units. None checks
            the monochromatic source wavelength or the complete broadband
            source grid.

        Returns
        -------
        ValidationReport
            Independent errors, cautions, solver assumptions and limits, and
            sampled source and detector areas. Validation does not solve
            transport or establish mesh convergence.
        """
        if self.source.spectrum is not None:
            report = self._reference_simulation().validate(
                wavelengths=(
                    self.source.spectrum.wavelengths
                    if wavelengths is None
                    else wavelengths
                )
            )
            if isinstance(self.solver, MomentRT.transport.transient_p1.TransientP1):
                report.errors = (
                    *report.errors,
                    "TransientP1 requires a monochromatic source.",
                )
            return report
        return self._validate_configuration(wavelengths=wavelengths)

    def _validate_configuration(
        self,
        *,
        wavelengths: Quantity | None,
    ) -> MomentRT.validation.ValidationReport:
        """Inspect supported simulation configurations without solving a PDE."""
        errors: list[str] = []

        warnings: list[str] = []

        detector_areas: dict[str, Quantity] = {}

        source_area = None

        medium = self.medium

        source = self.source

        solver = self.solver

        mesh = self.mesh

        finite = isinstance(
            solver,
            (
                MomentRT.transport.diffusion.Diffusion,
                MomentRT.transport.sp3.SP3,
                MomentRT.transport.p3.P3,
                MomentRT.transport.transient_p1.TransientP1,
            ),
        )

        infinite = isinstance(solver, MomentRT.transport.diffusion.InfiniteDiffusion)

        if not finite and not infinite:
            errors.append(
                "Preflight validation supports Diffusion, SP3, P3, TransientP1, and InfiniteDiffusion solvers."
            )

        if finite and not isinstance(medium, MomentRT.medium.StructuredMedium):
            errors.append("Finite solvers require a StructuredMedium.")

        if finite and not isinstance(source, MomentRT.source.UniformSurfaceSource):
            errors.append("Finite solvers require a surface source.")

        if infinite and not isinstance(medium, MomentRT.medium.HomogeneousMedium):
            errors.append("InfiniteDiffusion requires a HomogeneousMedium.")

        if infinite and not isinstance(source, MomentRT.source.PointSource):
            errors.append("InfiniteDiffusion requires a PointSource.")

        if infinite and self.detectors:
            errors.append("Contact detectors require finite-sample boundary fields.")

        if not isinstance(
            mesh, (MomentRT.mesh.CartesianMesh, MomentRT.mesh.TetrahedralMesh)
        ):
            errors.append("mesh must be a CartesianMesh or TetrahedralMesh.")

        if infinite and not isinstance(mesh, MomentRT.mesh.CartesianMesh):
            errors.append("InfiniteDiffusion requires a Cartesian observation mesh.")

        if isinstance(
            source, (MomentRT.source.PointSource, MomentRT.source.UniformSurfaceSource)
        ):
            for name, unit in (("wavelength", "m"), ("optical_power", "W")):
                try:
                    TypedUnit.validate_quantity(
                        registry=MomentRT.units.ureg,
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
                TypedUnit.validate_quantity(
                    registry=MomentRT.units.ureg,
                    value=wavelengths,
                    unit="m",
                    name="wavelengths",
                    positive=True,
                )

                if wavelengths.ndim != 1 or not len(wavelengths):
                    raise ValueError("wavelengths must be a nonempty vector.")

                selected_wavelengths = list(wavelengths)
            except (ValueError, TypeError) as error:
                errors.append(f"Wavelengths: {error}")

        for name, detector in self.detectors.items():
            if detector.responsivity is not None:
                for wavelength in selected_wavelengths:
                    try:
                        detector.responsivity.at(wavelength=wavelength)
                    except (ValueError, TypeError) as error:
                        errors.append(f"Detector {name!r} responsivity: {error}")

        if isinstance(medium, MomentRT.medium.StructuredMedium):
            materials = [
                ("background", medium.background),
                *((shape.name, shape.material) for shape in medium.structures),
            ]

            if medium.background is None:
                errors.append("Sample background is not configured.")

            if isinstance(
                solver, (MomentRT.transport.sp3.SP3, MomentRT.transport.p3.P3)
            ):
                for name, material in materials:
                    if (
                        material is not None
                        and material.anisotropy is None
                        and isinstance(
                            material.phase_function,
                            MomentRT.phase_function.HenyeyGreenstein,
                        )
                    ):
                        errors.append(
                            f"{type(solver).__name__} requires anisotropy for material {name!r}."
                        )

            for wavelength in selected_wavelengths:
                indices = []

                for name, material in materials:
                    if material is None:
                        continue

                    try:
                        resolved = material.at(wavelength=wavelength)

                        indices.append((name, resolved.refractive_index))

                        if isinstance(
                            solver,
                            (MomentRT.transport.sp3.SP3, MomentRT.transport.p3.P3),
                        ) and (
                            resolved.anisotropy is not None
                            or not isinstance(
                                resolved.phase_function,
                                MomentRT.phase_function.HenyeyGreenstein,
                            )
                        ):
                            resolved._phase_moments()
                    except (ValueError, TypeError) as error:
                        errors.append(f"Material {name!r} at {wavelength:~P}: {error}")

                if indices and any(index != indices[0][1] for _, index in indices):
                    details = ", ".join(
                        f"{name!r}={index:g}" for name, index in indices
                    )

                    errors.append(
                        f"Internal refractive-index jumps are unsupported at "
                        f"{wavelength:~P}: {details}."
                    )
        elif isinstance(medium, MomentRT.medium.HomogeneousMedium):
            for wavelength in selected_wavelengths:
                try:
                    medium.at(wavelength=wavelength)
                except (ValueError, TypeError) as error:
                    errors.append(f"Material at {wavelength:~P}: {error}")

        if isinstance(
            solver,
            (
                MomentRT.transport.diffusion.Diffusion,
                MomentRT.transport.sp3.SP3,
                MomentRT.transport.p3.P3,
            ),
        ):
            try:
                MomentRT.transport.diffusion.Diffusion(boundaries=solver.boundaries)

                if isinstance(solver, MomentRT.transport.sp3.SP3):
                    MomentRT.transport.sp3.SP3(boundaries=solver.boundaries)

                if isinstance(solver, MomentRT.transport.p3.P3):
                    MomentRT.transport.p3.P3(boundaries=solver.boundaries)
            except (ValueError, TypeError, AttributeError) as error:
                errors.append(f"Boundaries: {error}")

            if isinstance(
                source, MomentRT.source.UniformSurfaceSource
            ) and not isinstance(
                solver.boundaries.get(source.face), MomentRT.boundary.RobinBoundary
            ):
                errors.append(
                    f"The illuminated face {source.face!r} requires a RobinBoundary."
                )

        if (
            finite
            and isinstance(medium, MomentRT.medium.StructuredMedium)
            and isinstance(
                mesh, (MomentRT.mesh.CartesianMesh, MomentRT.mesh.TetrahedralMesh)
            )
        ):
            mesh_lower, mesh_upper = mesh.bounds.to("m").magnitude

            lower = medium.origin.to("m").magnitude

            upper = (medium.origin + medium.size).to("m").magnitude

            mesh_valid = np.allclose(
                mesh_lower, lower, rtol=1e-12, atol=1e-15
            ) and np.allclose(mesh_upper, upper, rtol=1e-12, atol=1e-15)

            if not mesh_valid:
                errors.append("The mesh must span the finite sample.")

            if isinstance(mesh, MomentRT.mesh.TetrahedralMesh):
                try:
                    mesh._validate_sample(sample=medium)
                except ValueError as error:
                    errors.append(str(error))
                    mesh_valid = False

            valid_detectors = {}

            for name, detector in self.detectors.items():
                try:
                    detector._validate_sample(
                        sample_origin=medium.origin, sample_size=medium.size
                    )

                    valid_detectors[name] = detector
                except (ValueError, TypeError) as error:
                    errors.append(f"Detector {name!r}: {error}")

            source_valid = isinstance(source, MomentRT.source.UniformSurfaceSource)

            if isinstance(source, MomentRT.source.SurfaceSourceGroup):
                for member in source.sources:
                    try:
                        axis = "xyz".index(member.face[0])

                        face_coordinate = (
                            lower[axis] if member.face.endswith("min") else upper[axis]
                        )

                        if not np.isclose(
                            member.position.to("m").magnitude[axis],
                            face_coordinate,
                            rtol=1e-12,
                            atol=1e-15,
                        ):
                            raise ValueError(
                                "Source position must lie on its configured face."
                            )

                        MomentRT._face_geometry.validate_face_disk(
                            sample=medium,
                            face=member.face,
                            position=member.position,
                            radius=member.radius,
                        )
                    except (ValueError, TypeError) as error:
                        source_valid = False

                        errors.append(f"Source {member.name!r}: {error}")

            elif isinstance(source, MomentRT.source.SurfaceSource):
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
                        raise ValueError(
                            "Source position must lie on its configured face."
                        )

                    MomentRT._face_geometry.validate_face_disk(
                        sample=medium,
                        face=source.face,
                        position=source.position,
                        radius=source.radius,
                    )
                except (ValueError, TypeError) as error:
                    source_valid = False

                    errors.append(f"Source: {error}")

            if mesh_valid:
                geometry = MomentRT.fem.FiniteElementGeometry(mesh=mesh)

                if source_valid and isinstance(
                    source, MomentRT.source.UniformSurfaceSource
                ):
                    try:
                        _, sampled = geometry.source_irradiance(source=source)

                        source_area = sampled * MomentRT.units.ureg.meter**2

                        axis = "xyz".index(source.face[0])

                        if isinstance(source, MomentRT.source.SurfaceSourceGroup):
                            geometric = sum(
                                np.pi * member.radius**2 for member in source.sources
                            )
                        elif isinstance(source, MomentRT.source.SurfaceSource):
                            geometric = np.pi * source.radius**2
                        else:
                            geometric = np.prod(
                                medium.size[
                                    [index for index in range(3) if index != axis]
                                ]
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
                        errors.append(
                            f"Detector {name!r} is not sampled; refine the mesh."
                        )
                    elif (
                        abs(
                            sampled / (np.pi * detector.radius.to("m").magnitude ** 2)
                            - 1
                        )
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

        return MomentRT.validation.ValidationReport(
            errors=tuple(errors),
            warnings=tuple(warnings),
            source_sampled_area=source_area,
            detector_sampled_areas=detector_areas,
            model_notes=tuple(getattr(solver, "model_notes", ())),
        )

    def run(
        self,
    ) -> (
        MomentRT.results.SimulationResult
        | MomentRT.results.BroadbandSimulationResult
        | MomentRT.results.TransientResult
    ):
        """Solve the configured transport problem and evaluate all named detectors.

        Returns
        -------
        SimulationResult, BroadbandSimulationResult, or TransientResult
            Monochromatic fields and measurements, or broadband detector spectra,
            trapezoidal optical and electrical integrals, and reference solves.
            Broadband integration uses exactly the source spectrum's grid and
            endpoints; refine that grid to assess spectral quadrature error.
        """

        if isinstance(self.solver, MomentRT.transport.transient_p1.TransientP1):
            if not isinstance(self.medium, MomentRT.medium.StructuredMedium):
                raise TypeError("TransientP1 requires a finite StructuredMedium.")
            if not isinstance(self.source, MomentRT.source.UniformSurfaceSource):
                raise TypeError("TransientP1 requires a surface source.")
            if self.source.spectrum is not None:
                raise ValueError("TransientP1 requires a monochromatic source.")
            for detector in self.detectors.values():
                if detector.responsivity is not None:
                    detector.responsivity.at(wavelength=self.source.wavelength)
            result = self.solver.solve(
                medium=self.medium,
                source=self.source,
                mesh=self.mesh,
            )
            for frame in result.frames:
                frame.detectors.update(
                    {
                        name: detector.measure(result=frame)
                        for name, detector in self.detectors.items()
                    }
                )
            return result

        if self.source.spectrum is not None:
            started = time.perf_counter()
            spectrum = copy.deepcopy(self.source.spectrum)
            reference = self._reference_simulation()
            broadband = MomentRT.results.BroadbandSimulationResult(
                spectrum=spectrum,
                reference_results=reference.sweep(wavelengths=spectrum.wavelengths),
                reference_optical_power=reference.source.optical_power,
            )
            broadband.diagnostics["total_time"] = (
                time.perf_counter() - started
            ) * MomentRT.units.ureg.second
            return broadband
        return self._run_monochromatic()

    def _reference_simulation(self) -> MomentRT.simulation.Simulation:
        """Copy broadband geometry with explicit one-watt monochromatic emission."""
        spectrum = self.source.spectrum
        assert spectrum is not None
        reference = copy.copy(self)
        reference.source = copy.copy(self.source)
        reference.source.spectrum = None
        reference.source.wavelength = spectrum.wavelengths[0].copy()
        reference.source.optical_power = 1 * MomentRT.units.ureg.watt
        return reference

    def _run_monochromatic(self) -> MomentRT.results.SimulationResult:
        """Solve a monochromatic configuration and measure every detector."""
        for detector in self.detectors.values():
            if detector.responsivity is not None:
                detector.responsivity.at(wavelength=self.source.wavelength)

        solution = self.solver.solve(
            medium=self.medium, source=self.source, mesh=self.mesh
        )

        if isinstance(solution, MomentRT.results.TransientResult):
            raise TypeError("Transient solvers must be run through Simulation.run().")
        if not isinstance(solution, MomentRT.results.SimulationResult):
            solution = MomentRT.results.SimulationResult(fluence=solution)

        measurements = {
            name: detector.measure(result=solution)
            for name, detector in self.detectors.items()
        }

        solution.detectors.update(measurements)

        return solution

    def sweep(
        self, *, wavelengths: Quantity, optical_power: Quantity | None = None
    ) -> MomentRT.results.SpectralSimulationResult:
        """Solve wavelengths in supplied order without mutating this configuration.

        optical_power may be a scalar or a same-length vector of powers. When
        omitted, every solve uses the configured source's optical_power.
        This is a set of monochromatic solves, not broadband spectral integration.

        Parameters
        ----------
        wavelengths : pint.Quantity
            Nonempty one-dimensional vector of positive vacuum wavelengths.
        optical_power : pint.Quantity or None, optional
            Positive scalar power or vector matching wavelengths. None retains
            the configured source power for every monochromatic solve.
            Default is None.

        Returns
        -------
        SpectralSimulationResult
            Ordered independent monochromatic solutions. Finite sweeps reuse
            geometry but rebuild the wavelength-dependent operator and factorization.

        Notes
        -----
        These are independent monochromatic powers, not a broadband spectral
        integral. Finite sweeps reuse geometry, material labels, source footprint,
        and boundary quadrature once per call. Each wavelength gets a new bulk
        operator and sparse factorization.
        """

        started = time.perf_counter()

        if isinstance(self.solver, MomentRT.transport.transient_p1.TransientP1):
            raise ValueError("TransientP1 does not support wavelength sweeps.")

        if self.source.spectrum is not None:
            raise ValueError(
                "sweep uses independent monochromatic powers; use run() for a source spectrum."
            )

        TypedUnit.validate_quantity(
            registry=MomentRT.units.ureg,
            value=wavelengths,
            unit="m",
            name="wavelengths",
            positive=True,
        )

        if wavelengths.ndim != 1 or not len(wavelengths):
            raise ValueError("wavelengths must be a nonempty vector.")

        if optical_power is not None:
            TypedUnit.validate_quantity(
                registry=MomentRT.units.ureg,
                value=optical_power,
                unit="W",
                name="optical_power",
                positive=True,
            )

            if optical_power.shape not in ((), wavelengths.shape):
                raise ValueError(
                    "optical_power must be scalar or shaped like wavelengths."
                )

        # Resolve every table first, so an out-of-range sweep fails before solving.
        wavelength: Quantity

        for wavelength in wavelengths:
            self.medium.at(wavelength=wavelength)
            for detector in self.detectors.values():
                if detector.responsivity is not None:
                    detector.responsivity.at(wavelength=wavelength)

        prepare = getattr(self.solver, "_prepare", None)

        preparation_started = time.perf_counter()

        solver = (
            self.solver
            if prepare is None
            else prepare(medium=self.medium, source=self.source, mesh=self.mesh)
        )

        preparation_time = time.perf_counter() - preparation_started

        results = []

        solve_times = []

        for index, wavelength in enumerate(wavelengths):
            source = copy.copy(self.source)

            source.wavelength = wavelength.copy()

            if isinstance(source, MomentRT.source.SurfaceSourceGroup):
                original_power = source.optical_power.copy()

                source.sources = tuple(copy.copy(member) for member in source.sources)

                for member in source.sources:
                    member.wavelength = wavelength.copy()

            if optical_power is not None:
                source.optical_power = (
                    optical_power if optical_power.shape == () else optical_power[index]
                ).copy()

                if isinstance(source, MomentRT.source.SurfaceSourceGroup):
                    ratio = source.optical_power / original_power

                    for member in source.sources:
                        member.optical_power *= ratio

            simulation = Simulation(
                medium=self.medium,
                source=source,
                mesh=self.mesh,
                solver=solver,
                detectors=self.detectors,
            )

            solve_started = time.perf_counter()

            results.append(simulation._run_monochromatic())

            solve_times.append(time.perf_counter() - solve_started)

        return MomentRT.results.SpectralSimulationResult(
            wavelengths=wavelengths,
            results=results,
            diagnostics={
                "geometry_preparations": int(prepare is not None),
                "preparation_time": preparation_time * MomentRT.units.ureg.second,
                "solve_times": np.array(solve_times) * MomentRT.units.ureg.second,
                "total_time": (time.perf_counter() - started)
                * MomentRT.units.ureg.second,
            },
        )

    def mesh_convergence(
        self,
        *,
        mesh_shapes: Iterable[tuple[int, int, int]] | None = None,
        mesh_resolutions: Iterable[Quantity] | None = None,
        contact_resolution: Quantity | None = None,
        contact_distance: Quantity | None = None,
        power_rtol: float = 0.01,
        area_rtol: float = 0.02,
        power_atol: Quantity | None = None,
        required_consecutive: int = 2,
    ) -> MomentRT.results.DetectorConvergenceResult:
        """Compare detector signals and footprint sampling across explicit meshes.

        Parameters
        ----------
        mesh_shapes : iterable of tuple of int or None, optional
            At least two coarse-to-fine Cartesian node shapes. Each shape must
            refine at least one axis and coarsen none. Supply this or
            ``mesh_resolutions``, but not both. Default is None.
        mesh_resolutions : iterable of pint.Quantity or None, optional
            At least two strictly decreasing global tetrahedral target sizes.
            Supply this or ``mesh_shapes``, but not both. Default is None.
        contact_resolution : pint.Quantity or None, optional
            Local target element size near disk sources and detectors. Required
            with ``contact_distance``; applies only to tetrahedral studies.
            Default is None.
        contact_distance : pint.Quantity or None, optional
            Positive distance defining each local contact refinement zone.
            Default is None.
        power_rtol : float, optional
            Nonnegative relative power tolerance against adjacent and finest meshes.
            Default is 0.01.
        area_rtol : float, optional
            Nonnegative relative tolerance of sampled source and detector areas
            against their geometric areas.
            Default is 0.02.
        power_atol : pint.Quantity or None, optional
            Nonnegative absolute power tolerance. None means zero watts.
            Default is None.
        required_consecutive : int, optional
            Positive number of consecutive passing refinements required at the
            fine end of the mesh sequence.
            Default is 2.

        Returns
        -------
        DetectorConvergenceResult
            Per-level fields, powers, areas, failures, and a mesh recommendation
            only when the final refinement sequence meets all requested criteria.
        """

        if isinstance(self.solver, MomentRT.transport.transient_p1.TransientP1):
            raise ValueError(
                "mesh_convergence currently supports steady-state detector results only."
            )

        return MomentRT.convergence._DetectorConvergenceStudy(simulation=self).run(
            mesh_shapes=mesh_shapes,
            mesh_resolutions=mesh_resolutions,
            contact_resolution=contact_resolution,
            contact_distance=contact_distance,
            power_rtol=power_rtol,
            area_rtol=area_rtol,
            power_atol=power_atol,
            required_consecutive=required_consecutive,
        )

    # Existing callers retain the same signature and behavior.
    convergence_study = mesh_convergence

    def spectral_convergence(
        self,
        *,
        wavelength_grids: Iterable[Quantity],
        include_table_knots: bool = True,
        power_rtol: float = 0.01,
        photocurrent_rtol: float = 0.01,
        power_atol: Quantity | None = None,
        photocurrent_atol: Quantity | None = None,
        required_consecutive: int = 2,
    ) -> MomentRT.results.SpectralConvergenceResult:
        """Check broadband power and current stability across wavelength grids.

        Parameters
        ----------
        wavelength_grids : iterable of pint.Quantity
            At least two nested coarse-to-fine wavelength vectors spanning the
            original source interval. Source knots are always retained.
        include_table_knots : bool, optional
            Include material and responsivity knots inside the interval. The
            resulting inspectable grids must still add samples. Default is True.
        power_rtol : float, optional
            Nonnegative relative optical-power tolerance. Default is 0.01.
        photocurrent_rtol : float, optional
            Nonnegative relative current tolerance for detectors with response.
            Default is 0.01.
        power_atol : pint.Quantity or None, optional
            Nonnegative absolute optical-power tolerance. None means zero watts.
            Default is None.
        photocurrent_atol : pint.Quantity or None, optional
            Nonnegative absolute current tolerance. None means zero amperes.
            Default is None.
        required_consecutive : int, optional
            Consecutive passing refinements required at the fine end. Default is 2.

        Returns
        -------
        SpectralConvergenceResult
            Per-level spectra, measurements, changes, and failures. A wavelength
            grid is recommended only after a stable final refinement sequence.

        Notes
        -----
        Geometry, source data, and detector response remain fixed. This checks
        wavelength integration; assess spatial mesh and model accuracy separately.
        """
        return MomentRT.convergence._SpectralConvergenceStudy(simulation=self).run(
            wavelength_grids=wavelength_grids,
            include_table_knots=include_table_knots,
            power_rtol=power_rtol,
            photocurrent_rtol=photocurrent_rtol,
            power_atol=power_atol,
            photocurrent_atol=photocurrent_atol,
            required_consecutive=required_consecutive,
        )

    def plot_system(
        self,
        *,
        wavelength: Quantity | None = None,
        length_unit: str = "mm",
        show: bool = True,
    ) -> MomentRT.plotting.Figure:
        """Display the physical sample, source, and detector geometry in 3D.

        Parameters
        ----------
        wavelength : pint.Quantity or None, optional
            Display wavelength for material properties. None uses the source
            wavelength for monochromatic sources. Broadband sources with
            spectral materials require an explicit display wavelength.
        length_unit : str, optional
            Length unit used for displayed or exported coordinates.
            Default is 'mm'.
        show : bool, optional
            Whether to display the figure immediately. The figure is always returned.
            Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Sample construction geometry with source and detector overlays.
        """

        if wavelength is None and self.source.spectrum is None:
            wavelength = self.source.wavelength

        if wavelength is None:
            medium = self.medium
        else:
            medium = self.medium.at(wavelength=wavelength)

        if isinstance(medium, MomentRT.medium.HomogeneousMedium) and not isinstance(
            self.mesh, MomentRT.mesh.CartesianMesh
        ):
            raise TypeError("Infinite observation windows require a Cartesian mesh.")

        figure = (
            medium.plot(mesh=self.mesh, length_unit=length_unit, show=False)
            if isinstance(medium, MomentRT.medium.StructuredMedium)
            else medium.plot(
                mesh=cast(MomentRT.mesh.CartesianMesh, self.mesh),
                length_unit=length_unit,
                show=False,
            )
        )

        context = {
            "length_unit": length_unit,
            "sample_origin": getattr(medium, "origin", None),
            "sample_size": getattr(medium, "size", None),
        }

        MomentRT.plotting._add_overlays(
            figure=figure,
            overlays=[self.source],
            context=context,
        )

        MomentRT.plotting._add_overlays(
            figure=figure,
            overlays=list(self.detectors.values()),
            context=context,
        )

        if show:
            figure.show()

        return figure
