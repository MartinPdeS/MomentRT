"""Deterministic light transport using angular-moment methods."""

from __future__ import annotations
import MomentRT._version
import MomentRT.boundary
import MomentRT.convergence as convergence  # noqa: F401 -- public submodule
import MomentRT.detector
import MomentRT.fem as fem  # noqa: F401 -- public submodule
import MomentRT.fields
import MomentRT.geometry
import MomentRT.material
import MomentRT.materials
import MomentRT.medium
import MomentRT.mesh
import MomentRT.phase_function
import MomentRT.plotting as plotting  # noqa: F401 -- public submodule
import MomentRT.results
import MomentRT.sections
import MomentRT.simulation
import MomentRT.source
import MomentRT.transport
import MomentRT.transport.diffusion
import MomentRT.transport.p3
import MomentRT.transport.sp3
import MomentRT.units
import MomentRT.validation

__version__ = MomentRT._version.__version__
ReflectingBoundary = MomentRT.boundary.ReflectingBoundary
RobinBoundary = MomentRT.boundary.RobinBoundary
Boundaries = MomentRT.boundary.Boundaries
DetectorResolutionError = MomentRT.detector.DetectorResolutionError
SurfaceDetector = MomentRT.detector.SurfaceDetector
Diffusion = MomentRT.transport.diffusion.Diffusion
InfiniteDiffusion = MomentRT.transport.diffusion.InfiniteDiffusion
BoundaryField = MomentRT.fields.BoundaryField
ScalarField = MomentRT.fields.ScalarField
VectorField = MomentRT.fields.VectorField
Box = MomentRT.geometry.Box
Layer = MomentRT.geometry.Layer
Sphere = MomentRT.geometry.Sphere
Material = MomentRT.material.Material
SpectralMaterial = MomentRT.material.SpectralMaterial
Water = MomentRT.materials.liquids.Water
Intralipid = MomentRT.materials.liquids.Intralipid
Epidermis = MomentRT.materials.skin.Epidermis
Dermis = MomentRT.materials.skin.Dermis
HomogeneousMedium = MomentRT.medium.HomogeneousMedium
StructuredMedium = MomentRT.medium.StructuredMedium
CartesianMesh = MomentRT.mesh.CartesianMesh
DetectorConvergenceResult = MomentRT.results.DetectorConvergenceResult
DetectorResult = MomentRT.results.DetectorResult
SimulationResult = MomentRT.results.SimulationResult
SpectralSimulationResult = MomentRT.results.SpectralSimulationResult
SampleSection = MomentRT.sections.SampleSection
TransportSolver = MomentRT.simulation.TransportSolver
Simulation = MomentRT.simulation.Simulation
SP3 = MomentRT.transport.sp3.SP3
P3 = MomentRT.transport.p3.P3
PhaseFunction = MomentRT.phase_function.PhaseFunction
HenyeyGreenstein = MomentRT.phase_function.HenyeyGreenstein
Isotropic = MomentRT.phase_function.Isotropic
Rayleigh = MomentRT.phase_function.Rayleigh
TabulatedPhaseFunction = MomentRT.phase_function.Tabulated
ValidationReport = MomentRT.validation.ValidationReport
PointSource = MomentRT.source.PointSource
SourceResolutionError = MomentRT.source.SourceResolutionError
SurfaceSource = MomentRT.source.SurfaceSource
UniformSurfaceSource = MomentRT.source.UniformSurfaceSource
ureg = MomentRT.units.ureg

__all__ = [
    "PhaseFunction",
    "HenyeyGreenstein",
    "Isotropic",
    "Rayleigh",
    "TabulatedPhaseFunction",
    "P3",
    "ValidationReport",
    "Boundaries",
    "SP3",
    "TransportSolver",
    "DetectorConvergenceResult",
    "DetectorResolutionError",
    "SourceResolutionError",
    "BoundaryField",
    "DetectorResult",
    "SpectralMaterial",
    "SpectralSimulationResult",
    "SurfaceDetector",
    "Diffusion",
    "RobinBoundary",
    "ReflectingBoundary",
    "SurfaceSource",
    "UniformSurfaceSource",
    "VectorField",
    "Box",
    "Layer",
    "Sphere",
    "StructuredMedium",
    "CartesianMesh",
    "HomogeneousMedium",
    "InfiniteDiffusion",
    "Material",
    "Water",
    "Intralipid",
    "Epidermis",
    "Dermis",
    "PointSource",
    "SampleSection",
    "ScalarField",
    "Simulation",
    "SimulationResult",
    "ureg",
    "__version__",
]
