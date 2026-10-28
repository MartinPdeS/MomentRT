"""Deterministic light transport using angular-moment methods."""

import MomentRT._version
import MomentRT.boundary
import MomentRT.convergence
import MomentRT.detector
import MomentRT.diffusion
import MomentRT.fem
import MomentRT.fields
import MomentRT.geometry
import MomentRT.material
import MomentRT.medium
import MomentRT.mesh
import MomentRT.plotting
import MomentRT.results
import MomentRT.sections
import MomentRT.simulation
import MomentRT.source
import MomentRT.units

# Expose tool modules alongside the convenience class aliases.
convergence = MomentRT.convergence
fem = MomentRT.fem
plotting = MomentRT.plotting

__version__ = MomentRT._version.__version__
ReflectingBoundary = MomentRT.boundary.ReflectingBoundary
RobinBoundary = MomentRT.boundary.RobinBoundary
DetectorResolutionError = MomentRT.detector.DetectorResolutionError
SurfaceDetector = MomentRT.detector.SurfaceDetector
Diffusion = MomentRT.diffusion.Diffusion
InfiniteDiffusion = MomentRT.diffusion.InfiniteDiffusion
BoundaryField = MomentRT.fields.BoundaryField
ScalarField = MomentRT.fields.ScalarField
VectorField = MomentRT.fields.VectorField
Box = MomentRT.geometry.Box
Layer = MomentRT.geometry.Layer
Sphere = MomentRT.geometry.Sphere
Material = MomentRT.material.Material
SpectralMaterial = MomentRT.material.SpectralMaterial
HomogeneousMedium = MomentRT.medium.HomogeneousMedium
StructuredMedium = MomentRT.medium.StructuredMedium
CartesianMesh = MomentRT.mesh.CartesianMesh
DetectorConvergenceResult = MomentRT.results.DetectorConvergenceResult
DetectorResult = MomentRT.results.DetectorResult
SimulationResult = MomentRT.results.SimulationResult
SpectralSimulationResult = MomentRT.results.SpectralSimulationResult
SampleSection = MomentRT.sections.SampleSection
Simulation = MomentRT.simulation.Simulation
PointSource = MomentRT.source.PointSource
SourceResolutionError = MomentRT.source.SourceResolutionError
SurfaceSource = MomentRT.source.SurfaceSource
UniformSurfaceSource = MomentRT.source.UniformSurfaceSource
ureg = MomentRT.units.ureg

__all__ = [
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
    "PointSource",
    "SampleSection",
    "ScalarField",
    "Simulation",
    "SimulationResult",
    "ureg",
    "__version__",
]
