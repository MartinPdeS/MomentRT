Python class reference
======================

All dimensional arguments require quantities from ``MomentRT.ureg``.
Array shapes and physical conventions are described for each class.
Objects are grouped by domain module; top-level convenience imports remain
available.

Materials and sample geometry
-----------------------------

.. autoclass:: MomentRT.phase_function.PhaseFunction
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.phase_function.HenyeyGreenstein
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.phase_function.Isotropic
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.phase_function.Rayleigh
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.phase_function.Tabulated
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.material.Material
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.material.SpectralMaterial
   :members:
   :member-order: bysource

Premade materials
~~~~~~~~~~~~~~~~~

.. autoclass:: MomentRT.materials.Epidermis
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.materials.Dermis
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.materials.Intralipid
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.materials.Water
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.medium.HomogeneousMedium
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.medium.StructuredMedium
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.geometry.Box
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.geometry.Sphere
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.geometry.Layer
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.sections.SampleSection
   :members:
   :member-order: bysource

Sources, detectors, and boundaries
----------------------------------

.. autoclass:: MomentRT.source.SourceSpectrum
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.detector.DetectorResponsivity
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.source.PointSource
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.source.UniformSurfaceSource
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.source.SurfaceSource
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.detector.SurfaceDetector
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.boundary.RobinBoundary
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.boundary.ReflectingBoundary
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.boundary.Boundaries
   :members:
   :member-order: bysource

Meshes and transport
--------------------

.. autoclass:: MomentRT.mesh.CartesianMesh
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.mesh.TetrahedralMesh
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.transport.InfiniteP1
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.transport.P1
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.transport.TransientP1
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.results.TransientResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.transport.SP3
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.transport.P3
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.simulation.TransportSolver
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.simulation.Simulation
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.validation.ValidationReport
   :members:
   :member-order: bysource

Fields and results
------------------

.. autoclass:: MomentRT.fields.ScalarField
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.fields.VectorField
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.fields.BoundaryField
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.results.SimulationResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.results.DetectorResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.convergence.MeshConvergence
   :members:

.. autoclass:: MomentRT.results.DetectorConvergenceResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.results.SpectralSimulationResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.results.BroadbandSimulationResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.results.BroadbandDetectorResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.results.SpectralConvergenceResult
   :members:
   :member-order: bysource
