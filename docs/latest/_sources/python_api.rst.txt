Python class reference
======================

All dimensional arguments require quantities from ``MomentRT.ureg``.
Array shapes and physical conventions are described for each class.
Plotting methods require the relevant optional visualization extra.

Materials and sample geometry
-----------------------------

.. autoclass:: MomentRT.Material
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.SpectralMaterial
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.HomogeneousMedium
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.StructuredMedium
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.Box
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.Sphere
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.Layer
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.SampleSection
   :members:
   :member-order: bysource

Sources, detectors, and boundaries
----------------------------------

.. autoclass:: MomentRT.PointSource
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.UniformSurfaceSource
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.SurfaceSource
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.SurfaceDetector
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.RobinBoundary
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.ReflectingBoundary
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.Boundaries
   :members:
   :member-order: bysource

Meshes and transport
--------------------

.. autoclass:: MomentRT.CartesianMesh
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.InfiniteDiffusion
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.Diffusion
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.SP3
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.TransportSolver
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.Simulation
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.ValidationReport
   :members:
   :member-order: bysource

Fields and results
------------------

.. autoclass:: MomentRT.ScalarField
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.VectorField
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.BoundaryField
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.SimulationResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.DetectorResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.DetectorConvergenceResult
   :members:
   :member-order: bysource

.. autoclass:: MomentRT.SpectralSimulationResult
   :members:
   :member-order: bysource
