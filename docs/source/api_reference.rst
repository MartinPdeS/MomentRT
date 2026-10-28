Public API reference
====================

All public constructors and methods use keyword arguments, apart from
``add_structures(*structures)`` and Python's implicit ``self`` and ``cls``.
Dimensional values require quantities from ``MomentRT.ureg``. See
:doc:`getting_started` for examples and :doc:`theory` and
:doc:`finite_element` for their physical interpretation.

Materials and media
-------------------

``Material`` requires ``absorption``, ``reduced_scattering``, and
``refractive_index``. The first two are nonnegative/positive inverse-length
quantities; refractive index is a positive plain number.

``HomogeneousMedium`` wraps one material for the infinite-medium model.
``StructuredMedium`` requires unit-bearing ``origin`` and positive ``size``
vectors. Set its background with ``add_background(material=...)`` and add
ordered ``Layer``, ``Box``, or ``Sphere`` objects with
``add_structures(*structures)``. ``add_layer(...)`` stacks below the deepest
existing layer.

Sources
-------

``PointSource`` requires ``position``, ``wavelength``, and positive
``optical_power``. It is supported by ``InfiniteDiffusion`` only.

``SurfaceSource`` requires a boundary ``face``, absolute 3D ``position``,
positive ``radius``, ``wavelength``, and ``optical_power``. The spot is a
circular top-hat and must fit fully on the face. ``UniformSurfaceSource``
requires a face, wavelength, and optical power, and illuminates the complete
face. Both finite-source classes specify power already coupled into the
sample.

Meshes and solvers
------------------

``CartesianMesh`` requires ``shape``, unit-bearing vector ``spacing``, and
unit-bearing vector ``origin``. ``CartesianMesh.from_sample(sample=...,
shape=...)`` creates a grid spanning every sample face.

``InfiniteDiffusion`` evaluates the analytical P1 Green function for a
``HomogeneousMedium`` and ``PointSource``. ``Diffusion`` requires a mapping
that configures exactly six faces with ``RobinBoundary`` or
``ReflectingBoundary`` objects. The finite solver accepts ``SurfaceSource`` or
``UniformSurfaceSource``.

``RobinBoundary`` requires a plain ``effective_reflectance`` in [0, 1).
``ReflectingBoundary`` imposes zero normal current and does not accept
illumination.

Simulation and results
----------------------

``Simulation`` composes ``medium``, ``source``, ``mesh``, and ``solver``.
Call ``simulation.run()`` to receive a ``SimulationResult`` containing
``fluence`` and, for the finite solver, ``flux``, ``energy_balance``,
``outgoing_power_by_face``, and ``diagnostics``.

``ScalarField.plot(...)`` returns a Plotly figure when Plotly is installed.
Its ``view`` may be ``volume`` (default), ``slices``, or ``isosurfaces``;
``scale`` may be ``linear`` or ``log``; ``length_unit`` controls displayed
coordinates; ``opacity`` chooses the initial opacity for volume and
isosurface views; and ``show`` controls whether the figure is immediately
opened. ``VectorField.plot(...)`` returns a Plotly cone plot and accepts
``max_vectors``, ``length_unit``, and ``show``.

The public classes are exported from :mod:`MomentRT`; plotting dependencies
are loaded only when a plotting method is called.
