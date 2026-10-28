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

``SpectralMaterial`` requires a strictly increasing ``wavelengths`` vector,
same-length ``absorption`` and ``reduced_scattering`` quantities, and a
scalar or same-length plain ``refractive_index``. Its
``at(wavelength=...)`` method linearly interpolates a monochromatic
``Material`` and rejects extrapolation. Fixed ``Material`` objects also
support ``at(wavelength=...)`` and retain their fixed coefficients.

``HomogeneousMedium`` wraps one material for the infinite-medium model.
``StructuredMedium`` requires unit-bearing ``origin`` and positive ``size``
vectors. Set its background with ``add_background(material=...)`` and add
ordered ``Layer``, ``Box``, or ``Sphere`` objects with
``add_structures(*structures)``. ``add_layer(...)`` stacks below the deepest
existing layer.

Both media support ``at(wavelength=...)``. Finite sample geometry and
insertion order are preserved on the resolved copy. Property sampling and
display methods accept an explicit ``wavelength`` for spectral samples.
See :doc:`spectral`.

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

Detectors
---------

``SurfaceDetector`` requires a sample boundary ``face``, unit-bearing
``position`` and positive ``radius``, a nonzero plain ``collection_axis``
pointing outward, and a unit-bearing ``acceptance_angle`` from zero to
90 degrees. The circular contact aperture lies in the sample face. The angle
is a cone half-angle around the collection axis.

``measure(result=...)`` returns a ``DetectorResult`` with
``collected_power``, ``sampled_area``, ``geometric_area``,
``angular_fraction``, ``minimum_escaping_irradiance``, wavelength,
and a snapshot of the detector configuration. This uses a Lambertian escape
assumption. See :doc:`detectors` for the integral, conventions, and limits.

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

``Simulation`` composes ``medium``, ``source``, ``mesh``, and ``solver``,
plus an optional ``detectors`` mapping from names to ``SurfaceDetector``
objects.
Call ``simulation.run()`` to receive a ``SimulationResult`` containing
``fluence`` and, for the finite solver, ``flux``, ``energy_balance``,
``outgoing_power_by_face``, ``boundary_fields``, and ``diagnostics``.
Configured measurements are in ``result.detectors``; the
``result.detector_power`` property maps their names to collected powers.
Boundary fields are ``BoundaryField`` snapshots of triangle-quadrature
positions, area weights, escaping irradiance, geometry, and wavelength.

``simulation.sweep(wavelengths=..., optical_power=...)`` returns a
``SpectralSimulationResult`` with ordered ``wavelengths``, per-wavelength
``results``, and detector-power vectors in ``detector_power``.
Optional optical power is a scalar or a vector matching the wavelength
shape; omitting it retains the configured source power.
Its ``plot(wavelength_unit=..., power_unit=..., show=...)`` method plots
collected power for configured detectors. Each solve evaluates material
tables at the current wavelength without mutating the original
configuration.

``ScalarField.plot(...)`` returns a Plotly figure when Plotly is installed.
Its ``view`` may be ``volume`` (default), ``slices``, or ``isosurfaces``;
``scale`` may be ``linear`` or ``log``; ``length_unit`` controls displayed
coordinates; ``opacity`` chooses the initial opacity for volume and
isosurface views; ``isosurface_count`` sets the number of iso levels in volume
and isosurface views while leaving their default counts intact when omitted
(it only applies to those views); and ``show`` controls whether the figure is
immediately opened.
``VectorField.plot(...)`` returns a Plotly cone plot and accepts
``max_vectors``, ``length_unit``, and ``show``.

The public classes are exported from :mod:`MomentRT`; plotting dependencies
are loaded only when a plotting method is called.
