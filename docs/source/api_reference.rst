Public API reference
====================

All public constructors and methods use keyword arguments, apart from
``add_structures(*structures)`` and Python's implicit ``self`` and ``cls``.
Dimensional values require quantities from ``MomentRT.ureg``. See
:doc:`getting_started` for examples and :doc:`theory` and
:doc:`finite_element` for their physical interpretation.

Materials and media
-------------------

``Material`` requires ``absorption``, ``reduced_scattering``,
``refractive_index``, and a ``phase_function`` object. The first two are
nonnegative/positive inverse-length quantities; refractive index is a positive
plain number. Henyey–Greenstein also requires ``anisotropy`` for SP3/P3.

``SpectralMaterial`` requires a strictly increasing ``wavelengths`` vector,
same-length ``absorption`` and ``reduced_scattering`` quantities, and a
scalar or same-length plain ``refractive_index``, and an explicit
``phase_function`` object. Its
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

Parameterized materials are available as ``MomentRT.Intralipid(concentration=,
wavelength=, refractive_index=)``, ``MomentRT.Epidermis(melanin_fraction=,
wavelength=, refractive_index=)``, ``MomentRT.Dermis(hydration=, wavelength=,
refractive_index=)``, and ``MomentRT.Water(wavelength=, refractive_index=)``.
Intralipid concentration is particle volume fraction in percent: 22.7 is
undiluted 20% stock. Its empirical equation has a wavelength-dependent validity
limit. The simple skin models put melanin in the epidermis; the dermis model
mixes water and a baseline absorber and omits blood-vessel structure. These
literature parameterizations are starting points, not patient-specific
properties. Pure water's very low scattering also makes it outside the useful
regime of MomentRT's diffusion solvers.

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

``SP3`` solves two coupled moment equations on the same finite sample and
uses each material's explicit :doc:`phase function <phase_functions>`. The
Henyey–Greenstein function requires material ``anisotropy``. It supports reflecting
faces and zero-reflectance Robin faces with isotropic incoming diffuse
radiance. See :doc:`sp3` for the equations and limitations.

``P3`` solves the full three-dimensional angular expansion through degree
three with six even-parity nodal fields. It requires material phase moments,
isotropic incoming diffuse sources, and zero-reflectance Robin conditions on
all faces. Reflecting and refracting boundaries are not supported. See
:doc:`p3` for the formulation and validation limits.

``RobinBoundary`` requires a plain ``effective_reflectance`` in [0, 1).
``ReflectingBoundary`` imposes zero normal current and does not accept
illumination.

``Boundaries.uniform(condition=...)`` applies an explicitly chosen condition
to all six faces. ``boundaries.with_overrides(conditions={...})`` returns a
new mapping with selected faces replaced. Both solvers also accept ordinary
six-face mappings.

For circular surface objects, ``SurfaceSource.on_face(...)`` and
``SurfaceDetector.on_face(...)`` accept the sample, face, and two absolute
in-plane coordinates. Coordinates follow Cartesian axis order: ``(y, z)``
for x faces, ``(x, z)`` for y faces, and ``(x, y)`` for z faces. A detector
requires ``orientation="outward_normal"``; use its full constructor for a
tilted collection axis. The resulting objects expose their full 3D positions
and detector axis.

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

``simulation.validate(wavelengths=...)`` checks configuration and, for finite
samples, measures source and detector footprint areas using boundary
quadrature. Its ``ValidationReport`` contains ``valid``, ``errors``,
``warnings``, and sampled areas. Omit wavelengths to check the source
wavelength. Call ``report.summary()`` or ``result.summary()`` for readable
diagnostics. A valid configuration does not establish mesh convergence or
accuracy of the physical model.

``simulation.sweep(wavelengths=..., optical_power=...)`` returns a
``SpectralSimulationResult`` with ordered ``wavelengths``, per-wavelength
``results``, and detector-power vectors in ``detector_power``.
Optional optical power is a scalar or a vector matching the wavelength
shape; omitting it retains the configured source power.
Its ``plot(wavelength_unit=..., power_unit=..., show=...)`` method plots
collected power for configured detectors. Each solve evaluates material
tables at the current wavelength without mutating the original
configuration.

Finite sweeps prepare geometry, material labels, region matrices, and
boundary quadrature once per call. Spectral ``diagnostics`` include
preparation and per-wavelength timing quantities. See :doc:`spectral`.

``simulation.convergence_study(mesh_shapes=..., power_rtol=...,
area_rtol=..., power_atol=..., required_consecutive=...)`` returns a
``DetectorConvergenceResult`` for configured finite detectors. It retains
per-level powers, footprint areas, errors, fields, and convergence flags.
The ``recommended_mesh`` and ``recommended_shape`` are None when no
supported recommendation can be made. Its ``plot(power_unit=...,
area_unit=..., show=...)`` method compares powers and areas.
See :doc:`convergence` for the criteria and unresolved levels.

``ScalarField.plot(...)`` returns a Plotly figure when Plotly is installed.
Its ``view`` may be ``volume`` (default), ``slices``, or ``isosurfaces``;
``scale`` may be ``linear`` or ``log``; ``length_unit`` controls displayed
coordinates; ``opacity`` chooses the initial opacity for volume and
isosurface views; ``isosurface_count`` sets the number of iso levels in volume
and isosurface views while leaving their default counts intact when omitted
(it only applies to those views); and ``show`` controls whether the figure is
immediately opened.
The optional ``overlays`` argument adds sources and detectors to any scalar
field's 3D view. Supply a sequence of objects or a mapping from legend labels
to objects. Their spatial coordinates are converted to ``length_unit``.
Full-face sources require the physical sample bounds saved by the finite
solver. See :doc:`results` for examples.
``VectorField.plot(...)`` returns a Plotly cone plot and accepts
``max_vectors``, ``length_unit``, and ``show``.

The public classes are exported from :mod:`MomentRT`; plotting dependencies
are optional and guarded at module scope.

Detailed Python API
-------------------

The reference below is generated from NumPy-style docstrings and type hints.
Dimensional annotations use ``pint.Quantity``; the docstrings specify the
required physical units and array shapes. Runtime unit validation remains
in force. Third-party solvers can implement the ``TransportSolver`` protocol
without inheriting from it.

.. toctree::
   :maxdepth: 1

   python_api
