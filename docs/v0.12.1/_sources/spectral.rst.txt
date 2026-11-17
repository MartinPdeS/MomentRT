Wavelength-dependent materials and sweeps
=========================================

``SpectralMaterial`` describes tabulated absorption, reduced scattering,
and refractive index. A solver evaluates these properties at the source
wavelength before calculating a monochromatic field. ``Simulation.sweep``
repeats that workflow while preserving geometry and instrument configuration.

See the :doc:`spectral detector example
<gallery/spectral/spectral_detectors>` for a layered phantom, two
contact apertures, collected-power spectra, and a 3D fluence volume.

Supplying a material table
---------------------------

.. code-block:: python

   import numpy as np
   import MomentRT
   from MomentRT import SpectralMaterial, ureg

   material = SpectralMaterial(
       wavelengths=np.array([700, 800, 900]) * ureg.nm,
       absorption=np.array([0.01, 0.02, 0.04]) / ureg.mm,
       reduced_scattering=np.array([1.2, 1.0, 0.8]) / ureg.mm,
       refractive_index=1.4,
       phase_function=MomentRT.phase_function.HenyeyGreenstein(),
   )

The wavelength vector must contain at least two positive, strictly increasing
entries. Absorption is nonnegative; reduced scattering and refractive index
are positive. The coefficient vectors must have the same shape as the
wavelength vector. Refractive index is a plain scalar shared across the
table or a plain vector of the same shape. All dimensional inputs require
quantities, and the supplied units are retained.

The example values are illustrative. The package does not ship a validated
tissue spectrum or derive a material spectrum from a wavelength label.
The caller is responsible for optical data appropriate to the material,
temperature, and experimental conditions.

Evaluation and interpolation
------------------------------

.. code-block:: python

   monochromatic = material.at(
       wavelength=750 * ureg.nm,
   )

   print(f"Absorption at 750 nm: {monochromatic.absorption:~P}")

For an optical property :math:`f` tabulated at adjacent wavelengths
:math:`\lambda_i,\lambda_{i+1}`, the interpolation rule is

.. math::

   t=\frac{\lambda-\lambda_i}{\lambda_{i+1}-\lambda_i},\qquad
   f(\lambda)=(1-t)f_i+tf_{i+1}.

Interpolation is linear in wavelength and in the supplied property, not in
its logarithm or optical frequency. Query wavelengths are converted to the
table's wavelength units. Out-of-range wavelengths are rejected; endpoint
conversion roundoff is allowed, but physical extrapolation is not performed.
Interpolating a positive table preserves positivity between its entries.
It does not establish that sparse tabulation resolves a real spectral peak.

The returned object is an ordinary ``Material``. A wavelength-independent
``Material.at(wavelength=...)`` returns the same fixed properties. Fixed and
tabulated materials can coexist in one sample.

Sample composition and plotting
--------------------------------

Assign a spectral material to a background, layer, sphere, or box in exactly
the same way as a fixed material. Both implemented solvers resolve it at the
source wavelength, including direct calls to ``solver.solve(...)``.
For a finite sample, every internal material must still have the same
refractive index at the wavelength being solved. Dispersion does not remove
the present equal-index interface restriction.

``sample.at(wavelength=...)`` returns a monochromatic sample with the same
extent, structure geometry, names, and insertion order. It does not replace
the original material tables. Property views on the original spectral sample
require an explicit wavelength:

.. code-block:: python

   figure = sample.plot(
       wavelength=800 * ureg.nm,
       show=False,
   )

   fields = sample.optical_fields(
       mesh=mesh,
       wavelength=800 * ureg.nm,
   )

   section = sample.section(
       normal="y",
       position=0 * ureg.mm,
       wavelength=800 * ureg.nm,
   )

``plot_section`` also accepts ``wavelength``. A simulation's system plot
uses its configured source wavelength. No arbitrary display wavelength is
chosen for a spectral material.

Sweeping an instrument configuration
--------------------------------------

Construct a ``Simulation`` once with a sample, source, mesh, solver, and
optional named :doc:`detectors`. Then run

.. code-block:: python

   wavelengths = np.array([700, 750, 800, 850, 900]) * ureg.nm

   spectrum = simulation.sweep(
       wavelengths=wavelengths,
   )

   powers = spectrum.detector_power["reflection"]

   print(f"Collected powers: {powers.to('mW'):~P}")

   figure = spectrum.plot(
       wavelength_unit="nm",
       power_unit="mW",
       show=False,
   )

The wavelength vector must be nonempty and positive. Results preserve
supplied order, including repeated wavelengths. Every requested wavelength
is checked against the material tables before solving. At each wavelength,
the solver evaluates

.. math::

   D(\mathbf r,\lambda)=
   \frac{1}{3[\mu_a(\mathbf r,\lambda)+\mu_s'(\mathbf r,\lambda)]}

and solves the same spatial boundary-value problem with those coefficients.
The sample geometry, source face and footprint, mesh, boundary configuration,
and detector aperture and acceptance remain fixed.

By default, every solve uses the configured incoming optical power.
For explicitly wavelength-dependent source powers, supply a positive scalar
or a power vector matching the requested wavelengths:

.. code-block:: python

   spectrum = simulation.sweep(
       wavelengths=wavelengths,
       optical_power=np.array([0.8, 1.0, 1.1, 1.0, 0.9]) * ureg.mW,
   )

This vector contains total monochromatic powers in watts, not a spectral
density in watts per nanometre. A scalar power applies to every wavelength.
The original source wavelength and power are retained after the sweep.
Boundary effective reflectance is held fixed; it is not automatically
recomputed from the material's refractive index.

Results and interpretation
----------------------------

``SpectralSimulationResult.wavelengths`` records the wavelength vector.
Its ``results`` tuple contains one ordinary ``SimulationResult`` per
wavelength, with fluence, flux, power balance, boundary fields, and configured
detector measurements. For example, ``spectrum.results[2].fluence.plot()``
shows the third monochromatic field.

``spectrum.detector_power`` maps detector names to unit-bearing power
vectors in the same order. Each constituent detector result retains its
sampled aperture area and angular fraction. Check spatial convergence for
the most demanding wavelengths, rather than assuming a mesh adequate at one
wavelength is adequate for the whole table.

The sweep remains a collection of independent monochromatic experiments.
Use a source configured with ``spectrum=...`` and ``Simulation.run()`` for
the explicit broadband integration described below. Chromophore fitting and
fluorescence that couples wavelengths remain outside this model.

Broadband sources and detector response
---------------------------------------

The same source and detector classes accept spectral emission and response.
Supply either ``wavelength`` and ``optical_power``, or a ``SourceSpectrum``;
the two emission forms are mutually exclusive. For example:

.. code-block:: python

   source_spectrum = MomentRT.SourceSpectrum(
       wavelengths=np.array([750, 800, 850, 900]) * ureg.nm,
       spectral_power=np.array([1, 4, 3, 1]) * ureg.uW / ureg.nm,
   )
   response = MomentRT.DetectorResponsivity(
       wavelengths=np.array([750, 900]) * ureg.nm,
       responsivity=np.array([0.4, 0.6]) * ureg.A / ureg.W,
   )
   source = MomentRT.SurfaceSource.on_face(
       sample=sample,
       face="z_min",
       coordinates=np.array([1.5, 1.5]) * ureg.mm,
       radius=0.5 * ureg.mm,
       spectrum=source_spectrum,
   )
   detector = MomentRT.SurfaceDetector.on_face(
       sample=sample,
       face="z_max",
       coordinates=np.array([1.5, 1.5]) * ureg.mm,
       radius=0.7 * ureg.mm,
       acceptance="full_hemisphere",
       responsivity=response,
   )
   simulation = MomentRT.Simulation(
       medium=sample,
       source=source,
       mesh=mesh,
       solver=solver,
       detectors={"sensor": detector},
   )
   result = simulation.run()
   measurement = result.detectors["sensor"]
   print(f"Optical power: {measurement.collected_power:~P}")
   print(f"Photocurrent: {measurement.photocurrent:~P}")

See the :doc:`complete executable example
<gallery/spectral/broadband_detector>` for the sample, solver, mesh,
and spectral plots. Its values are illustrative, not a calibrated instrument.

``SourceSpectrum`` requires increasing positive wavelengths and a nonnegative
power density in units such as W/nm. Zero-valued endpoints are supported;
the integrated source power must be positive. ``source.optical_power`` gives
the trapezoidal total over the stated interval. A broadband source has no
single ``source.wavelength``. Grouped surface sources currently require
monochromatic members.

``DetectorResponsivity`` linearly interpolates nonnegative A/W values without
physical extrapolation. Its table and every material table must cover the
entire source wavelength grid; coverage is checked before transport solves.
A monochromatic detector result also has ``photocurrent`` when a response
is supplied, equal to collected optical power times the interpolated
responsivity. Without a response, photocurrent is None.

For broadband emission, ``run()`` returns a ``BroadbandSimulationResult``.
Each named detector provides ``spectral_power`` in W per wavelength unit,
``collected_power`` in W, and optional ``spectral_photocurrent`` and
``photocurrent`` in A per wavelength unit and A. If the source density is
:math:`S_\lambda`, the computed dimensionless collection response is
:math:`T_\lambda`, and electrical responsivity is :math:`R_\lambda`, then

.. math::

   P_{\mathrm{det}} = \int S_\lambda T_\lambda\,d\lambda,
   \qquad I_{\mathrm{det}} = \int R_\lambda S_\lambda T_\lambda\,d\lambda.

Both integrals use trapezoidal quadrature on exactly the source grid and its
endpoints. Refine that grid to resolve source, material, and responsivity
features; intermediate response-table features are not automatically added
to the integration grid. Spectral quadrature convergence is separate from
spatial mesh and transport-model validation. Negative numerical detector
values are retained through integration.

Monochromatic and broadband detector measurements both expose
``sampled_area``, ``geometric_area``, ``angular_fraction``, and the detector
snapshot directly. Broadband ``reference_measurements`` retain the individual
reference-power measurements on the same aperture and quadrature.

``result.reference_results`` holds independent one-watt monochromatic solves
used to determine the collection response. Their fluence, flux, energy
diagnostics, and source and detector areas retain their reference excitation
and units; they are not spectral densities. Geometry is prepared once per
broadband call, while each wavelength gets its own operator and factorization.
The original source, detector, material, and mesh configuration is preserved.

The response is a linear electrical conversion after the existing aperture
and Lambertian angular collection model. This extension does not supply
Fresnel-resolved optics, detector saturation, electronic gain, dark current,
or acquisition noise. Analytical slab and wavelength-quadrature checks verify
the implementation; experimental instrument calibration remains necessary
for measured-device predictions.

Inspecting and refining spectral inputs
---------------------------------------

The tabulated objects provide composable interactive plots and explicit
source-density resampling:

.. code-block:: python

   source_figure = source.spectrum.plot(show=False)
   response_figure = detector.responsivity.plot(show=False)
   refined = source.spectrum.resample(
       wavelengths=np.linspace(750, 900, 61) * ureg.nm,
   )

Resampling linearly interpolates densities with their original units and
rejects extrapolation. Include original source knots when preserving the
measured piecewise linear shape and integrated source power. Coarsening away
a measured peak changes the emission configuration.

``simulation.spectral_convergence(...)`` holds the source data, geometry,
spatial mesh, and detector response fixed while refining wavelength quadrature:

.. code-block:: python

   study = simulation.spectral_convergence(
       wavelength_grids=[
           np.linspace(750, 900, count) * ureg.nm
           for count in (31, 61, 121)
       ],
       include_table_knots=True,
       power_rtol=0.01,
       photocurrent_rtol=0.01,
       required_consecutive=2,
   )
   figure = study.plot(show=False)
   print(study.summary())
   print(f"Recommended wavelengths: {study.recommended_wavelengths}")

Every requested grid must span the original source interval. Source knots are
always retained, and material and responsivity knots inside that interval are
included by default. The resulting grids must be nested and add wavelengths
at each level. ``study.wavelength_grids`` exposes the grids actually used.
This grid merging applies to the study; ordinary ``run()`` uses exactly the
source's configured grid.

``SpectralConvergenceResult`` provides per-level powers, photocurrents,
relative changes against adjacent and finest levels, and explicit failures
through ``records``. A recommendation requires the configured number of
passing refinements at the fine end for every detector's optical power and,
where supplied, photocurrent. Failed levels, unresolved footprints, and
negative outputs remain visible; a later failure invalidates an earlier
plateau. These are stability estimates, not rigorous error bounds or spatial
mesh validation. Use :doc:`mesh convergence <convergence>` separately at
representative monochromatic wavelengths.

Reusing finite-element preparation
-----------------------------------

For SP3 or P3 with the Henyey–Greenstein model, each
``SpectralMaterial`` must additionally specify ``anisotropy`` as a plain
scalar or a vector on the wavelength grid. A vector is interpolated
linearly along with the other material coefficients; separate materials can
have different anisotropies. Fixed isotropic, Rayleigh, and tabulated phase
functions supply their own first moment and may omit anisotropy. P1 does not
require this property.

For finite ``P1``, ``SP3``, and ``P3`` sweeps, the package prepares mesh geometry,
tetrahedra, element volumes and gradients, centroid material labels, boundary
triangles, and boundary quadrature once per sweep. It also prepares the
fixed boundary matrix, source footprint, and normalized source load.

Centroid labels partition the elements into material regions. Unit
diffusion and absorption matrices are assembled once per occupied region:

.. math::

   K^{(m)}_{ij}=\int_{\Omega_m}\nabla N_i\cdot\nabla N_j\,dV,
   \qquad
   M^{(m)}_{ij}=\int_{\Omega_m}N_iN_j\,dV.

At a wavelength, the global matrix is formed as

.. math::

   A(\lambda)=B+
   \sum_m\left[D_m(\lambda)K^{(m)}+
   \mu_{a,m}(\lambda)M^{(m)}\right],
   \qquad
   b(\lambda)=P_{\mathrm{in}}(\lambda)b_{\mathrm{unit}}.

Here :math:`B` contains the fixed Robin boundary terms, and
:math:`b_{\mathrm{unit}}` includes the boundary source factor.
Source normalization is preserved when optical power varies with wavelength.

Material interpolation, internal-index checks, matrix combination, sparse
factorization, and field evaluation still occur at each wavelength.
Factorization is not reused for a changed matrix. Boundary snapshots remain
independent outputs, and detector readouts use the same quadrature as
standalone solves.

Preparation is scoped to the current call; there is no persistent cache
requiring manual invalidation. Changing mesh, structures, source footprint,
material tables, or boundaries before another call creates fresh
preparation. Infinite-medium sweeps retain their analytical path.

``spectrum.diagnostics`` reports ``geometry_preparations``, unit-bearing
``preparation_time``, per-wavelength ``solve_times``, and ``total_time``.
Per-wavelength times include field and detector evaluation. Total time
also includes validation and preparation. Individual finite results record
``diagnostics["preparation_reused"]``: False for the first solve and True
for later wavelengths in the same prepared sweep.

Time saved depends on geometry size, material regions, wavelength count,
and sparse solve cost. Use these timings for the actual configuration
rather than assuming a fixed speedup.
