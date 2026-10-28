Wavelength-dependent materials and sweeps
=========================================

``SpectralMaterial`` describes tabulated absorption, reduced scattering,
and refractive index. A solver evaluates these properties at the source
wavelength before calculating a monochromatic field. ``Simulation.sweep``
repeats that workflow while preserving geometry and instrument configuration.

See the :doc:`spectral detector example
<gallery/finite_samples/spectral_detectors>` for a layered phantom, two
contact apertures, collected-power spectra, and a 3D fluence volume.

Supplying a material table
---------------------------

.. code-block:: python

   material = SpectralMaterial(
       wavelengths=np.array([700, 800, 900]) * ureg.nm,
       absorption=np.array([0.01, 0.02, 0.04]) / ureg.mm,
       reduced_scattering=np.array([1.2, 1.0, 0.8]) / ureg.mm,
       refractive_index=1.4,
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

The sweep is a collection of independent monochromatic experiments.
It does not integrate an LED emission spectrum, interpolate detector
responsivity, calculate an electrical current, fit chromophore
concentrations, or simulate fluorescence that couples wavelengths.
Those operations require additional physical data and normalization.

Reusing finite-element preparation
-----------------------------------

For finite ``Diffusion`` sweeps, the package prepares mesh geometry,
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
