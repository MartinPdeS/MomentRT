Choosing a mesh for detector collection
============================================================

``simulation.convergence_study(...)`` compares detector power and footprint
sampling across explicit finite meshes. It keeps the sample, source,
detector apertures, wavelength, and boundaries fixed, returning a
``DetectorConvergenceResult`` without changing the original configuration.
See the :doc:`worked refinement example
<gallery/finite_samples/detector_convergence>` for comparison plots and a
fluence volume with instrument overlays.

Running a study
-----------------

Configure a finite ``Diffusion`` or ``SP3`` simulation with at least one named detector,
then supply meshes from coarse to fine:

.. code-block:: python

   study = simulation.convergence_study(
       mesh_shapes=[
           (9, 7, 7),
           (13, 10, 10),
           (17, 13, 13),
           (25, 19, 19),
       ],
       power_rtol=0.01,
       area_rtol=0.02,
       power_atol=0 * ureg.watt,
       required_consecutive=2,
   )

   figure = study.plot(
       power_unit="mW",
       area_unit="mm^2",
       show=False,
   )

   print(f"Recommended mesh shape: {study.recommended_shape}")

Every shape contains three node counts of at least two. Each successive mesh
must refine at least one axis and coarsen none. At least two levels are
required, but two consecutive passing refinements need at least three
resolved levels. The study evaluates all supplied levels rather than
stopping at the first apparent plateau.

These tolerances are numerical choices, not claims about transport-model
error. A study uses the configured source wavelength. Repeat it at demanding
wavelengths when assessing a spectral experiment.

Power and area criteria
------------------------

Let :math:`P_{d,k}` be detector power at mesh level :math:`k` and
:math:`A_{d,k}` the sampled detector area. Its geometric disk area is
:math:`A_d=\pi r_d^2`. The source has sampled area :math:`A_{s,k}` and
geometric area :math:`A_s`: a disk area for ``SurfaceSource`` or the
complete face area for ``UniformSurfaceSource``.

.. math::

   e_{d,k}=\left|\frac{A_{d,k}}{A_d}-1\right|,
   \qquad
   e_{s,k}=\left|\frac{A_{s,k}}{A_s}-1\right|.

Both errors compare sampled areas with geometry, rather than only with the
previous mesh. Two coarse meshes could sample the same incorrect area.
Source normalization preserves incoming power even when a footprint is
poorly sampled, so energy conservation cannot replace this area check.

Successive powers pass when

.. math::

   |P_{d,k}-P_{d,k-1}|
   \le P_{\mathrm{atol}}+r_P|P_{d,k}|.

The report also records relative change
:math:`|P_{d,k}-P_{d,k-1}|/|P_{d,k}|`. The first value is NaN because no
previous estimate exists. Two zero powers have relative change zero. A
change to zero from nonzero power has infinite relative change, but may
still pass an explicitly supplied absolute tolerance. Missing observations
never pass.

Defaults are ``power_rtol=0.01``, ``area_rtol=0.02``, zero absolute
power tolerance, and ``required_consecutive=2``. A nonzero
``power_atol`` requires an explicit power quantity; choose it according to
the numerical accuracy needed for a small signal.

Recommendation rule
---------------------

A level passes only when every detector passes its successive power check,
all source and detector area errors are within ``area_rtol``, and there
are no unresolved footprints or reported negative fields at that level.

The finest end of the sequence must have at least
``required_consecutive`` uninterrupted passing refinements. A later failure
discards an earlier apparent plateau. Among candidates in that final stable
sequence, the first acceptable mesh must also agree with the final level's
detector powers:

.. math::

   |P_{d,k}-P_{d,K}|
   \le P_{\mathrm{atol}}+r_P|P_{d,K}|.

This prevents small successive changes from accumulating into a large
discrepancy with the finest supplied observation. The recommendation is
evidence from a finite sequence, not a rigorous error bound relative to an
infinitely refined solution.

``recommended_mesh`` and ``recommended_shape`` are None and ``converged``
is False when the criteria are not met. The final mesh is not automatically
declared sufficient.

Inspecting the report
----------------------

The report retains:

* ``meshes``, ``mesh_shapes``, ``node_counts``, and ``tetrahedron_counts``.
* Unit-bearing ``detector_power[name]``, ``detector_sampled_area[name]``,
  and ``detector_geometric_area[name]``.
* ``source_sampled_area`` and ``source_geometric_area``.
* ``power_relative_change[name]`` and ``power_relative_to_finest[name]``.
* ``detector_area_relative_error[name]`` and ``source_area_relative_error``.
* ``meets_tolerances``, ``stable_refinements``, and ``level_errors``.
* ``results``, containing the available finite solution at each level.

``records`` groups quantities and diagnostics by mesh level:

.. code-block:: python

   for record in study.records:
       measurement = record["detectors"]["reflection"]

       print(f"Mesh: {record['mesh_shape']}")

       print(f"Collected power: {measurement['collected_power'].to('mW'):~P}")

       print(f"Aperture area error: {measurement['area_relative_error']:.2%}")

The plot compares powers and sampled footprint areas against mesh node
count; dashed area curves are the geometric references. Missing observations
remain gaps.

Unresolved levels and configuration errors
-------------------------------------------

If no quadrature point samples a source disk, the level records
``SourceResolutionError``, source area zero, and ``results[k]=None``.
Detector powers and areas for that unsolved level are unknown and NaN.
The study continues on finer meshes.

If a solve succeeds but an aperture contains no quadrature points,
``DetectorResolutionError`` is recorded for that detector. Its sampled area
is zero and its power is unknown. Other measurements remain available.
An unresolved level breaks a stable sequence. Negative nodal fluence or
escaping irradiance in an aperture is preserved and flagged rather than
clipped or accepted for a recommendation.

Both resolution exceptions inherit from ValueError, preserving existing
exception handling for single solves. Other errors, such as an aperture
outside the sample or unsupported index mismatch, propagate as configuration
errors rather than being disguised as coarse-mesh failures.

This study assesses resolution of the tested observables and footprints.
It does not establish P1 or SP3 accuracy against RTE, validate material spectra or
detector optics, or guarantee convergence of unrelated local field values.
See :doc:`finite_element` and :doc:`detectors` for those distinctions.
