Choosing a mesh for detector collection
============================================================

``simulation.mesh_convergence(...)`` compares detector power and footprint
sampling across explicit finite meshes. It keeps the sample, source,
detector apertures, wavelength, and boundaries fixed, returning a
``DetectorConvergenceResult`` without changing the original configuration.
``simulation.convergence_study(...)`` remains a compatibility alias for
``simulation.mesh_convergence(...)`` with the same arguments and behavior.
See the :doc:`worked refinement example
<gallery/sources_detectors/detector_convergence>` for comparison plots and a
fluence volume with instrument overlays.

Running a study
-----------------

Configure a finite ``Diffusion``, ``SP3``, or ``P3`` simulation with at least one named detector,
then supply meshes from coarse to fine:

.. code-block:: python

   study = simulation.mesh_convergence(
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

For interface-fitted tetrahedra, give global target sizes from coarse to fine.
Optional contact refinement applies the same local element target near the
source disks and all detector apertures at every level:

.. code-block:: python

   study = simulation.mesh_convergence(
       mesh_resolutions=[
           1.5 * ureg.mm,
           1.0 * ureg.mm,
           0.7 * ureg.mm,
       ],
       contact_resolution=0.25 * ureg.mm,
       contact_distance=0.5 * ureg.mm,
       power_rtol=0.02,
       area_rtol=0.03,
   )

   print(f"Recommended target resolution: {study.recommended_resolution}")

The global sizes must decrease strictly. A contact target must be no larger
than any global target and is applied within ``contact_distance`` of each disk,
with a gradual transition to the global size. For a uniform full-face source,
the local rule covers detectors; the source already spans the face. Each level
is generated from the original sample and its own contact geometry. The report
continues to compare sampled source and aperture areas with their geometric
areas and detector power across levels. Target element size is not a strict
bound on every tetrahedron edge, and satisfying these criteria does not
establish transport accuracy.

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

``recommended_mesh`` is None and ``converged`` is False when the criteria are
not met. ``recommended_shape`` applies to Cartesian meshes;
``recommended_resolution`` reports the target size for a generated tetrahedral
mesh. The final mesh is not automatically declared sufficient.

Inspecting the report
----------------------

The report retains:

* ``meshes``, storage ``mesh_shapes``, ``mesh_resolutions``, ``node_counts``,
  and actual ``tetrahedron_counts``.
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
It does not establish P1, SP3, or P3 accuracy against RTE, validate material spectra or
detector optics, or guarantee convergence of unrelated local field values.
See :doc:`finite_element` and :doc:`detectors` for those distinctions.

Performance measurements
------------------------

The repository includes a repeatable finite spectral sweep benchmark. For example::

   .venv/bin/python tools/benchmark_sweeps.py --solver p1 --shape 9 9 9 --wavelengths 3
   .venv/bin/python tools/benchmark_sweeps.py --solver sp3 --shape 15 15 11 --wavelengths 5

The JSON output records the mesh, solver, software versions, geometry preparation
and individual solve times, and the peak memory tracked by Python's
``tracemalloc`` during the sweep and, where available, peak process RSS since
process start (including imports and setup). Traced memory may exclude native
solver allocations; peak RSS is unavailable on Windows. Compare runs on the same machine and software
configuration; these illustrative cases do not establish numerical convergence.

For isolated tetrahedral mesh and direct-solver scaling, use:

.. code-block:: console

   .venv/bin/python tools/benchmark_scaling.py \
       --solvers p1 sp3 p3 \
       --resolutions-mm 2.0 1.5 1.0

Each solver and target-size pair runs in a fresh process. The JSON report
includes mesh generation time, solve time, node and tetrahedron counts, and
process peak RSS. RSS includes imports, mesh generation, assembly, and solve;
compare runs on the same machine and software configuration. Target size is
only a mesher request, so use the reported element counts when comparing
meshes. The results characterize current direct-solver scaling; they do not
estimate accuracy or recommend a mesh.
