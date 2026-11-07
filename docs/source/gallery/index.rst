Examples
========

Choose a task below to find a complete, executable example. New to MomentRT?
Start with an analytical calculation, then try a finite sample. The values are
illustrative; each example describes its physical assumptions and limits.

Getting started
---------------

Run a first simulation and construct a structured sample.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Model a skin contact sensor
      :link: getting_started/skin_contact
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_skin_contact_thumb.png
      :img-alt: Fluence in a three-layer skin sample with contact source and detector

      Build a layered skin sample, refine a tetrahedral mesh, and solve for fluence.

   .. grid-item-card:: Infinite-medium diffusion
      :link: getting_started/infinite_diffusion
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_infinite_diffusion_thumb.png
      :img-alt: Infinite-medium MomentRT fluence field

      Configure a point source and observation mesh, solve the analytical P1
      Green function, and explore its 3D fluence field.

   .. grid-item-card:: Circular surface illumination
      :link: getting_started/surface_diffusion
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_surface_diffusion_thumb.png
      :img-alt: Fluence produced by a diffuse circular surface source

      Inspect fluence, flux, and energy accounting for a finite-element solve
      with a disk source.

   .. grid-item-card:: Build a heterogeneous sample
      :link: getting_started/structured_sample
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_structured_sample_thumb.png
      :img-alt: Fluence volume in a layered sample with inclusions

      Combine layers, spheres, and a box, solve the illuminated sample, and
      explore its 3D fluence volume.

Samples and meshes
------------------

Use named sample models, then create Cartesian or tetrahedral meshes.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Build samples with named models
      :link: samples_and_meshes/sample_models
      :link-type: doc
      :img-top: samples_and_meshes/images/thumb/sphx_glr_sample_models_thumb.png
      :img-alt: Fluence in an explicitly configured two-layer skin sample

      Build a two-layer skin sample and a layered phantom from explicit materials.

   .. grid-item-card:: Cartesian mesh from node counts
      :link: samples_and_meshes/cartesian_shape
      :link-type: doc
      :img-top: meshes/images/thumb/sphx_glr_cartesian_shape_thumb.png
      :img-alt: Regular Cartesian grid shown inside its sample bounds

      Set explicit node counts for a regular grid and plot its cell edges.

   .. grid-item-card:: Cartesian mesh from physical resolution
      :link: samples_and_meshes/cartesian_resolution
      :link-type: doc
      :img-top: meshes/images/thumb/sphx_glr_cartesian_resolution_thumb.png
      :img-alt: Cartesian grid with different spacing by axis

      Request spacing in physical units, including anisotropic spacing by axis.

   .. grid-item-card:: Uniform tetrahedral mesh
      :link: samples_and_meshes/tetrahedral_global
      :link-type: doc
      :img-top: meshes/images/thumb/sphx_glr_tetrahedral_global_thumb.png
      :img-alt: Tetrahedral elements inside the sample boundary

      Generate and plot an unstructured tetrahedral mesh at a global target size.

   .. grid-item-card:: Locally refined tetrahedral mesh
      :link: samples_and_meshes/tetrahedral_local
      :link-type: doc
      :img-top: meshes/images/thumb/sphx_glr_tetrahedral_local_thumb.png
      :img-alt: Tetrahedral mesh refined around a skin layer and fiber contacts

      Refine a named layer and the neighborhoods around a source and detector.

Sources and detectors
---------------------

Configure diffuse illumination and contact apertures, and check detector mesh convergence.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Uniform face illumination
      :link: sources_and_detectors/uniform_surface_source
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_uniform_surface_source_thumb.png
      :img-alt: Finite-sample fluence under uniform face illumination

      Use ``UniformSurfaceSource`` to illuminate an entire sample face and
      render the resulting field in 3D.

   .. grid-item-card:: Collect light with finite apertures
      :link: sources_and_detectors/surface_detector
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_surface_detector_thumb.png
      :img-alt: Finite-sample fluence with source and detector apertures

      Define detector positions and acceptance cones, inspect the system, and
      calculate reflection and transmission signals.

   .. grid-item-card:: Seven-fiber skin probe
      :link: sources_and_detectors/seven_fiber_skin
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_seven_fiber_skin_thumb.png
      :img-alt: Six illumination fibers and a central detector on layered skin

      Simulate simultaneous illumination from six fibers around a central
      collection fiber in a layered skin model.

   .. grid-item-card:: Check detector mesh convergence
      :link: sources_and_detectors/detector_convergence
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_detector_convergence_thumb.png
      :img-alt: Finest-mesh fluence with source and detector overlays

      Compare detector signals and sampled source and aperture areas over
      successively finer meshes.

Spectral calculations
---------------------

Sweep wavelength-dependent materials and integrate source spectra and detector responsivity.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Build analytical source spectra
      :link: spectral/source_spectra
      :link-type: doc
      :img-top: spectral/images/thumb/sphx_glr_source_spectra_thumb.png
      :img-alt: Infinite-medium fluence for an analytical source spectrum

      Compare Gaussian and Lorentzian wavelength profiles and attach one to a source.

   .. grid-item-card:: Sweep a spectral sensor
      :link: spectral/spectral_detectors
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_spectral_detectors_thumb.png
      :img-alt: Monochromatic fluence from a spectral sensor sweep

      Interpolate material spectra and plot collected power across wavelengths
      with fixed sample and instrument geometry.

   .. grid-item-card:: Integrate source and detector spectra
      :link: spectral/broadband_detector
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_broadband_detector_thumb.png
      :img-alt: Incoming and collected spectral power and detector photocurrent

      Configure the existing source with a power density and the detector with
      responsivity, then integrate optical power and electrical current.

   .. grid-item-card:: Sweep spectral anisotropy
      :link: spectral/sp3_spectral_anisotropy
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_sp3_spectral_anisotropy_thumb.png
      :img-alt: SP3 fluence from wavelength-dependent material tables

      Interpolate layer-specific anisotropy and detector power across
      wavelengths.

Transport models
----------------

Compare angular approximations and explore phase functions and layer anisotropy.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Compare P1 and SP3 transport
      :link: transport_models/compare_sp3
      :link-type: doc
      :img-top: transport_models/images/thumb/sphx_glr_compare_sp3_thumb.png
      :img-alt: SP3 fluence volume with contact source and detector

      Compare detector signals from P1 and SP3 on the same sample.

   .. grid-item-card:: Compare P1, SP3, and full P3
      :link: transport_models/compare_p3
      :link-type: doc
      :img-top: transport_models/images/thumb/sphx_glr_compare_p3_thumb.png
      :img-alt: Full P3 fluence with source and contact detector

      Solve the same sensor geometry with three angular approximations.

   .. grid-item-card:: Compare scattering phase functions
      :link: transport_models/phase_functions
      :link-type: doc
      :img-top: transport_models/images/thumb/sphx_glr_phase_functions_thumb.png
      :img-alt: SP3 fluence for a Rayleigh scattering example

      Compare phase functions with matching anisotropy but different higher
      moments.

   .. grid-item-card:: Assign anisotropy by layer
      :link: transport_models/sp3_layered_anisotropy
      :link-type: doc
      :img-top: transport_models/images/thumb/sphx_glr_sp3_layered_anisotropy_thumb.png
      :img-alt: SP3 fluence in a layered sample with a source and detector

      Give each material its own scattering anisotropy and solve with SP3.

Parameter studies
-----------------

Vary optical properties and boundary reflectance, and examine computational scaling.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: How absorption changes fluence
      :link: parameter_studies/absorption_sweep
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_absorption_sweep_thumb.png
      :img-alt: Infinite-medium fluence at several absorption values

      Compare fluence profiles while holding reduced scattering and source
      configuration fixed.

   .. grid-item-card:: Change boundary reflectance
      :link: parameter_studies/reflectance_sweep
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_reflectance_sweep_thumb.png
      :img-alt: Fluence volume for an intermediate boundary reflectance

      Compare centerline fluence and escaping power as diffuse reflectance
      changes.

   .. grid-item-card:: Benchmark P1, SP3, and P3 scaling
      :link: parameter_studies/solver_scaling
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_solver_scaling_thumb.png
      :img-alt: Solve time and peak process memory versus tetrahedron count

      Measure solve time and peak process memory across tetrahedral mesh sizes,
      and use the results to decide whether iterative P3 solvers warrant study.

Validation
----------

Assess mesh refinement and compare calculations with independent references.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Mesh refinement and FEM timings
      :link: validation/mesh_refinement
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_mesh_refinement_thumb.png
      :img-alt: Finite-sample fluence on the finest mesh

      Compare fluence over three grids and inspect recorded FEM preparation,
      assembly, factorization, and solve timings.

   .. grid-item-card:: Compare against planar slab references
      :link: validation/planar_slab_references
      :link-type: doc
      :img-alt: P1 and SP3 finite-element fluence compared with analytical slab solutions

      Compare P1 diffusion and SP3 FEM profiles with independent closed-form
      and matrix-exponential slab solutions.

   .. grid-item-card:: Vary layered sensor validation cases
      :link: validation/sensor_validation_cases
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_sensor_validation_cases_thumb.png
      :img-alt: P1, SP3, and P3 detector-power errors across eight sensor cases

      Explore published comparisons while varying layer thickness, optical
      properties, source size, detector radius, and source-detector separation.

Inverse problems
----------------

Explore synthetic inverse calculations and their numerical and model limitations.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Check a mesh before fitting
      :link: inverse_problems/mesh_convergence
      :link-type: doc
      :img-top: inverse_problems/images/thumb/sphx_glr_mesh_convergence_thumb.png
      :img-alt: Fluence field in the layered sample used for mesh convergence

      Check how refinement affects detector powers for a layered sample before
      attempting a synthetic fit.

   .. grid-item-card:: Explore an objective landscape
      :link: inverse_problems/objective_landscape
      :link-type: doc
      :img-top: inverse_problems/images/thumb/sphx_glr_objective_landscape_thumb.png
      :img-alt: Infinite-medium fluence for a synthetic absorption example

      Plot the mismatch between synthetic fluence observations and candidate
      absorption values under the same analytical model.

   .. grid-item-card:: Fit a simulated epidermis change
      :link: inverse_problems/epidermis_change
      :link-type: doc
      :img-top: inverse_problems/images/thumb/sphx_glr_epidermis_change_thumb.png
      :img-alt: Fluence field in a sample used for a synthetic epidermis fit

      Recover a model parameter from ideal simulated measurements using SciPy;
      the same forward model generates and fits the data.

.. only:: html

   Plotly figures embedded in the example pages are interactive: rotate, zoom,
   and hover to inspect the generated MomentRT data.

.. toctree::
   :hidden:

   getting_started/index
   samples_and_meshes/index
   sources_and_detectors/index
   spectral/index
   transport_models/index
   parameter_studies/index
   validation/index
   inverse_problems/index
