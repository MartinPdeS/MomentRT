Examples
========

Browse worked examples, from the analytical infinite-medium solution through
finite-sample models, instruments, and numerical studies. Each card opens an
executable example with its complete configuration and generated plots. Input
values are illustrative.

Getting started
---------------

Learn the basic workflow, from a point source in an infinite medium to layered
finite samples and surface illumination.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Infinite-medium diffusion
      :link: getting_started/infinite_diffusion
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_infinite_diffusion_thumb.png
      :img-alt: Infinite-medium MomentRT fluence field

      Configure a point source and observation mesh, solve the analytical P1
      Green function, and explore its 3D fluence field.

   .. grid-item-card:: Build a heterogeneous sample
      :link: getting_started/structured_sample
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_structured_sample_thumb.png
      :img-alt: Fluence volume in a layered sample with inclusions

      Combine layers, spheres, and a box, solve the illuminated sample, and
      explore its 3D fluence volume.

   .. grid-item-card:: Circular surface illumination
      :link: getting_started/surface_diffusion
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_surface_diffusion_thumb.png
      :img-alt: Fluence produced by a diffuse circular surface source

      Inspect fluence, flux, and energy accounting for a finite-element solve
      with a disk source.

   .. grid-item-card:: Uniform face illumination
      :link: getting_started/uniform_surface_source
      :link-type: doc
      :img-top: getting_started/images/thumb/sphx_glr_uniform_surface_source_thumb.png
      :img-alt: Finite-sample fluence under uniform face illumination

      Use ``UniformSurfaceSource`` to illuminate an entire sample face and
      render the resulting field in 3D.

Sources and detectors
---------------------

Configure source footprints and detector apertures, including spectral
collection and a multi-fiber contact probe.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Collect light with finite apertures
      :link: sources_detectors/surface_detector
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_surface_detector_thumb.png
      :img-alt: Finite-sample fluence with source and detector apertures

      Define detector positions and acceptance cones, inspect the system, and
      calculate reflection and transmission signals.

   .. grid-item-card:: Sweep a spectral sensor
      :link: sources_detectors/spectral_detectors
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_spectral_detectors_thumb.png
      :img-alt: Monochromatic fluence from a spectral sensor sweep

      Interpolate material spectra and plot collected power across wavelengths
      with fixed sample and instrument geometry.

   .. grid-item-card:: Seven-fiber skin probe
      :link: sources_detectors/seven_fiber_skin
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_seven_fiber_skin_thumb.png
      :img-alt: Six illumination fibers and a central detector on layered skin

      Simulate simultaneous illumination from six fibers around a central
      collection fiber in a layered skin model.

   .. grid-item-card:: Check detector mesh convergence
      :link: sources_detectors/detector_convergence
      :link-type: doc
      :img-top: sources_detectors/images/thumb/sphx_glr_detector_convergence_thumb.png
      :img-alt: Finest-mesh fluence with source and detector overlays

      Compare detector signals and sampled source and aperture areas over
      successively finer meshes.

Transport models
----------------

Compare angular approximations and inspect how scattering phase functions and
anisotropy influence model outputs.

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

Sweeps and convergence
----------------------

Vary optical properties and wavelength, and assess how numerical resolution
affects calculated fields and observables.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: How absorption changes fluence
      :link: sweeps_convergence/absorption_sweep
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_absorption_sweep_thumb.png
      :img-alt: Infinite-medium fluence at several absorption values

      Compare fluence profiles while holding reduced scattering and source
      configuration fixed.

   .. grid-item-card:: Check finite-element mesh refinement
      :link: sweeps_convergence/mesh_refinement
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_mesh_refinement_thumb.png
      :img-alt: Finite-sample fluence on the finest mesh

      Compare a fixed sample and source over three grids to assess sensitivity
      to discretization.

   .. grid-item-card:: Change boundary reflectance
      :link: sweeps_convergence/reflectance_sweep
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_reflectance_sweep_thumb.png
      :img-alt: Fluence volume for an intermediate boundary reflectance

      Compare centerline fluence and escaping power as diffuse reflectance
      changes.

   .. grid-item-card:: Sweep spectral anisotropy
      :link: sweeps_convergence/sp3_spectral_anisotropy
      :link-type: doc
      :img-top: sweeps_convergence/images/thumb/sphx_glr_sp3_spectral_anisotropy_thumb.png
      :img-alt: SP3 fluence from wavelength-dependent material tables

      Interpolate layer-specific anisotropy and detector power across
      wavelengths.

.. only:: html

   Plotly figures embedded in the example pages are interactive: rotate, zoom,
   and hover to inspect the generated MomentRT data.

.. toctree::
   :hidden:

   getting_started/index
   sources_detectors/index
   transport_models/index
   sweeps_convergence/index
