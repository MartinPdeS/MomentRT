Examples
========

Browse worked examples, from the analytical infinite-medium solution to
finite-sample FEM. Each card opens an executable example with its complete
configuration and generated plots. Input values are illustrative.

Fundamentals
------------

Start with the infinite-medium solver, then change one optical property to see
how it affects the predicted field.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Infinite-medium diffusion
      :link: fundamentals/infinite_diffusion
      :link-type: doc
      :img-top: fundamentals/images/thumb/sphx_glr_infinite_diffusion_thumb.png
      :img-alt: Infinite-medium MomentRT fluence field

      Configure a point source and observation mesh, solve the analytical P1
      Green function, and explore its 3D fluence field.

   .. grid-item-card:: How absorption changes fluence
      :link: fundamentals/absorption_sweep
      :link-type: doc
      :img-top: fundamentals/images/thumb/sphx_glr_absorption_sweep_thumb.png
      :img-alt: Three-dimensional infinite-medium fluence for the highest absorption

      Compare fluence profiles while holding reduced scattering and the source
      fixed.

Finite samples
--------------

Construct geometry, illuminate a finite sample, and use the examples to
understand the influence of source profile, boundary reflectance, and mesh
resolution.

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Build a heterogeneous sample
      :link: finite_samples/structured_sample
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_structured_sample_thumb.png
      :img-alt: Fluence volume in the layered sample with inclusions

      Combine three layers, two spheres, and a box, solve the illuminated
      sample, and explore its 3D fluence volume.

   .. grid-item-card:: Circular surface illumination
      :link: finite_samples/surface_diffusion
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_surface_diffusion_thumb.png
      :img-alt: Fluence produced by a diffuse circular surface source

      Run the finite-element workflow with a disk source and inspect fluence,
      flux, and energy accounting.

   .. grid-item-card:: Uniform face illumination
      :link: finite_samples/uniform_surface_source
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_uniform_surface_source_thumb.png
      :img-alt: Finite-sample fluence under uniform face illumination

      Use ``UniformSurfaceSource`` to illuminate an entire sample face and
      render the resulting field in 3D.

   .. grid-item-card:: Change boundary reflectance
      :link: finite_samples/reflectance_sweep
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_reflectance_sweep_thumb.png
      :img-alt: Fluence volume for an intermediate effective boundary reflectance

      Compare centerline fluence and escaping power as the effective diffuse
      reflectance changes.

   .. grid-item-card:: Check mesh refinement
      :link: finite_samples/mesh_refinement
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_mesh_refinement_thumb.png
      :img-alt: Finite-sample fluence volume on the finest mesh

      Compare a fixed sample and source over three grids to assess this
      observable's sensitivity to discretization.

   .. grid-item-card:: Collect light with finite apertures
      :link: finite_samples/surface_detector
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_surface_detector_thumb.png
      :img-alt: Finite-sample fluence volume with a source and detector apertures

      Define detector positions and acceptance cones, inspect the 3D system,
      and calculate reflection and transmission signals.

   .. grid-item-card:: Sweep a spectral sensor
      :link: finite_samples/spectral_detectors
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_spectral_detectors_thumb.png
      :img-alt: Monochromatic fluence volume from the spectral sensor sweep

      Interpolate two material spectra and plot collected power across five
      wavelengths with a fixed geometry.

   .. grid-item-card:: Compare P1 and SP3
      :link: finite_samples/compare_sp3
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_compare_sp3_thumb.png
      :img-alt: SP3 fluence volume with contact source and detector

      Compare detector signals from the two moment models on the same sample.

   .. grid-item-card:: Compare scattering phase functions
      :link: finite_samples/phase_functions
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_phase_functions_thumb.png
      :img-alt: SP3 fluence for the Rayleigh scattering example

      Compare phase functions with the same anisotropy but different higher moments.

   .. grid-item-card:: Compare P1, SP3, and full P3
      :link: finite_samples/compare_p3
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_compare_p3_thumb.png
      :img-alt: Full P3 fluence with a source and contact detector

      Solve the same sensor geometry with three angular approximations.

   .. grid-item-card:: Seven-fiber skin probe
      :link: finite_samples/seven_fiber_skin
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_seven_fiber_skin_thumb.png
      :img-alt: Six illumination fibers and a central detector on layered skin

      Explore a two-layer sample and a six-around-one contact fiber arrangement.

   .. grid-item-card:: Assign anisotropy by layer
      :link: finite_samples/sp3_layered_anisotropy
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_sp3_layered_anisotropy_thumb.png
      :img-alt: SP3 fluence for a layered sample with a source and detector

      Give each material its own scattering anisotropy and solve with SP3.

   .. grid-item-card:: Sweep spectral anisotropy
      :link: finite_samples/sp3_spectral_anisotropy
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_sp3_spectral_anisotropy_thumb.png
      :img-alt: SP3 fluence from wavelength-dependent material tables

      Interpolate layer-specific anisotropy and detector power across wavelengths.

   .. grid-item-card:: Choose a detector mesh
      :link: finite_samples/detector_convergence
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_detector_convergence_thumb.png
      :img-alt: Finest-mesh fluence volume with a source and detector

      Compare detector signals and sampled source and aperture areas, then
      inspect the finest fluence volume with instrument overlays.

.. only:: html

   Plotly figures embedded in the example pages are interactive: rotate, zoom,
   and hover to inspect the generated MomentRT data.

.. toctree::
   :hidden:

   fundamentals/index
   finite_samples/index
