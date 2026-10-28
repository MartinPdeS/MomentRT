Examples
========

Browse seven worked examples, from the analytical infinite-medium solution to
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
      :img-alt: Radial fluence profiles for three absorption values

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
      :img-alt: Layered MomentRT sample with overlapping spherical and box regions

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
      :img-alt: Centerline fluence profiles under three effective reflectances

      Compare centerline fluence and escaping power as the effective diffuse
      reflectance changes.

   .. grid-item-card:: Check mesh refinement
      :link: finite_samples/mesh_refinement
      :link-type: doc
      :img-top: finite_samples/images/thumb/sphx_glr_mesh_refinement_thumb.png
      :img-alt: Fluence profiles computed on coarse, medium, and fine meshes

      Compare a fixed sample and source over three grids to assess this
      observable's sensitivity to discretization.

.. only:: html

   Plotly figures embedded in the example pages are interactive: rotate, zoom,
   and hover to inspect the generated MomentRT data.

.. toctree::
   :hidden:

   fundamentals/index
   finite_samples/index
