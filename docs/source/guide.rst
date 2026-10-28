Choose a workflow
=================

MomentRT currently exposes two steady-state P1 workflows. Choose based on
whether the physical domain is unbounded and homogeneous or a finite sample
with an explicitly modeled illuminated boundary.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Infinite homogeneous medium
      :link: theory
      :link-type: doc

      Evaluate the analytical point-source Green function on a Cartesian
      observation grid. The grid is a viewing window, not a physical boundary.

   .. grid-item-card:: Finite rectangular sample
      :link: finite_element
      :link-type: doc

      Solve P1 diffusion with linear tetrahedral finite elements and explicit
      conditions on all six sample faces.

   .. grid-item-card:: Build the sample
      :link: /gallery/finite_samples/structured_sample
      :link-type: doc

      Define layers and inclusions, inspect geometry, and resolve material
      cross-sections.

   .. grid-item-card:: Inspect outputs
      :link: results
      :link-type: doc

      Explore fluence, flux, energy accounting, diagnostics, and 3D plotting.

   .. grid-item-card:: Collect detector power
      :link: detectors
      :link-type: doc

      Define contact apertures, collection axes, and angular acceptance.

   .. grid-item-card:: Sweep wavelengths
      :link: spectral
      :link-type: doc

      Interpolate optical tables and reuse a fixed instrument geometry.

.. note::

   Both models use the P1 diffusion approximation. The finite solver is not a
   full radiative-transfer model and does not infer boundary parameters from
   refractive index.

.. toctree::
   :hidden:

   Finite-aperture detectors <detectors>
   Materials and wavelength sweeps <spectral>
