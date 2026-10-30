Choose a workflow
=================

MomentRT exposes steady-state P1, SP3, and full P3 workflows. Choose based on
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

   .. grid-item-card:: Compare P1 and SP3
      :link: sp3
      :link-type: doc

      Solve coupled SP3 moments with explicit scattering anisotropy and
      matched-index boundary conditions.

   .. grid-item-card:: Full three-dimensional P3
      :link: p3
      :link-type: doc

      Retain all angular moments through order three using six coupled
      even-parity fields and matched-index faces.

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

   .. grid-item-card:: Choose a detector mesh
      :link: convergence
      :link-type: doc

      Compare collected power, aperture sampling, and source area under refinement.

.. note::

   P1, SP3, and P3 are approximations to transport. No model infers boundary
   parameters from refractive index or replaces validation against transport
   references for the geometry of interest.

.. toctree::
   :hidden:

   Finite-aperture detectors <detectors>
   SP3 transport <sp3>
   Full P3 transport <p3>
   Materials and wavelength sweeps <spectral>
   Detector mesh convergence <convergence>
