Development roadmap
===================

1. Initial baseline (implemented)
   Explicit-unit material, ideal point source, Cartesian observation mesh,
   infinite-medium diffusion, nodal fluence, optional 3D visualization,
   VTK export, analytical energy-balance and rendering checks.

2. Physical geometry and fields
   Layered materials and box/sphere construction are implemented for geometry
   and nodal optical-property sampling. Layer stacking and resolved planar
   material/property sections, surface-source geometry, and element-centered
   vector flux are also implemented. Circular contact detector collection
   includes explicit angular acceptance under a Lambertian escape assumption.
   Next: configurable field slices, shared color scales, and external optics.

3. Boundary-aware transport
   Linear tetrahedral FEM diffusion, explicit Robin/reflecting face conditions,
   and diffuse disk/full-face sources are implemented. Independent analytical
   slab comparisons cover homogeneous and equal-index layered media.
   SP3 and full 3D P3 provide coupled moments with Henyey–Greenstein scattering
   and matched-index boundaries. Next: refracting interfaces and comparisons against independent
   Monte Carlo results across sensor wavelengths and source-detector distances.
   Do not substitute infinite-medium fields for surface reflectance.

4. Spectral sensor model
   Tabulated wavelength-dependent materials and independent monochromatic
   sweeps with fixed source and detector geometry are implemented.
   Next: provenance-tracked chromophore data, measured LED spectra, detector
   responsivity, broadband spectral integration, calibration and acquisition
   noise.

5. Inference
   Explicit parameter references, unit-aware sensitivities, nuisance parameters,
   whitened Jacobian analysis and independent synthetic recovery studies.

API principles
--------------

Keep the package flat, use keyword-only calls and explicit quantities, and
return structured results. Keep numerical models separate from experimental
configuration. Use type hints and NumPy-style API docstrings. Transport fields retain physical
units and mesh association; rendering is optional and must not change values.
