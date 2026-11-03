# Changelog

User-visible changes and scientific corrections are recorded here. Historical
entries start at 0.9.0; earlier changes remain available in the repository tags.
The compatibility and deprecation policy is in [CONTRIBUTING.rst](CONTRIBUTING.rst).

## [Unreleased]

## [0.10.0] - 2026-11-03

### Fixed

- Correct the standalone Monte Carlo sensor reference: stop at every layer
  interface and sample scattering events using the scattering coefficient while
  applying continuous absorption along traveled segments.
- Regenerate the sensor reference JSON and error plots using reference revision
  2. Earlier reference data distributed with 0.9.0 and 0.9.1 are superseded and
  must not be used to judge P1, SP3, or P3 accuracy. The correction concerns the
  independent reference generator; it does not change the moment solvers.
- Remove obsolete `viz` extras from the README's setup commands and use the
  development environment's interpreter consistently.

### Added

- `MeshConvergence` studies with incremental progress tables, timing,
  diagnostic assessments, and shaded power and footprint tolerances.
- Executable inverse and mesh-convergence examples for a layered sample,
  including a documented synthetic fit and numerical convergence limits.
- Public `P1` and `InfiniteP1` solver names and dedicated transport modules.

- Independent reference checks for absorbing-layer transmission, unscattered
  transmission, a layered isotropic RTE solution, artificial-interface
  invariance, Henyey–Greenstein moments, and all-face energy accounting.
- Reference diagnostics for face escape, unscattered power, absorption, and
  the uncertainty of the correlated energy-balance residual.
- A contributor guide covering setup, verification, documentation, packaging,
  releases, and the API compatibility policy.
- Changelog integration with the shared maintenance release tools.


### Changed

- Organize result containers into a `results` package while retaining public
  imports. Preserve `Diffusion` and `InfiniteDiffusion` as P1 aliases.
- Draw all tetrahedral and Cartesian mesh overlay edges instead of sampling
  a fixed number; very fine mesh overlays can take longer to render.

## [0.9.1]

### Added

- Source spectral densities and detector responsivity on the existing source and
  detector classes, with explicit broadband optical power and photocurrent
  integration and retained monochromatic reference solves.
- Spectrum plotting and resampling, shared face-placement conventions, and
  direct aperture-area access on broadband measurements.
- `mesh_convergence` and `spectral_convergence`, including explicit refinement
  grids, failure records, and recommendations from stable final refinements.
- Executable solver-scaling and expanded sensor-validation examples.

### Compatibility

- Keep `convergence_study` and `full_hemisphere_on_face` as supported entry points.
- Preserve monochromatic source configuration and independent wavelength sweeps.
- The Monte Carlo reference dataset in this release is superseded by the
  correction under Unreleased above.

## [0.9.0]

### Added

- Detector convergence studies using tetrahedral target resolutions and contact
  refinement, with sampled footprint and detector-power diagnostics.
- A P1, SP3, and P3 solver-scaling benchmark utility.
- An initial layered-sensor Monte Carlo comparison and published reference
  report. Its reference data are superseded by the correction under Unreleased.

[Unreleased]: https://github.com/MartinPdeS/MomentRT/compare/v0.10.0...HEAD
[0.10.0]: https://github.com/MartinPdeS/MomentRT/releases/tag/v0.10.0
[0.9.1]: https://github.com/MartinPdeS/MomentRT/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/MartinPdeS/MomentRT/compare/v0.8.4...v0.9.0
