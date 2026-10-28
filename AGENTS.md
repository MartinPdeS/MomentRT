# Working on MomentRT

MomentRT is a Python package for deterministic light transport using angular-moment methods. Its first application is a contact optical sensor, but the transport API should remain useful for other media and instruments.

## Package and API

- Keep `MomentRT/` flat. Define physical objects, numerical models, results, and visualization in their existing modules; expose the main public classes through `MomentRT/__init__.py`.
- Keep material properties in `material.py`, sample composition in `medium.py`, geometry in `geometry.py`, sources in `source.py`, boundary conditions in `boundary.py`, meshes in `mesh.py`, transport models in `diffusion.py`, and FEM kernels in `fem.py`. Keep orchestration in `simulation.py`, result containers in `results.py`, fields in `fields.py`, sections in `sections.py`, and plotting in `plotting.py`.
- Preserve a clear workflow: construct named physical objects, configure the solver separately, construct a `Simulation`, and then call `result = simulation.run()`. Keep sample, source, boundary, and solver choices explicit and independently replaceable.
- Require keyword-only arguments for constructors, functions, and methods, apart from implicit `self` and `cls`. `add_structures(*structures)` is the explicit exception and accepts positional shape objects.
- Do not add Python type annotations, stub files, or `from __future__ import annotations`.
- Use four spaces and Ruff formatting. Format calls with multiple arguments over multiple lines, one argument per line, with a trailing comma. Separate statements and calculation steps with blank lines; keep imports and related comments grouped.
- Use specific physical names such as `refractive_index`, `effective_reflectance`, `reduced_scattering`, and `optical_power`. Use `index` for array positions. Print calculated values with labelled f-strings.
- Avoid arbitrary defaults for physical configurations. Wavelengths, optical properties, sample dimensions, sources, and boundary conditions require explicit values. Numerical controls and neutral coordinate conventions may have defaults.

## Units and sample geometry

- Require unit-bearing quantities from `MomentRT.ureg` for dimensional inputs. Validate supplied quantities directly with `MomentRT.units.validate_units` and preserve them; do not reassign fields merely to validate them. Construct quantities by multiplying values by units and use `.to(...)` for conversions.
- Refractive indices and effective reflectance are plain dimensionless numbers, not quantities. Counts and other dimensionless controls may be bare numbers. Extract magnitudes only at numerical, plotting, or serialization boundaries and preserve units on results.
- `StructuredMedium` is a finite rectangular sample with an explicit background. `add_background` replaces the background and keeps existing structures; `add_structures` appends a batch in order and validates it before mutation. Later structures replace earlier materials, including at shared interfaces.
- Layers span the sample's lateral extent and use absolute z depth. `add_layer` stacks below the deepest existing layer. The sample's positive z direction points into the sample. Keep geometry, construction plots, and material sampling consistent with these conventions.
- Distinguish transparent 3D construction surfaces from the resolved material map. Sections sample cell centers at an explicit resolution; sampled fields are snapshots. Do not present visual blending or a coarse sampling as exact interface geometry.

## Transport and scientific meaning

- `InfiniteDiffusion` is the analytical steady-state P1 solution for an embedded point source in an infinite homogeneous medium. Its Cartesian mesh is an observation window, not a physical boundary. Keep the point-source singularity explicit.
- `Diffusion` is a steady-state, continuous-fluence P1 FEM solver on a finite `StructuredMedium`. It uses linear tetrahedra, centroid-sampled material coefficients, and explicitly configured conditions on all six faces. Unaligned interfaces and finite source footprints need mesh-convergence checks. Internal refractive-index jumps are currently rejected.
- `RobinBoundary.effective_reflectance` is an explicitly supplied diffuse partial-current coefficient, not normal-incidence Fresnel reflectance. `ReflectingBoundary` has zero normal current. Preserve the documented source and boundary sign convention and do not silently infer a refractive-index mismatch model.
- Surface-source `optical_power` is incoming diffuse power already coupled into the sample. It does not include collimated LED propagation or specular entry losses. The disk quadrature is normalized to that power; report its sampled area and do not treat energy conservation as proof that the source footprint is resolved.
- Keep incoming power, absorption, net current, and escaping transmitted power distinct. Per-face escaping power is not a finite-aperture detector measurement. Preserve energy-balance diagnostics, residuals, physical units, and field provenance.
- Fluence is nodal; FEM flux is element-centered and is not smoothed across interfaces. Do not silently clip negative values or conceal numerical artifacts. A logarithmic display must reject nonpositive data unless an explicit, documented policy is introduced.
- Clearly state each approximation and its validity limits. P1/diffusion, SPN, full PN, and radiative transfer are different models. Do not present planned SP3, detector, spectral, or inverse features as implemented, or claim that agreement with diffusion references proves accuracy for skin or short source-detector separations.

## Visualization and documentation

- Use Plotly for optional interactive 3D and section plots. Import it lazily and return composable figures. Keep PyVista optional and lazily imported for VTK export; numerical calculations must work without either visualization extra.
- Label plot axes, coordinates, color scales, and physical units. Preserve physical values on hover and distinguish transformed display values from stored data. Show source and sample geometry consistently with the modeled configuration.
- `README.rst` and `docs/source` use reStructuredText. Keep examples in `docs/examples` executable Sphinx-Gallery scripts: start with an RST title and introduction, explain meaningful steps in `# %%` comment cells, and execute at top level.
- Gallery examples construct figures with `show=False`; they must not open browsers or write HTML as side effects. The gallery scraper embeds each figure. Keep example values clearly illustrative, and avoid unrelated institution names.
- Keep the MIT license consistent across `LICENSE`, package metadata, the conda recipe, `CITATION.cff`, and `.zenodo.json`; include the license in distributable archives. Do not invent release dates, DOI values, or repository URLs.
- Keep version metadata aligned between `pyproject.toml`, `MomentRT/_version.py`, `conda.recipe/meta.yaml`, `CITATION.cff`, and `.zenodo.json` using the release tools. Keep the conda recipe's minimum Python and runtime requirements aligned with `pyproject.toml`; do not turn optional visualization dependencies into core runtime requirements. Use the release tools and do not commit, tag, or push unless explicitly requested. Preserve unrelated work.

## Verification

- Reuse `.venv`. Choose checks appropriate to the change; run focused validation for minor changes and `make check PYTHON=.venv/bin/python` when broader code changes warrant it. Run `git diff --check` when applicable. Tests live under `tests/`, including `analytical/` and `numerical/`.
- For numerical or physical changes, add meaningful independent references and validation: test units, boundary/source conventions, energy accounting, mesh and source resolution, and convergence as applicable. Energy balance alone does not validate a transport approximation.
- For documentation or packaging changes, run the relevant build or metadata checks when the change warrants them: `make docs PYTHON=.venv/bin/python`, `make package-check PYTHON=.venv/bin/python`, or `make release-check PYTHON=.venv/bin/python`. Inspect generated gallery output when examples or the scraper change.
- For visualization changes, check the rendered interactive figures in a real browser in addition to relevant automated tests. Verify coordinates, labels, units, and interaction without requiring visualization dependencies for core imports.
- Keep GitHub workflows aligned with the `master` branch and the MPSActions reusable workflow contracts. The conda recipe is noarch: publish a single artifact rather than a matrix of duplicate filenames. Validate workflow syntax and check wheel installation outside the checkout when changing CI.
