Contributing to MomentRT
========================

MomentRT combines numerical transport models with explicit physical inputs.
Contributions should make behavior reproducible and document its scientific
meaning. Follow the repository's
`AGENTS.md <https://github.com/MartinPdeS/MomentRT/blob/master/AGENTS.md>`_
for module ownership, units, Python style, and model conventions.

Development environment
-----------------------

Use Python 3.11 or newer. From a checkout, reuse ``.venv`` if it already exists;
otherwise create it with ``python3 -m venv .venv``. Install the declared extras
with that environment's interpreter:

.. code-block:: console

   .venv/bin/python -m pip install --upgrade pip
   .venv/bin/python -m pip install -e ".[testing,documentation,dev,release]"

``testing``, ``dev``, and ``release`` use the shared
``martinpdes-maintenance`` tooling; ``documentation`` supplies the Sphinx
configuration. Plotly is a runtime dependency, so no ``viz`` extra is needed.
On Debian or Ubuntu, the Gmsh wheel also needs the system runtime library
``libglu1-mesa``, installed with ``sudo apt-get install libglu1-mesa``.

On Windows, use ``py -m venv .venv`` when creating the environment, then replace
``.venv/bin/python`` with ``.venv\Scripts\python.exe`` in these commands.
GNU Make is optional: the module commands below can be run directly.

Checks before submitting a change
---------------------------------

.. code-block:: console

   make check typecheck docstring-check PYTHON=.venv/bin/python
   git diff --check

``make check`` runs Ruff linting, formatting checks, and pytest. The separate
``typecheck`` and ``docstring-check`` targets run mypy and public API
annotation/docstring validation. Their direct equivalents are:

.. code-block:: console

   .venv/bin/python -m ruff check MomentRT tests docs/examples tools
   .venv/bin/python -m ruff format --check MomentRT tests docs/examples tools
   .venv/bin/python -m pytest
   .venv/bin/python -m mypy MomentRT
   .venv/bin/python tools/check_docstrings.py

Tests under ``tests/analytical`` compare against mathematical or independent
references. ``tests/numerical`` covers discretization, geometry, and workflow
behavior; ``tests/packaging`` checks metadata. Select focused checks during
development, then run the appropriate broader checks before submitting.

Numerical changes need independent evidence for the affected observable.
For photon transport, cover layer crossings and collision/absorption
conventions as well as energy accounting. A repeatable seeded answer or a
small energy-balance residual alone does not establish correctness. Preserve
units, unresolved levels, negative outputs, and uncertainty diagnostics.

Documentation and reference data
--------------------------------

.. code-block:: console

   make docs PYTHON=.venv/bin/python

The documentation build executes changed Sphinx-Gallery scripts and treats
Sphinx warnings as errors. Examples belong in ``docs/examples`` and start
with an RST title and introduction. Use ``# %%`` cells to explain meaningful
steps and ``show=False`` for figures; examples must not open browsers or write
HTML as side effects. Inspect changed figures, including interactive behavior.

For a changed 3D gallery preview, install Playwright and Chromium if needed,
build the docs, run the preview generator, and rebuild:

.. code-block:: console

   .venv/bin/python -m pip install playwright
   .venv/bin/python -m playwright install chromium
   make docs PYTHON=.venv/bin/python
   .venv/bin/python tools/generate_gallery_thumbnails.py
   make docs PYTHON=.venv/bin/python

Include the generated preview with the change. Ordinary documentation builds
then use that image without requiring a browser.

The full sensor reference study is an explicit maintenance operation:

.. code-block:: console

   .venv/bin/python tools/validate_sensor_models.py

It regenerates the committed JSON and plots under
``docs/source/_static/images/sensor_validation``. Review changed reference
results and uncertainties, rerun their independent checks, and explain any
algorithm correction in the changelog and validation guide. Documentation
examples read the committed report rather than rerunning Monte Carlo.

API and compatibility policy
----------------------------

This policy applies to future releases; it does not relabel existing tags.

* Patch releases preserve documented public import paths and calling
  conventions. They contain fixes and documentation updates. Correcting
  numerical results or rejecting previously accepted invalid input can be a
  patch fix, but must be described in the changelog with its scientific impact.
* Minor releases add capabilities. During the 0.x series, an incompatible
  change requires a minor release and migration instructions describing the
  previous and replacement APIs, units, and result meaning.
* A deprecation must identify a replacement and a target minor release for
  removal. Keep the old API through at least one subsequent minor release;
  do not remove it in a patch release. Existing compatibility aliases remain
  supported until a deprecation is explicitly announced.
* Physical units, source/boundary signs, normalization, array locations, and
  result provenance are part of the public contract. Document changes to them
  as carefully as renamed parameters. Keep unsupported physical cases explicit.
* Numerical comparisons use justified tolerances and statistical uncertainty.
  Bit-for-bit equality across meshes, dependency versions, platforms, or
  Monte Carlo generator revisions is not guaranteed.

Add a concise entry under ``Unreleased`` in
`CHANGELOG.md <https://github.com/MartinPdeS/MomentRT/blob/master/CHANGELOG.md>`_
for user-visible changes. For scientific corrections, identify the affected
results and whether published data must be regenerated. Keep existing aliases
such as ``convergence_study`` and ``full_hemisphere_on_face`` working unless a
documented deprecation is part of the change.

Packaging and releases
----------------------

.. code-block:: console

   make package-check release-check PYTHON=.venv/bin/python

``package-check`` builds the source archive and wheel. For package changes,
also run ``tools/check_wheel.py`` with the newly built wheel's exact path;
it installs and imports the wheel in a temporary environment outside the
checkout. ``release-check`` validates version metadata and the changelog's
section for the current version.

Maintainers create releases from a clean, checked commit. To inspect the next
version without writing files, run:

.. code-block:: console

   .venv/bin/python -m martinpdes_maintenance.cli next patch

The CLI's ``tag vX.Y.Z --dry-run`` command previews metadata and changelog
changes for a chosen tag, also requiring a clean tree. Use the release tools
to align ``pyproject.toml``, ``_version.py``, the conda recipe, citation, and
Zenodo metadata. The configured changelog rolls ``Unreleased`` into the new
version and records the actual release date.

``make tag VERSION=vX.Y.Z`` creates a local release commit and annotated tag.
``make release patch``, ``make release minor``, and ``make release major``
also push the commit and tag and trigger publication. Run those publishing
commands only when a release is explicitly intended, then check the publishing
workflow and package artifacts.
