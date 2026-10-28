.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Badge
     - Status
   * - Python versions
     - |python|
   * - Documentation
     - |docs|
   * - Tests
     - |tests|
   * - Coverage workflow
     - |coverage|
   * - PyPI package
     - |pypi|
   * - PyPI downloads
     - |pypi_downloads|
   * - Anaconda package
     - |anaconda|
   * - Anaconda downloads
     - |anaconda_downloads|
   * - Latest Anaconda release
     - |anaconda_date|

MomentRT
========

**MomentRT** is an open-source Python package for deterministic light
transport using angular-moment methods. It models steady-state P1 diffusion
in homogeneous infinite media and finite structured samples. Its API keeps the
sample, illumination, mesh, and solver explicit, and its results include
unit-bearing optical fields and power diagnostics.

Features
--------

* Analytical P1 diffusion from an isotropic point source in an infinite,
  homogeneous medium.
* Boundary-aware P1 diffusion for finite samples, solved with linear
  tetrahedral finite elements.
* Layered samples with box and sphere structures, plus resolved material
  cross-sections.
* Explicit diffuse surface illumination and Robin or reflecting conditions
  on all six sample faces.
* Interactive Plotly views for sample geometry, fluence, sections, and
  element-centered flux. Volume plots include an opacity control.
* Unit-bearing fields, per-face escaping power, and solver diagnostics.

MomentRT uses the P1 diffusion approximation. It does not currently solve the
full radiative transfer equation or implement P3/SP3, collimated LED
propagation, detector apertures, spectral integration, time dependence, or
inverse estimation. Read the `scientific guide`_ for the equations,
approximations, validation limits, and complete workflow.

Installation
------------

Install the released Python package with Plotly visualization support:

.. code-block:: console

   python -m pip install "MomentRT[viz]"

Or install the conda package and its dependencies:

.. code-block:: console

   conda install -c martinpdes -c conda-forge momentrt

Verify the installation with the same Python interpreter you will use for your
calculations:

.. code-block:: console

   python -c "import MomentRT; print(MomentRT.__version__)"

The base package includes the numerical dependencies. Plotly is optional for
interactive plots; PyVista is optional for VTK export. Install both plotting
extras with ``python -m pip install "MomentRT[viz,vtk]"``.

First simulation
----------------

The infinite-medium solver is a quick way to calculate and inspect a fluence
field. This example uses explicit units and places the source between grid
nodes, avoiding the singular point-source location.

.. code-block:: python

   import numpy as np

   from MomentRT import (
       CartesianMesh,
       HomogeneousMedium,
       InfiniteDiffusion,
       Material,
       PointSource,
       Simulation,
       ureg,
   )

   material = Material(
       absorption=0.02 / ureg.millimeter,
       reduced_scattering=1.0 / ureg.millimeter,
       refractive_index=1.4,
   )

   medium = HomogeneousMedium(material=material)

   mesh = CartesianMesh(
       shape=(40, 40, 40),
       spacing=np.full(3, 0.4) * ureg.millimeter,
       origin=np.full(3, -7.8) * ureg.millimeter,
   )

   source = PointSource(
       position=np.zeros(3) * ureg.millimeter,
       wavelength=970 * ureg.nanometer,
       optical_power=1 * ureg.milliwatt,
   )

   solver = InfiniteDiffusion()

   simulation = Simulation(
       medium=medium,
       source=source,
       mesh=mesh,
       solver=solver,
   )

   result = simulation.run()

   figure = result.fluence.plot(
       scale="log",
   )

``result.fluence`` is a nodal fluence-rate field with units of W/m². The
infinite solver evaluates an analytical Green function; the observation grid
does not create physical sample boundaries.

Finite-sample FEM workflow
--------------------------

For a finite sample, define its geometry and materials, create a mesh, specify
surface illumination and all six boundary faces, then run ``Diffusion``. This
example uses the same explicit effective reflectance on each face and a
circular diffuse source on the upper face.

.. code-block:: python

   import numpy as np

   from MomentRT import (
       CartesianMesh,
       Diffusion,
       Material,
       RobinBoundary,
       Simulation,
       StructuredMedium,
       SurfaceSource,
       ureg,
   )

   background = Material(
       absorption=0.02 / ureg.millimeter,
       reduced_scattering=1.0 / ureg.millimeter,
       refractive_index=1.4,
   )

   upper_material = Material(
       absorption=0.04 / ureg.millimeter,
       reduced_scattering=1.5 / ureg.millimeter,
       refractive_index=1.4,
   )

   sample = StructuredMedium(
       origin=np.array([-4, -4, 0]) * ureg.millimeter,
       size=np.array([8, 8, 4]) * ureg.millimeter,
   )

   sample.add_background(material=background)

   sample.add_layer(
       name="Upper layer",
       thickness=1 * ureg.millimeter,
       material=upper_material,
   )

   mesh = CartesianMesh.from_sample(
       sample=sample,
       shape=(25, 25, 17),
   )

   source = SurfaceSource(
       face="z_min",
       position=np.array([0, 0, 0]) * ureg.millimeter,
       radius=1 * ureg.millimeter,
       wavelength=970 * ureg.nanometer,
       optical_power=1 * ureg.milliwatt,
   )

   boundary = RobinBoundary(effective_reflectance=0.5)

   solver = Diffusion(
       boundaries={
           "x_min": boundary,
           "x_max": boundary,
           "y_min": boundary,
           "y_max": boundary,
           "z_min": boundary,
           "z_max": boundary,
       },
   )

   simulation = Simulation(
       medium=sample,
       source=source,
       mesh=mesh,
       solver=solver,
   )

   geometry = simulation.plot_system(
       show=False,
   )

   result = simulation.run()

   figure = result.fluence.plot(
       view="volume",
       scale="log",
   )

   print(f"Incoming power: {result.energy_balance['incoming_power']}")
   print(f"Absorbed power: {result.energy_balance['absorbed_power']}")
   print(f"Escaping power: {result.energy_balance['escaping_power']}")

The finite ``Diffusion`` solver uses P1 diffusion with linear tetrahedral FEM.
Its sample mesh must span the sample exactly. Internal refractive indices
must be equal because the solver models continuous fluence and does not
resolve refracting interfaces. ``SurfaceSource.optical_power`` is diffuse
power already coupled into the sample; incident beam propagation and
specular entry losses are not included. The supplied effective reflectance
is a diffuse partial-current coefficient, not a Fresnel value inferred from
refractive index.

Visualization
-------------

The default fluence plot is a Plotly volume view. Adjust opacity from 0 to 1
using its browser slider, or choose a starting value with ``opacity``:

.. code-block:: python

   figure = result.fluence.plot(
       view="volume",
       scale="log",
       opacity=0.12,
       show=False,
   )

   figure.write_html("fluence.html")

Other scalar-field views include central mesh ``slices`` and
``isosurfaces``. ``result.flux.plot()`` displays element-centered diffusive
flux vectors. ``sample.plot()`` shows the transparent construction surfaces;
``sample.plot_section(...)`` displays a sampled material or optical-property
map. Log color scales change the displayed colors only; hover values retain
their physical units.

The following figures were generated from MomentRT examples and results.
The live documentation contains interactive Plotly versions.

.. image:: https://raw.githubusercontent.com/MartinPdeS/MomentRT/master/docs/source/_static/images/sample-geometry.png
   :width: 900
   :align: center
   :alt: MomentRT view of a layered sample with a spherical inclusion.

.. image:: https://raw.githubusercontent.com/MartinPdeS/MomentRT/master/docs/source/_static/images/infinite-fluence-volume.png
   :width: 900
   :align: center
   :alt: Logarithmic volume rendering of fluence around an isotropic point source.

.. image:: https://raw.githubusercontent.com/MartinPdeS/MomentRT/master/docs/source/_static/images/finite-fluence-slices.png
   :width: 900
   :align: center
   :alt: Logarithmic fluence on three central planes in a finite sample.

Documentation
-------------

The `online documentation <https://martinpdes.github.io/MomentRT/>`_ includes:

* The P1 moment equations, infinite-medium Green function, and approximation
  limits.
* The finite-sample Robin boundary condition, weak form, and FEM assembly.
* Complete workflows, visualization guidance, results, and API reference.
* An executable `example gallery
  <https://martinpdes.github.io/MomentRT/docs/latest/auto_examples/index.html>`_.

Build the docs and run the examples locally with:

.. code-block:: console

   python -m pip install -e ".[documentation,viz]"
   make docs PYTHON=.venv/bin/python

Common issues
-------------

* **Point-source singularity:** do not place an observation node exactly at
  the point source. MomentRT reports this instead of silently regularizing it.
* **No visible internal interface:** construction plots use transparent
  geometry surfaces. Use a sampled section to inspect the resolved material
  precedence.
* **Finite solver rejects refractive indices:** the current FEM formulation
  requires equal internal refractive indices; it does not model internal
  refraction.
* **Unexpected surface-source power:** source power is already coupled diffuse
  power in the sample. Collimated illumination and specular entry losses are
  outside the model.
* **Slow or changing fields under refinement:** element-centroid materials,
  source quadrature, interface alignment, and P1 validity can all affect the
  result. Refine meshes and compare independent references; energy balance
  alone does not validate the physical approximation.
* **Log plot rejects values:** logarithmic rendering requires strictly
  positive values and does not clip zeros or negative values.

Developer setup
---------------

Clone the repository, install the development extras, and build the docs with
the same Python interpreter used for development:

.. code-block:: console

   git clone https://github.com/MartinPdeS/MomentRT.git
   cd MomentRT
   python3 -m venv .venv
   .venv/bin/python -m pip install -e ".[testing,documentation,dev,viz]"
   make check PYTHON=.venv/bin/python
   make docs PYTHON=.venv/bin/python

``make check`` runs quality checks and tests. The codebase intentionally keeps
the package flat and does not use Python type annotations.

Citation and license
--------------------

MomentRT is distributed under the MIT license. If you use it in research,
cite the version recorded in ``CITATION.cff``. For bibliographic metadata,
see the `CITATION.cff file
<https://github.com/MartinPdeS/MomentRT/blob/master/CITATION.cff>`_.

.. _scientific guide: https://martinpdes.github.io/MomentRT/docs/latest/theory.html

.. |python| image:: https://img.shields.io/pypi/pyversions/MomentRT.svg
   :alt: Python versions
   :target: https://pypi.org/project/MomentRT/
.. |docs| image:: https://github.com/MartinPdeS/MomentRT/actions/workflows/deploy_documentation.yml/badge.svg?branch=master
   :alt: Documentation build status
   :target: https://martinpdes.github.io/MomentRT/
.. |tests| image:: https://github.com/MartinPdeS/MomentRT/actions/workflows/tests.yml/badge.svg?branch=master
   :alt: Test workflow status
   :target: https://github.com/MartinPdeS/MomentRT/actions/workflows/tests.yml
.. |coverage| image:: https://github.com/MartinPdeS/MomentRT/actions/workflows/deploy_coverage.yml/badge.svg?branch=master
   :alt: Coverage workflow status
   :target: https://github.com/MartinPdeS/MomentRT/actions/workflows/deploy_coverage.yml
.. |pypi| image:: https://img.shields.io/pypi/v/MomentRT.svg
   :alt: PyPI version
   :target: https://pypi.org/project/MomentRT/
.. |pypi_downloads| image:: https://img.shields.io/pypi/dm/MomentRT.svg
   :alt: Monthly PyPI downloads
   :target: https://pypi.org/project/MomentRT/
.. |anaconda| image:: https://anaconda.org/martinpdes/momentrt/badges/version.svg
   :alt: Anaconda version
   :target: https://anaconda.org/martinpdes/momentrt
.. |anaconda_downloads| image:: https://anaconda.org/martinpdes/momentrt/badges/downloads.svg
   :alt: Anaconda downloads
   :target: https://anaconda.org/martinpdes/momentrt
.. |anaconda_date| image:: https://anaconda.org/martinpdes/momentrt/badges/latest_release_relative_date.svg
   :alt: Latest Anaconda release date
   :target: https://anaconda.org/martinpdes/momentrt
