Getting started
===============

This page takes a calculation from optical properties to a plotted field.
MomentRT separates the sample, source, observation mesh, and transport model,
then combines them in a :class:`~MomentRT.Simulation`.

Install
-------

For numerical calculations and interactive Plotly figures, install the
editable project with its visualization extra:

.. code-block:: console

   python -m pip install -e '.[viz]'

For the executable documentation gallery, install the documentation extra and
build the docs:

.. code-block:: console

   python -m pip install -e '.[documentation]'
   make docs PYTHON=.venv/bin/python

The base package requires NumPy, SciPy, Pint, and Matplotlib. Plotly is an
optional dependency for interactive figures. PyVista is an optional extra for
VTK image export; it is not needed to solve transport problems.

Choose the physical model first
-------------------------------

Use :class:`~MomentRT.InfiniteDiffusion` for a homogeneous, unbounded medium
with an embedded isotropic point source. Its Cartesian grid only selects where
the analytical field is evaluated; the grid edges are not physical boundaries.

Use :class:`~MomentRT.Diffusion` for a finite :class:`~MomentRT.StructuredMedium`
illuminated on one of its faces. This model solves the steady P1 diffusion
equation with linear tetrahedral finite elements and explicitly configured
conditions on all six faces. It currently requires a common refractive index
throughout the sample.

Both models are diffusion approximations. See :doc:`theory` for the moment
equations and their validity limits, and :doc:`finite_element` for the
finite-domain boundary model and its weak form.

Units and optical properties
-----------------------------

Create dimensional values using :data:`~MomentRT.ureg`. Absorption and reduced
scattering have units of inverse length; sample geometry, source positions,
wavelengths, and optical power carry their corresponding units. Refractive
index is a plain dimensionless number.

.. code-block:: python

   from MomentRT import Material, ureg

   material = Material(
       absorption=0.02 / ureg.millimeter,
       reduced_scattering=1.0 / ureg.millimeter,
       refractive_index=1.4,
   )

MomentRT converts quantities at numerical boundaries as required and keeps
units attached to returned fields and power values. Material values are
assumed to apply at the source wavelength. No wavelength-dependent material
law or chromophore model is inferred.

Infinite-medium workflow
------------------------

The compact example below evaluates the analytical point-source solution.
Choose an even number of grid points along each axis if the source lies at the
center, so that no mesh node coincides with the point-source singularity.

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

The example is also available as an executable gallery page:
:doc:`auto_examples/infinite_diffusion`.

Finite-sample workflow
----------------------

For a bounded sample, first define its full extent and background, then add
named structures in precedence order. The following outline shows the order
of operations; the gallery example contains a complete runnable configuration.

.. code-block:: python

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

   result = simulation.run()

The finite sample, source footprint, mesh, and boundary coefficients must all
be compatible. The source face must use a Robin boundary, the source disk must
fit on the face, and the mesh must span the sample exactly. The accepted power
is diffuse power already coupled into the sample; MomentRT does not model a
collimated LED or its specular entry loss. Review :doc:`finite_element` before
interpreting the predicted field.

The full surface-source example is :doc:`auto_examples/surface_diffusion`.
For geometry construction and resolved material cross-sections, see
:doc:`auto_examples/structured_sample`.

Explore fields and diagnostics
------------------------------

The returned :class:`~MomentRT.SimulationResult` exposes the nodal fluence,
optional element-centered diffusive flux, energy accounting, per-face escaping
power, and solver diagnostics. For finite samples, inspect incoming,
absorbed, and escaping power together. A small residual checks numerical
accounting; it does not validate P1 against the radiative transfer equation.

The fluence plot defaults to a Plotly volume rendering with an in-figure
opacity slider. Use ``view="slices"`` to inspect central grid planes or
``view="isosurfaces"`` to show equal-fluence surfaces. See :doc:`results` for
the meaning of each field and view.
