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

For finite samples, :class:`~MomentRT.SP3` and :class:`~MomentRT.P3` provide
higher-order moment approximations with different boundary and phase-function
requirements. See :doc:`theory`, :doc:`sp3`, and :doc:`p3` for their model
limits, and :doc:`finite_element` for the P1 finite-domain formulation.

Units and optical properties
-----------------------------

Create dimensional values using :data:`~MomentRT.ureg`. Absorption and reduced
scattering have units of inverse length; sample geometry, source positions,
wavelengths, and optical power carry their corresponding units. Refractive
index is a plain dimensionless number.

.. code-block:: python

   import MomentRT
   from MomentRT import Material, ureg

   material = Material(
       absorption=0.02 / ureg.millimeter,
       reduced_scattering=1.0 / ureg.millimeter,
       refractive_index=1.4,
       phase_function=MomentRT.phase_function.HenyeyGreenstein(),
   )

MomentRT converts quantities at numerical boundaries as required and keeps
units attached to returned fields and power values. Material values are
assumed to apply at the source wavelength. For tabulated wavelength-dependent
properties, use ``SpectralMaterial`` and the :doc:`spectral` workflow.
No tissue spectrum or chromophore model is inferred from a wavelength alone.

Infinite-medium workflow
------------------------

The compact example below evaluates the analytical point-source solution.
Choose an even number of grid points along each axis if the source lies at the
center, so that no mesh node coincides with the point-source singularity.

.. code-block:: python

   import numpy as np
   import MomentRT

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
       phase_function=MomentRT.phase_function.HenyeyGreenstein(),
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

The complete executable example is in the :doc:`fundamentals gallery
<gallery/index>`.

Finite-sample workflow
----------------------

For a bounded sample, first define its full extent and background, then add
named structures in precedence order. The following outline shows the order
of operations; the gallery example contains a complete runnable configuration.

.. code-block:: python

   sample = StructuredMedium(
       origin=(-4, -4, 0) * ureg.millimeter,
       size=(8, 8, 4) * ureg.millimeter,
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

   source = SurfaceSource.on_face(
       sample=sample,
       face="z_min",
       coordinates=np.array([-1, 0]) * ureg.millimeter,
       radius=0.5 * ureg.millimeter,
       wavelength=970 * ureg.nanometer,
       optical_power=1 * ureg.milliwatt,
   )

   detector = SurfaceDetector.full_hemisphere_on_face(
       sample=sample,
       face="z_min",
       coordinates=np.array([1, 0]) * ureg.millimeter,
       radius=0.5 * ureg.millimeter,
   )

   boundaries = Boundaries.uniform(
       condition=RobinBoundary(effective_reflectance=0.5),
   )

   simulation = Simulation(
       medium=sample,
       source=source,
       mesh=mesh,
       solver=Diffusion(boundaries=boundaries),
       detectors={"reflection": detector},
   )

   report = simulation.validate()

   print(report.summary())

   if not report.valid:
       raise ValueError("Correct the configuration errors before solving.")

   result = simulation.run()

   print(result.summary())

   measurement = result.detectors["reflection"]

   print(f"Collected power: {measurement.collected_power:~P}")
   print(f"Sampled aperture area: {measurement.sampled_area:~P}")
   print(f"Geometric aperture area: {measurement.geometric_area:~P}")

The finite sample, source footprint, mesh, and boundary coefficients must all
be compatible. The source face must use a Robin boundary, the source disk must
fit on the face, and the mesh must span the sample exactly. The accepted power
is diffuse power already coupled into the sample; MomentRT does not model a
collimated LED or its specular entry loss. The detector collects the modeled
escaping irradiance over the entire outward hemisphere; it does not model
external optics. Preflight validation reports configuration errors and footprint
sampling but does not establish mesh convergence or transport accuracy. Review
:doc:`finite_element` before interpreting the predicted field.

The full surface-source example and the sample-construction example are in
the :doc:`finite-sample gallery <gallery/index>`.

Explore fields and diagnostics
------------------------------

The returned :class:`~MomentRT.SimulationResult` exposes the nodal fluence,
optional element-centered diffusive flux, energy accounting, per-face escaping
power, and solver diagnostics. For finite samples, inspect incoming,
absorbed, and escaping power together. A small residual checks numerical
accounting; it does not validate P1 against the radiative transfer equation.

The fluence plot defaults to a Plotly volume rendering with an in-figure
opacity slider. Use ``view="slices"`` to inspect central grid planes or
``view="isosurfaces"`` to show equal-fluence surfaces. Set
``isosurface_count`` to control how many levels these 3D views render. See
:doc:`results` for the meaning of each field and view.
