MomentRT
========

Deterministic light transport using angular-moment methods.

Early development
-----------------

The analytical baseline is steady-state P1 diffusion from an
isotropic point source in an infinite homogeneous medium. It produces a
unit-bearing nodal fluence-rate field and interactive 3D visualization.
The Cartesian grid is an observation window, not a physical sample boundary.
The source singularity is rejected rather than regularized silently.
Diffusion is an approximation; near-source fields are not reliable transport
predictions. Refractive index is recorded but has no effect without interfaces.

Finite sample construction and visualization support layers, boxes, and spheres.
Boundary-aware tetrahedral FEM diffusion supports finite structured samples
and diffuse surface illumination. SP3, detectors, spectral integration, and
inverse estimation are planned and are not implemented yet.

Installation
------------

.. code-block:: console

   python -m pip install -e ".[viz,dev]"
   python docs/examples/infinite_diffusion.py

Conda package from a local checkout
-----------------------------------

The conda recipe builds the checked-out source as a pure-Python package.
With ``conda-build`` installed, build and install it locally with:

.. code-block:: console

   conda build -c conda-forge conda.recipe
   conda install --use-local -c conda-forge momentrt

The base conda package contains the numerical dependencies. Install Plotly
from conda-forge to use interactive plots; PyVista is needed only for VTK
export. The recipe does not publish a package.

Visualization
-------------

``simulation.plot_system()`` displays the observation window and source.
``result.fluence.plot(view="slices", scale="log")`` displays orthogonal slices.
Other views are ``volume`` and ``isosurfaces``. Plots return a Plotly figure;
use ``show=False`` to compose a scene and ``figure.write_html("fluence.html")``
to save a standalone interactive document.
Field export uses ``result.fluence.export(path="fluence.vti")`` and records
value units, coordinate units, and transport metadata for ParaView.
Plotly is optional and imported only when plotting. Install ``MomentRT[vtk]``
for optional PyVista-based VTK export. Plotting does not require VTK.
Log color scales transform display values only; hover labels retain physical
values. Slices select the central grid plane on each axis.

Development
-----------

.. code-block:: console

   make check PYTHON=.venv/bin/python
   make package-check PYTHON=.venv/bin/python

The package is flat. Python type annotations are intentionally omitted.

Reference
---------

The implemented infinite-medium Green function is
``Phi(r) = P exp(-sqrt(mu_a / D) r) / (4 pi D r)``, with
``D = 1 / (3 (mu_a + mu_s_prime))``. It solves
``-D Laplacian(Phi) + mu_a Phi = P delta(r)``.

Plotly volume rendering:
https://plotly.com/python/3d-volume-plots/

Defining a sample
-----------------

Build a finite ``StructuredMedium`` with explicit ``origin`` and ``size``
quantities (three Cartesian components). Configure its background with
``sample.add_background(material=material)``. Add named ``Layer``, ``Box``,
and ``Sphere`` objects using ``sample.add_structures(epidermis, dermis, inclusion)``.
These calls mutate the builder in place and return None.

A ``Layer`` uses an absolute z ``depth`` and positive ``thickness`` and spans
the sample's full lateral extent. A ``Box`` uses its lower-corner ``origin``
and full ``size``; a ``Sphere`` uses ``center`` and ``radius``. Positive z
points into the sample. Dimensional inputs require explicit units.
All structures must fit inside the sample; duplicate names are rejected.
Later structures replace earlier materials, including at shared interfaces.
Replacing the background preserves structures. Uncovered regions use the
background; the exterior is intentionally unspecified.

Run the complete layered-sample example:

.. code-block:: console

   python docs/examples/structured_sample.py

``figure = sample.plot(show=False)`` creates a Plotly figure without a mesh.
Rotate the sample, toggle regions through the legend, and hover over regions
to inspect their optical properties. Use ``opacity`` to expose internal shapes
and ``length_unit`` to change displayed coordinates. Save standalone HTML with
``figure.write_html("sample.html")``.

The plot shows transparent construction surfaces, not Boolean-subtracted
interfaces: insertion order defines the actual sampled material map, not
visual blending. ``sample.material_labels(points=points)`` returns 0 for the
background and 1..N for structures in insertion order. Query points must lie
inside the sample. ``sample.optical_fields(mesh=mesh)`` returns nodal
absorption, reduced-scattering, and refractive-index ``ScalarField`` objects.
Generated fields are independent snapshots.

Use ``Diffusion`` for finite-sample FEM transport. ``InfiniteDiffusion`` accepts
only ``HomogeneousMedium`` and does not simulate finite boundaries.

Stacking layers and inspecting cross-sections
---------------------------------------------

``sample.add_layer(name="Epidermis", thickness=200 * ureg.micrometer,
material=epidermal_material)`` appends a layer without a manually calculated
depth. The first starts at ``sample.origin[2]``; later layers start at the
deepest existing layer's lower interface. Explicit ``Layer`` regions also
participate; spheres and boxes do not change the stack depth. Overflow and
duplicate names are rejected before changing the sample. Layers added after
inclusions follow ordinary insertion precedence and can replace them.
Add the layer stack before inclusions when inclusions should take precedence.

.. code-block:: python

   figure = sample.plot_section(
       normal="y",
       position=0 * ureg.millimeter,
       quantity="material",
       resolution=(300, 150),
       show=False,
   )

   figure.write_html("sample_section.html")

``normal`` selects the plane normal; ``position`` is an absolute coordinate
with length units. The remaining axes occur in Cartesian order, with depth
increasing downward when the vertical axis is z. The plot has equal physical
axis scales. ``quantity`` may be ``material``, ``absorption``,
``reduced_scattering``, or ``refractive_index``. Optical properties support
``value_unit`` conversions and ``scale="log"`` for strictly positive values.
Material maps use discrete colors; no smoothing is applied.

For numerical access, use ``section = sample.section(...)``. It retains
``labels``, region ``names``, physical ``horizontal`` and ``vertical``
coordinates, and unit-bearing optical-property ``fields``. Arrays have
(vertical, horizontal) order. Resolution specifies (horizontal, vertical)
cell counts; fields are evaluated at cell centers and remain snapshots if
the sample changes. Thin structures require sufficient resolution. The
cross-section resolves replacement precedence, unlike the transparent
construction surfaces in the 3D view.

Boundary-aware diffusion
------------------------

Run a finite, layered-sample calculation with a surface spot:

.. code-block:: console

   python docs/examples/surface_diffusion.py

``CartesianMesh.from_sample(sample=sample, shape=(25, 25, 17))`` spans the finite
sample. Each Cartesian cell is split into six conforming tetrahedra; fluence
uses continuous piecewise-linear basis functions. The solver assembles
consistent volume and Robin surface mass matrices and solves the sparse system
with SciPy. Material coefficients are sampled at element centroids. Layer
interfaces should align with mesh planes; curved and unaligned boundaries are
approximated and require geometry and mesh-convergence checks. This is not an
interface-conforming automatic mesher.

Configure ``Diffusion(boundaries=...)`` with all six faces: ``x_min``,
``x_max``, ``y_min``, ``y_max``, ``z_min``, ``z_max``. Each requires either
``ReflectingBoundary()`` (zero normal current) or
``RobinBoundary(effective_reflectance=...)``. Effective reflectance is the
specified diffuse partial-current reflection coefficient, in [0, 1), not the
normal-incidence Fresnel reflectance. No exterior refractive index or Fresnel
approximation is silently assumed.

``SurfaceSource`` specifies ``face``, absolute ``position``, disk ``radius``,
``wavelength``, and ``optical_power``. ``UniformSurfaceSource`` illuminates a
whole face and requires ``face``, ``wavelength``, and ``optical_power``.
Source power is incoming diffuse partial-current power already transmitted
into the sample; incident LED power and specular entry losses are not modeled.
Illumination requires a Robin face. A disk must fit fully on that face.
Three-point triangle quadrature samples the top-hat spot and is normalized to
the requested power. Inspect ``result.diagnostics['source_quadrature_area']``
and refine the mesh to resolve the footprint; power normalization alone does
not establish spatial accuracy.

With outward current J.n, effective reflectance R, incoming partial current
I_in, and beta = (1-R)/(2(1+R)), the boundary convention is
``J.n = beta Phi - 2 I_in/(1+R)``. The weak equation therefore has Robin
term ``integral(beta v Phi)`` and illumination load
``integral(2 v I_in/(1+R))``. Escaping transmitted partial current is
``beta Phi - 2 beta I_in``. It is distinct from the gross Robin term and
from the net outward current. ``result.outgoing_power_by_face`` integrates
this escaping quantity; it is not a finite-aperture detector measurement.

``result.energy_balance`` contains incoming, absorbed, escaping, and residual
power. ``result.flux`` contains unsmoothed element-centered values of
``J = -D grad(Phi)`` in W/m^2, and ``result.flux.plot()`` displays relative-length
arrows. Fluence retains the standard 3D slice, volume, and isosurface views.
Numerical positivity is not enforced by clipping; minimum fluence and linear
relative residual are reported in ``result.diagnostics``. Log plots reject
nonpositive values.

The solver currently requires equal internal refractive indices. Its continuous
fluence treatment is not a model of refracting internal index jumps, which are
rejected explicitly. Optical coefficients apply at the source wavelength;
spectral models are not implemented. Diffusion is approximate near boundaries,
strong absorbers, and short source-detector separations. Energy conservation
and agreement with diffusion references do not establish agreement with RTE
or measured tissue data.

Boundary-condition and FEM reference:
https://pmc.ncbi.nlm.nih.gov/articles/PMC3195519/


Building the example gallery
----------------------------

Examples in ``docs/examples`` are executable Sphinx-Gallery scripts with an
RST introduction and explanatory ``# %%`` cells. They construct figures at
top level without opening browsers or writing export files. The documentation
scraper embeds each figure immediately after its corresponding code cell.

.. code-block:: console

   python -m pip install -e ".[documentation]"
   make docs PYTHON=.venv/bin/python

Open ``docs/_build/html/auto_examples/index.html`` to browse the interactive
examples. A saved figure can still be exported explicitly using its
``write_html`` method in your own scripts.

Citation and license
--------------------

MomentRT is distributed under the MIT license; see ``LICENSE``. If you use
MomentRT in research, please cite the software version using the metadata in
``CITATION.cff``. The ``.zenodo.json`` file provides metadata for a future
Zenodo archive. No DOI or publication date has been assigned yet.

GitHub workflows
----------------

The workflows follow PyMieSim's reusable MPSActions setup. Changes to ``master``
and pull requests run tests and quality checks; documentation and coverage are
also built for branch pushes. A wheel smoke-test matrix verifies installation
outside the checkout on Linux, macOS, and Windows.

A pushed ``v*`` tag triggers the release workflow, which validates version
metadata and publishes a pure Python distribution, a single noarch conda
package, and versioned documentation. The repository needs ``PYPI_API_TOKEN``
and ``ANACONDA_API_TOKEN`` secrets and GitHub Pages configured for publication.
No tag or publication is created by the local workflow files themselves.
