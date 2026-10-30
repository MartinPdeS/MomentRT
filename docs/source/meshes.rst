Meshes and local refinement
===========================

``CartesianMesh`` provides a uniform grid with an explicit physical resolution.
``TetrahedralMesh`` provides unstructured linear tetrahedra with geometric
refinement rules. P1, SP3, and full P3 accept either finite-sample mesh.
``InfiniteDiffusion`` continues to use a Cartesian observation window.

Uniform grids
-------------

Request a maximum spacing in length units. A scalar applies to all axes; a
three-component quantity sets x/y/z spacing independently. Actual spacing is
adjusted downward so the mesh spans the sample exactly.

.. code-block:: python

   mesh = MomentRT.CartesianMesh.from_sample(
       sample=sample,
       resolution=(50, 50, 200) * MomentRT.ureg.um,
   )

``mesh.resolution`` preserves the request, ``mesh.spacing`` reports actual
spacing, and ``mesh.shape`` gives the node counts. Existing ``shape=...``
configurations remain supported; supply exactly one of shape or resolution.

Tetrahedral generation
----------------------

Gmsh is installed as a core MomentRT dependency. Its native runtime libraries
must also be available; on Ubuntu these can require ``libglu1-mesa``. Numerical
calculations with an existing tetrahedral mesh still require the Gmsh package
to import MomentRT, but do not invoke the mesher. Plotly and PyVista remain
optional.

The sample, its layers, boxes, and spheres are meshed together so internal
interfaces share nodes. Curved interfaces are represented by planar triangles.
Material assignment retains the sample's rule that later structures replace
earlier ones. All current transport restrictions, including refractive-index
requirements, still apply.

.. code-block:: python

   mesh = MomentRT.TetrahedralMesh.from_sample(
       sample=sample,
       resolution=0.5 * MomentRT.ureg.mm,
   )

   mesh.refine_region(
       name="Epidermis",
       resolution=50 * MomentRT.ureg.um,
   )

   mesh.refine_near(
       objects=[*sources, detector],
       resolution=25 * MomentRT.ureg.um,
       distance=300 * MomentRT.ureg.um,
   )

For tetrahedra, resolution is a **scalar target element size**, not a strict
maximum edge length. A region rule refines that shape's geometric support,
including portions overridden by later regions. A contact rule measures distance
from the complete circular disk, including its radius and interior. The fine
target applies within ``distance``; the size field transitions toward the global
target over another ``distance``. Source groups are accepted and expanded into
their member disks. Contacts must fit on their configured sample faces.

Each refinement call regenerates the mesh and applies all accumulated rules.
Only a successful generation replaces the previous mesh. The sample and contact
geometry are snapshotted: changing their original objects does not modify the
mesh. Regenerate after changing sample geometry; stale fitted interfaces are
rejected by preflight and finite solvers. Existing simulation results keep their
own independent mesh snapshots.

The simulation interface remains the same:

.. code-block:: python

   simulation = MomentRT.Simulation(
       medium=sample,
       source=sources,
       mesh=mesh,
       solver=solver,
       detectors={"collection": detector},
   )

   report = simulation.validate()
   result = simulation.run()

These rules are geometric refinement, not an automatic error estimator.
Compare fluence and detector power across refinements and inspect sampled
source and detector areas. Fitting a contact outline or satisfying energy
balance does not establish transport accuracy or mesh convergence.

Inspecting and exporting results
--------------------------------

Tetrahedral ``mesh.points`` has shape ``(N, 3)``, connectivity
``mesh.tetrahedra`` has shape ``(T, 4)``, and ``mesh.faces`` identifies the six
sample boundaries. ``mesh.shape`` is the nodal storage shape ``(N,)``.
Fluence and even moments have one value per node; flux remains element-centered
and is not smoothed across interfaces.

Plot the sample's translucent boundary and the tetrahedral edges to inspect
the mesh and its local refinement together. The returned Plotly figure is
composable, and ``max_edges`` limits large meshes to a representative subset
of edges. Install Plotly with the ``viz`` extra if it is not already available:

.. code-block:: python

   mesh_figure = mesh.plot(
       length_unit="mm",
       max_edges=20000,
       show=False,
   )

Volume plots explicitly request a display grid:

.. code-block:: python

   figure = result.fluence.plot(
       display_resolution=100 * MomentRT.ureg.um,
       overlays=[mesh, *sources, detector],
       show=False,
   )

   display_field = result.fluence.resample(
       resolution=100 * MomentRT.ureg.um,
   )

Interpolation uses the linear basis functions and connectivity of the original
tetrahedra. The returned display field records its interpolation provenance and
does not alter the computed FEM field. A finer display grid cannot improve FEM
accuracy. Mesh overlays draw a representative subset of tetrahedral or
Cartesian edges over the field. Log plots retain the requirement for strictly
positive values.

With the ``vtk`` extra, ``result.fluence.export(path="fluence.vtu")`` exports
the original tetrahedra and nodal values. Cartesian fields continue to use
``.vti``. External tetrahedral connectivity can be supplied directly through
``TetrahedralMesh(points=..., tetrahedra=...)``; it must fill the finite
rectangular bounds. Imported meshes use centroid-sampled coefficients and
cannot use the generation refinement methods without an original sample recipe.
