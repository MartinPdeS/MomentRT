MomentRT documentation
======================

Deterministic light transport with angular-moment diffusion models.

MomentRT provides an analytical P1 solution for an infinite homogeneous
medium and a finite-sample P1 solver using linear tetrahedral finite elements.
The guides explain the equations, boundary model, assumptions, and returned
fields so that numerical output can be interpreted in context.

**Version:** |version|

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Get started
      :link: getting_started
      :link-type: doc

      Set up a sample, source, mesh, and solver.

   .. grid-item-card:: Choose a workflow
      :link: guide
      :link-type: doc

      Select the solver that matches your geometry.

   .. grid-item-card:: Explore examples
      :link: gallery/index
      :link-type: doc

      Run complete examples with interactive visualizations.

   .. grid-item-card:: Theory and finite elements
      :link: theory/index
      :link-type: doc

      Read the P1 model, assumptions, and FEM formulation.

   .. grid-item-card:: Results and API
      :link: reference/index
      :link-type: doc

      Interpret outputs and find the public objects.

Before interpreting a field
---------------------------

MomentRT solves steady-state P1 diffusion, SP3, and full P3 moment models,
not the full radiative transfer equation. Their fields are model predictions;
they are not measurements or automatically validated tissue properties. Read
the assumptions in :doc:`theory`, :doc:`sp3`, :doc:`p3`, and
:doc:`finite_element` before interpreting fine structure, near-source values,
or boundary behavior.

.. toctree::
   :maxdepth: 2
   :hidden:

   Start <getting_started>
   Guide <guide>
   Choose a transport model <model_selection>
   Examples <gallery/index>
   Theory <theory/index>
   Reference <reference/index>
