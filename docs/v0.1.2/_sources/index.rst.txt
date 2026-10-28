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

   .. grid-item-card:: Understand the model
      :link: theory
      :link-type: doc

      Read the P1 equations, approximations, and limits.

   .. grid-item-card:: Finite element details
      :link: finite_element
      :link-type: doc

      Follow the weak form, tetrahedral discretization, and boundaries.

   .. grid-item-card:: API reference
      :link: api_reference
      :link-type: doc

      Find the public objects and their responsibilities.

Before interpreting a field
---------------------------

MomentRT solves the diffusion approximation, not the radiative transfer
equation. Its fields are model predictions; they are not measurements or
automatically validated tissue properties. Read the assumptions in
:doc:`theory` and :doc:`finite_element` before interpreting fine structure,
near-source values, or boundary behavior.

.. toctree::
   :maxdepth: 2
   :hidden:

   Start <getting_started>
   Guide <guide>
   Examples <gallery/index>
   Theory <theory>
   FEM <finite_element>
   Results <results>
   API <api_reference>
   Roadmap <roadmap>
