Theory and finite element method
================================

MomentRT's two solvers share the P1 diffusion approximation but apply it to
different domains. Start with the transport model, then follow its
finite-sample weak form and numerical discretization.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: P1 transport model
      :link: ../theory
      :link-type: doc

      Angular moments, diffusion closure, infinite-medium Green function, and
      approximation limits.

   .. grid-item-card:: Finite-sample FEM
      :link: ../finite_element
      :link-type: doc

      Robin boundaries, weak form, tetrahedral assembly, and source loading.

.. toctree::
   :maxdepth: 1
   :hidden:

   P1 transport model <../theory>
   Finite-sample FEM <../finite_element>
