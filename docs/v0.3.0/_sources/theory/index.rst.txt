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

      Radiance and units, angular moments, the P1 closure, conservation,
      infinite-medium Green function, and approximation limits.

   .. grid-item-card:: Finite-sample FEM
      :link: ../finite_element
      :link-type: doc

      Robin boundaries, weak form, tetrahedral assembly, and source loading.

   .. grid-item-card:: P1 and the telegrapher equation
      :link: ../time_dependent
      :link-type: doc

      Current relaxation, transient diffusion, finite propagation speed,
      initial conditions, and the steady-state limit.

.. toctree::
   :maxdepth: 1
   :hidden:

   P1 transport model <../theory>
   Finite-sample FEM <../finite_element>
   Time-dependent P1 and telegrapher equation <../time_dependent>
