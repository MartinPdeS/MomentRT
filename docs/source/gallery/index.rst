Examples
========

The examples use MomentRT to configure physical inputs, run supported P1
calculations, and inspect the output. Input values are illustrative; the
examples are executable documentation, not calibrated tissue measurements.

Fundamentals
------------

Start with the analytical infinite-medium solution. It is the smallest
complete solver workflow and introduces units, sources, meshes, results, and
interactive 3D fluence plotting.

.. toctree::
   :maxdepth: 1

   fundamentals/index

Finite samples
--------------

Build and inspect a structured sample before solving light transport through
it. The surface illumination example then shows the finite element workflow,
including boundary conditions and energy accounting.

.. toctree::
   :maxdepth: 1

   finite_samples/index

.. only:: html

   The Plotly figures embedded in these pages are interactive: rotate, zoom,
   and hover to inspect the generated MomentRT data.
