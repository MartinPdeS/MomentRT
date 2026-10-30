Choose a transport model
========================

MomentRT provides steady-state, monochromatic moment approximations. Select
the model from the physical domain and desired angular approximation, then
check the configuration and its assumptions before solving:

.. code-block:: python

   simulation = MomentRT.Simulation(...)
   report = simulation.validate()
   print(report.summary())

The report lists errors that prevent a solve, mesh and footprint cautions, and
``model_notes`` containing the configured solver's assumptions and known
limits. Validation does not solve the equation, establish mesh convergence,
or determine whether an approximation is accurate for a particular sample.

.. list-table:: Solver choice and principal limits
   :header-rows: 1
   :widths: 18 34 32 16

   * - Solver
     - Use it for
     - Key assumptions and limits
     - Example
   * - ``InfiniteDiffusion``
     - A homogeneous, unbounded medium with an embedded isotropic point source.
     - Analytical P1 Green function; no physical boundaries or finite source. Singular at the source; diffusion is least reliable nearby.
     - :doc:`Infinite medium <gallery/getting_started/infinite_diffusion>`
   * - ``Diffusion``
     - A finite sample illuminated by diffuse light on a face.
     - P1 diffusion with continuous fluence, centroid-sampled properties, and explicit six-face conditions. Internal refractive-index jumps are rejected.
     - :doc:`Finite sample <gallery/getting_started/surface_diffusion>`
   * - ``SP3``
     - A finite sample when the P1 closure is too restrictive and a higher angular moment model is useful.
     - Simplified P3; requires phase-function moments, zero-reflectance Robin conditions on illuminated faces, and equal internal indices. Unilluminated faces may be reflecting. It is not full P3 or the radiative transfer equation.
     - :doc:`P1 and SP3 <gallery/transport_models/compare_sp3>`
   * - ``P3``
     - A finite sample when retaining all three-dimensional harmonics through degree three is useful.
     - Full P3 moment approximation; requires matched-index zero-reflectance Robin faces. More angular detail does not guarantee accuracy for a given geometry.
     - :doc:`P1, SP3, and P3 <gallery/transport_models/compare_p3>`

P1 diffusion is generally appropriate when multiple scattering supports a
diffusive description and the observable is not dominated by the source
boundary layer. Higher-order moment models retain more angular information,
but remain approximations and have their own boundary and phase-function
limits. For optically thin, weakly scattering, strongly absorbing, or very
short source-detector configurations, compare against a suitable radiative
transfer or Monte Carlo reference rather than assuming a higher order is
accurate.

Model notes are available programmatically as ``report.model_notes``. The
complete report also identifies unsupported internal refractive-index jumps
with region names and resolved indices. Fixing configuration errors does not
remove numerical or physical model error.

Independent references
----------------------

The :doc:`planar slab example
<gallery/transport_models/planar_slab_references>` compares P1 and SP3 finite-element
fields with independently calculated one-dimensional analytical solutions.
This checks the FEM implementation and spatial discretization for a planar
homogeneous case. It does not validate the approximations against the full
radiative transfer equation.

For a broader view of model breakdown in scattering slabs, see Akdemir et al.,
`Breakdown of light transport models in photonic scattering slabs with strong
absorption and anisotropy <https://doi.org/10.1103/PhysRevA.105.033517>`_. Its
Monte Carlo comparisons show why P1 and P3 applicability depends on optical
thickness, absorption, anisotropy, and refractive-index contrast. Those results
are not a direct validation of MomentRT's finite-element implementation.
