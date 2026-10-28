Finite-sample P1 solver
=======================

``Diffusion`` solves for steady fluence rate :math:`\Phi` in a finite
rectangular sample. The geometry is represented by a structured Cartesian
mesh split into linear tetrahedra. All six outer faces receive explicit
boundary conditions; the solver does not infer a boundary model from a
material's refractive index.

Continuous equation and face conditions
----------------------------------------

Inside a source-free sample the P1 equation is

.. math::

   -\nabla\cdot(D\nabla\Phi)+\mu_a\Phi=0,
   \qquad \mathbf J=-D\nabla\Phi.

For a Robin face, the user supplies an effective diffuse reflectance
:math:`R\in[0,1)`. MomentRT defines

.. math::

   \beta=\frac{1-R}{2(1+R)},
   \qquad
   c_R=\frac{2}{1+R},
   \qquad
   \mathbf J\cdot\mathbf n=\beta\Phi-c_R I_{\mathrm{in}},

where :math:`\mathbf n` points outward and :math:`I_{\mathrm{in}}` is incoming
diffuse partial-current irradiance already coupled into the sample. This
effective coefficient is not normal-incidence Fresnel reflectance. The solver
does not choose an index-mismatch formula. A reflecting face uses
:math:`\mathbf J\cdot\mathbf n=0` and cannot be illuminated.

To derive the finite-element equation, multiply the diffusion equation by a
test function :math:`v`, integrate over the sample :math:`\Omega`, and
integrate the diffusion term by parts. Substituting the Robin condition yields

.. math::

   \int_\Omega D\nabla v\cdot\nabla\Phi\,dV
   +\int_\Omega \mu_a v\Phi\,dV
   +\int_{\partial\Omega_R}\beta v\Phi\,dA
   =\int_{\partial\Omega_R}c_R v I_{\mathrm{in}}\,dA.

Reflecting faces contribute no boundary matrix or load. The illuminated face
must be Robin. For a uniform full-face source,
:math:`I_{\mathrm{in}}=P/A_{\mathrm{face}}`. For a disk source it is a
top-hat over the sampled disk area, normalized so its integrated power equals
the requested source power.

.. figure:: _static/images/finite-fluence-slices.png
   :align: center
   :alt: Three intersecting central planes showing logarithmic FEM fluence in a finite sample.
   :width: 85%

   Finite-sample result from the P1 FEM solve. The horizontal plane lies at the
   middle depth; the vertical planes pass through the center of the sample.
   These are sampled mesh planes, not smoothed reconstructions.

Finite-element assembly
-----------------------

``CartesianMesh.from_sample(sample=sample, shape=(n_x, n_y, n_z))`` places
nodes on the sample faces. Each Cartesian cell is split into six tetrahedra
that share a consistent body diagonal. Therefore the total element count is

.. math::

   N_{\mathrm{tet}}=6(n_x-1)(n_y-1)(n_z-1).

The unknowns are nodal fluence values. On each tetrahedron, the basis functions
:math:`N_i` are linear and their gradients are constant. Absorption and
reduced-scattering properties are sampled at each element centroid, giving
constant element coefficients

.. math::

   D_e=\frac{1}{3(\mu_{a,e}+\mu'_{s,e})}.

The element stiffness and absorption matrices are

.. math::

   K^{(e)}_{ij}=D_e V_e\,\nabla N_i\cdot\nabla N_j,
   \qquad
   M^{(e)}_{ij}=\mu_{a,e}\frac{V_e}{20}(1+\delta_{ij}).

On a boundary triangle with area :math:`A_f` and constant :math:`\beta`, the
Robin matrix is

.. math::

   B^{(f)}_{ij}=\beta\frac{A_f}{12}(1+\delta_{ij}).

These local matrices are assembled into a sparse global system

.. math::

   \mathbf A\boldsymbol\Phi=\mathbf b.

The source load is integrated with a three-point triangle quadrature rule.
For a circular spot, quadrature points outside the disk receive zero
irradiance. MomentRT reports the total sampled area. Power normalization
ensures the requested integrated power, but it cannot recover a poorly
resolved spot edge; refine the mesh and compare results to assess that error.

Material map and interface approximation
-----------------------------------------

``StructuredMedium`` first assigns the background, then applies structures in
insertion order. Later structures replace earlier ones wherever they overlap.
This precedence is used when the solver labels element centroids and samples
their material. It does not create an interface-conforming mesh. A layer
interface aligned to mesh planes is represented cleanly; an unaligned planar
interface or curved sphere is represented only through the centroids of nearby
elements. Refine the grid and align important planar boundaries where
possible.

The formulation uses one continuous fluence field. MomentRT therefore rejects
internal refractive-index jumps rather than implying a refracting interface
model. The optical-property map may still vary in absorption and reduced
scattering between structures, subject to the geometry and mesh limitations
above.

Power accounting and validation
-------------------------------

The solution is differentiated within each tetrahedron to form element
centered current

.. math::

   \mathbf J_e=-D_e\sum_{i=1}^{4}\Phi_i\nabla N_i.

The returned ``result.flux`` stores these vectors without smoothing across
interfaces. Absorbed power is integrated from the element absorption and mean
nodal fluence. Escaping transmitted partial current is computed per face from
the Robin solution, with the incoming contribution on the illuminated face
removed. It differs from both gross Robin leakage and net outward current.
It is an integral across an entire face, not a detector measurement over an
aperture.

The reported balance is

.. math::

   P_{\mathrm{in}}=P_{\mathrm{abs}}+P_{\mathrm{escape}}+P_{\mathrm{residual}}.

Inspect ``result.diagnostics`` for the linear relative residual, minimum
fluence, tetrahedron count, and sampled source area. Inspect
``result.energy_balance`` for unit-bearing powers. A small residual and good
power balance show that the discrete linear system and bookkeeping are
consistent. They do not show that P1 approximates RTE well, that a source disk
is spatially resolved, or that a mesh is converged. For physical confidence,
compare independent references and repeat with refined and interface-aligned
meshes.

Limitations to keep visible
---------------------------

The finite solver is steady-state and monochromatic. It has no internal
Fresnel/refraction law, collimated illumination, detector aperture, time
response, spectral integration, or inverse model. P1 diffusion is
least reliable near the illumination footprint, near boundaries, in
low-scattering areas, around unresolved interfaces, and at short source-to-
detector distances. These are limitations of the physical approximation and
discretization; changing the plot scale cannot remove them.
