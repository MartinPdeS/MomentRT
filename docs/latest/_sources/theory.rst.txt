Transport theory and approximations
===================================

MomentRT implements steady-state P1 diffusion, simplified P3 (SP3), and
full three-dimensional P3. These are distinct angular-moment approximations
to the radiative transfer equation (RTE); MomentRT does not solve the full
RTE. P1 is available as an analytical Green function in an infinite
homogeneous medium and as a finite-domain FEM solver. SP3 and P3 use finite
sample FEM solvers, with their own boundary and material restrictions.

This page derives the P1 diffusion model. For the higher-order models and
their limits, see :doc:`sp3` and :doc:`p3`.

How to read this theory
-----------------------

The derivation follows a sequence of modeling choices:

#. Describe incoherent radiant power by a directional radiance.
#. Integrate over directions to obtain balance equations for angular moments.
#. Close those equations with the first-order spherical-harmonic expansion.
#. In steady state, eliminate current to obtain a scalar diffusion equation.
#. Select boundary conditions and discretize the spatial equation.

The angular closure determines which directional information the model
retains. Boundary conditions determine how light enters and leaves a finite
sample. The finite element method determines how the chosen continuous
equation is approximated on a mesh. Refining the mesh improves the spatial
approximation; it does not restore angular information discarded by P1.

For the boundary and spatial discretization, continue with
:doc:`finite_element`. For time dependence and the telegrapher equation,
continue with :doc:`time_dependent`. The time-dependent equations are
background theory, not implemented solver options.

Quantities and units
---------------------

The discussion is monochromatic: optical coefficients describe one specified
wavelength. MomentRT does not calculate a material spectrum from a wavelength.
The caller supplies appropriate coefficients for the chosen wavelength.

.. list-table:: Symbols used in the derivation
   :header-rows: 1
   :widths: 18 52 30

   * - Symbol
     - Meaning
     - Units
   * - :math:`\mathbf r`, :math:`\mathbf s`
     - Position and unit travel direction
     - :math:`\mathrm m`, dimensionless
   * - :math:`L`
     - Radiance per unit solid angle
     - :math:`\mathrm{W\,m^{-2}\,sr^{-1}}`
   * - :math:`\Phi`
     - Fluence rate: radiance integrated over all directions
     - :math:`\mathrm{W\,m^{-2}}`
   * - :math:`\mathbf J`
     - Net radiant current per unit area
     - :math:`\mathrm{W\,m^{-2}}`
   * - :math:`\mu_a`, :math:`\mu_s`, :math:`\mu_s'`
     - Absorption, scattering, and reduced scattering coefficients
     - :math:`\mathrm{m^{-1}}`
   * - :math:`g`
     - Mean cosine of the scattering angle
     - Dimensionless
   * - :math:`D`
     - Diffusion coefficient in MomentRT's steady transport convention
     - :math:`\mathrm m`
   * - :math:`Q`, :math:`q_0`
     - Directional volume source and its angular integral
     - :math:`\mathrm{W\,m^{-3}\,sr^{-1}}`, :math:`\mathrm{W\,m^{-3}}`
   * - :math:`I_{\mathrm{in}}`
     - Prescribed incoming diffuse partial-current irradiance at a face
     - :math:`\mathrm{W\,m^{-2}}`

Fluence rate is not the flux through a particular plane. For perfectly
isotropic radiance :math:`L=L_0`, one has :math:`\Phi=4\pi L_0` but
:math:`\mathbf J=\mathbf 0`: light crosses a plane equally in both
directions. A high fluence therefore need not imply a large net current.
Absorbed power density is :math:`\mu_a\Phi`, and integrating it over a
volume gives absorbed power in watts.

Radiative transfer and angular moments
---------------------------------------

Let :math:`L(\mathbf r,\mathbf s)` denote radiance at position
:math:`\mathbf r` in direction :math:`\mathbf s`. A steady, monochromatic
transport equation with absorption, scattering, and source :math:`Q` can be
written

.. math::

   \mathbf s\cdot\nabla L + (\mu_a+\mu_s)L
   = \mu_s \int_{4\pi} p(\mathbf s,\mathbf s')
     L(\mathbf r,\mathbf s')\,d\Omega' + Q(\mathbf r,\mathbf s).

Here :math:`\mu_a` is the absorption coefficient, :math:`\mu_s` the scattering
coefficient, and :math:`p` a normalized phase function. The angular moments
describe how radiance is distributed in direction at each point. The zeroth
moment is the fluence rate and the first moment is the vector current:

.. math::

   \Phi(\mathbf r)=\int_{4\pi}L(\mathbf r,\mathbf s)\,d\Omega,
   \qquad
   \mathbf J(\mathbf r)=\int_{4\pi}\mathbf sL(\mathbf r,\mathbf s)\,d\Omega.

The streaming term :math:`\mathbf s\cdot\nabla L` describes propagation along
a ray. Absorption removes radiant power. Scattering removes radiance from
one direction and adds it to other directions through the integral term.
The source :math:`Q` adds radiance inside the volume.

This transport model describes incoherent intensity. It does not evolve the
electromagnetic field phase, interference, or polarization. A telegrapher
equation derived from these moments remains an intensity-transport
approximation; its wave-like mathematical form does not turn it into
Maxwell's equations.

The phase function is normalized over outgoing directions, so scattering
redistributes radiance without creating or destroying radiant power. Its
first angular moment is summarized by the anisotropy :math:`g`, the mean
cosine of the scattering angle. Positive :math:`g` describes predominantly
forward scattering; :math:`g=0` describes isotropic scattering. In the RTE
above, the phase normalization is
:math:`\int_{4\pi}p(\mathbf s,\mathbf s')\,d\Omega=1` for each incoming
direction :math:`\mathbf s'`. For a rotationally symmetric scattering
medium, its first directional moment is
:math:`\int_{4\pi}\mathbf s\,p(\mathbf s,\mathbf s')\,d\Omega
=g\mathbf s'`. This identity is why the first moment of the scattering gain
is :math:`g\mathbf J`.

The transport equation couples these moments to higher angular moments.
:math:`\Phi` is the radiance integrated over all directions, while
:math:`\mathbf J` weights each direction by its unit vector to record the
local net flow. With the convention used here, both quantities have units of
:math:`\mathrm{W\,m^{-2}}`.

Taking angular moments produces a hierarchy: the balance of :math:`\Phi`
depends on :math:`\mathbf J`, and the balance of :math:`\mathbf J` depends on
the second moment tensor

.. math::

   \mathbf P(\mathbf r)=\int_{4\pi}\mathbf s\mathbf s
   L(\mathbf r,\mathbf s)\,d\Omega.

The dyadic :math:`\mathbf s\mathbf s` records how strongly radiance points
along each pair of axes. Equations for :math:`\mathbf P` involve still higher
moments, so the hierarchy does not close by itself. A finite model needs a
rule that estimates the next moment from the moments it retains.

How P1 closes the moment equations
----------------------------------

P1 expands angular radiance in spherical harmonics and retains the isotropic
term and the three first-order directional terms. Written in terms of the
moments above, this approximation is

.. math::

   L(\mathbf r,\mathbf s)\approx
   \frac{1}{4\pi}\left[\Phi(\mathbf r)+3\,\mathbf s\cdot\mathbf J(\mathbf r)\right].

In a full expansion, the angular basis is :math:`Y_{\ell m}(\mathbf s)`,
with degree :math:`\ell` and azimuthal index :math:`m`. P1 keeps
:math:`\ell=0,1`, so there are :math:`1+3=4` real angular coefficients in
three dimensions. These can be represented by one scalar :math:`\Phi`
and the three components of :math:`\mathbf J`. The name P1 refers to angular
order; it does not mean that only one unknown is retained.

The factor :math:`1/(4\pi)` normalizes the isotropic part: integrating that
term over the sphere returns :math:`\Phi`. The factor 3 normalizes the first
order term, since
:math:`\int_{4\pi}\mathbf s\mathbf s\,d\Omega=(4\pi/3)\mathbf I`;
integrating :math:`\mathbf sL` therefore returns :math:`\mathbf J`.

Substitute this approximation into the second moment. The isotropic term
contributes :math:`\Phi\mathbf I/3`; the first-order term integrates to zero
because every odd power of direction cancels between opposite directions.
P1 therefore uses the closure

.. math::

   \mathbf P\approx\frac{\Phi}{3}\mathbf I.

This is the precise step that removes the unknown second moment and closes the
first two equations. It says that, to leading order, radiation is spread
equally over directions. The first-order correction allows a modest net flow,
but P1 cannot preserve a narrow beam, a sharply directional source, or strong
angular anisotropy.

The closure can also be checked component by component. The necessary
angular integrals are

.. math::

   \int_{4\pi}s_i\,d\Omega=0,\qquad
   \int_{4\pi}s_i s_j\,d\Omega=\frac{4\pi}{3}\delta_{ij},\qquad
   \int_{4\pi}s_i s_j s_k\,d\Omega=0.

The first and third integrals vanish by symmetry under
:math:`\mathbf s\mapsto-\mathbf s`. The second is isotropic and has trace
:math:`\int|\mathbf s|^2d\Omega=4\pi`, which fixes the factor :math:`1/3`.
Consequently :math:`P_{ij}=\Phi\delta_{ij}/3`.

P1 does not automatically preserve nonnegative radiance. Its reconstructed
minimum over directions is

.. math::

   L_{\min}=\frac{\Phi-3|\mathbf J|}{4\pi}.

A nonnegative P1 angular reconstruction requires
:math:`|\mathbf J|\le\Phi/3`. General nonnegative radiance only requires
:math:`|\mathbf J|\le\Phi`, so a physically possible strongly directional
field may already be outside the useful range of the linear P1 ansatz.
MomentRT solves for fluence and current rather than reconstructing a
guaranteed nonnegative angular distribution.

Diffusion equation
------------------

To derive the first two moment equations, define the zeroth and first moments
of the source by

.. math::

   q_0=\int_{4\pi}Q\,d\Omega,
   \qquad
   \mathbf q_1=\int_{4\pi}\mathbf sQ\,d\Omega.

Integrating the transport equation over all angles gives the zeroth balance.
Scattering gain and loss cancel because scattering conserves power. For the
first balance, multiply the transport equation by :math:`\mathbf s` before
integrating. A fraction :math:`g` of the directed scattering current remains
in the same direction on average, leaving the net current loss
:math:`\mu_s(1-g)\mathbf J=\mu_s'\mathbf J`. Before making the P1 closure,
these balances are

.. math::

   \nabla\cdot\mathbf J+\mu_a\Phi=q_0,
   \qquad
   \nabla\cdot\mathbf P+(\mu_a+\mu_s')\mathbf J=\mathbf q_1.

The cancellation in the zeroth balance follows explicitly from phase
normalization:

.. math::

   \mu_s\int_{4\pi}\int_{4\pi}
   p(\mathbf s,\mathbf s')L(\mathbf s')\,d\Omega'\,d\Omega
   =\mu_s\Phi.

That gain cancels the scattering loss :math:`\mu_s\Phi`. The first moment
instead retains directional memory :math:`g`, so its gain is
:math:`\mu_s g\mathbf J`; subtracting it from
:math:`(\mu_a+\mu_s)\mathbf J` leaves
:math:`[\mu_a+\mu_s(1-g)]\mathbf J`.

Now substitute :math:`\mathbf P\approx\Phi\mathbf I/3` from the P1
approximation above. An isotropic volume source, including the point source
used by ``InfiniteDiffusion``, has no preferred direction and therefore
:math:`\mathbf q_1=\mathbf 0`. The resulting first two P1 equations are

.. math::

   \nabla\cdot\mathbf J + \mu_a\Phi = q_0,

   \frac{1}{3}\nabla\Phi + (\mu_a+\mu_s')\mathbf J = 0,

Here :math:`q_0` is source power per volume and
:math:`\mu_s'=\mu_s(1-g)` is reduced scattering. Absorption and effective
scattering remove directed current, which is why they appear together in the
current balance.

The second equation can now be solved algebraically for :math:`\mathbf J`:

.. math::

   \mathbf J=-D\nabla\Phi,
   \qquad
   D=\frac{1}{3(\mu_a+\mu_s')},

.. math::

   -\nabla\cdot(D\nabla\Phi)+\mu_a\Phi=q_0.

This is the steady P1 diffusion equation used by MomentRT. It is a scalar
conservation equation for fluence rate. In a source-free region,
:math:`q_0=0`; an isotropic volume source contributes its power density to
:math:`q_0`. MomentRT's finite solver has no volume source term: it supplies
illumination on a sample face, which enters through the boundary condition
derived in :doc:`finite_element` rather than through :math:`q_0` or
:math:`\mathbf q_1`.

MomentRT asks for :math:`\mu_a` and :math:`\mu_s'` directly, both in inverse
length units. It does not infer :math:`\mu_s'` from :math:`\mu_s` and :math:`g`.
In this convention :math:`D` has units of length, :math:`\Phi` and
:math:`\mathbf J` have units of :math:`\mathrm{W\,m^{-2}}`, and integrated
source and power terms have units of watts. With spatially varying material
properties, :math:`D` varies in space and remains inside the divergence
operator; replacing :math:`\nabla\cdot(D\nabla\Phi)` by
:math:`D\nabla^2\Phi` would incorrectly drop gradients of :math:`D`.

The transport mean free path
:math:`\ell_{\mathrm{tr}}=1/(\mu_a+\mu_s')` is the distance over which the
first angular moment loses its direction. The diffusion approximation is
most useful when fluence changes over distances long compared with this
scale, so that local angular distributions are close to the P1 form.

For example, illustrative coefficients
:math:`\mu_a=0.02\,\mathrm{mm^{-1}}` and
:math:`\mu_s'=1.0\,\mathrm{mm^{-1}}` give
:math:`D\approx0.327\,\mathrm{mm}` and
:math:`\ell_{\mathrm{tr}}\approx0.980\,\mathrm{mm}`. These numbers describe
the chosen model. They do not establish that P1 is accurate at a particular
distance from a source or interface.

Conservation over a sample
--------------------------

Integrating the zeroth moment over a volume and using the divergence theorem
gives

.. math::

   \int_{\partial\Omega}\mathbf J\cdot\mathbf n\,dA
   +\int_\Omega\mu_a\Phi\,dV
   =\int_\Omega q_0\,dV.

For the finite solver, :math:`q_0=0` inside the sample. Incoming illumination
makes the net outward current negative on the illuminated part of the
boundary; absorption and escape balance that incoming power. A source-free
interior therefore does not imply a dark sample.

The local equation and the global power balance are two forms of the same
conservation law. This is why a numerical solution can conserve total power
very accurately while retaining a mesh error or a transport-model error.

Infinite homogeneous point source
---------------------------------

For constant optical coefficients in unbounded three-dimensional space, an
isotropic point source with power :math:`P` gives

.. math::

   -D\nabla^2\Phi+\mu_a\Phi=P\,\delta(\mathbf r-\mathbf r_s).

Away from the point source, the radial field satisfies
:math:`D\nabla^2\Phi=\mu_a\Phi`. For a radial function,
:math:`\nabla^2\Phi=r^{-2}\partial_r(r^2\partial_r\Phi)`. Setting
:math:`u=r\Phi` reduces the equation to
:math:`u''=\kappa^2u`, where :math:`\kappa=\sqrt{\mu_a/D}`. Only the decaying
solution :math:`u\propto e^{-\kappa r}` is finite at infinity. The point
source fixes the amplitude by conserving power through a small sphere around
it. If :math:`\Phi=Ae^{-\kappa r}/r`, the outward radial current is
:math:`J_r=-D\partial_r\Phi`. Its integral over a sphere of radius :math:`r`
is

.. math::

   4\pi r^2J_r=4\pi DA(1+\kappa r)e^{-\kappa r}
   \longrightarrow 4\pi DA\quad (r\to0).

Setting that limit equal to the source power :math:`P` gives
:math:`A=P/(4\pi D)`, and hence the Green function

.. math::

   \Phi(r)=\frac{P}{4\pi D r}\exp(-\kappa r),
   \qquad
   \kappa=\sqrt{\frac{\mu_a}{D}},
   \qquad
   r=\lVert\mathbf r-\mathbf r_s\rVert.

``InfiniteDiffusion`` evaluates this expression at every node in a
``CartesianMesh``. The mesh is only a window into the infinite solution. There
is no finite-sample boundary condition, interface, or external reflectance in
this result. A node exactly at :math:`r=0` is rejected because the ideal point
source is singular; the model does not silently soften or spread the source.

The exponential sets the attenuation length :math:`1/\kappa`; the factor
:math:`1/r` is the geometric spreading of power over a sphere of area
:math:`4\pi r^2`. Together they show the model's two length scales:
directional randomization through :math:`\ell_{\mathrm{tr}}` and absorption
through :math:`1/\kappa`. The :doc:`absorption sweep example
<gallery/sweeps_convergence/absorption_sweep>` compares the resulting radial profiles
while holding reduced scattering and source power fixed.

.. figure:: _static/images/infinite-fluence-volume.png
   :align: center
   :alt: MomentRT logarithmic volume rendering of fluence from an embedded point source.
   :width: 85%

   An infinite-medium P1 field evaluated on a finite observation grid. The
   observation-window box is not a physical boundary, and the near-source
   values remain subject to the diffusion approximation.

Diffusion is an approximation
-----------------------------

P1 is most plausible after scattering has made radiance nearly diffuse and
when the field varies over distances larger than the transport mean free path

.. math::

   \ell_{\mathrm{tr}}=\frac{1}{\mu_a+\mu_s'}.

It becomes less reliable near a point source, near an illuminated boundary,
in weakly scattering or optically thin regions, around sharp optical
transitions, and when source-detector separation is short. The diffusion law
replaces angular detail by a local gradient relation; a smooth-looking field
does not establish that the physical field is accurate.

MomentRT provides passive contact :doc:`detectors` with an additional
Lambertian escape assumption and independent monochromatic :doc:`spectral`
sweeps. It does not provide collimated-source propagation, time dependence,
broadband spectral integration,
chromophore fitting, or inverse estimation. Refractive index is stored with
materials, but no automatic Fresnel rule is derived from it. In the finite
solver all internal refractive indices must be equal; an external diffuse
boundary coefficient must be supplied explicitly.

What the order labels mean
--------------------------

:math:`P_1` is the first-order spherical-harmonic angular expansion used to
derive diffusion. :math:`P_N` denotes a general finite spherical-harmonic
expansion through order :math:`N`; it contains multiple coupled angular
moments. :math:`SP_N` is a simplified transport approximation, often arranged
as coupled diffusion-like equations for scalar combinations of moments.
Asymptotic derivations treat it as a correction to diffusion under specified
transport scalings. In multiple dimensions it does not retain the complete
angular information of full :math:`P_N`.

In three dimensions, a full P3 expansion retains every angular degree
:math:`\ell=0,1,2,3`, giving :math:`(3+1)^2=16` angular coefficients.
It does not run separate P1, P2, and P3 solves. P2 also exists: even angular
orders are legitimate truncations. The common sequence P1, P3, P5 reflects
practical parity and boundary-condition choices in many formulations; it is
not a rule forbidding P2.

For the same reason, SP3 and P3 are not interchangeable names. A common
steady SP3 formulation has two coupled scalar equations after eliminating
auxiliary moments, but it uses its own closure and boundary conditions.
Replacing a diffusion coefficient or increasing the number of iso-surfaces
in a plot cannot turn the present solver into SP3.

For the asymptotic basis of simplified moment systems, see the original
`SPN derivation <https://www.osti.gov/servlets/purl/6808825>`_ and
`the work on time-dependent SPN equations
<https://doi.org/10.1016/j.jcp.2007.07.009>`_.
Neither label means that a P1 solve secretly contains P2 or P3.
MomentRT offers the P1 diffusion equation above, the separate :doc:`sp3`
solver, and the full three-dimensional :doc:`p3` moment system.

Angular order and finite element order
--------------------------------------

The transport closure and spatial approximation can both carry a P1 label:

* Transport P1 truncates the angular spherical-harmonic expansion.
* Linear tetrahedral finite elements use polynomials of degree one in
  position inside each tetrahedron.

MomentRT currently uses both. Raising the spatial polynomial degree would
improve how the diffusion field is represented within an element, while
raising the transport angular order would change the physical equations and
their unknowns. They are independent changes.

From theory to the public API
-----------------------------

``Material.absorption`` supplies :math:`\mu_a`, and
``Material.reduced_scattering`` supplies :math:`\mu_s'`. The material's
``refractive_index`` does not enter the implemented bulk steady diffusion
coefficient. Internal index changes are rejected by the finite solver, and
external diffuse reflectance is configured separately.

``InfiniteDiffusion`` evaluates the Green function in this chapter.
``Diffusion`` assembles the finite weak problem in :doc:`finite_element`.
``RobinBoundary`` selects a diffuse partial-current boundary law.
``result.fluence`` stores nodal :math:`\Phi`, while the finite solver's
``result.flux`` stores element-centered :math:`\mathbf J`. The source
wavelength selects coefficients from a supplied ``SpectralMaterial`` table.
A fixed ``Material`` retains its configured coefficients. Neither choice
infers an optical spectrum from the wavelength label.
