Transport theory and approximations
===================================

MomentRT currently implements steady-state P1 diffusion. It does not solve
the full radiative transfer equation (RTE), a higher-order :math:`P_N` system,
or an :math:`SP_N` approximation. Those models have different angular
closures and capabilities. The two implemented solvers use the same diffusion
physics in different geometries: an analytical Green function for an infinite
homogeneous medium and a finite-domain FEM solve.

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
<gallery/fundamentals/absorption_sweep>` compares the resulting radial profiles
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

MomentRT does not currently provide a finite-aperture detector model,
collimated-source propagation, time dependence, spectral integration,
chromophore fitting, or inverse estimation. Refractive index is stored with
materials, but no automatic Fresnel rule is derived from it. In the finite
solver all internal refractive indices must be equal; an external diffuse
boundary coefficient must be supplied explicitly.

What the order labels mean
--------------------------

:math:`P_1` is the first-order spherical-harmonic angular expansion used to
derive diffusion. :math:`P_N` denotes a general finite spherical-harmonic
expansion through order :math:`N`; it contains multiple coupled angular
moments. :math:`SP_N` is a different, simplified system derived from selected
even-order transport equations, often arranged as coupled diffusion-like
equations. Neither label means that a P1 solve secretly contains P2 or P3.
MomentRT currently solves the P1 diffusion equation above, not a P3 or SP3
system.
