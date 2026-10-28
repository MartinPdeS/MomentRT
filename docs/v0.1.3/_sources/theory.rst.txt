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
used in diffusion are the fluence rate and vector current:

.. math::

   \Phi(\mathbf r)=\int_{4\pi}L(\mathbf r,\mathbf s)\,d\Omega,
   \qquad
   \mathbf J(\mathbf r)=\int_{4\pi}\mathbf sL(\mathbf r,\mathbf s)\,d\Omega.

The transport equation couples these moments to higher angular moments. The
P1 closure truncates the angular expansion after the first moment. Equivalently,
it approximates radiance by

.. math::

   L(\mathbf r,\mathbf s)\approx
   \frac{1}{4\pi}\left[\Phi(\mathbf r)+3\,\mathbf s\cdot\mathbf J(\mathbf r)\right].

This closure treats the angular distribution as nearly isotropic, with a
first-order directional correction. It does not retain a narrow beam or
resolve strongly directional propagation.

Diffusion equation
------------------

For a source whose first angular moment is zero, the first two moment
equations under P1 closure are

.. math::

   \nabla\cdot\mathbf J + \mu_a\Phi = q_0,

   \frac{1}{3}\nabla\Phi + (\mu_a+\mu_s')\mathbf J = 0,

where :math:`q_0` is the isotropic source power per volume and
:math:`\mu_s'=\mu_s(1-g)` is reduced scattering, with :math:`g` the mean
cosine of the scattering angle. Eliminating :math:`\mathbf J` gives Fick's
law and the scalar diffusion equation:

.. math::

   \mathbf J=-D\nabla\Phi,
   \qquad
   D=\frac{1}{3(\mu_a+\mu_s')},

.. math::

   -\nabla\cdot(D\nabla\Phi)+\mu_a\Phi=q_0.

MomentRT asks for :math:`\mu_a` and :math:`\mu_s'` directly, both in inverse
length units. It does not infer :math:`\mu_s'` from :math:`\mu_s` and :math:`g`.
In this convention :math:`D` has units of length, :math:`\Phi` and
:math:`\mathbf J` have units of :math:`\mathrm{W\,m^{-2}}`, and integrated
source and power terms have units of watts.

Infinite homogeneous point source
---------------------------------

For constant optical coefficients in unbounded three-dimensional space, an
isotropic point source with power :math:`P` gives

.. math::

   -D\nabla^2\Phi+\mu_a\Phi=P\,\delta(\mathbf r-\mathbf r_s).

The decaying Green function is

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
