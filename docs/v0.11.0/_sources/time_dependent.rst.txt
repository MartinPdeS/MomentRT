Time-dependent P1 and the telegrapher equation
============================================================

P1 is an angular approximation. Steady state is a separate assumption.
Keeping time derivatives in the first two P1 moment equations produces a
system for fluence and current that can be written as a telegrapher equation.
Neglecting the time derivative of current produces the usual time-dependent
diffusion equation.

.. important::

   ``TransientP1`` solves the coupled time-dependent P1 equations on finite
   samples. It supports monochromatic surface sources, a sampled pulse, and
   time-resolved contact-detector output. Frequency-domain modulation and a
   time-dependent diffusion solver are not implemented.

Adding time to radiative transfer
---------------------------------

Let :math:`v` be the energy propagation speed in a homogeneous,
nondispersive medium. The common approximation :math:`v=c_0/n` uses the
vacuum light speed :math:`c_0` and refractive index :math:`n`. For dispersive
media, the appropriate energy or group speed needs a separate physical
choice; a steady refractive-index parameter does not settle that question.

The time-dependent transport equation is

.. math::

   \frac{1}{v}\frac{\partial L}{\partial t}
   +\mathbf s\cdot\nabla L+(\mu_a+\mu_s)L
   =\mu_s\int_{4\pi}p(\mathbf s,\mathbf s')L(\mathbf s',t)\,d\Omega'
   +Q.

The angular moments use the same definitions as in :doc:`theory`, now with
time as an additional argument. For constant :math:`v`, the radiation energy
density is :math:`u=\Phi/v`, with units :math:`\mathrm{J\,m^{-3}}`.
Consequently :math:`(1/v)\partial_t\Phi` represents a local rate of energy
storage rather than a loss mechanism.

For an isotropic volume source, :math:`\mathbf q_1=\mathbf 0`. Define
:math:`\mu_{\mathrm{tr}}=\mu_a+\mu_s'`. Applying the P1 closure gives

.. math::

   \frac{1}{v}\partial_t\Phi+\nabla\cdot\mathbf J+\mu_a\Phi=q_0,
   \qquad
   \frac{1}{v}\partial_t\mathbf J+\frac{1}{3}\nabla\Phi
   +\mu_{\mathrm{tr}}\mathbf J=\mathbf 0.

The first equation conserves radiant energy. The second describes how the
directional current changes: a fluence gradient drives flow, while absorption
and transport scattering relax its directional memory. Derivations retaining
these two time derivatives are presented in the
`radiative transfer P1 study hosted by OSTI
<https://www.osti.gov/servlets/purl/878230>`_.

Current relaxation and the diffusion approximation
--------------------------------------------------

Dividing the current equation by :math:`\mu_{\mathrm{tr}}` and defining
:math:`D=1/(3\mu_{\mathrm{tr}})` gives

.. math::

   \tau_J\partial_t\mathbf J+\mathbf J=-D\nabla\Phi,
   \qquad
   \tau_J=\frac{1}{v\mu_{\mathrm{tr}}}.

The current relaxation time :math:`\tau_J` is the time to travel one
transport mean free path. Current does not adjust instantly to a new
gradient in this system. For an approximately fixed local gradient, the
homogeneous part of its relaxation decays as
:math:`\exp(-t/\tau_J)`.

If variations occur over times :math:`T` for which
:math:`\tau_J/T\ll1`, dropping the current time derivative gives
:math:`\mathbf J\approx-D\nabla\Phi`. Substituting into energy conservation
then yields the parabolic diffusion equation

.. math::

   \frac{1}{v}\partial_t\Phi
   -\nabla\cdot(D\nabla\Phi)+\mu_a\Phi=q_0.

There are two common diffusion-coefficient conventions:

.. math::

   D=\frac{1}{3\mu_{\mathrm{tr}}}\quad[\mathrm m],
   \qquad
   \mathcal D=vD\quad[\mathrm{m^2\,s^{-1}}].

Multiplying the whole diffusion equation by :math:`v` gives
:math:`\partial_t\Phi-\nabla\cdot(\mathcal D\nabla\Phi)
+v\mu_a\Phi=vq_0` when :math:`v` is spatially constant.
MomentRT uses the first convention, so its steady :math:`D` must not be
inserted as a diffusivity in :math:`\mathrm{m^2\,s^{-1}}`.

Eliminating current without dropping its time derivative
--------------------------------------------------------

The scalar derivation below assumes spatially and temporally constant
:math:`v`, :math:`\mu_a`, and :math:`\mu_{\mathrm{tr}}`, and an isotropic
source. With heterogeneous coefficients, spatial derivatives of those
coefficients must be accounted for; the coupled moment equations are the
safer starting point.

Differentiate energy conservation with respect to time:

.. math::

   \frac{1}{v}\partial_t^2\Phi
   +\nabla\cdot\partial_t\mathbf J
   +\mu_a\partial_t\Phi=\partial_t q_0.

Take the divergence of the current equation:

.. math::

   \nabla\cdot\partial_t\mathbf J
   =-\frac{v}{3}\nabla^2\Phi
   -v\mu_{\mathrm{tr}}\nabla\cdot\mathbf J.

Use the original conservation law to replace
:math:`\nabla\cdot\mathbf J=q_0-\mu_a\Phi-(1/v)\partial_t\Phi`.
After collecting terms and multiplying by :math:`v`, the result is

.. math::

   \boxed{
   \partial_t^2\Phi
   +v(\mu_{\mathrm{tr}}+\mu_a)\partial_t\Phi
   -\frac{v^2}{3}\nabla^2\Phi
   +v^2\mu_{\mathrm{tr}}\mu_a\Phi
   =v\partial_tq_0+v^2\mu_{\mathrm{tr}}q_0.
   }

Both absorption contributions must be retained: absorption appears in energy
conservation and in the loss of directional current. This is why the damping
coefficient contains :math:`\mu_{\mathrm{tr}}+\mu_a`.
The derivative of the source also follows from the elimination; it cannot
be omitted for an arbitrary time-varying source.

In a nonabsorbing source-free homogeneous region, this simplifies to

.. math::

   \tau_J\partial_t^2\Phi+\partial_t\Phi
   =\mathcal D\nabla^2\Phi.

This is the classical damped-wave, or telegrapher, form. It is related to
the equations for signals on electrical transmission lines through its
mathematical structure. It describes propagation and relaxation of an
intensity distribution, not electromagnetic phase oscillations.

Propagation speed and validity
-------------------------------

For the homogeneous P1 system, the characteristic propagation speed is

.. math::

   v_{\mathrm{P1}}=\sqrt{\frac{\mathcal D}{\tau_J}}
   =\frac{v}{\sqrt{3}}.

The finite front speed corrects the mathematical instantaneous tails of a
parabolic diffusion solution. It is still smaller than the transport speed
:math:`v` of unscattered photons. Retaining current inertia therefore does
not make P1 an accurate ballistic or early-arrival model. The angular closure
still suppresses strongly directional structure. The difference between
diffusion, conventional P1, and alternative transient closures is examined
in `Diffusion, P1, and other approximate forms of radiation transport
<https://doi.org/10.1016/S0022-4073(99)00150-8>`_.

For times much longer than current relaxation, diffusion can describe the
slowly varying part of the solution. With nonzero absorption, it is clearer
to take that approximation in the coupled current equation: merely deleting
the second derivative from the general scalar equation can leave an extra
absorption-dependent factor on :math:`\partial_t\Phi`.

Initial and boundary data
--------------------------

The coupled transient system requires initial fluence and current. A scalar
second-order equation instead requires initial fluence and its initial time
derivative, constrained by conservation:

.. math::

   \partial_t\Phi(\mathbf r,0)
   =v\left[q_0(\mathbf r,0)-\mu_a\Phi(\mathbf r,0)
   -\nabla\cdot\mathbf J(\mathbf r,0)\right].

An arbitrary choice of both scalar initial data may not correspond to the
intended initial angular field. A transient boundary model must also be
specified. The partial-current reasoning in :doc:`finite_element` can be
applied to instantaneous P1 moments, but extending the solver requires a
consistent treatment of time-dependent sources, interfaces, and boundary
conditions rather than adding one time derivative to the steady matrix.

Finite-element implementation
-----------------------------

``TransientP1`` solves the coupled equations above rather than applying the
homogeneous scalar telegrapher equation across a heterogeneous sample. The
fluence uses continuous linear tetrahedral nodal basis functions, while the
three current components are constant within each tetrahedron. Optical
coefficients are sampled at element centroids; they may vary between elements.
The discrete current is not smoothed across interfaces. Internal refractive
indices must match and the energy speed is an explicit uniform input.

The instantaneous outward P1 partial-current condition on a Robin face is

.. math::

   \mathbf J\cdot\mathbf n = \beta\Phi - \frac{2I_\mathrm{in}}{1+R},
   \qquad \beta=\frac{1-R}{2(1+R)}.

Here :math:`R` is the supplied diffuse effective reflectance and
:math:`I_\mathrm{in}` is diffuse incoming partial current already coupled into
the sample. Reflecting faces have zero normal current. The face term is
inserted in the weak conservation equation at each time step; it is not a
normal-incidence Fresnel model. Boundary escaping irradiance and contact
detectors retain the existing Lambertian angular approximation.

The method starts with zero fluence and current. A dimensionless nonnegative
``power_profile`` multiplies the monochromatic source power at the supplied
``times`` and must start at zero. It is linearly interpolated over each time
interval. Implicit midpoint time integration advances both P1 moments together;
the element current can be eliminated locally for each step, leaving one
sparse nodal solve. The time step may vary. Midpoint is second-order for smooth
solutions but does not suppress all rapidly oscillating numerical modes; resolve
the pulse and current-relaxation timescale, and inspect signed results.

The returned ``TransientResult`` contains ``times`` and per-time ``frames``
with nodal fluence, element-centered current, and boundary fields. A simulation
with named detectors also provides ``detector_power`` arrays. It reports
incoming, absorbed, and escaping powers and stored energy at each time; the
``interval_energy_residual`` checks midpoint energy conservation. Negative
values are not clipped. Storing every frame costs memory proportional to the
number of time samples and mesh elements.

For each interval, the reported residual is

.. math::

   \frac{U_{n+1}-U_n}{t_{n+1}-t_n}
   +\overline{P_\mathrm{abs}}+\overline{P_\mathrm{escape}}
   -\overline{P_\mathrm{in}},
   \qquad U=\int_\Omega \frac{\Phi}{v}\,dV,

where overbars denote endpoint averages. The difference between escaping and
incoming power is the signed net outward current. Detector-aperture power is
reported separately and is not part of this whole-boundary energy balance.

Use the `pulsed layered sample example
<gallery/transport_models/transient_p1.html>`_ for a complete calculation.
The solver handles one monochromatic source configuration at a time; spectral
sweeps, broadband integration, nonzero initial fields, and the steady-state
detector-convergence helper are not implemented for transient results. Refine
the mesh, source and detector footprints, and time grid independently. Compare
against an appropriate transport reference for the intended geometry: energy
conservation or convergence of the discrete P1 equations alone does not prove
P1 is accurate for early-time or near-source light transport.

Steady state and modulation
----------------------------

Setting all time derivatives to zero recovers

.. math::

   \mathbf J=-D\nabla\Phi,\qquad
   -\nabla\cdot(D\nabla\Phi)+\mu_a\Phi=q_0.

This is the implemented bulk model. In steady state there is no additional
telegrapher time scale for the solver to resolve.

For reference, with the harmonic convention :math:`e^{i\omega t}`, the
full homogeneous P1 current law becomes

.. math::

   \widehat{\mathbf J}=-D_\omega\nabla\widehat\Phi,\qquad
   D_\omega=\frac{1}{3(\mu_{\mathrm{tr}}+i\omega/v)}.

Its fluence equation is

.. math::

   -\nabla\cdot(D_\omega\nabla\widehat\Phi)
   +(\mu_a+i\omega/v)\widehat\Phi=\widehat q_0.

The approximation :math:`|\omega|\tau_J\ll1` replaces :math:`D_\omega`
by the real steady coefficient :math:`D`. This is the distinction between
frequency-domain P1 and frequency-domain diffusion. These complex equations
are explanatory; the current real-valued steady solver does not implement
them.
