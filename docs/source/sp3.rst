SP3 transport model
===================

The ``SP3`` solver solves a steady-state simplified third-order moment model
on the same finite tetrahedral mesh as ``Diffusion``. Its unknowns are two
coupled composite moments. In planar one-dimensional geometry these equations
coincide with the corresponding P3 moment equations; in three dimensions SP3
is a simplification of P3, not the full spherical-harmonics system.

Scattering and bulk equations
-----------------------------

The reduced scattering coefficient alone does not determine the higher
angular moments needed by SP3. Each material supplies a scalar phase function
whose normalized Legendre moments are :math:`g_0=1,g_1,g_2,g_3`. For
absorption :math:`\mu_a` and reduced scattering
:math:`\mu'_s=\mu_s(1-g_1)`, define

.. math::

   \mu_l=\mu_a+\mu'_s\frac{1-g_l}{1-g_1},\quad l=1,2,3.

The Henyey–Greenstein model requires material ``anisotropy``
:math:`g=g_1` and has :math:`g_l=g^l`. For that specific model,

.. math::

   \mu_1 &= \mu_a + \mu'_s,\\
   \mu_2 &= \mu_a + \mu'_s(1+g),\\
   \mu_3 &= \mu_a + \mu'_s(1+g+g^2),\\
   D_1 &= (3\mu_1)^{-1},\qquad D_2 = (7\mu_3)^{-1}.

The two composite nodal fields :math:`u` and :math:`v` obey

.. math::

   -\nabla\!\cdot(D_1\nabla u) + \mu_a u - \tfrac23\mu_a v &= 0,\\
   -\nabla\!\cdot(D_2\nabla v) - \tfrac23\mu_a u
      + \tfrac19(4\mu_a+5\mu_2)v &= 0.

Each coefficient is evaluated with the local material's phase-function
moments at the source wavelength. See :doc:`phase_functions` for the available
models. The physical fluence and first angular current are

.. math::

   \Phi = u - \tfrac23 v,\qquad \mathbf J=-D_1\nabla u.

Both composite fields are available as ``result.moments['composite_1']`` and
``result.moments['composite_2']``. ``result.fluence`` contains :math:`\Phi`,
and ``result.flux`` contains element-centered :math:`\mathbf J`.

Marshak face conditions
-----------------------

Let :math:`j_1=-D_1\mathbf n\cdot\nabla u` and
:math:`j_2=-D_2\mathbf n\cdot\nabla v`, with :math:`\mathbf n` pointing
outward. On a nonreflecting, matched-index face, the incoming radiance is
assumed constant over the inward hemisphere. Its partial current is
:math:`I_{\mathrm{in}}`. To derive the face conditions, use the local planar
angular expansion in :math:`\mu=\mathbf\Omega\cdot\mathbf n`:

.. math::

   L(\mu)=\frac{1}{4\pi}\left[\Phi+3\mu j_1
   +\frac53 v P_2(\mu)+7j_2P_3(\mu)\right].

The second angular moment is :math:`v/3`. For a constant prescribed incoming
radiance, the two half-range conditions are

.. math::

   \int_{\mu<0}(-\mu)L\,d\Omega
       &= \tfrac14u-\tfrac1{16}v-\tfrac12j_1=I_{\mathrm{in}},\\
   \int_{\mu<0}P_3(\mu)L\,d\Omega
       &= \tfrac1{16}u-\tfrac7{48}v+\tfrac12j_2
       =\tfrac14I_{\mathrm{in}}.

Solving these two equations for the outward currents gives

.. math::

   j_1 &= \tfrac12 u-\tfrac18 v-2I_{\mathrm{in}},\\
   j_2 &= -\tfrac18 u+\tfrac7{24}v+\tfrac12 I_{\mathrm{in}}.

The second source term follows from integrating
:math:`(5|\mu|^3-3|\mu|)` over the inward hemisphere. It cannot be inferred
from incoming power for an arbitrary angular distribution; the isotropic
incoming assumption fixes it. ``UniformSurfaceSource`` and ``SurfaceSource``
use their usual spatial footprints and coupled optical power, now with this
explicit angular assumption. A reflecting face has :math:`j_1=j_2=0` and
cannot be illuminated.

Only ``RobinBoundary(effective_reflectance=0)`` is supported by SP3. A single
P1 effective reflectance does not specify the higher angular reflection
moments required by an SP3 refractive boundary. Nonzero values are rejected,
as are internal refractive-index jumps. No index matching is inferred from
the material's numerical refractive index: the user explicitly chooses the
matched-index boundary model.

On a nonreflecting face the outgoing partial current is
:math:`j_1+I_{\mathrm{in}} = u/2-v/8-I_{\mathrm{in}}`. These values are
integrated at boundary quadrature points for face powers and detector
apertures. Detector angular acceptance continues to use MomentRT's separate
Lambertian escape approximation; SP3 does not provide Fresnel-resolved
detector optics. Negative fields are retained and reported rather than clipped.

Finite-element workflow
-----------------------

Multiplying each bulk equation by a nodal test function and integrating by
parts gives two diffusion stiffness matrices, three bulk mass terms, and a
coupled boundary mass matrix with coefficients :math:`1/2`, :math:`-1/8`,
and :math:`7/24`. The boundary load is :math:`(2I_{\mathrm{in}},
-I_{\mathrm{in}}/2)`. The resulting symmetric sparse system has twice as many
nodal unknowns as P1. Material coefficients are sampled at tetrahedron
centroids; source and detector footprints use boundary triangle quadrature.
Spectral sweeps reuse the geometry and quadrature, then assemble and solve a
new coupled operator at each wavelength.

The test suite independently solves homogeneous and layered planar SP3 systems using a
matrix exponential and checks mesh refinement, outgoing power, and energy
balance. This validates the implementation against that SP3 problem; it does
not validate SP3 against radiative transfer or biological measurements.
Check mesh convergence, and compare against Monte Carlo or other transport
references for the intended sensor distances and tissue properties.

For runnable examples, see :doc:`gallery/finite_samples/sp3_layered_anisotropy`
for material-specific anisotropy and
:doc:`gallery/finite_samples/sp3_spectral_anisotropy` for a wavelength sweep.

Reference: H. K. Kim et al., *Frequency-domain optical tomographic image
reconstruction algorithm with the simplified spherical harmonics (SP3) light
propagation model*, International Journal of Thermal Sciences 116 (2017),
265–277. `Author-hosted paper <https://bpb-us-e1.wpmucdn.com/wp.nyu.edu/dist/6/19597/files/2020/11/2017_IJTS_Hyun_FD_DOT_SP3.pdf>`_.
