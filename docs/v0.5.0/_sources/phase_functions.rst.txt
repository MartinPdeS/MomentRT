Scattering phase functions
==========================

A phase function describes how a single scattering event redistributes light
with scattering angle. MomentRT's scalar, axisymmetric models use the cosine
:math:`\mu=\cos\theta` and return probability density per steradian,
normalized so that :math:`2\pi\int_{-1}^{1}p(\mu)\,d\mu=1`. They do not
represent polarization or azimuth-dependent scattering.

Available models
----------------

``MomentRT.phase_function.HenyeyGreenstein()`` uses the material's supplied
``anisotropy`` :math:`g` and has normalized Legendre moments
:math:`(1,g,g^2,g^3)`. Pass the phase-function object explicitly to each
material. A spectral material can provide a scalar or wavelength-dependent
anisotropy table.

Every ``Material`` and ``SpectralMaterial`` requires an explicit
``phase_function`` object. P1 does not use its higher moments, but the object
still records the user's physical choice.

``MomentRT.phase_function.Isotropic()`` has moments :math:`(1,0,0,0)`.
``MomentRT.phase_function.Rayleigh()`` is the unpolarized scalar Rayleigh
phase function, with moments :math:`(1,0,1/10,0)`. Their first moments are
zero, so material ``anisotropy`` may be omitted; if supplied, it must be zero.
Rayleigh and isotropic scattering thus have the same first moment but different
second moments.

``MomentRT.phase_function.Tabulated(cosine=..., values=...)`` takes a
strictly increasing plain cosine grid spanning exactly -1 to 1 and nonnegative
relative density values. It linearly interpolates in cosine, normalizes the
result over 4 pi, and integrates its first three Legendre moments. The table
is fixed across a spectral material's wavelength range; it does not automatically
compute wavelength-dependent Mie scattering. If an anisotropy is supplied, it
must agree with the table's calculated first moment.

.. code-block:: python

   import numpy as np
   import MomentRT

   ureg = MomentRT.ureg

   phase = MomentRT.phase_function.Tabulated(
       cosine=np.array([-1.0, 0.0, 1.0]),
       values=np.array([0.4, 1.0, 1.6]),
   )

   material = MomentRT.Material(
       absorption=0.03 / ureg.millimeter,
       reduced_scattering=1.0 / ureg.millimeter,
       refractive_index=1.4,
       phase_function=phase,
   )

   moments = phase.moments()

   print(f"Phase moments: {moments}")

The input values are relative densities per unit solid angle sampled at cosine
values, not a probability histogram per angle. Converting measured or Mie
data requires the correct angular measure and normalization. Tabulation is
only as accurate as the supplied resolution, especially near sharp forward
peaks. A user-defined ``PhaseFunction`` subclass can provide its own density
and moments; moments used by the solvers are checked for normalization and
compatibility with a nonnegative distribution on [-1, 1].

How the solvers use moments
---------------------------

The material's ``reduced_scattering`` is :math:`\mu_s'=\mu_s(1-g_1)`, where
:math:`g_1` is the selected phase function's first moment. P1 only uses that
coefficient. SP3 and P3 additionally use :math:`g_2` and :math:`g_3` through

.. math::

   \sigma_l = \mu_a + \mu'_s\frac{1-g_l}{1-g_1},\qquad l=1,2,3.

Consequently, changing the phase function while holding :math:`\mu_s'` fixed
can change SP3/P3 even when the first moment is unchanged. Phase-function
selection does not change the solver's source, boundary, or detector angular
assumptions. See :doc:`gallery/finite_samples/phase_functions` for a runnable
comparison and `NIST SCATMECH's phase-function reference
<https://pages.nist.gov/SCATMECH/docs/phase.htm>`_ for related scalar models.
