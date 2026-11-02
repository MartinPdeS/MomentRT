Sensor-focused model validation
===============================

This benchmark compares P1, SP3, and P3 detector-power predictions with an
independent Monte Carlo photon-transport reference for two layered samples. It
varies the source-detector separation and records detector power, Monte Carlo
uncertainty, mesh size, source and aperture quadrature areas, and runtime. The
plot reports each moment model's relative detector-power error against the
Monte Carlo estimate; the shaded band shows the reference's approximate 95%
sampling uncertainty.

At the finest tested mesh (8,844 tetrahedra, with 96.2% of the exact aperture
area sampled), relative detector-power errors were:

.. list-table:: Finest-mesh detector-power error against Monte Carlo
   :header-rows: 1
   :widths: 24 18 16 16 16 16

   * - Optical case
     - Separation (mm)
     - MC 95% uncertainty
     - P1 error
     - SP3 error
     - P3 error
   * - Moderate absorption
     - 0.75
     - ±2.2%
     - +105.0%
     - +92.6%
     - +24.9%
   * - Moderate absorption
     - 1.5
     - ±4.7%
     - +59.8%
     - +55.0%
     - −1.7%
   * - Moderate absorption
     - 2.5
     - ±9.5%
     - +47.7%
     - +47.8%
     - −3.1%
   * - Higher absorption
     - 0.75
     - ±2.4%
     - +152.6%
     - +136.9%
     - +47.0%
   * - Higher absorption
     - 1.5
     - ±4.7%
     - +90.5%
     - +86.9%
     - −2.6%
   * - Higher absorption
     - 2.5
     - ±9.3%
     - +77.4%
     - +82.4%
     - +1.4%

The errors are signed ``(model − Monte Carlo) / Monte Carlo`` percentages.
The table is one mesh level, not an extrapolated continuum result; use the
plotted sequence to assess mesh behavior. The complete data, including each
mesh level's power and area diagnostics, are available as
:download:`JSON <_static/images/sensor_validation/sensor_validation.json>`.

.. image:: _static/images/sensor_validation/sensor_validation.png
   :alt: Relative detector-power error for P1, SP3, and P3 across mesh sizes,
         source-detector separations, and two layered optical-property cases.
   :align: center
   :width: 100%

Reproduce the benchmark from the repository root with::

   .venv/bin/python tools/validate_sensor_models.py

The command writes ``sensor_validation.json`` and the plot next to each other
under ``docs/source/_static/images/sensor_validation``. It uses one million
photon histories per optical-property case and a fixed random seed. Use
``--photons`` and ``--seed`` to choose another run; use ``--output-dir`` to
write elsewhere. The JSON file contains the full detector-by-detector mesh
study and sampled-area diagnostics.

Reference definition and limits
-------------------------------

The independent reference launches cosine-weighted diffuse photons uniformly
over the source disk. It samples Henyey-Greenstein scattering using each
layer's anisotropy and converts reduced scattering to total scattering. It
uses continuous absorption weighting and tallies the weighted photons that
leave through each circular top-face aperture. All six sample faces transmit
without Fresnel reflection, matching the benchmark's zero-reflectance Robin
boundaries. The detector accepts the full outward hemisphere, matching the
Lambertian angular acceptance fraction of one used by these MomentRT
detectors.

These are deliberately matched boundary and detector conventions, not a
complete model of an optical instrument. The reference uses independent
photon histories and reports their sampling error; it does not include
Fresnel reflection, a collimated LED, external optics, or detector
responsivity. P1, SP3, and P3 all solve the same finite layered geometry, but
their transport approximations can differ from the angularly resolved
reference. Mesh refinement cannot remove model-form error. Read the plotted
errors together with the Monte Carlo uncertainty and area fractions in the
JSON; agreement in energy balance or sampled area alone is not evidence of
transport accuracy.

The benchmark covers illustrative layered cases and a small set of lateral
separations. It is a diagnostic baseline, not a claim of accuracy for skin,
all source-detector distances, or every contact-sensor design.
