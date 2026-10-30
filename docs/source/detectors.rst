Finite-aperture surface detectors
=================================

``SurfaceDetector`` collects diffuse escaping power from a circular contact
aperture on one sample face. Position and radius define the spatial aperture;
an outward collection axis and an acceptance half-angle define the accepted
directions. Detectors are passive measurements of the solved field and do
not modify the boundary conditions or remove power from the simulation.

This distinguishes a detector signal from ``outgoing_power_by_face``:
the latter integrates all modeled escaping power across an entire face.
For a complete runnable configuration and 3D geometry, see the
:doc:`detector example <gallery/finite_samples/surface_detector>`.

Defining a detector
--------------------

.. code-block:: python

   detector = SurfaceDetector(
       face="z_min",
       position=np.array([1.5, 0, 0]) * ureg.mm,
       radius=0.75 * ureg.mm,
       collection_axis=np.array([0, 0, -1]),
       acceptance_angle=30 * ureg.degree,
   )

The position is an absolute 3D coordinate on the selected face, and the
entire disk must fit within that face. The aperture plane coincides with the
sample face. ``collection_axis`` is a finite nonzero plain vector, normalized
by the constructor. It points along outgoing ray propagation and must have
a positive dot product with the outward face normal. On ``z_min`` this
normal is :math:`-\mathbf e_z`.

The collection axis may tilt within the outward hemisphere. Tilting it
changes angular acceptance; it does not rotate the contact aperture plane.
This models an angular acceptance cone attached to a surface patch, not an
arbitrary tilted photodiode in free space. The acceptance half-angle must be
an explicit quantity from zero to 90 degrees. Zero accepts no power.

For an ideal detector that collects all escaping light through the aperture,
use the full outward hemisphere:

.. code-block:: python

   detector = SurfaceDetector.full_hemisphere_on_face(
       sample=sample,
       face="z_min",
       coordinates=np.array([1.5, 0]) * ureg.mm,
       radius=0.75 * ureg.mm,
   )

This sets the axis to the outward normal and the acceptance half-angle to
90 degrees; it does not model external detector optics.

Attach named detectors to the configuration:

.. code-block:: python

   simulation = Simulation(
       medium=sample,
       source=source,
       mesh=mesh,
       solver=solver,
       detectors={"reflection": detector},
   )

   result = simulation.run()

   measurement = result.detectors["reflection"]

   print(f"Collected power: {measurement.collected_power.to('mW'):~P}")

   print(f"Sampled aperture area: {measurement.sampled_area.to('mm^2'):~P}")

The shorter ``result.detector_power["reflection"]`` returns the same power
quantity. A detector can also read an existing finite solution with
``detector.measure(result=result)``. This uses stored boundary snapshots;
changing a source object afterwards does not change the solved outgoing
field. ``simulation.plot_system()`` draws the aperture and collection axis
alongside the source and sample.

The same detector can be drawn over a solved field without constructing a
separate system figure:

.. code-block:: python

   figure = result.fluence.plot(
       overlays=[source, detector],
       scale="log",
       show=False,
   )

Set ``name`` and ``color`` on the detector to choose its legend label and
plot color. Adding an overlay is visual composition and does not
attach a new detector measurement to the result.

Escaping flux and angular acceptance
-------------------------------------

The Robin closure in :doc:`finite_element` gives escaping irradiance

.. math::

   E_{\mathrm{esc}}(\mathbf r)
   =(1-R)J_{\mathrm{out}}(\mathbf r)
   =\beta\Phi(\mathbf r)-2\beta I_{\mathrm{in}}(\mathbf r).

On unilluminated faces, the incoming term vanishes. On the illuminated face,
it must be subtracted so that a detector reads the transmitted outward
partial current rather than the gross Robin leakage or the net current.

An effective diffuse reflectance alone does not specify the transmitted
angular distribution. ``SurfaceDetector`` therefore makes an additional,
explicit Lambertian escape assumption:

.. math::

   L_{\mathrm{esc}}(\mathbf r,\mathbf s)=\frac{E_{\mathrm{esc}}(\mathbf r)}{\pi},
   \qquad \mathbf s\cdot\mathbf n>0.

The factor :math:`\pi` follows from
:math:`\int_{\mathbf s\cdot\mathbf n>0}\mathbf s\cdot\mathbf n\,d\Omega=\pi`.
This is a constant outgoing radiance whose projected integral matches the
computed escaping irradiance. The radiometric cosine factor is explained in
`Working with Radiometric Integrals
<https://www.pbr-book.org/3ed-2018/Color_and_Radiometry/Working_with_Radiometric_Integrals>`_.
The Lambertian reconstruction is an additional collection model, not an
angularly resolved RTE result.

For collection axis :math:`\mathbf a` and half-angle :math:`\theta_a`,
define the accepted fraction

.. math::

   F=\frac{1}{\pi}
   \int_{\substack{\mathbf s\cdot\mathbf n>0\\
                   \mathbf s\cdot\mathbf a\ge\cos\theta_a}}
   \mathbf s\cdot\mathbf n\,d\Omega.

The collected power over aperture :math:`A_d` is

.. math::

   P_d=F\int_{A_d}E_{\mathrm{esc}}(\mathbf r)\,dA.

For an axis aligned with the outward normal,
:math:`F=\sin^2\theta_a`. A 30-degree half-angle accepts one quarter of
Lambertian escaping power, while 90 degrees accepts the complete outward
hemisphere. With tilt :math:`\alpha` and a cone wholly inside the outward
hemisphere, :math:`F=\cos\alpha\sin^2\theta_a`. For a cone that crosses
the hemisphere edge, MomentRT integrates the clipped projected solid angle
numerically. There is no extra multiplication by aperture tilt: the spatial
aperture remains in the sample face.

Quadrature and resolution
--------------------------

Finite solutions retain ``result.boundary_fields[face]`` as
``BoundaryField`` snapshots. Each contains triangle quadrature positions,
area weights, and escaping irradiance, plus sample extent and wavelength.
The property ``escaping_power`` integrates that entire boundary field and
agrees with the corresponding face power.

For a circular detector, only quadrature points inside the aperture are
included:

.. math::

   A_{d,h}=\sum_{f,q}w_{fq}\chi_d(\mathbf r_{fq}),\qquad
   P_{d,h}=F\sum_{f,q}w_{fq}\chi_d(\mathbf r_{fq})E_{\mathrm{esc},fq}.

The aperture is not normalized to :math:`\pi r_d^2`. Both sampled and
geometric areas are reported in ``DetectorResult``, allowing direct
resolution checks. A disk with no quadrature points is rejected. Refine the
mesh around small detectors and compare both areas and collected power
across meshes.

The outgoing field uses the same three-point boundary quadrature and
incoming-source normalization as the FEM loading. This also keeps the
incident correction consistent when a detector overlaps a source footprint.
Negative escaping irradiance caused by a deficient P1 solution or unresolved
configuration is not clipped: the result reports
``minimum_escaping_irradiance`` and preserves the resulting signed power.
Such a value requires investigation before physical interpretation.

Scope of the measurement
-------------------------

Use :doc:`convergence` to compare collected power and sampled source and
detector areas on explicit successively refined meshes. The study retains
unresolved levels and reports a recommendation only when the stated
criteria hold at the fine end of the sequence.

This detector model applies to finite ``Diffusion``, ``SP3``, and ``P3`` results. An observation
window from ``InfiniteDiffusion`` is not a sample surface and cannot supply
these measurements. Reflecting faces have zero modeled escape.

The detector does not model angle-dependent Fresnel transmission, refraction,
lenses, a stand-off distance, quantum efficiency, electrical responsivity, or
noise. Acceptance is specified as an angular cone, not inferred from a
numerical aperture. Detector results retain a copy of the aperture
configuration, wavelength, and angular model for provenance.

Several detectors read the same undisturbed field independently. Summing
overlapping apertures can count the same escaping power more than once;
these results are separate passive instrument configurations rather than a
model of competing absorbers on the boundary.

Spectral measurements
---------------------

For a spectral sweep, ``spectrum.detectors["reflection"]`` returns the
individual :class:`~MomentRT.DetectorResult` measurements in the order of
``spectrum.wavelengths``. Each provides collected power and both sampled and
geometric aperture areas. These are independent monochromatic measurements,
not a broadband integral. ``spectrum.detector_power["reflection"]`` remains
available as a unit-bearing power array.
