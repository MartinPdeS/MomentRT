"""
Seven-fiber bundle on a three-layer skin sample
================================================

Build a contact fiber bundle with six circular illumination footprints
surrounding a central collection aperture. The three layer optical properties
and fiber dimensions are illustrative, not measured probe parameters.

The layer boundaries follow a published three-layer skin model. Optical
coefficients are linearly interpolated to 850 nm from ex vivo measurements at
800 and 980 nm; they represent one Asian skin dataset, not universal skin
properties. The model uses a shared refractive index to satisfy the current
SP3 solver's matched-index interface assumption.

References
----------
Shimojo, Y. et al. (2020), *Journal of Biomedical Optics*, 25(4), 045002.
https://doi.org/10.1117/1.JBO.25.4.045002
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/seven_fiber_skin.png'

import os
import numpy as np
import MomentRT

ureg = MomentRT.ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Define the three tissue layers
# ------------------------------
# Positive z points into the sample, so add the epidermis, dermis, then
# hypodermis. The 4 mm sample contains 0.2 mm epidermis, 1.5 mm dermis, and
# 2.3 mm of the requested 2.35 mm hypodermis. Optical coefficients are
# approximate 850 nm values.
# SP3 currently requires equal refractive indices across internal interfaces.

epidermis = MomentRT.Material(
    absorption=0.329 / ureg.millimeter,
    reduced_scattering=3.75 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.9,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

dermis = MomentRT.Material(
    absorption=0.055 / ureg.millimeter,
    reduced_scattering=1.85 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.9,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

hypodermis = MomentRT.Material(
    absorption=0.04 / ureg.millimeter,
    reduced_scattering=1.33 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.9,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

sample = MomentRT.StructuredMedium(
    origin=(-6, -6, 0) * ureg.millimeter / 10,
    size=(12, 12, 40) * ureg.millimeter / 10,
)

sample.add_background(material=hypodermis)

sample.add_layer(
    name="Epidermis",
    thickness=0.20 * ureg.millimeter,
    material=epidermis,
    clip_to_sample=True,
)

sample.add_layer(
    name="Dermis",
    thickness=1.5 * ureg.millimeter,
    material=dermis,
    clip_to_sample=True,
)

sample.add_layer(
    name="Hypodermis (subcutaneous fat)",
    thickness=2.35 * ureg.millimeter,
    material=hypodermis,
    clip_to_sample=True,
)

# %%
# Arrange the seven fiber contacts
# --------------------------------
# Six illumination fibers are equally spaced on a ring around the collection
# fiber. Each source supplies one sixth of the total coupled optical power.
# Disk footprints represent diffuse light already coupled into the tissue;
# fiber propagation, entry losses, and coherent interference are not modeled.

fiber_radius = 100 * ureg.micrometer

ring_radius = 250 * ureg.micrometer

total_power = 1 * ureg.milliwatt

wavelength = 850 * ureg.nanometer

sources = []

for fiber_number in range(6):
    angle = fiber_number * np.pi / 3

    coordinates = np.array([np.cos(angle), np.sin(angle)]) * ring_radius

    source = MomentRT.SurfaceSource.on_face(
        sample=sample,
        face="z_min",
        coordinates=coordinates,
        radius=fiber_radius,
        wavelength=wavelength,
        optical_power=total_power / 6,
        name=f"Source {fiber_number + 1}",
    )

    sources.append(source)

detector = MomentRT.SurfaceDetector.on_face(
    sample=sample,
    face="z_min",
    coordinates=np.array([0, 0]) * ureg.millimeter,
    radius=fiber_radius,
    orientation="outward_normal",
    acceptance_angle=30 * ureg.degree,
    name="Collection",
)

mesh = MomentRT.CartesianMesh.from_sample(
    sample=sample,
    shape=(31 * 1, 21 * 1, 21 * 1),
)

boundary = MomentRT.RobinBoundary(effective_reflectance=0)

boundaries = MomentRT.Boundaries.uniform(condition=boundary)

# %%
# Inspect the tissue geometry
# ---------------------------
# The transparent surfaces show the epidermis, dermis, and hypodermis in
# surface-to-depth order; the resolved fluence and seven fiber contacts follow.

sample_figure = sample.plot(show=show_figures)

# %%
# Solve all six illumination fibers together with SP3
# ---------------------------------------------------
# Simulation converts the source list to a SurfaceSourceGroup and adds the six
# disk profiles into one boundary load, so the transport system is solved once.
# This is simultaneous diffuse power illumination, not coherent interference.
# SP3 is not the full 3D P3 system.

solver = MomentRT.SP3(boundaries=boundaries)

simulation = MomentRT.Simulation(
    medium=sample,
    source=sources,
    mesh=mesh,
    solver=solver,
    detectors={"collection": detector},
)

report = simulation.validate()

if not report.valid:
    raise ValueError(f"Seven-fiber configuration: {report.summary()}")

result = simulation.run()

print(
    f"Combined sampled source area: {result.source_sampled_area.to('mm^2'):~P}; "
    f"geometric area: {result.source_geometric_area.to('mm^2'):~P}"
)

collected_power = result.detector_power["collection"]

detector_result = result.detectors["collection"]

print(f"Total coupled illumination power: {total_power:~P}")

print(f"Central fiber collected power: {collected_power.to('mW'):~P}")

print(f"Collection sampled area: {detector_result.sampled_area.to('mm^2'):~P}")

print(f"Collection geometric area: {detector_result.geometric_area.to('mm^2'):~P}")

# %%
# Visualize the combined fluence and seven contacts
# -------------------------------------------------
# These curves show the six source footprints and the detector aperture on
# the tissue surface. The sampled footprint areas above differ from the ideal
# circular areas. Refine the mesh and compare the combined detector signal
# before interpreting a result; the approximations themselves also require
# validation against a suitable radiative-transfer or Monte Carlo reference.

overlays = [*sources, detector, *sample.structures]

fluence_figure = result.fluence.plot(
    view="volume",
    scale="linear",
    overlays=overlays,
    show=show_figures,
)


# %%
# Inspect the modeled geometry
# ----------------------------
# This view shows the configured geometry and instrument footprints.

gallery_system_figure = sample.plot(overlays=overlays, show=False)
