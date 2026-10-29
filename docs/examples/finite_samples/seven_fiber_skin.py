"""
Seven-fiber bundle on a two-layer skin sample
================================================

Build a contact fiber bundle with six circular illumination footprints
surrounding a central collection aperture. The epidermis and dermis optical
properties and fiber dimensions are illustrative, not measured skin or probe
parameters.
"""

# sphinx_gallery_thumbnail_path = '_static/images/seven-fiber-skin.png'

import os

import numpy as np

import MomentRT

ureg = MomentRT.ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Define the two tissue layers
# ----------------------------
# The lower dermis fills the background; the epidermis spans the top of the
# sample. SP3 currently requires equal refractive indices across internal
# interfaces and a material-specific Henyey–Greenstein anisotropy.

dermis = MomentRT.Material(
    absorption=0.025 / ureg.millimeter,
    reduced_scattering=1.0 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.8,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

epidermis = MomentRT.Material(
    absorption=0.06 / ureg.millimeter,
    reduced_scattering=1.6 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.7,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

sample = MomentRT.StructuredMedium(
    origin=np.array([-4, -4, 0]) * ureg.millimeter,
    size=np.array([8, 8, 4]) * ureg.millimeter,
)

sample.add_background(material=dermis)


sample.add_layer(
    name="Dermis",
    thickness=3.2 * ureg.millimeter,
    material=dermis,
)

sample.add_layer(
    name="Epidermis",
    thickness=0.8 * ureg.millimeter,
    material=epidermis,
)

# %%
# Arrange the seven fiber contacts
# --------------------------------
# Six illumination fibers are equally spaced on a ring around the collection
# fiber. Each source supplies one sixth of the total coupled optical power.
# Disk footprints represent diffuse light already coupled into the tissue;
# fiber propagation, entry losses, and coherent interference are not modeled.

fiber_radius = 0.45 * ureg.millimeter

ring_radius = 1.5 * ureg.millimeter

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
    )

    sources.append(source)

detector = MomentRT.SurfaceDetector.on_face(
    sample=sample,
    face="z_min",
    coordinates=np.array([0, 0]) * ureg.millimeter,
    radius=fiber_radius,
    orientation="outward_normal",
    acceptance_angle=30 * ureg.degree,
)

mesh = MomentRT.CartesianMesh.from_sample(
    sample=sample,
    shape=(21, 21, 11),
)

boundary = MomentRT.RobinBoundary(effective_reflectance=0)

boundaries = MomentRT.Boundaries.uniform(condition=boundary)

# %%
# Inspect the tissue geometry
# ---------------------------
# The transparent surfaces show the epidermis layer and dermis background;
# the resolved fluence and seven fiber contacts are shown together below.

sample_figure = sample.plot(show=show_figures)

# %%
# Solve the six illumination fibers with SP3
# -------------------------------------------
# MomentRT currently accepts one surface source per simulation. SP3 is
# linear for this fixed sample and boundary configuration, so six independent
# monochromatic solves can be summed to obtain the simultaneous diffuse
# fluence and collected power. This is incoherent power superposition, not a
# coherent optical field calculation. SP3 is not the full 3D P3 system.

solver = MomentRT.SP3(boundaries=boundaries)

results = []

for fiber_number, source in enumerate(sources, start=1):
    simulation = MomentRT.Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=solver,
        detectors={"collection": detector},
    )

    report = simulation.validate()

    if not report.valid:
        raise ValueError(f"Fiber {fiber_number} configuration: {report.summary()}")

    result = simulation.run()

    results.append(result)

    print(
        f"Fiber {fiber_number}: sampled source area "
        f"{result.source_sampled_area.to('mm^2'):~P}; "
        f"geometric area {result.source_geometric_area.to('mm^2'):~P}"
    )

# %%
# Combine the fluence and collection signals
# ------------------------------------------
# The detector is passive and is evaluated for each source independently.
# Its accepted power is therefore additive under this diffuse linear model.
# The angular acceptance uses a Lambertian escape assumption; this is not a
# model of guided fiber modes, numerical aperture coupling, or Fresnel loss.

fluence_values = sum(result.fluence.values for result in results)

combined_fluence = MomentRT.ScalarField(
    mesh=results[0].fluence.mesh,
    values=fluence_values,
    name="Combined SP3 fluence rate",
    metadata={
        **results[0].fluence.metadata,
        "source_count": len(sources),
        "superposition": "six independent diffuse sources",
    },
)

collected_power = sum(result.detector_power["collection"] for result in results)

detector_result = results[0].detectors["collection"]

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

overlays = {
    **{
        f"Illumination fiber {number}": source
        for number, source in enumerate(sources, start=1)
    },
    "Central collection fiber": detector,
}

fluence_figure = combined_fluence.plot(
    view="volume",
    scale="linear",
    overlays=overlays,
    show=show_figures,
)
