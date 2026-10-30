"""
Use different anisotropies in a layered sample
================================================

Assign a scattering anisotropy to each material and solve a layered sample
with SP3. The values are illustrative, not measured skin properties. SP3
uses the Henyey–Greenstein phase function for each material.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/sp3_layered_anisotropy.png'

import os
import numpy as np
import MomentRT

ureg = MomentRT.ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Set optical properties for each layer
# -------------------------------------
# Anisotropy belongs to the material, not the solver. SP3 requires a value
# for every region, including the background. This model assumes the same
# refractive index across internal interfaces.

dermis = MomentRT.Material(
    absorption=0.03 / ureg.millimeter,
    reduced_scattering=1.0 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.8,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

epidermis = MomentRT.Material(
    absorption=0.06 / ureg.millimeter,
    reduced_scattering=1.6 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.65,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

sample = MomentRT.StructuredMedium(
    origin=np.array([-4, -4, 0]) * ureg.millimeter,
    size=np.array([8, 8, 4]) * ureg.millimeter,
)

sample.add_background(material=dermis)

sample.add_layer(
    name="Epidermis",
    thickness=0.8 * ureg.millimeter,
    material=epidermis,
)

# %%
# Configure the source, detector, and solver
# ------------------------------------------
# The incoming source is diffuse and isotropic over the inward hemisphere.
# A zero-reflectance Robin face selects the matched-index SP3 boundary model.

source = MomentRT.SurfaceSource.on_face(
    sample=sample,
    face="z_min",
    coordinates=np.array([-1.5, 0]) * ureg.millimeter,
    radius=1.0 * ureg.millimeter,
    wavelength=850 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

detector = MomentRT.SurfaceDetector.on_face(
    sample=sample,
    face="z_min",
    coordinates=np.array([1.5, 0]) * ureg.millimeter,
    radius=0.8 * ureg.millimeter,
    orientation="outward_normal",
    acceptance_angle=90 * ureg.degree,
)

mesh = MomentRT.CartesianMesh.from_sample(
    sample=sample,
    shape=(15, 15, 11),
)

boundary = MomentRT.RobinBoundary(effective_reflectance=0)

boundaries = MomentRT.Boundaries.uniform(condition=boundary)

solver = MomentRT.SP3(boundaries=boundaries)

simulation = MomentRT.Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
    detectors={"contact": detector},
)

# %%
# Inspect geometry and solve the sample
# -------------------------------------
# Regions have different absorption, reduced scattering, and anisotropy.
# The displayed surfaces show the geometry used to construct the sample.

report = simulation.validate()

print(report.summary())

if not report.valid:
    raise ValueError("Fix the reported configuration errors before solving.")

system_figure = simulation.plot_system(show=show_figures)

result = simulation.run()

print(f"Dermis anisotropy: {dermis.anisotropy:.2f}")

print(f"Epidermis anisotropy: {epidermis.anisotropy:.2f}")

print(result.summary())

# %%
# View the physical fluence in three dimensions
# ----------------------------------------------
# The source and detector are overlaid on the solved fluence volume. Refine
# the mesh to assess both footprint sampling and the detector signal; SP3 is
# an approximation and this plot does not establish transport accuracy.

fluence_figure = result.fluence.plot(
    view="volume",
    scale="linear",
    overlays={"Source": source, "Contact detector": detector},
    show=show_figures,
)


# %%
# Inspect the modeled geometry
# ----------------------------
# This view shows the configured geometry and instrument footprints.

gallery_system_figure = simulation.plot_system(show=False)
