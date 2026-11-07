"""
Simulate a contact sensor on layered skin
=========================================

Build an illustrative three-layer skin sample, generate a tetrahedral mesh,
and calculate diffuse fluence and the power reaching a contact detector. The
material coefficients describe one ex vivo dataset near 850 nm, not a
universal or patient-specific skin model. P1 and the Lambertian detector model
are approximations and require independent validation for the intended sensor.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/skin_contact.png'
import MomentRT
from MomentRT.models.skin import ShimojoThreeLayer850nm

ureg = MomentRT.ureg

# %%
# Build the sample and contact sensor
# -----------------------------------
# The preset uses 0.2 mm epidermis, 1.5 mm dermis, and a hypodermis background.
# This 4 mm sample clips the requested hypodermis layer to 2.3 mm. Source
# optical power is diffuse power already coupled into the tissue.

sample = ShimojoThreeLayer850nm(
    origin=(-3, -3, 0) * ureg.mm,
    size=(6, 6, 4) * ureg.mm,
)
source = MomentRT.SurfaceSource.on_face(
    sample=sample,
    face="z_min",
    coordinates=(-1, 0) * ureg.mm,
    radius=0.8 * ureg.mm,
    wavelength=850 * ureg.nm,
    optical_power=1 * ureg.mW,
)
detector = MomentRT.SurfaceDetector.on_face(
    sample=sample,
    face="z_min",
    coordinates=(1, 0) * ureg.mm,
    radius=0.8 * ureg.mm,
    orientation="outward_normal",
    acceptance_angle=90 * ureg.degree,
)

# %%
# Generate and inspect the mesh
# -----------------------------
# Refine near the contacts to sample their circular footprints. The plot
# displays a subset of tetrahedral edges, not every element. Resolution and
# source/detector sampling need a separate mesh-convergence study.

mesh = MomentRT.TetrahedralMesh.from_sample(
    sample=sample,
    resolution=0.4 * ureg.mm,
)
mesh.refine_near(
    objects=[source, detector],
    resolution=0.2 * ureg.mm,
    distance=0.3 * ureg.mm,
)
solver = MomentRT.P1(
    boundaries=MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0),
    ),
)
simulation = MomentRT.Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
    detectors={"contact": detector},
)

system_figure = simulation.plot_system(show=False)
mesh_figure = mesh.plot(max_edges=4000, show=False)
print(f"Tetrahedra: {len(mesh.tetrahedra)}")

# %%
# Solve and inspect the results
# -----------------------------
# The detector integrates sampled escaping irradiance over its aperture using
# a Lambertian angular model. The displayed volume is interpolated from nodal
# FEM fluence onto a separate grid. Neither visualization nor energy balance
# establishes spatial convergence or P1 accuracy near the source.

result = simulation.run()
measurement = result.detectors["contact"]
print(f"Collected power: {measurement.collected_power.to('uW'):~P}")
print(f"Sampled source area: {result.source_sampled_area.to('mm^2'):~P}")
print(f"Geometric source area: {result.source_geometric_area.to('mm^2'):~P}")
print(f"Sampled detector area: {measurement.sampled_area.to('mm^2'):~P}")
print(f"Geometric detector area: {measurement.geometric_area.to('mm^2'):~P}")
print(f"Minimum nodal fluence: {result.fluence.values.min().to('W/m^2'):~P}")

gallery_system_figure = result.fluence.plot(
    view="volume",
    scale="log",
    display_resolution=0.3 * ureg.mm,
    opacity=0.3,
    overlays=[source, detector],
    show=False,
)
