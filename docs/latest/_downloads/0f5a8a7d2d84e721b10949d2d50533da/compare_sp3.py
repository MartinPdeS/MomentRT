"""
Compare P1 and SP3 transport
============================

Compare the coupled SP3 moment model and P1 diffusion in the same finite
sample with matched-index boundary conditions and incoming diffuse light.
This comparison illustrates model differences; it does not establish which
model matches radiative transfer for a particular sensor.
"""

# sphinx_gallery_thumbnail_path = '_static/images/sp3-fluence.png'

import os

import numpy as np

import MomentRT

ureg = MomentRT.ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Define the sample and instrument
# --------------------------------
# Both models receive the same material, source, detector, mesh, and faces.
# SP3 additionally uses each material's Henyey–Greenstein anisotropy.

material = MomentRT.Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    absorption=0.03 / ureg.millimeter,
    reduced_scattering=1 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.8,
)

sample = MomentRT.StructuredMedium(
    origin=np.array([-4, -4, 0]) * ureg.millimeter,
    size=np.array([8, 8, 4]) * ureg.millimeter,
)

sample.add_background(material=material)

mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=(13, 13, 9))

source = MomentRT.SurfaceSource(
    face="z_min",
    position=np.array([0, 0, 0]) * ureg.millimeter,
    radius=1.5 * ureg.millimeter,
    wavelength=970 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

detector = MomentRT.SurfaceDetector(
    face="z_min",
    position=np.array([2.5, 0, 0]) * ureg.millimeter,
    radius=0.8 * ureg.millimeter,
    collection_axis=[0, 0, -1],
    acceptance_angle=90 * ureg.degree,
)

boundary = MomentRT.RobinBoundary(effective_reflectance=0)

boundaries = {face: boundary for face in MomentRT.boundary.FACES}

# %%
# Run both approximations
# -----------------------
# A zero-reflectance Robin face means a matched-index boundary for this SP3
# formulation. The illumination is isotropic in the incoming hemisphere.

p1_solver = MomentRT.Diffusion(boundaries=boundaries)

sp3_solver = MomentRT.SP3(boundaries=boundaries)

p1_simulation = MomentRT.Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=p1_solver,
    detectors={"contact": detector},
)

sp3_simulation = MomentRT.Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=sp3_solver,
    detectors={"contact": detector},
)

p1_result = p1_simulation.run()

sp3_result = sp3_simulation.run()

print(f"P1 collected power: {p1_result.detector_power['contact'].to('mW'):~P}")

print(f"SP3 collected power: {sp3_result.detector_power['contact'].to('mW'):~P}")

# %%
# Inspect the SP3 fluence
# -----------------------
# The displayed fluence combines both SP3 moments. Refine the mesh and
# compare against an independent transport reference before interpreting a
# particular source–detector configuration.

figure = sp3_result.fluence.plot(
    view="volume",
    scale="linear",
    overlays=[source, detector],
    show=show_figures,
)
