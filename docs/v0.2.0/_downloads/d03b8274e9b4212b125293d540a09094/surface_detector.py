"""
Collect light with finite apertures
============================================================

Place two circular contact detectors on a finite optical phantom. Each
aperture samples escaping power over a disk and accepts outgoing rays within
an explicit cone. Collection uses a Lambertian escape assumption; it does not
model angle-dependent Fresnel transmission or detector back-action.
"""

import os

import numpy as np

from MomentRT import (
    CartesianMesh,
    Diffusion,
    Material,
    RobinBoundary,
    Simulation,
    StructuredMedium,
    SurfaceDetector,
    SurfaceSource,
    ureg,
)

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Define a monochromatic phantom and incoming diffuse spot
# ------------------------------------------------------------
# The optical coefficients are illustrative at 800 nm. Positive z points
# into the sample, so reflected light exits toward negative z.

material = Material(
    absorption=0.02 / ureg.mm,
    reduced_scattering=1.0 / ureg.mm,
    refractive_index=1.4,
)

sample = StructuredMedium(
    origin=np.array([-4, -3, 0]) * ureg.mm,
    size=np.array([8, 6, 4]) * ureg.mm,
)

sample.add_background(material=material)

source = SurfaceSource(
    face="z_min",
    position=np.array([-1.5, 0, 0]) * ureg.mm,
    radius=0.8 * ureg.mm,
    wavelength=800 * ureg.nm,
    optical_power=1 * ureg.mW,
)

mesh = CartesianMesh.from_sample(
    sample=sample,
    shape=(17, 13, 13),
)

boundary = RobinBoundary(effective_reflectance=0.4)

solver = Diffusion(
    boundaries={
        face: boundary
        for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
    },
)

# %%
# Specify aperture geometry and angular acceptance
# ------------------------------------------------
# collection_axis points along accepted ray propagation, out of the sample.
# The radius describes a disk in the face plane. acceptance_angle is the
# half-angle around that axis. A normal 30-degree cone collects one quarter
# of Lambertian escaping power from its sampled aperture.

reflection_detector = SurfaceDetector(
    face="z_min",
    position=np.array([1.5, 0, 0]) * ureg.mm,
    radius=0.75 * ureg.mm,
    collection_axis=np.array([0, 0, -1]),
    acceptance_angle=30 * ureg.degree,
)

transmission_detector = SurfaceDetector(
    face="z_max",
    position=np.array([0, 0, 4]) * ureg.mm,
    radius=0.9 * ureg.mm,
    collection_axis=np.array([0, 0, 1]),
    acceptance_angle=60 * ureg.degree,
)

simulation = Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
    detectors={
        "reflection": reflection_detector,
        "transmission": transmission_detector,
    },
)

system_figure = simulation.plot_system(
    show=show_figures,
)

# %%
# Solve, inspect power, and check aperture sampling
# -------------------------------------------------
# Detector power is a passive readout of the solved boundary field. It does
# not reduce the escaping power in the energy balance. The sampled aperture
# area should approach the geometric area under mesh refinement.

result = simulation.run()

for name, measurement in result.detectors.items():
    print(f"{name} collected power: {measurement.collected_power.to('mW'):~P}")

    print(f"{name} sampled area: {measurement.sampled_area.to('mm^2'):~P}")

    print(f"{name} geometric area: {measurement.geometric_area.to('mm^2'):~P}")

fluence_figure = result.fluence.plot(
    view="volume",
    scale="log",
    show=show_figures,
)
