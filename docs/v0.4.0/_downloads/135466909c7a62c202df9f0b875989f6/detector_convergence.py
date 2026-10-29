"""
Choose a mesh for detector collection
=====================================

Compare collected power, sampled detector area, and sampled illumination area
over four explicit meshes. The study preserves the physical sample and
instrument geometry. Its recommendation is a numerical resolution check for
these observables, not validation of P1 against radiative transfer.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_volumes/detector_convergence.png'

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
# Fix the sample and optical configuration
# ----------------------------------------
# These coefficients describe an illustrative homogeneous phantom at 800 nm.
# The disk source and the separated reflection aperture lie on the top face.

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

detector = SurfaceDetector(
    face="z_min",
    position=np.array([1.5, 0, 0]) * ureg.mm,
    radius=0.75 * ureg.mm,
    collection_axis=np.array([0, 0, -1]),
    acceptance_angle=45 * ureg.degree,
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

simulation = Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
    detectors={"reflection": detector},
)

# %%
# Compare explicit mesh levels
# ----------------------------
# Power must change by at most 5 percent between successive levels. Both
# sampled footprint areas must be within 5 percent of their geometric areas.
# Two consecutive passing refinements are required at the fine end of the
# study; a recommended level must also agree with the finest supplied power.
# Each detector legend entry controls both its power and area curves.
# Node-count ticks identify the exact tested mesh levels.
# A missing recommendation means more refinement or different tolerances
# need investigation, rather than automatic acceptance of the last mesh.

study = simulation.convergence_study(
    mesh_shapes=[
        (9, 7, 7),
        (13, 10, 10),
        (17, 13, 13),
        (25, 19, 19),
    ],
    power_rtol=0.05,
    area_rtol=0.05,
    required_consecutive=2,
)

convergence_figure = study.plot(
    power_unit="mW",
    area_unit="mm^2",
    show=show_figures,
)

for record in study.records:
    measurement = record["detectors"]["reflection"]

    print(f"Mesh shape: {record['mesh_shape']}")

    print(f"Reflection power: {measurement['collected_power'].to('mW'):~P}")

    print(f"Detector area error: {measurement['area_relative_error']:.2%}")

    print(f"Source area error: {record['source_area_relative_error']:.2%}")

print(f"Recommended shape: {study.recommended_shape}")

# %%
# Inspect the finest field with the same instrument overlays
# ----------------------------------------------------------
# The study retains each result. A source footprint with no quadrature points
# is recorded with result=None; this example's meshes resolve the source.
# The original simulation.mesh remains unchanged.

finest_result = study.results[-1]

fluence_figure = finest_result.fluence.plot(
    view="volume",
    scale="log",
    overlays={"Illumination": source, "Reflection detector": detector},
    show=show_figures,
)
