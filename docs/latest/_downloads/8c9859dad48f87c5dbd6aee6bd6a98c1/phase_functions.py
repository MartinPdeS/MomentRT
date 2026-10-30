"""
Compare scattering phase functions
==================================

Compare isotropic and unpolarized Rayleigh scattering in the same finite
sample. Both have zero mean scattering cosine, but their second angular
moments differ. All optical properties are illustrative.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/phase_functions.png'

import os
import numpy as np
import MomentRT

ureg = MomentRT.ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Choose phase functions
# ----------------------
# Each model returns density per steradian and normalized Legendre moments.
# Isotropic and Rayleigh both have g₁ = 0, while g₂ is 0 and 0.1 respectively.
# The existing Henyey–Greenstein model instead takes anisotropy from its
# material; tabulated functions can be supplied on a cosine grid.

phases = {
    "Isotropic": MomentRT.phase_function.Isotropic(),
    "Rayleigh": MomentRT.phase_function.Rayleigh(),
}

for name, phase in phases.items():
    print(f"{name} moments: {phase.moments()}")

# %%
# Configure the common sample and instrument
# ------------------------------------------
# The sample is solved twice with equal absorption and reduced scattering.
# SP3 uses the selected phase function's higher moments. P1 only uses the
# reduced scattering coefficient and would therefore give the same result.

sample = MomentRT.StructuredMedium(
    origin=(-4, -4, 0) * ureg.millimeter,
    size=(8, 8, 4) * ureg.millimeter,
)

source = MomentRT.SurfaceSource.on_face(
    sample=sample,
    face="z_min",
    coordinates=np.array([0, 0]) * ureg.millimeter,
    radius=1.5 * ureg.millimeter,
    wavelength=850 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

detector = MomentRT.SurfaceDetector.on_face(
    sample=sample,
    face="z_min",
    coordinates=np.array([2.5, 0]) * ureg.millimeter,
    radius=0.8 * ureg.millimeter,
    acceptance_angle=90 * ureg.degree,
    orientation="outward_normal",
)

sample.add_background(
    material=MomentRT.Material(
        absorption=0.03 / ureg.millimeter,
        reduced_scattering=1.0 / ureg.millimeter,
        refractive_index=1.4,
        phase_function=phases["Isotropic"],
    )
)

mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=(11, 11, 9))

boundary = MomentRT.RobinBoundary(effective_reflectance=0)

boundaries = MomentRT.Boundaries.uniform(condition=boundary)

solver = MomentRT.SP3(boundaries=boundaries)

results = {}

for name, phase in phases.items():
    material = MomentRT.Material(
        absorption=0.03 / ureg.millimeter,
        reduced_scattering=1.0 / ureg.millimeter,
        refractive_index=1.4,
        phase_function=phase,
    )

    sample.add_background(material=material)

    simulation = MomentRT.Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=solver,
        detectors={"contact": detector},
    )

    result = simulation.run()

    results[name] = result

    print(f"{name} collected power: {result.detector_power['contact'].to('mW'):~P}")

# %%
# Visualize the Rayleigh solution
# -------------------------------
# This volume displays the computed fluence, not the single-scattering phase
# function. Source and detector overlays use their physical contact positions.
# Verify mesh convergence and compare with a transport reference for the
# intended application; the contact detector assumes Lambertian escape.

fluence_figure = results["Rayleigh"].fluence.plot(
    view="volume",
    scale="linear",
    overlays=[source, detector],
    show=show_figures,
)

# %%
# Preview the sample and instrument geometry
# ------------------------------------------
# Keep the physical geometry as the final figure so the gallery thumbnail
# shows the sample, source, and detector rather than the phase-function result.

gallery_system_figure = sample.plot(
    overlays=[source, detector],
    show=show_figures,
)
