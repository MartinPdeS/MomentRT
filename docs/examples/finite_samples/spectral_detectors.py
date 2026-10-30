"""
Sweep wavelengths with a fixed sensor geometry
==============================================

Assign tabulated spectra to two regions, keep source and detector geometry
fixed, and compare collected power at five wavelengths. All coefficients are
illustrative. The sweep runs independent monochromatic P1 solves; it is not
a broadband source integration or a validated tissue spectrum.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/spectral_detectors.png'
import os
import numpy as np
import MomentRT
from MomentRT import (
    CartesianMesh,
    Diffusion,
    RobinBoundary,
    Simulation,
    SpectralMaterial,
    StructuredMedium,
    SurfaceDetector,
    SurfaceSource,
    ureg,
)

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Define tabulated optical properties
# ------------------------------------
# Tables are linearly interpolated. Out-of-range wavelengths are rejected.
# The refractive index is constant across the internal interface at every
# wavelength, as required by the current continuous-fluence solver.

table_wavelengths = np.array([650, 750, 850, 950]) * ureg.nm

background = SpectralMaterial(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    wavelengths=table_wavelengths,
    absorption=np.array([0.025, 0.014, 0.018, 0.035]) / ureg.mm,
    reduced_scattering=np.array([1.5, 1.3, 1.1, 0.9]) / ureg.mm,
    refractive_index=1.4,
)

upper_material = SpectralMaterial(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    wavelengths=table_wavelengths,
    absorption=np.array([0.06, 0.03, 0.025, 0.05]) / ureg.mm,
    reduced_scattering=np.array([1.8, 1.6, 1.4, 1.2]) / ureg.mm,
    refractive_index=1.4,
)

sample = StructuredMedium(
    origin=np.array([-4, -3, 0]) * ureg.mm,
    size=np.array([8, 6, 4]) * ureg.mm,
)

sample.add_background(material=background)

sample.add_layer(
    name="Superficial region",
    thickness=1 * ureg.mm,
    material=upper_material,
)

# %%
# Configure the source, mesh, boundaries, and detectors once
# ----------------------------------------------------------
# Power is already coupled into the sample. The sweep retains that same
# optical power at every wavelength unless a power vector is supplied.
# Boundary effective reflectance and detector acceptance also stay fixed.

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

reflection_detector = SurfaceDetector(
    face="z_min",
    position=np.array([1.5, 0, 0]) * ureg.mm,
    radius=0.75 * ureg.mm,
    collection_axis=np.array([0, 0, -1]),
    acceptance_angle=35 * ureg.degree,
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
# Run the ordered wavelength sweep
# ----------------------------------
# Each result keeps its own fields, wavelength, and detector measurements.
# The original sample, source wavelength, and detector geometry are retained.

wavelengths = np.array([700, 750, 800, 850, 900]) * ureg.nm

spectrum = simulation.sweep(
    wavelengths=wavelengths,
)

print(f"Geometry preparations: {spectrum.diagnostics['geometry_preparations']}")

print(f"Preparation time: {spectrum.diagnostics['preparation_time'].to('s'):~P}")

print(f"Per-wavelength run times: {spectrum.diagnostics['solve_times'].to('s'):~P}")

spectral_figure = spectrum.plot(
    wavelength_unit="nm",
    power_unit="mW",
    show=show_figures,
)

for wavelength, power in zip(wavelengths, spectrum.detector_power["reflection"]):
    print(f"Reflection at {wavelength:~P}: {power.to('mW'):~P}")

# %%
# Inspect one monochromatic fluence volume
# -----------------------------------------
# The third solve corresponds to 800 nm. Its geometry is the same as in the
# other solves; its optical coefficients are interpolated at this wavelength.
# Overlays place the instrument on the spatial field; the spectrum above keeps
# its wavelength and power axes.
# Named overlays give the two apertures distinct legend entries.

fluence_figure = spectrum.results[2].fluence.plot(
    view="volume",
    scale="log",
    overlays={
        "Illumination": source,
        "Reflection detector": reflection_detector,
        "Transmission detector": transmission_detector,
    },
    show=show_figures,
)


# %%
# Inspect the modeled geometry
# ----------------------------
# This view shows the configured geometry and instrument footprints.

gallery_system_figure = simulation.plot_system(show=False)
