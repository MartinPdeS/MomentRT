"""
Sweep wavelength-dependent scattering anisotropy
==================================================

Specify different spectral anisotropy tables for two layers and evaluate an
SP3 sensor at several wavelengths. The numbers are illustrative. Each result
is monochromatic; the sweep is not a broadband integration.
"""

# sphinx_gallery_thumbnail_path = '_static/images/sp3_spectral_anisotropy.png'

import os

import numpy as np

import MomentRT

ureg = MomentRT.ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Give each material its own spectral table
# -----------------------------------------
# Anisotropy can be a plain scalar or a vector matching the wavelength grid.
# MomentRT linearly interpolates the vector without extrapolation. SP3 uses
# the interpolated value with a Henyey–Greenstein phase function.

table_wavelengths = np.array([700, 800, 900]) * ureg.nanometer

dermis = MomentRT.SpectralMaterial(
    wavelengths=table_wavelengths,
    absorption=np.array([0.035, 0.025, 0.03]) / ureg.millimeter,
    reduced_scattering=np.array([1.3, 1.1, 0.9]) / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=np.array([0.75, 0.80, 0.84]),
)

epidermis = MomentRT.SpectralMaterial(
    wavelengths=table_wavelengths,
    absorption=np.array([0.06, 0.045, 0.05]) / ureg.millimeter,
    reduced_scattering=np.array([1.8, 1.6, 1.4]) / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=np.array([0.60, 0.65, 0.70]),
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
# Keep the instrument geometry fixed
# ----------------------------------
# Every wavelength uses the same source footprint and detector aperture. The
# zero-reflectance faces select the matched-index SP3 boundary condition.

source = MomentRT.SurfaceSource(
    face="z_min",
    position=np.array([-1.5, 0, 0]) * ureg.millimeter,
    radius=1.0 * ureg.millimeter,
    wavelength=800 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

detector = MomentRT.SurfaceDetector(
    face="z_min",
    position=np.array([1.5, 0, 0]) * ureg.millimeter,
    radius=0.8 * ureg.millimeter,
    collection_axis=[0, 0, -1],
    acceptance_angle=90 * ureg.degree,
)

mesh = MomentRT.CartesianMesh.from_sample(
    sample=sample,
    shape=(15, 15, 11),
)

boundary = MomentRT.RobinBoundary(effective_reflectance=0)

boundaries = {face: boundary for face in MomentRT.boundary.FACES}

solver = MomentRT.SP3(boundaries=boundaries)

simulation = MomentRT.Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
    detectors={"contact": detector},
)

# %%
# Interpolate anisotropy and solve each wavelength
# -------------------------------------------------
# The 750 nm value lies halfway between the first two tabulated values.
# Geometry and boundary quadrature are prepared once for the sweep.

wavelengths = np.array([700, 750, 800, 850, 900]) * ureg.nanometer

spectrum = simulation.sweep(wavelengths=wavelengths)

for wavelength, power in zip(wavelengths, spectrum.detector_power["contact"]):
    dermis_at_wavelength = dermis.at(wavelength=wavelength)

    epidermis_at_wavelength = epidermis.at(wavelength=wavelength)

    print(
        f"{wavelength:~P}: dermis g={dermis_at_wavelength.anisotropy:.3f}, "
        f"epidermis g={epidermis_at_wavelength.anisotropy:.3f}, "
        f"collected power={power.to('mW'):~P}"
    )

spectral_figure = spectrum.plot(
    wavelength_unit="nm",
    power_unit="mW",
    show=show_figures,
)

# %%
# View a monochromatic fluence volume
# ------------------------------------
# The 800 nm solution uses the third entry in the sweep. The detector uses a
# separate Lambertian angular-acceptance approximation. Check mesh convergence
# and compare with a transport reference for the intended application.

fluence_figure = spectrum.results[2].fluence.plot(
    view="volume",
    scale="linear",
    overlays={"Source": source, "Contact detector": detector},
    show=show_figures,
)
