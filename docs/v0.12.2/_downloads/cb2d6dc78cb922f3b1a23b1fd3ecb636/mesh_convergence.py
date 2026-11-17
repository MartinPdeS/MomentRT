"""
Check detector mesh convergence
===============================

Study how refining a tetrahedral mesh changes detector powers for a layered
sample at one wavelength. This is a numerical check to perform before an
inverse fit; it does not establish the accuracy of the P1 model for skin.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/mesh_convergence.png'
import os
import numpy as np
import MomentRT

ureg = MomentRT.ureg

# %%
# Configure the sample and instrument
# -----------------------------------
# The layer thickness and optical properties are illustrative. The epidermis
# preset's melanin fraction is a model parameter, not a clinical measurement.
# Both materials have the same refractive index because the solver does not
# model internal refractive-index jumps.

wavelength = 650 * ureg.nm

background = MomentRT.Material(
    absorption=0.025 / ureg.mm,
    reduced_scattering=1.2 / ureg.mm,
    refractive_index=1.4,
    phase_function=MomentRT.Isotropic(),
)

epidermis = MomentRT.Epidermis(
    melanin_fraction=0.01,
    wavelength=wavelength,
    refractive_index=1.4,
)

sample = MomentRT.StructuredMedium(
    origin=(-3, -3, 0) * ureg.mm,
    size=(6, 6, 3) * ureg.mm,
)

sample.add_background(material=background)

sample.add_layer(name="Epidermis", thickness=0.5 * ureg.mm, material=epidermis)

source = MomentRT.SurfaceSource(
    face="z_min",
    position=np.array([-1.5, 0, 0]) * ureg.mm,
    radius=0.8 * ureg.mm,
    wavelength=wavelength,
    optical_power=1 * ureg.mW,
)

detectors = {
    name: MomentRT.SurfaceDetector(
        face="z_min",
        position=np.array([position, 0, 0]) * ureg.mm,
        radius=0.75 * ureg.mm,
        collection_axis=[0, 0, -1],
        acceptance_angle=90 * ureg.degree,
    )
    for name, position in (("near", 0.5), ("far", 1.5))
}

boundary = MomentRT.RobinBoundary(effective_reflectance=0.4)

solver = MomentRT.P1(
    boundaries={face: boundary for face in MomentRT.boundary.FACES},
)

mesh = MomentRT.TetrahedralMesh.from_sample(
    sample=sample,
    resolution=0.8 * ureg.mm,
)

simulation = MomentRT.Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
    detectors=detectors,
)

# %%
# Compare successively refined meshes
# ------------------------------------
# The initial mesh belongs to the simulation; the study creates independent
# meshes at the explicit target resolutions below. It compares detector powers
# and sampled source and aperture areas. Two consecutive refinements must meet
# the tolerances at the fine end, and negative outputs remain visible.
# The plot shades the power tolerance around the finest detector values and
# area tolerances around the geometric footprint areas. These bands are not
# sufficient on their own: the study also checks adjacent refinements.
# The source power is incoming diffuse light already coupled into the sample;
# detector power uses the modeled Lambertian angular approximation. The
# progress table prints each completed mesh level as the study runs.

study = MomentRT.MeshConvergence(
    simulation=simulation,
    mesh_resolutions=[
        0.8 * ureg.mm,
        0.4 * ureg.mm,
        0.3 * ureg.mm,
        0.2 * ureg.mm,
        0.15 * ureg.mm,
        0.10 * ureg.mm,
    ],
    power_rtol=0.05,
    area_rtol=0.05,
    required_consecutive=2,
    progress=True,
)

result = study.run()

convergence_figure = result.plot(
    power_unit="uW",
    area_unit="mm^2",
    show=False,
)

# %%
# Inspect the finest solution
# ---------------------------
# Only recommend a mesh when the final sequence meets the stated criteria.
# The fluence volume below is
# interpolated onto a separate display grid and is not a convergence test.
# Transparent layer surfaces show construction geometry, not a material map.
# A convergence result at this wavelength and baseline does not establish
# convergence for other wavelengths or tissue properties, or P1 accuracy.

finest = result.results[-1]

if finest is None:
    raise RuntimeError("The finest source footprint was not resolved.")

background_layer = MomentRT.Layer(
    name="Underlying background",
    depth=0.5 * ureg.mm,
    thickness=2.5 * ureg.mm,
    material=background,
)

gallery_system_figure = finest.fluence.plot(
    view="volume",
    scale="linear",
    display_resolution=0.5 * ureg.mm,
    overlays=[*sample.structures, background_layer, source, *detectors.values()],
    show=False,
)

# %%
# Show the figures when running locally
# -------------------------------------
# Sphinx-Gallery embeds both figures without opening browser windows.

if os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1":
    convergence_figure.show()
    gallery_system_figure.show()
