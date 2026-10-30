"""
Create a Cartesian mesh from physical resolution
================================================

Choose grid spacing in length units instead of selecting node counts. A
three-component resolution can make the spacing different along each axis.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/cartesian_resolution.png'
import MomentRT
from MomentRT import CartesianMesh, Material, StructuredMedium, ureg

# %%
# Define the sample and optical material
# --------------------------------------

material = Material(
    absorption=0.02 / ureg.millimeter,
    reduced_scattering=1.0 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.8,
    phase_function=MomentRT.HenyeyGreenstein(),
)

sample = StructuredMedium(
    origin=(0, 0, 0) * ureg.millimeter,
    size=(4, 4, 3) * ureg.millimeter,
)
sample.add_background(material=material)

# %%
# Request anisotropic physical spacing
# ------------------------------------
# The z spacing is finer to resolve depth variations. MomentRT adjusts the
# actual spacing downward as needed so that the grid reaches every sample face.

mesh = CartesianMesh.from_sample(
    sample=sample,
    resolution=(0.5, 0.5, 0.25) * ureg.millimeter,
)

figure = mesh.plot(
    length_unit="mm",
    show=False,
)

print(f"Requested resolution: {mesh.resolution.to('mm')}")
print(f"Node counts: {mesh.shape}")
print(f"Actual spacing: {mesh.spacing.to('mm')}")
