"""
Create a Cartesian mesh from node counts
========================================

Build a regular finite-element grid by choosing the node count along each
sample axis. The plot shows the grid inside the sample bounds.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/cartesian_shape.png'
import MomentRT
from MomentRT import CartesianMesh, Material, StructuredMedium, ureg

# %%
# Define a finite sample
# ----------------------

material = Material(
    absorption=0.02 / ureg.millimeter,
    reduced_scattering=1.0 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.8,
    phase_function=MomentRT.HenyeyGreenstein(),
)

sample = StructuredMedium(
    origin=(-2, -2, 0) * ureg.millimeter,
    size=(4, 4, 3) * ureg.millimeter,
)
sample.add_background(material=material)

# %%
# Create and plot the regular grid
# --------------------------------

mesh = CartesianMesh.from_sample(
    sample=sample,
    shape=(9, 9, 7),
)

figure = mesh.plot(
    length_unit="mm",
    show=False,
)

print(f"Cartesian mesh node counts: {mesh.shape}")
print(f"Cartesian mesh spacing: {mesh.spacing.to('mm')}")
