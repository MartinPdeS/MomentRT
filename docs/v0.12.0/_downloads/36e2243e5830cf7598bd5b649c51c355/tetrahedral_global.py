"""
Create a uniform tetrahedral mesh
=================================

Generate an unstructured tetrahedral discretization with a global target
element size, then plot the sample boundary and tetrahedral edges.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/tetrahedral_global.png'
import MomentRT
from MomentRT import Material, StructuredMedium, TetrahedralMesh, ureg

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
    origin=(0, 0, 0) * ureg.millimeter,
    size=(2, 2, 1) * ureg.millimeter,
)
sample.add_background(material=material)

# %%
# Generate and plot tetrahedra
# ----------------------------

mesh = TetrahedralMesh.from_sample(
    sample=sample,
    resolution=0.5 * ureg.millimeter,
)

figure = mesh.plot(
    length_unit="mm",
    show=False,
)

print(f"Tetrahedral mesh nodes: {len(mesh.points)}")
print(f"Tetrahedral elements: {len(mesh.tetrahedra)}")
