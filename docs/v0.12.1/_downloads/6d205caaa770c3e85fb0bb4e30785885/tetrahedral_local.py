"""
Refine a tetrahedral mesh near layers and contacts
==================================================

Create a layered sample, refine the epidermis, and add finer elements near a
source and detector on the surface. Plot the resulting unstructured mesh.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/tetrahedral_local.png'
# Limit displayed edges to keep the local size transition visible.
import MomentRT
from MomentRT import (
    Material,
    StructuredMedium,
    SurfaceDetector,
    SurfaceSource,
    TetrahedralMesh,
    ureg,
)

# %%
# Configure a layered optical sample
# ----------------------------------

epidermis_material = Material(
    absorption=0.1 / ureg.millimeter,
    reduced_scattering=1.5 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.85,
    phase_function=MomentRT.HenyeyGreenstein(),
)

dermis_material = Material(
    absorption=0.03 / ureg.millimeter,
    reduced_scattering=1.0 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.8,
    phase_function=MomentRT.HenyeyGreenstein(),
)

sample = StructuredMedium(
    origin=(0, 0, 0) * ureg.millimeter,
    size=(2, 2, 1) * ureg.millimeter,
)
sample.add_background(material=dermis_material)
sample.add_layer(
    name="Epidermis",
    material=epidermis_material,
    thickness=0.2 * ureg.millimeter,
)

# %%
# Place surface contacts and generate a base mesh
# -----------------------------------------------

source = SurfaceSource.on_face(
    sample=sample,
    face="z_min",
    coordinates=(0.7, 1.0) * ureg.millimeter,
    radius=0.15 * ureg.millimeter,
    wavelength=850 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

detector = SurfaceDetector.on_face(
    sample=sample,
    face="z_min",
    coordinates=(1.3, 1.0) * ureg.millimeter,
    radius=0.15 * ureg.millimeter,
    orientation="outward_normal",
    acceptance_angle=30 * ureg.degree,
)

mesh = TetrahedralMesh.from_sample(
    sample=sample,
    resolution=0.5 * ureg.millimeter,
)

# %%
# Add layer and contact refinement
# --------------------------------

mesh.refine_region(
    name="Epidermis",
    resolution=0.1 * ureg.millimeter,
)

mesh.refine_near(
    objects=[source, detector],
    resolution=0.08 * ureg.millimeter,
    distance=0.2 * ureg.millimeter,
)

figure = mesh.plot(
    length_unit="mm",
    max_edges=4000,
    show=False,
)

print(f"Locally refined mesh nodes: {len(mesh.points)}")
print(f"Locally refined tetrahedra: {len(mesh.tetrahedra)}")
