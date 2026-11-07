"""
Build samples with named models
================================

Construct a two-layer skin sample and a general layered phantom with explicit
materials. The models return ordinary ``StructuredMedium`` objects, so they
work with the same meshes, sources, and solvers as a manually built sample.
Optical values are illustrative and do not represent an individual.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/sample_models.png'
import numpy as np
import MomentRT
from MomentRT.models.phantom import LayeredPhantom
from MomentRT.models.skin import TwoLayerSkin

ureg = MomentRT.ureg

# %%
# Choose materials and geometry
# -----------------------------
# The parameterized epidermis and dermis materials use approximate literature
# fits at one wavelength. Their refractive indices match because the finite
# solvers do not currently model refractive-index jumps within the sample.

epidermis = MomentRT.Epidermis(
    melanin_fraction=0.05,
    wavelength=850 * ureg.nm,
    refractive_index=1.4,
)
dermis = MomentRT.Dermis(
    hydration=0.7,
    wavelength=850 * ureg.nm,
    refractive_index=1.4,
)
origin = (-3, -3, 0) * ureg.mm
size = (6, 6, 3) * ureg.mm

# %%
# Build a two-layer sample
# ------------------------
# The epidermis starts at the sample's upper face. The dermis fills the
# remaining depth and is also the background. Positive z points into the sample.

sample = TwoLayerSkin(
    origin=origin,
    size=size,
    epidermis=epidermis,
    dermis=dermis,
    epidermis_thickness=0.5 * ureg.mm,
)
print(f"Skin layers: {[layer.name for layer in sample.structures]}")

# %%
# Build a general layered phantom
# -------------------------------
# A phantom accepts a sequence of named layers, thicknesses, and materials.
# The explicit background fills the space below the stack. Materials and
# thicknesses here are chosen to match the skin geometry above; they do not
# describe a manufactured phantom.

phantom = LayeredPhantom(
    origin=origin,
    size=size,
    background=dermis,
    layers=(("Upper layer", 0.5 * ureg.mm, epidermis),),
)
print(f"Phantom layers: {[layer.name for layer in phantom.structures]}")

# %%
# Solve the skin sample
# ---------------------
# The source is diffuse power already coupled into the sample. This coarse
# mesh is illustrative: resolving the layer and source footprint, and the
# validity of P1 at the relevant distances, require separate checks.

mesh = MomentRT.CartesianMesh.from_sample(
    sample=sample,
    shape=(13, 13, 13),
)
source = MomentRT.SurfaceSource(
    face="z_min",
    position=np.zeros(3) * ureg.mm,
    radius=1 * ureg.mm,
    wavelength=850 * ureg.nm,
    optical_power=1 * ureg.mW,
)
solver = MomentRT.P1(
    boundaries=MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0),
    ),
)
simulation = MomentRT.Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
)
result = simulation.run()
print(f"Sampled source area: {result.source_sampled_area.to('mm^2'):~P}")
print(f"Geometric source area: {result.source_geometric_area.to('mm^2'):~P}")

gallery_system_figure = result.fluence.plot(
    view="volume",
    scale="linear",
    overlays=[*sample.structures, source],
    show=False,
)
