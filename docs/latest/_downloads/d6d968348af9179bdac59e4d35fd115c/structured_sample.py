"""
A layered sample with an inclusion
==================================

Illustrative properties only, assumed at a single wavelength. Positive z points
into the sample. Geometry is ready to inspect and sample; a boundary-aware
transport solver is still required to propagate light through this sample.
"""

import numpy as np

from MomentRT import Material, Sphere, StructuredMedium, ureg

# %%
# Optical materials
# -----------------
# Use explicit unit-bearing optical coefficients. These values are illustrative.

background = Material(
    absorption=0.01 / ureg.millimeter,
    reduced_scattering=0.8 / ureg.millimeter,
    refractive_index=1.4,
)

epidermal_material = Material(
    absorption=0.1 / ureg.millimeter,
    reduced_scattering=1.5 / ureg.millimeter,
    refractive_index=1.42,
)

dermal_material = Material(
    absorption=0.03 / ureg.millimeter,
    reduced_scattering=1.0 / ureg.millimeter,
    refractive_index=1.4,
)

inclusion_material = Material(
    absorption=0.2 / ureg.millimeter,
    reduced_scattering=1.2 / ureg.millimeter,
    refractive_index=1.4,
)

# %%
# Finite sample extent
# --------------------
# Set the lower corner and full dimensions. Positive z points into the sample.

sample = StructuredMedium(
    origin=np.array([-5, -5, 0]) * ureg.millimeter,
    size=np.array([10, 10, 4]) * ureg.millimeter,
)

sample.add_background(
    material=background,
)

# %%
# Stack tissue layers
# -------------------
# Each layer starts below the preceding one, without manually calculating depth.

sample.add_layer(
    name="Epidermis",
    thickness=0.2 * ureg.millimeter,
    material=epidermal_material,
)

sample.add_layer(
    name="Dermis",
    thickness=1.8 * ureg.millimeter,
    material=dermal_material,
)

# %%
# Add an inclusion
# ----------------
# Add the inclusion after the layers so it replaces their material in overlaps.

inclusion = Sphere(
    name="Absorbing inclusion",
    center=np.array([1, 0, 1.2]) * ureg.millimeter,
    radius=0.5 * ureg.millimeter,
    material=inclusion_material,
)

sample.add_structures(inclusion)

# %%
# Explore the sample in 3D
# ------------------------
# Rotate the geometry and toggle regions in the legend. Transparent surfaces
# show construction geometry; hover labels report the optical properties.

figure = sample.plot(
    show=False,
)


# %%
# Inspect the resolved cross-section
# ----------------------------------
# Unlike transparent construction surfaces, this map resolves material
# replacement. The y=0 plane intersects the inclusion.

section_figure = sample.plot_section(
    normal="y",
    position=0 * ureg.millimeter,
    show=False,
)
