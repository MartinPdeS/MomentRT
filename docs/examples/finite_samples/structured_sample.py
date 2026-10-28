"""
A heterogeneous layered sample
==============================

Build a synthetic optical phantom with three layers, two spherical inclusions,
and a box-shaped heterogeneity. Its materials are illustrative and assumed at
one wavelength. Positive z points into the sample. This example explores
geometry only; a boundary-aware transport solve is not performed here.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from MomentRT import Box, Material, Sphere, StructuredMedium, ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

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

papillary_material = Material(
    absorption=0.04 / ureg.millimeter,
    reduced_scattering=1.3 / ureg.millimeter,
    refractive_index=1.4,
)

reticular_material = Material(
    absorption=0.02 / ureg.millimeter,
    reduced_scattering=0.9 / ureg.millimeter,
    refractive_index=1.4,
)

inclusion_material = Material(
    absorption=0.18 / ureg.millimeter,
    reduced_scattering=1.6 / ureg.millimeter,
    refractive_index=1.4,
)

second_inclusion_material = Material(
    absorption=0.08 / ureg.millimeter,
    reduced_scattering=2.0 / ureg.millimeter,
    refractive_index=1.4,
)

box_material = Material(
    absorption=0.06 / ureg.millimeter,
    reduced_scattering=0.6 / ureg.millimeter,
    refractive_index=1.4,
)

# %%
# Finite sample extent
# --------------------
# Set the lower corner and full dimensions. Positive z points into the sample.

sample = StructuredMedium(
    origin=np.array([-5, -5, 0]) * ureg.millimeter,
    size=np.array([10, 10, 5]) * ureg.millimeter,
)

sample.add_background(
    material=background,
)

# %%
# Stack background layers
# -----------------------
# Each layer starts below the preceding one. The remaining depth keeps the
# explicitly configured background material.

sample.add_layer(
    name="Epidermis",
    thickness=0.2 * ureg.millimeter,
    material=epidermal_material,
)

sample.add_layer(
    name="Papillary layer",
    thickness=0.8 * ureg.millimeter,
    material=papillary_material,
)

# %%
# Add a deeper layer and three localized regions
# ----------------------------------------------
# The two spheres overlap layer interfaces. The box intersects the y=0 section.
# Batch insertion makes the region replacement order explicit.

sample.add_layer(
    name="Reticular layer",
    thickness=2.0 * ureg.millimeter,
    material=reticular_material,
)

absorbing_inclusion = Sphere(
    name="Absorbing inclusion",
    center=np.array([-1.3, 0, 1.8]) * ureg.millimeter,
    radius=0.65 * ureg.millimeter,
    material=inclusion_material,
)

scattering_inclusion = Sphere(
    name="Scattering inclusion",
    center=np.array([1.5, 0, 2.5]) * ureg.millimeter,
    radius=0.5 * ureg.millimeter,
    material=second_inclusion_material,
)

rectangular_region = Box(
    name="Rectangular heterogeneity",
    origin=np.array([-4, -1.25, 3.3]) * ureg.millimeter,
    size=np.array([1.5, 2.5, 0.9]) * ureg.millimeter,
    material=box_material,
)

sample.add_structures(
    absorbing_inclusion,
    scattering_inclusion,
    rectangular_region,
)

# %%
# Explore the sample in 3D
# ------------------------
# Rotate the geometry and toggle regions in the legend. Transparent surfaces
# show primitive boundaries; the section below resolves overlaps.

figure = sample.plot(
    show=show_figures,
)


# %%
# Inspect the resolved cross-section
# ----------------------------------
# Unlike transparent construction surfaces, this map resolves material
# replacement. The y=0 plane intersects both spheres and the box.

section_figure = sample.plot_section(
    normal="y",
    position=0 * ureg.millimeter,
    show=show_figures,
)

# %%
# Show a compact resolved-material overview
# -----------------------------------------
# The static section highlights insertion precedence. The interactive sample
# and section views above remain available for geometry and property queries.

section = sample.section(
    normal="y",
    position=0 * ureg.millimeter,
    resolution=(160, 120),
)

extent = section.extent.to("millimeter").magnitude

figure_static, axes = plt.subplots(figsize=(7, 3.5))

image = axes.imshow(
    section.labels,
    origin="lower",
    extent=(extent[0, 0], extent[0, 1], extent[1, 0], extent[1, 1]),
    interpolation="nearest",
    cmap=plt.get_cmap("Set2", len(section.names)),
    vmin=-0.5,
    vmax=len(section.names) - 0.5,
    aspect="equal",
)

colorbar = figure_static.colorbar(image, ax=axes, ticks=range(len(section.names)))

colorbar.ax.set_yticklabels(section.names)

axes.set_xlabel("x (mm)")
axes.set_ylabel("z (mm)")
axes.set_title("Resolved materials on the y = 0 plane")
figure_static.tight_layout()

if show_figures:
    plt.show()
