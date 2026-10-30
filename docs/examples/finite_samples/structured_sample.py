"""
A heterogeneous layered sample
==============================

Build a synthetic optical phantom with three layers, two spherical inclusions,
and a box-shaped heterogeneity. Its materials are illustrative and assumed at
one wavelength. Positive z points into the sample. Then illuminate the sample
and inspect the finite-sample P1 fluence volume.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/structured_sample.png'
import os
import numpy as np
import MomentRT
from MomentRT import (
    Box,
    CartesianMesh,
    Diffusion,
    Material,
    RobinBoundary,
    Simulation,
    Sphere,
    StructuredMedium,
    UniformSurfaceSource,
    ureg,
)

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Optical materials
# -----------------
# Use explicit unit-bearing optical coefficients. These values are illustrative.

background = Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    absorption=0.01 / ureg.millimeter,
    reduced_scattering=0.8 / ureg.millimeter,
    refractive_index=1.4,
)

epidermal_material = Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    absorption=0.1 / ureg.millimeter,
    reduced_scattering=1.5 / ureg.millimeter,
    refractive_index=1.4,
)

papillary_material = Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    absorption=0.04 / ureg.millimeter,
    reduced_scattering=1.3 / ureg.millimeter,
    refractive_index=1.4,
)

reticular_material = Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    absorption=0.02 / ureg.millimeter,
    reduced_scattering=0.9 / ureg.millimeter,
    refractive_index=1.4,
)

inclusion_material = Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    absorption=0.18 / ureg.millimeter,
    reduced_scattering=1.6 / ureg.millimeter,
    refractive_index=1.4,
)

second_inclusion_material = Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
    absorption=0.08 / ureg.millimeter,
    reduced_scattering=2.0 / ureg.millimeter,
    refractive_index=1.4,
)

box_material = Material(
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
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
# Configure a finite-sample solve
# ---------------------------------
# Use equal refractive indices across all regions because internal refraction
# at refractive-index jumps is outside the current continuous-fluence FEM model.

mesh = CartesianMesh.from_sample(
    sample=sample,
    shape=(21, 21, 21),
)

source = UniformSurfaceSource(
    face="z_min",
    wavelength=970 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

boundary = RobinBoundary(
    effective_reflectance=0.5,
)

solver = Diffusion(
    boundaries={
        "x_min": boundary,
        "x_max": boundary,
        "y_min": boundary,
        "y_max": boundary,
        "z_min": boundary,
        "z_max": boundary,
    },
)

simulation = Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
)

# %%
# Solve and inspect fluence in 3D
# ----------------------------------
# The logarithmic volume displays the nodal fluence rate throughout the sample.
# The mesh gives a useful visualization; refine it to assess interface accuracy.

result = simulation.run()

fluence_figure = result.fluence.plot(
    view="volume",
    scale="log",
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
# Inspect the modeled geometry
# ----------------------------
# This view shows the configured geometry and instrument footprints.

gallery_system_figure = simulation.plot_system(show=False)
