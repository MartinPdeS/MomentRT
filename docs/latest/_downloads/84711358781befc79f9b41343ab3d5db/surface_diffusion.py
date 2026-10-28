"""
Surface illumination and FEM diffusion
======================================

Illustrative single-wavelength properties. Effective reflectance is supplied
explicitly, not inferred from refractive index. All internal refractive indices
are equal: refracting internal interfaces are outside this solver's scope.
"""

import numpy as np

from MomentRT import (
    CartesianMesh,
    Diffusion,
    Material,
    RobinBoundary,
    Simulation,
    StructuredMedium,
    SurfaceSource,
    ureg,
)

# %%
# Single-wavelength materials
# ---------------------------
# Use equal internal refractive indices: this continuous-fluence FEM model
# does not support refracting internal interfaces.

background = Material(
    absorption=0.02 / ureg.millimeter,
    reduced_scattering=1 / ureg.millimeter,
    refractive_index=1.4,
)

upper_material = Material(
    absorption=0.04 / ureg.millimeter,
    reduced_scattering=1.5 / ureg.millimeter,
    refractive_index=1.4,
)

# %%
# Build the layered sample
# ------------------------
# Create a finite rectangular sample and add a one-millimeter upper layer.

sample = StructuredMedium(
    origin=np.array([-4, -4, 0]) * ureg.millimeter,
    size=np.array([8, 8, 4]) * ureg.millimeter,
)

sample.add_background(
    material=background,
)

sample.add_layer(
    name="Upper layer",
    thickness=1 * ureg.millimeter,
    material=upper_material,
)

# %%
# Discretize the sample
# ---------------------
# Align the layer interface with mesh planes. Each Cartesian cell is split
# into six tetrahedra; mesh refinement is necessary to assess spatial accuracy.

mesh = CartesianMesh.from_sample(
    sample=sample,
    shape=(25, 25, 17),
)

# %%
# Illuminate the surface
# ----------------------
# Define a finite circular diffuse spot on the upper face. Optical power is
# already coupled into the sample, excluding specular entry losses.

source = SurfaceSource(
    face="z_min",
    position=np.array([0, 0, 0]) * ureg.millimeter,
    radius=1 * ureg.millimeter,
    wavelength=970 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

# %%
# Configure boundary conditions
# -----------------------------
# Supply the effective diffuse reflectance explicitly on every face. It is
# not the normal-incidence Fresnel reflectance.

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

# %%
# Compose the simulation
# ----------------------
# Keep sample, mesh, illumination, and solver configuration distinct.

simulation = Simulation(
    medium=sample,
    source=source,
    mesh=mesh,
    solver=solver,
)

# %%
# Inspect the illuminated sample
# ------------------------------
# The red circle shows the illumination footprint on the sample surface.

geometry = simulation.plot_system(
    show=False,
)


# %%
# Solve and check energy balance
# ------------------------------
# Compare incoming, absorbed, and escaping power. Conservation verifies
# numerical accounting, not the accuracy of diffusion against radiative transfer.

result = simulation.run()

print(f"Incoming power: {result.energy_balance['incoming_power']}")

print(f"Absorbed power: {result.energy_balance['absorbed_power'].to('milliwatt')}")

print(f"Escaping power: {result.energy_balance['escaping_power'].to('milliwatt')}")

print(f"Energy residual: {result.energy_balance['residual_power'].to('watt')}")

# %%
# Explore fluence in 3D
# ---------------------
# Inspect orthogonal slices. The logarithmic colors retain physical hover
# values; neither negative values nor zeros are silently clipped.

figure = result.fluence.plot(
    view="slices",
    scale="log",
    show=False,
)


# %%
# Inspect diffusive flux
# ----------------------
# Arrows show element-centered flux without smoothing across interfaces.

flux_figure = result.flux.plot(
    show=False,
)
