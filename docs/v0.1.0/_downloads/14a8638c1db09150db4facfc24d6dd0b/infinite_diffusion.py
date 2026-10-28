"""
Infinite-medium diffusion
=========================

All physical values are illustrative. The plotted box is an observation
window, not a tissue boundary. Optical properties apply at 970 nm.
"""

import numpy as np

from MomentRT import (
    CartesianMesh,
    HomogeneousMedium,
    InfiniteDiffusion,
    Material,
    PointSource,
    Simulation,
    ureg,
)

# %%
# Optical properties
# ------------------
# Define illustrative absorption and reduced scattering at 970 nm.

material = Material(
    absorption=0.02 / ureg.millimeter,
    reduced_scattering=1 / ureg.millimeter,
    refractive_index=1.4,
)

medium = HomogeneousMedium(
    material=material,
)

# %%
# Observation window
# ------------------
# This grid samples an infinite medium; its edges are not physical boundaries.
# The even node counts avoid sampling the singular point source.

mesh = CartesianMesh(
    shape=(50, 50, 50),
    spacing=np.array([0.4, 0.4, 0.4]) * ureg.millimeter,
    origin=np.array([-9.8, -9.8, -9.8]) * ureg.millimeter,
)

# %%
# Embedded illumination
# ---------------------
# Place an isotropic point source at the center of the observation window.

source = PointSource(
    position=np.zeros(3) * ureg.millimeter,
    wavelength=970 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

# %%
# Configure transport
# -------------------
# Keep the numerical model separate from the physical configuration.

solver = InfiniteDiffusion()

# %%
# Compose the experiment
# ----------------------
# Combine the medium, illumination, observation mesh, and solver.

simulation = Simulation(
    medium=medium,
    source=source,
    mesh=mesh,
    solver=solver,
)

# %%
# Inspect the configuration
# -------------------------
# Rotate the three-dimensional observation window and locate the source.

geometry = simulation.plot_system(
    show=False,
)

# %%
# Calculate fluence
# -----------------
# Solve the analytical infinite-medium diffusion model.

result = simulation.run()

# %%
# Inspect the fluence field
# -------------------------
# The orthogonal slices use a logarithmic color scale with physical values
# retained on hover. Values close to the source remain diffusion approximations.

figure = result.fluence.plot(
    view="slices",
    scale="log",
    show=False,
)
