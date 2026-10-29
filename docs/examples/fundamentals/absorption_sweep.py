"""
How absorption changes fluence
==============================

Compare analytical P1 fields with the same source and scattering but three
absorption coefficients. Higher absorption shortens the distance over which
fluence reaches a given value. Inputs are illustrative.
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/absorption_sweep.png'
import os

import matplotlib.pyplot as plt
import numpy as np

import MomentRT
from MomentRT import (
    CartesianMesh,
    HomogeneousMedium,
    InfiniteDiffusion,
    Material,
    PointSource,
    Simulation,
    ureg,
)

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Set up the common source and observation window
# -----------------------------------------------
# Use an even node count so that no mesh node lands on the point-source
# singularity at the origin.

mesh = CartesianMesh(
    shape=(24, 24, 24),
    spacing=np.full(3, 0.4) * ureg.millimeter,
    origin=np.full(3, -4.6) * ureg.millimeter,
)

source = PointSource(
    position=np.zeros(3) * ureg.millimeter,
    wavelength=970 * ureg.nanometer,
    optical_power=1 * ureg.milliwatt,
)

solver = InfiniteDiffusion()

# %%
# Solve at several absorption levels
# ----------------------------------
# Keep reduced scattering fixed so that the change between curves is due to
# absorption. Sample an off-axis line so it does not pass through the source.

line_index = mesh.shape[1] // 2 - 1

coordinates = mesh.points.to("millimeter").magnitude

line_positions = coordinates[:, line_index, line_index, :]

radial_distance = np.linalg.norm(line_positions, axis=-1)

absorption_values = (0.005, 0.02, 0.05)

fluence_profiles = []

for absorption_value in absorption_values:
    material = Material(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        absorption=absorption_value / ureg.millimeter,
        reduced_scattering=0.8 / ureg.millimeter,
        refractive_index=1.4,
    )

    medium = HomogeneousMedium(
        material=material,
    )

    simulation = Simulation(
        medium=medium,
        source=source,
        mesh=mesh,
        solver=solver,
    )

    result = simulation.run()

    line_fluence = result.fluence.values.to("watt / meter**2").magnitude

    fluence_profiles.append(line_fluence[:, line_index, line_index])

# %%
# Compare the radial profiles
# ---------------------------
# The nearest plotted sample remains offset from the ideal point source.

order = np.argsort(radial_distance)

figure, axes = plt.subplots(figsize=(8, 5))

for absorption_value, profile in zip(absorption_values, fluence_profiles):
    axes.plot(
        radial_distance[order],
        profile[order],
        label=f"Absorption = {absorption_value:g} mm$^{{-1}}$",
    )

axes.set_yscale("log")
axes.set_xlabel("Distance from source (mm)")
axes.set_ylabel("Fluence rate (W m$^{-2}$)")
axes.set_title("Absorption controls attenuation in infinite-medium P1")
axes.grid(True, which="both", alpha=0.25)
axes.legend()
figure.tight_layout()

if show_figures:
    plt.show()

# %%
# Explore a three-dimensional fluence volume
# ------------------------------------------
# This is the highest-absorption solution above. The mesh is an observation
# window within the infinite medium, not the boundary of a physical sample.

volume_figure = result.fluence.plot(
    view="volume",
    scale="log",
    overlays=[source],
    show=show_figures,
)


# %%
# Inspect the modeled geometry
# ----------------------------
# This view shows the configured geometry and instrument footprints.

gallery_system_figure = simulation.plot_system(show=False)
