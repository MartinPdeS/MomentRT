"""
How dermis thickness changes a two-fiber scan
================================================================================

Keep the outer sample and contact probe fixed while varying dermis thickness
from 0.5 to 2.5 mm. Use an absorbing dermis with 3x the preset absorption to
explore how layer extent changes the collected signal. These geometries and
optical coefficients are illustrative, not a model of a diagnosed skin condition.

Compare all three thicknesses with the 1.5 mm reference at the SAME absorption.
Moving the dermis-hypodermis interface changes the volume occupied by both
materials, including their different scattering coefficients. This is a
geometrical sensitivity experiment; it does not isolate absorption alone or
establish that thickness can be uniquely recovered from a measurement.

Baseline optical properties follow MomentRT's 850 nm Shimojo preset. P1,
matched-index boundaries, and Lambertian collection approximate the instrument.

References
--------------------------------------------------------------------------------
Shimojo, Y. et al. (2020), *Journal of Biomedical Optics*, 25(4), 045002.
https://doi.org/10.1117/1.JBO.25.4.045002
"""

# sphinx_gallery_thumbnail_path = '_static/images/gallery_systems/skin_dermis_thickness.png'
import matplotlib.pyplot as plt
import numpy as np
import MomentRT
from MomentRT.models.skin import ShimojoThreeLayer850nm

ureg = MomentRT.ureg

# %%
# Configure one source and a scanning collection fiber
# ------------------------------------------------------------------------------
# The surface is z = 0 and positive z points into the sample. The preset uses
# 0.2 mm epidermis, 1.5 mm dermis, and hypodermis below. Source power is diffuse
# power already coupled into tissue. All layers share refractive index 1.4,
# with explicitly matched-index, zero-reflectance Robin conditions on all faces.
# Collection circles are alternative positions of ONE passive fiber, not a
# physical array. Each uses a 30 degree acceptance half-angle and a Lambertian
# angular model; fiber propagation, entry losses, and specular return are omitted.

baseline_sample = ShimojoThreeLayer850nm(
    origin=(-5, -4, 0) * ureg.mm,
    size=(12, 8, 6) * ureg.mm,
)

source = MomentRT.SurfaceSource.on_face(
    sample=baseline_sample,
    face="z_min",
    coordinates=(-2, 0) * ureg.mm,
    radius=0.3 * ureg.mm,
    wavelength=850 * ureg.nm,
    optical_power=1 * ureg.mW,
    name="Illumination fiber",
)

separations = np.linspace(0.8, 4.0, 9) * ureg.mm

detectors = {
    f"position_{index}": MomentRT.SurfaceDetector.on_face(
        sample=baseline_sample,
        face="z_min",
        coordinates=np.array([-2 + separation.to("mm").magnitude, 0]) * ureg.mm,
        radius=0.3 * ureg.mm,
        orientation="outward_normal",
        acceptance_angle=30 * ureg.degree,
        name=f"Collection at {separation.to('mm').magnitude:.1f} mm",
    )
    for index, separation in enumerate(separations)
}

solver = MomentRT.P1(
    boundaries=MomentRT.Boundaries.uniform(
        condition=MomentRT.RobinBoundary(effective_reflectance=0),
    ),
)

mesh_resolutions_mm = (0.55, 0.4)

contact_resolutions_mm = (0.12, 0.08)

# %%
# Move the lower dermal interface and remesh each geometry
# ------------------------------------------------------------------------------
# Keep epidermis thickness and total sample depth fixed. Hypodermis fills the
# remainder. Build a new conforming mesh for EACH geometry; a mesh fitted to one
# interface position would not resolve the other positions correctly.
# At each refinement level the same target sizes and contact refinements apply.

materials = {layer.name: layer.material for layer in baseline_sample.structures}

epidermis = materials["Epidermis"]

baseline_dermis = materials["Dermis"]

hypodermis = materials["Hypodermis (subcutaneous fat)"]

dermis = MomentRT.Material(
    absorption=3 * baseline_dermis.absorption,
    reduced_scattering=baseline_dermis.reduced_scattering,
    refractive_index=baseline_dermis.refractive_index,
    anisotropy=baseline_dermis.anisotropy,
    phase_function=baseline_dermis.phase_function,
)

dermis_thicknesses = np.array([0.5, 1.5, 2.5]) * ureg.mm

case_names = tuple(
    f"Dermis {thickness.to('mm').magnitude:g} mm" for thickness in dermis_thicknesses
)

reference_case = 1

collected_power_uw = np.empty(
    (len(mesh_resolutions_mm), len(case_names), len(detectors))
)

area_errors = np.empty((len(mesh_resolutions_mm), len(case_names)))

for level, (bulk_resolution, contact_resolution) in enumerate(
    zip(mesh_resolutions_mm, contact_resolutions_mm)
):
    for case, thickness in enumerate(dermis_thicknesses):
        sample = MomentRT.models.phantom.LayeredPhantom(
            origin=baseline_sample.origin,
            size=baseline_sample.size,
            background=hypodermis,
            layers=(
                ("Epidermis", 0.2 * ureg.mm, epidermis),
                ("Dermis", thickness, dermis),
                (
                    "Hypodermis",
                    baseline_sample.size[2] - 0.2 * ureg.mm - thickness,
                    hypodermis,
                ),
            ),
        )

        mesh = MomentRT.TetrahedralMesh.from_sample(
            sample=sample,
            resolution=bulk_resolution * ureg.mm,
        )

        mesh.refine_region(
            name="Epidermis",
            resolution=contact_resolution * ureg.mm,
        )

        mesh.refine_near(
            objects=[source, *detectors.values()],
            resolution=contact_resolution * ureg.mm,
            distance=0.3 * ureg.mm,
        )

        simulation = MomentRT.Simulation(
            medium=sample,
            source=source,
            mesh=mesh,
            solver=solver,
            detectors=detectors,
        )

        result = simulation.run()

        collected_power_uw[level, case] = [
            result.detector_power[name].to("uW").magnitude for name in detectors
        ]

        area_ratios = [
            (result.source_sampled_area / result.source_geometric_area)
            .to("")
            .magnitude,
            *[
                (measurement.sampled_area / measurement.geometric_area).to("").magnitude
                for measurement in result.detectors.values()
            ],
        ]

        area_errors[level, case] = np.max(np.abs(np.asarray(area_ratios) - 1))

        residual_fraction = (
            (result.energy_balance["residual_power"] / source.optical_power)
            .to("")
            .magnitude
        )

        print(
            f"Mesh {level + 1}, {case_names[case]}: "
            f"{len(mesh.tetrahedra)} tetrahedra, "
            f"maximum source/aperture area error = {area_errors[level, case]:.2%}, "
            f"energy residual / input = {residual_fraction:.2e}, "
            f"minimum fluence = {result.fluence.values.min().to('W/m^2'):~P}"
        )

# %%
# Report spatial sensitivity without hiding unresolved values
# ------------------------------------------------------------------------------
# These two meshes measure sensitivity to one refinement, not a converged
# sequence. Disk area checks do not establish field accuracy. Minimum fluence
# and energy residuals above retain numerical artifacts. Never clip negative
# powers before plotting or infer transport accuracy from energy conservation.
# Quantitative probe design also needs larger-domain checks and comparison with
# a radiative-transfer or Monte Carlo reference, especially near the source.

if not np.all(np.isfinite(collected_power_uw)) or np.any(collected_power_uw <= 0):
    raise ValueError("Unresolved or nonpositive detector powers; refine the study.")

relative_change = (
    np.abs(collected_power_uw[-1] - collected_power_uw[-2]) / collected_power_uw[-1]
)

print(f"Maximum collected-power change on refinement: {np.max(relative_change):.2%}")

for case, name in enumerate(case_names):
    print(
        f"{name}: finest-mesh power at 0.8 / 4 mm = "
        f"{collected_power_uw[-1, case, 0]:.4g} / "
        f"{collected_power_uw[-1, case, -1]:.4g} uW"
    )

# %%
# Compare thickness-dependent collected power
# ------------------------------------------------------------------------------
# Solid and dotted curves show finer and coarser meshes, respectively. Ratios
# compare with the 1.5 mm dermis at each mesh level, avoiding changes to source
# power or detector optics. An interface shift replaces hypodermis with dermis
# or vice versa; a change in this response cannot be assigned to a unique depth.

figure, axes = plt.subplots(
    1,
    2,
    figsize=(12, 4.5),
)

for case, name in enumerate(case_names):
    for level, linestyle in enumerate((":", "-")):
        axes[0].semilogy(
            separations.to("mm").magnitude,
            collected_power_uw[level, case],
            linestyle=linestyle,
            marker="o" if level else None,
            color=f"C{case}",
            label=name if level else None,
        )

        axes[1].plot(
            separations.to("mm").magnitude,
            collected_power_uw[level, case] / collected_power_uw[level, reference_case],
            linestyle=linestyle,
            marker="o" if level else None,
            color=f"C{case}",
        )

    print(
        f"{name}: power / 1.5 mm reference at 0.8 / 4 mm = "
        f"{collected_power_uw[-1, case, 0] / collected_power_uw[-1, reference_case, 0]:.3f} / "
        f"{collected_power_uw[-1, case, -1] / collected_power_uw[-1, reference_case, -1]:.3f}"
    )

axes[0].set_title("Vary dermis thickness at fixed optical coefficients")
axes[0].set_ylabel("Collected power (µW), 1 mW coupled input")
axes[0].legend()
axes[1].set_title("Contrast relative to a 1.5 mm dermis")
axes[1].set_ylabel("Collected power / reference collected power")

for axis in axes:
    axis.set_xlabel("Fiber center-to-center separation (mm)")
    axis.grid(
        True,
        which="both",
        alpha=0.25,
    )

figure.tight_layout()

# %%
# Inspect the final case in three dimensions
# ------------------------------------------------------------------------------
# Rings show alternative scan positions. Linear color preserves the sign of
# the interpolated field; sampled layer surfaces are construction overlays.

gallery_system_figure = result.fluence.plot(
    view="volume",
    scale="linear",
    display_resolution=0.3 * ureg.mm,
    overlays=[source, *detectors.values(), *sample.structures],
    show=False,
)
