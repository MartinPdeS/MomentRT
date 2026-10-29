"""
Compare full P3 with SP3 and P1
================================

Solve the same seven-fiber contact geometry on two skin layers with three
angular approximations.
Full P3 retains the three-dimensional spherical-harmonic moments through
order three. The supplied values are illustrative, and agreement between the
models does not establish accuracy against radiative transfer.
"""

# sphinx_gallery_thumbnail_path = '_static/images/p3-fluence.png'

import os

import numpy as np

import MomentRT

ureg = MomentRT.ureg

show_figures = os.environ.get("MOMENTRT_SPHINX_GALLERY") != "1"

# %%
# Define the common physical configuration
# ----------------------------------------
# The Henyey–Greenstein anisotropy belongs to the material. The three
# models share matched-index faces, six diffuse sources, mesh, and aperture.
# The finite mesh samples the layer interface at tetrahedron centroids, so
# the interface geometry and disk footprints also require mesh refinement.

dermis = MomentRT.Material(
    absorption=0.03 / ureg.millimeter,
    reduced_scattering=1.0 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.8,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

epidermis = MomentRT.Material(
    absorption=0.06 / ureg.millimeter,
    reduced_scattering=1.6 / ureg.millimeter,
    refractive_index=1.4,
    anisotropy=0.7,
    phase_function=MomentRT.phase_function.HenyeyGreenstein(),
)

sample = MomentRT.StructuredMedium(
    origin=np.array([-4, -4, 0]) * ureg.millimeter,
    size=np.array([8, 8, 4]) * ureg.millimeter,
)

sample.add_background(material=dermis)

sample.add_layer(
    name="Epidermis",
    thickness=0.8 * ureg.millimeter,
    material=epidermis,
)

sources = []

for fiber_number in range(6):
    angle = fiber_number * np.pi / 3

    coordinates = np.array([np.cos(angle), np.sin(angle)]) * 1.8 * ureg.millimeter

    source = MomentRT.SurfaceSource.on_face(
        sample=sample,
        face="z_min",
        coordinates=coordinates,
        radius=0.8 * ureg.millimeter,
        wavelength=850 * ureg.nanometer,
        optical_power=(1 / 6) * ureg.milliwatt,
    )

    sources.append(source)

detector = MomentRT.SurfaceDetector.on_face(
    sample=sample,
    face="z_min",
    coordinates=np.array([0, 0]) * ureg.millimeter,
    radius=0.8 * ureg.millimeter,
    acceptance_angle=90 * ureg.degree,
    orientation="outward_normal",
)

mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=(9, 9, 7))

boundary = MomentRT.RobinBoundary(effective_reflectance=0)

boundaries = MomentRT.Boundaries.uniform(condition=boundary)

# %%
# Solve and compare the detector signal
# -------------------------------------
# P3 solves six even angular moments and reconstructs ten odd moments from
# their gradients. SP3 instead solves two composite fields. The contact
# detector currently uses a Lambertian angular-acceptance approximation for
# all three models, including P3. Each linear diffuse solve uses one source;
# the six fluences and detector powers are added incoherently. This does not
# model optical interference, fiber modes, or Fresnel entry losses.

solvers = {
    "P1": MomentRT.Diffusion(boundaries=boundaries),
    "SP3": MomentRT.SP3(boundaries=boundaries),
    "P3": MomentRT.P3(boundaries=boundaries),
}

combined_fields = {}

for name, solver in solvers.items():
    results = []

    for fiber_number, source in enumerate(sources, start=1):
        simulation = MomentRT.Simulation(
            medium=sample,
            source=source,
            mesh=mesh,
            solver=solver,
            detectors={"contact": detector},
        )

        report = simulation.validate()

        if not report.valid:
            raise ValueError(report.summary())

        result = simulation.run()

        results.append(result)

        if name == "P1":
            print(
                f"Fiber {fiber_number} source area: sampled "
                f"{result.source_sampled_area.to('mm^2'):~P}, "
                f"geometric {result.source_geometric_area.to('mm^2'):~P}"
            )

    combined_values = sum(result.fluence.values for result in results)

    combined_fields[name] = MomentRT.ScalarField(
        mesh=results[0].fluence.mesh,
        values=combined_values,
        name=f"Combined {name} fluence rate",
        metadata={**results[0].fluence.metadata, "source_count": len(sources)},
    )

    collected_power = sum(result.detector_power["contact"] for result in results)

    residual_power = sum(result.energy_balance["residual_power"] for result in results)

    print(f"{name} collected power: {collected_power.to('mW'):~P}")

    print(f"{name} minimum fluence: {combined_values.min():~P}")

    print(f"{name} energy residual: {residual_power:~P}")

    print(
        f"{name} sampled detector area: {results[0].detectors['contact'].sampled_area.to('mm^2'):~P}"
    )

    print(
        f"{name} geometric detector area: {results[0].detectors['contact'].geometric_area.to('mm^2'):~P}"
    )

# %%
# Explore the full P3 fluence
# ----------------------------
# Linear scale retains any negative values from the angular or spatial
# approximation. Refine the mesh and compare against an independent transport
# reference before interpreting the result for a particular sensor. The
# gallery preview is captured from this same combined P3 volume.

overlays = {
    **{
        f"Illumination fiber {number}": source
        for number, source in enumerate(sources, start=1)
    },
    "Collection fiber": detector,
}

fluence_figure = combined_fields["P3"].plot(
    view="volume",
    scale="linear",
    overlays=overlays,
    show=show_figures,
)
