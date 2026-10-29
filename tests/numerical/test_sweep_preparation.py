"""Prepared sweeps must reuse geometry and agree with independent fresh solves."""

from copy import copy
from functools import wraps

import numpy as np
import pytest

import MomentRT
from MomentRT import (
    CartesianMesh,
    Diffusion,
    RobinBoundary,
    Simulation,
    SpectralMaterial,
    StructuredMedium,
    SurfaceDetector,
    SurfaceSource,
    fem,
    ureg,
)


def configuration():
    background = SpectralMaterial(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        wavelengths=np.array([700, 900]) * ureg.nm,
        absorption=np.array([0.01, 0.04]) / ureg.mm,
        reduced_scattering=np.array([1.2, 0.8]) / ureg.mm,
        refractive_index=1.4,
    )

    upper = SpectralMaterial(
        phase_function=MomentRT.phase_function.HenyeyGreenstein(),
        wavelengths=np.array([700, 900]) * ureg.nm,
        absorption=np.array([0.04, 0.08]) / ureg.mm,
        reduced_scattering=np.array([1.5, 1.0]) / ureg.mm,
        refractive_index=1.4,
    )

    sample = StructuredMedium(
        origin=np.zeros(3) * ureg.mm,
        size=np.array([4, 4, 4]) * ureg.mm,
    )

    sample.add_background(material=background)

    sample.add_layer(name="upper", thickness=1 * ureg.mm, material=upper)

    source = SurfaceSource(
        face="z_min",
        position=np.array([1, 2, 0]) * ureg.mm,
        radius=0.7 * ureg.mm,
        wavelength=800 * ureg.nm,
        optical_power=1 * ureg.mW,
    )

    detector = SurfaceDetector(
        face="z_max",
        position=np.array([2, 2, 4]) * ureg.mm,
        radius=0.9 * ureg.mm,
        collection_axis=[0, 0, 1],
        acceptance_angle=45 * ureg.degree,
    )

    mesh = CartesianMesh.from_sample(sample=sample, shape=(7, 7, 9))

    boundary = RobinBoundary(effective_reflectance=0.4)

    return Simulation(
        medium=sample,
        source=source,
        mesh=mesh,
        solver=Diffusion(
            boundaries={
                face: boundary
                for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
            },
        ),
        detectors={"transmission": detector},
    )


def test_geometry_and_labels_are_prepared_once_per_sweep(monkeypatch):
    calls = {}

    def count(function, name):
        @wraps(function)
        def wrapped(*args, **kwargs):
            calls[name] = calls.get(name, 0) + 1

            return function(*args, **kwargs)

        return wrapped

    for name in (
        "tetrahedralize",
        "element_geometry",
        "boundary_triangles",
        "triangle_quadrature",
    ):
        monkeypatch.setattr(fem, name, count(getattr(fem, name), name))

    monkeypatch.setattr(
        StructuredMedium,
        "material_labels",
        count(StructuredMedium.material_labels, "material_labels"),
    )

    simulation = configuration()

    spectrum = simulation.sweep(wavelengths=np.array([700, 800, 900]) * ureg.nm)

    assert calls == {
        "tetrahedralize": 1,
        "element_geometry": 1,
        "boundary_triangles": 1,
        "triangle_quadrature": 6,
        "material_labels": 1,
    }

    assert spectrum.diagnostics["geometry_preparations"] == 1

    assert spectrum.diagnostics["solve_times"].to("s").shape == (3,)

    assert [
        result.diagnostics["preparation_reused"] for result in spectrum.results
    ] == [False, True, True]

    # A later call creates fresh geometry rather than relying on stale state.
    simulation.mesh = CartesianMesh.from_sample(
        sample=simulation.medium, shape=(9, 9, 9)
    )

    simulation.source.radius = 0.9 * ureg.mm

    second = simulation.sweep(wavelengths=np.array([800, 900]) * ureg.nm)

    assert calls["material_labels"] == 2

    assert calls["triangle_quadrature"] == 12

    assert second.results[0].fluence.values.shape == (9, 9, 9)

    assert second.results[0].diagnostics["source_geometric_area"].to(
        "mm^2"
    ).magnitude == pytest.approx(np.pi * 0.9**2)


def test_prepared_fields_and_measurements_match_fresh_solves():
    simulation = configuration()

    wavelengths = np.array([900, 700, 800]) * ureg.nm

    powers = np.array([2, 1, 3]) * ureg.mW

    spectrum = simulation.sweep(wavelengths=wavelengths, optical_power=powers)

    for wavelength, power, reused in zip(wavelengths, powers, spectrum.results):
        source = copy(simulation.source)

        source.wavelength = wavelength

        source.optical_power = power

        independent = Simulation(
            medium=simulation.medium,
            source=source,
            mesh=simulation.mesh,
            solver=simulation.solver,
            detectors=simulation.detectors,
        )

        fresh = independent.run()

        np.testing.assert_allclose(
            reused.fluence.values.magnitude, fresh.fluence.values.magnitude, rtol=1e-12
        )

        np.testing.assert_allclose(
            reused.flux.values.magnitude, fresh.flux.values.magnitude, atol=1e-10
        )

        for face in reused.boundary_fields:
            np.testing.assert_allclose(
                reused.boundary_fields[face].escaping_irradiance.magnitude,
                fresh.boundary_fields[face].escaping_irradiance.magnitude,
                atol=1e-10,
            )

        assert reused.detector_power["transmission"].to("W").magnitude == pytest.approx(
            fresh.detector_power["transmission"].to("W").magnitude, rel=1e-12
        )

        assert abs(reused.energy_balance["residual_power"].to("W").magnitude) < 1e-13

    # Mutating one returned observation mesh cannot change another result or input.
    spectrum.results[0].fluence.mesh.origin[0] = 10 * ureg.mm

    assert spectrum.results[1].fluence.mesh.origin[0] == 0 * ureg.mm

    assert simulation.mesh.origin[0] == 0 * ureg.mm

    assert simulation.source.wavelength == 800 * ureg.nm

    assert simulation.source.optical_power == 1 * ureg.mW
