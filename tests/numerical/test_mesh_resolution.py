"""Physical mesh spacing, sample coverage, and finite-solver snapshots."""

import numpy as np
import pytest
from pint import DimensionalityError
import MomentRT


@pytest.fixture
def sample():
    unit = MomentRT.ureg

    sample = MomentRT.StructuredMedium(
        origin=np.array([-0.5, -1, -2]) * unit.mm,
        size=np.array([1, 2, 3]) * unit.mm,
    )

    sample.add_background(
        material=MomentRT.Material(
            absorption=0.02 / unit.mm,
            reduced_scattering=1 / unit.mm,
            refractive_index=1.4,
            anisotropy=0.8,
            phase_function=MomentRT.HenyeyGreenstein(),
        ),
    )

    return sample


@pytest.mark.parametrize(
    "resolution, expected_shape",
    [
        (300 * MomentRT.ureg.um, (5, 8, 11)),
        ((250, 500, 1000) * MomentRT.ureg.um, (5, 5, 4)),
        (10 * MomentRT.ureg.mm, (2, 2, 2)),
    ],
)
def test_resolution_fits_sample_with_bounded_spacing(
    sample, resolution, expected_shape
):
    mesh = MomentRT.CartesianMesh.from_sample(
        sample=sample,
        resolution=resolution,
    )

    assert mesh.shape == expected_shape

    np.testing.assert_allclose(
        mesh.points[0, 0, 0].to("mm").magnitude,
        sample.origin.to("mm").magnitude,
    )

    np.testing.assert_allclose(
        mesh.points[-1, -1, -1].to("mm").magnitude,
        (sample.origin + sample.size).to("mm").magnitude,
    )

    assert np.all(mesh.spacing <= resolution * (1 + 1e-14))

    np.testing.assert_array_equal(mesh.resolution.magnitude, resolution.magnitude)


@pytest.mark.parametrize(
    "resolution, exception",
    [
        (0.3, TypeError),
        (0 * MomentRT.ureg.mm, ValueError),
        (-1 * MomentRT.ureg.mm, ValueError),
        (np.nan * MomentRT.ureg.mm, ValueError),
        (np.inf * MomentRT.ureg.mm, ValueError),
        ((1, 2) * MomentRT.ureg.mm, ValueError),
        (np.ones((1, 3)) * MomentRT.ureg.mm, ValueError),
        (1 * MomentRT.ureg.second, DimensionalityError),
    ],
)
def test_invalid_resolution(sample, resolution, exception):
    with pytest.raises(exception):
        MomentRT.CartesianMesh.from_sample(
            sample=sample,
            resolution=resolution,
        )


def test_resolution_and_shape_are_exclusive(sample):
    with pytest.raises(ValueError, match="exactly one"):
        MomentRT.CartesianMesh.from_sample(sample=sample)

    with pytest.raises(ValueError, match="exactly one"):
        MomentRT.CartesianMesh.from_sample(
            sample=sample,
            shape=(3, 3, 3),
            resolution=1 * MomentRT.ureg.mm,
        )

    mesh = MomentRT.CartesianMesh.from_sample(sample=sample, shape=(3, 4, 5))

    assert mesh.shape == (3, 4, 5)

    assert mesh.resolution is None


def test_integral_interval_counts_tolerate_unit_conversion_roundoff():
    sample = MomentRT.StructuredMedium(
        origin=np.zeros(3) * MomentRT.ureg.mm,
        size=np.array([0.07, 0.14, 0.21]) * MomentRT.ureg.mm,
    )

    mesh = MomentRT.CartesianMesh.from_sample(
        sample=sample,
        resolution=70 * MomentRT.ureg.um,
    )

    assert mesh.shape == (2, 3, 4)


def test_mesh_copy_preserves_independent_resolution(sample):
    requested = np.array([250, 500, 1000]) * MomentRT.ureg.um

    mesh = MomentRT.CartesianMesh.from_sample(sample=sample, resolution=requested)

    snapshot = mesh.copy()

    requested[0] = 900 * MomentRT.ureg.um

    mesh.spacing[0] = 700 * MomentRT.ureg.um

    mesh.origin[0] = 0 * MomentRT.ureg.mm

    assert snapshot.spacing[0].to("um").magnitude == pytest.approx(250)

    assert snapshot.origin[0].to("mm").magnitude == pytest.approx(-0.5)

    assert snapshot.resolution[0].to("um").magnitude == pytest.approx(250)


@pytest.mark.parametrize(
    "solver_class", [MomentRT.Diffusion, MomentRT.SP3, MomentRT.P3]
)
def test_resolution_mesh_matches_explicit_counts_in_finite_solvers(
    sample, solver_class
):
    mesh = MomentRT.CartesianMesh.from_sample(
        sample=sample,
        resolution=0.5 * MomentRT.ureg.mm,
    )

    source = MomentRT.UniformSurfaceSource(
        face="z_min",
        wavelength=850 * MomentRT.ureg.nm,
        optical_power=1 * MomentRT.ureg.mW,
    )

    solver = solver_class(
        boundaries=MomentRT.Boundaries.uniform(
            condition=MomentRT.RobinBoundary(effective_reflectance=0),
        ),
    )

    result = solver.solve(medium=sample, source=source, mesh=mesh)

    reference = solver.solve(
        medium=sample,
        source=source,
        mesh=MomentRT.CartesianMesh.from_sample(sample=sample, shape=(3, 5, 7)),
    )

    np.testing.assert_allclose(
        result.fluence.values.magnitude,
        reference.fluence.values.magnitude,
        rtol=1e-12,
    )

    assert result.fluence.mesh.resolution == 0.5 * MomentRT.ureg.mm

    assert result.fluence.mesh is not mesh
