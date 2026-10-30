"""Unstructured connectivity, conforming refinement, and field provenance."""

import numpy as np
import pytest
import MomentRT

u = MomentRT.ureg
requires_gmsh = pytest.mark.skipif(
    MomentRT.mesh.gmsh is None,
    reason="Gmsh's native runtime library is unavailable",
)


@pytest.fixture
def sample():
    material = MomentRT.Material(
        absorption=0.02 / u.mm,
        reduced_scattering=1 / u.mm,
        refractive_index=1.4,
        anisotropy=0.8,
        phase_function=MomentRT.HenyeyGreenstein(),
    )
    sample = MomentRT.StructuredMedium(origin=(0, 0, 0) * u.mm, size=(1, 1, 1) * u.mm)
    sample.add_background(material=material)
    sample.add_layer(name="Epidermis", material=material, thickness=0.25 * u.mm)
    return sample


@pytest.fixture
def cartesian(sample):
    return MomentRT.CartesianMesh.from_sample(sample=sample, shape=(4, 4, 5))


@pytest.fixture
def imported(cartesian):
    geometry = MomentRT.fem._MeshGeometry(mesh=cartesian)
    permutation = np.random.default_rng(10).permutation(len(geometry.points))
    inverse = np.argsort(permutation)
    return MomentRT.TetrahedralMesh(
        points=geometry.points[permutation] * u.m,
        tetrahedra=inverse[geometry.tetrahedra][:, [1, 0, 2, 3]],
    ), permutation


@pytest.mark.parametrize(
    "solver_class", [MomentRT.Diffusion, MomentRT.SP3, MomentRT.P3]
)
def test_unstructured_node_order_matches_cartesian_solution(
    sample, cartesian, imported, solver_class
):
    mesh, permutation = imported
    source = MomentRT.UniformSurfaceSource(
        face="z_min", wavelength=850 * u.nm, optical_power=1 * u.mW
    )
    detector = MomentRT.SurfaceDetector.on_face(
        sample=sample,
        face="z_max",
        coordinates=(0.5, 0.5) * u.mm,
        radius=0.3 * u.mm,
        orientation="outward_normal",
        acceptance_angle=45 * u.degree,
    )
    solver = solver_class(
        boundaries=MomentRT.Boundaries.uniform(
            condition=MomentRT.RobinBoundary(effective_reflectance=0)
        )
    )
    results = []
    for candidate in (cartesian, mesh):
        simulation = MomentRT.Simulation(
            medium=sample,
            source=source,
            mesh=candidate,
            solver=solver,
            detectors={"collection": detector},
        )
        assert simulation.validate().valid
        results.append(simulation.run())
    reference, actual = results
    assert actual.fluence.values.shape == mesh.shape
    np.testing.assert_allclose(
        actual.fluence.values.magnitude,
        reference.fluence.values.magnitude.ravel()[permutation],
        rtol=1e-11,
    )
    np.testing.assert_allclose(
        actual.flux.values.magnitude,
        reference.flux.values.magnitude,
        rtol=1e-10,
        atol=1e-10,
    )
    assert actual.detector_power["collection"].to("W").magnitude == pytest.approx(
        reference.detector_power["collection"].to("W").magnitude, rel=1e-11
    )
    assert abs(actual.energy_balance["residual_power"].to("W").magnitude) < 1e-14


def test_interpolation_uses_original_tetrahedra_and_preserves_snapshots(imported):
    mesh, _ = imported
    points = mesh.points.to("mm").magnitude
    affine = 1 + points @ np.array([2, -3, 4])
    field = MomentRT.ScalarField(mesh=mesh, values=affine * u.W / u.m**2, name="Affine")
    sampled = field.resample(resolution=0.2 * u.mm)
    expected = 1 + sampled.mesh.points.to("mm").magnitude @ np.array([2, -3, 4])
    np.testing.assert_allclose(sampled.values.magnitude, expected, atol=1e-12)
    np.testing.assert_array_equal(field.values.magnitude, affine)
    assert "display_interpolation" in sampled.metadata
    assert "display_interpolation" not in field.metadata
    with pytest.raises(ValueError, match="inside"):
        mesh._interpolate(points=np.array([[2, 0, 0]]) * u.mm, values=affine)


def test_invalid_connectivity_is_rejected(imported):
    mesh, _ = imported
    cells = mesh.tetrahedra
    with pytest.raises(ValueError, match="Duplicate"):
        MomentRT.TetrahedralMesh(
            points=mesh.points, tetrahedra=np.vstack([cells, cells[:1]])
        )
    with pytest.raises(ValueError, match="integer indices"):
        MomentRT.TetrahedralMesh(points=mesh.points, tetrahedra=cells.astype(float))
    invalid = cells.copy()
    invalid[0, 0] = len(mesh.points)
    with pytest.raises(ValueError, match="integer indices"):
        MomentRT.TetrahedralMesh(points=mesh.points, tetrahedra=invalid)
    invalid = cells.copy()
    invalid[0, 0] = invalid[0, 1]
    with pytest.raises(ValueError, match="volumes"):
        MomentRT.TetrahedralMesh(points=mesh.points, tetrahedra=invalid)
    with pytest.raises(ValueError):
        MomentRT.TetrahedralMesh(points=mesh.points, tetrahedra=cells[1:])


def test_imported_mesh_cannot_refine(imported):
    mesh, _ = imported
    assert len(mesh.copy().tetrahedra) == len(mesh.tetrahedra)
    with pytest.raises(ValueError, match="from_sample"):
        mesh.refine_region(name="Epidermis", resolution=0.1 * u.mm)


def test_generation_reports_missing_native_runtime(sample, monkeypatch):
    monkeypatch.setattr(MomentRT.mesh, "gmsh", None)
    with pytest.raises(RuntimeError, match="native runtime could not be loaded"):
        MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.4 * u.mm)


@requires_gmsh
def test_region_interfaces_and_atomic_refinement(sample):
    mesh = MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.35 * u.mm)
    original = mesh.copy()
    mesh.refine_region(name="Epidermis", resolution=0.12 * u.mm)
    assert len(mesh.tetrahedra) > len(original.tetrahedra)
    depths = mesh.points[mesh.tetrahedra].to("mm").magnitude[..., 2]
    assert not np.any(
        (depths.min(axis=1) < 0.25 - 1e-8) & (depths.max(axis=1) > 0.25 + 1e-8)
    )
    before = mesh.tetrahedra.copy()
    with pytest.raises(ValueError, match="Unknown"):
        mesh.refine_region(name="missing", resolution=0.1 * u.mm)
    np.testing.assert_array_equal(mesh.tetrahedra, before)
    with pytest.raises(ValueError, match="global"):
        mesh.refine_region(name="Epidermis", resolution=1 * u.mm)
    np.testing.assert_array_equal(mesh.tetrahedra, before)
    assert original.shape != mesh.shape


@requires_gmsh
def test_refine_near_respects_contact_footprint(sample):
    mesh = MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.35 * u.mm)
    source = MomentRT.SurfaceSource.on_face(
        sample=sample,
        face="z_min",
        coordinates=(0.5, 0.5) * u.mm,
        radius=0.15 * u.mm,
        wavelength=850 * u.nm,
        optical_power=1 * u.mW,
    )
    mesh.refine_near(objects=[source], resolution=0.06 * u.mm, distance=0.1 * u.mm)
    vertices = mesh.points[mesh.tetrahedra].to("mm").magnitude
    centers = vertices.mean(axis=1)
    gap = np.sqrt(
        np.maximum(0, np.linalg.norm(centers[:, :2] - 0.5, axis=1) - 0.15) ** 2
        + centers[:, 2] ** 2
    )
    sizes = np.max(
        np.linalg.norm(vertices[:, :, None] - vertices[:, None, :], axis=-1),
        axis=(1, 2),
    )
    assert np.median(sizes[gap < 0.1]) < 0.6 * np.median(sizes[gap > 0.6])
    solver = MomentRT.Diffusion(
        boundaries=MomentRT.Boundaries.uniform(
            condition=MomentRT.RobinBoundary(effective_reflectance=0)
        )
    )
    simulation = MomentRT.Simulation(
        medium=sample, source=source, mesh=mesh, solver=solver
    )
    assert simulation.validate().valid
    result = simulation.run()
    assert result.diagnostics["linear_relative_residual"] < 1e-10
    before = mesh.points.copy()
    with pytest.raises(ValueError, match="positive"):
        mesh.refine_near(objects=[source], resolution=0.06 * u.mm, distance=0 * u.mm)
    np.testing.assert_array_equal(mesh.points.magnitude, before.magnitude)


@requires_gmsh
def test_external_gmsh_session_is_preserved(sample):
    gmsh = MomentRT.mesh.gmsh
    gmsh.initialize(interruptible=False)
    gmsh.model.add("external")
    terminal = gmsh.option.getNumber("General.Terminal")
    try:
        MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.4 * u.mm)
        assert gmsh.isInitialized()
        assert gmsh.model.getCurrent() == "external"
        assert gmsh.option.getNumber("General.Terminal") == terminal
    finally:
        gmsh.finalize()


def test_tetrahedral_vtk_export_and_explicit_display_resolution(imported, tmp_path):
    pytest.importorskip("plotly")
    pv = pytest.importorskip("pyvista")
    mesh, _ = imported
    field = MomentRT.ScalarField(
        mesh=mesh, values=np.ones(mesh.shape) * u.W / u.m**2, name="Fluence"
    )
    with pytest.raises(ValueError, match="display_resolution"):
        field.plot(show=False)
    figure = field.plot(display_resolution=0.25 * u.mm, show=False)
    assert figure.data[0].type == "volume"
    mesh_figure = mesh.plot(show=False)
    assert mesh_figure.data[0].type == "mesh3d"
    assert mesh_figure.data[1].type == "scatter3d"
    target = tmp_path / "fluence.vtu"
    field.export(path=target)
    loaded = pv.read(target)
    assert np.all(loaded.celltypes == 10)
    np.testing.assert_array_equal(loaded.point_data["Fluence"], field.values.magnitude)


@pytest.mark.parametrize("face", ["x_min", "x_max", "y_min", "y_max", "z_max"])
@requires_gmsh
def test_contact_refinement_on_every_face(sample, face):
    mesh = MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.4 * u.mm)
    detector = MomentRT.SurfaceDetector.on_face(
        sample=sample,
        face=face,
        coordinates=(0.5, 0.5) * u.mm,
        radius=0.15 * u.mm,
        orientation="outward_normal",
        acceptance_angle=30 * u.degree,
    )
    mesh.refine_near(objects=[detector], resolution=0.12 * u.mm, distance=0.1 * u.mm)
    assert len(mesh.faces[face]) > 10


@requires_gmsh
def test_fitted_mesh_rejects_changed_sample_geometry(sample):
    mesh = MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.4 * u.mm)
    source = MomentRT.UniformSurfaceSource(
        face="z_min", wavelength=850 * u.nm, optical_power=1 * u.mW
    )
    solver = MomentRT.Diffusion(
        boundaries=MomentRT.Boundaries.uniform(
            condition=MomentRT.RobinBoundary(effective_reflectance=0)
        )
    )
    simulation = MomentRT.Simulation(
        medium=sample, source=source, mesh=mesh, solver=solver
    )
    result = simulation.run()
    original = result.fluence.mesh.points.copy()
    mesh.refine_region(name="Epidermis", resolution=0.2 * u.mm)
    np.testing.assert_array_equal(
        result.fluence.mesh.points.magnitude, original.magnitude
    )
    sample.structures[0].thickness = 0.4 * u.mm
    assert not simulation.validate().valid
    with pytest.raises(ValueError, match="geometry changed"):
        simulation.run()


def test_unstructured_sweep_preserves_mesh_and_coefficients(sample, imported):
    mesh, _ = imported
    source = MomentRT.UniformSurfaceSource(
        face="z_min", wavelength=850 * u.nm, optical_power=1 * u.mW
    )
    solver = MomentRT.SP3(
        boundaries=MomentRT.Boundaries.uniform(
            condition=MomentRT.RobinBoundary(effective_reflectance=0)
        )
    )
    simulation = MomentRT.Simulation(
        medium=sample, source=source, mesh=mesh, solver=solver
    )
    original = mesh.points.copy()
    sweep = simulation.sweep(
        wavelengths=np.array([800, 850]) * u.nm, optical_power=np.array([1, 2]) * u.mW
    )
    assert sweep.diagnostics["geometry_preparations"] == 1
    np.testing.assert_allclose(
        sweep.results[1].fluence.values.magnitude,
        2 * sweep.results[0].fluence.values.magnitude,
        rtol=1e-12,
    )
    np.testing.assert_array_equal(mesh.points.magnitude, original.magnitude)


@requires_gmsh
def test_generation_failure_preserves_mesh_and_releases_session(sample, monkeypatch):
    mesh = MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.4 * u.mm)
    before = mesh.tetrahedra.copy()

    def fail(dimension):
        raise RuntimeError("mesher failed")

    monkeypatch.setattr(MomentRT.mesh.gmsh.model.mesh, "generate", fail)
    with pytest.raises(RuntimeError, match="mesher failed"):
        mesh.refine_region(name="Epidermis", resolution=0.2 * u.mm)
    np.testing.assert_array_equal(mesh.tetrahedra, before)
    assert not MomentRT.mesh.gmsh.isInitialized()


def test_log_display_checks_original_tetrahedral_values(imported):
    mesh, _ = imported
    values = np.ones(mesh.shape)
    values[0] = -1
    field = MomentRT.ScalarField(
        mesh=mesh, values=values * u.W / u.m**2, name="Fluence"
    )
    with pytest.raises(ValueError, match="original nodal"):
        field.plot(scale="log", display_resolution=2 * u.mm, show=False)


@requires_gmsh
def test_box_and_sphere_region_generation(sample):
    sample.add_structures(
        MomentRT.Box(
            name="Box",
            origin=(0.4, 0.4, 0.4) * u.mm,
            size=(0.3, 0.3, 0.3) * u.mm,
            material=sample.background,
        ),
        MomentRT.Sphere(
            name="Sphere",
            center=(0.72, 0.71, 0.68) * u.mm,
            radius=0.15 * u.mm,
            material=sample.background,
        ),
    )
    mesh = MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.3 * u.mm)
    mesh.refine_region(name="Sphere", resolution=0.12 * u.mm)
    mesh.refine_region(name="Box", resolution=0.15 * u.mm)
    centers = mesh.points[mesh.tetrahedra].mean(axis=1)
    assert set(np.unique(sample.material_labels(points=centers))) == {0, 1, 2, 3}


@requires_gmsh
def test_contact_groups_and_repeated_refinement(sample):
    sources = [
        MomentRT.SurfaceSource.on_face(
            sample=sample,
            face="z_min",
            coordinates=(x, 0.5) * u.mm,
            radius=0.1 * u.mm,
            wavelength=850 * u.nm,
            optical_power=0.5 * u.mW,
        )
        for x in (0.3, 0.7)
    ]
    group = MomentRT.SurfaceSourceGroup(sources=sources)
    detector = MomentRT.SurfaceDetector.on_face(
        sample=sample,
        face="z_min",
        coordinates=(0.3, 0.5) * u.mm,
        radius=0.1 * u.mm,
        orientation="outward_normal",
        acceptance_angle=30 * u.degree,
    )
    mesh = MomentRT.TetrahedralMesh.from_sample(sample=sample, resolution=0.4 * u.mm)
    mesh.refine_near(
        objects=[group, detector], resolution=0.15 * u.mm, distance=0.08 * u.mm
    )
    first = len(mesh.points)
    mesh.refine_near(objects=sources, resolution=0.04 * u.mm, distance=0.08 * u.mm)
    assert len(mesh.points) > first
