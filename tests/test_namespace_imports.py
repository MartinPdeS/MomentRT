"""Domain namespaces expose the same public classes as convenience imports."""

import MomentRT


def test_domain_namespace_exports() -> None:
    """Grouped API paths preserve the established top-level class identities."""
    modules = {
        "material": ("Material", "SpectralMaterial"),
        "materials": ("Epidermis", "Dermis", "Water", "Intralipid"),
        "medium": ("HomogeneousMedium", "StructuredMedium"),
        "geometry": ("Box", "Sphere", "Layer"),
        "source": (
            "SourceSpectrum",
            "PointSource",
            "SurfaceSource",
            "UniformSurfaceSource",
            "SurfaceSourceGroup",
        ),
        "detector": ("DetectorResponsivity", "SurfaceDetector"),
        "boundary": ("Boundaries", "RobinBoundary", "ReflectingBoundary"),
        "mesh": ("CartesianMesh", "TetrahedralMesh"),
        "simulation": ("Simulation", "TransportSolver"),
        "fields": ("ScalarField", "VectorField", "BoundaryField"),
        "validation": ("ValidationReport",),
        "convergence": ("MeshConvergence",),
        "sections": ("SampleSection",),
    }
    for module_name, names in modules.items():
        module = getattr(MomentRT, module_name)
        for name in names:
            assert getattr(MomentRT, name) is getattr(module, name)
