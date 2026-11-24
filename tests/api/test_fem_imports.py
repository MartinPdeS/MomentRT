"""FEM package exports preserve the standalone kernel import paths."""

import pickle
import MomentRT.fem
from MomentRT.fem.geometry import FiniteElementGeometry
from MomentRT.fem.kernels import (
    boundary_triangles,
    element_geometry,
    source_irradiance,
    source_load,
    tetrahedralize,
    triangle_areas,
    triangle_quadrature,
)


def test_fem_import_paths_share_implementations() -> None:
    """Package-level names and direct module imports refer to the same objects."""
    assert MomentRT.fem.FiniteElementGeometry is FiniteElementGeometry
    assert MomentRT.fem.boundary_triangles is boundary_triangles
    assert MomentRT.fem.element_geometry is element_geometry
    assert MomentRT.fem.source_irradiance is source_irradiance
    assert MomentRT.fem.source_load is source_load
    assert MomentRT.fem.tetrahedralize is tetrahedralize
    assert MomentRT.fem.triangle_areas is triangle_areas
    assert MomentRT.fem.triangle_quadrature is triangle_quadrature


def test_previous_class_reference_still_resolves() -> None:
    """A class reference serialized under the former module path resolves."""
    previous_reference = b"cMomentRT.fem\nFiniteElementGeometry\n."
    assert pickle.loads(previous_reference) is FiniteElementGeometry
