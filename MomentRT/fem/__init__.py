"""Linear tetrahedral FEM geometry and numerical kernels."""

from __future__ import annotations
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

__all__ = [
    "FiniteElementGeometry",
    "boundary_triangles",
    "element_geometry",
    "source_irradiance",
    "source_load",
    "tetrahedralize",
    "triangle_areas",
    "triangle_quadrature",
]
