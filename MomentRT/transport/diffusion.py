"""Compatibility imports for the P1 diffusion transport models."""

from __future__ import annotations
from MomentRT.transport.infinite_p1 import InfiniteP1 as InfiniteDiffusion
from MomentRT.transport.p1 import P1 as Diffusion
from MomentRT.transport.p1 import _finite_materials_at, _PreparedDiffusion

__all__ = [
    "Diffusion",
    "InfiniteDiffusion",
    "_PreparedDiffusion",
    "_finite_materials_at",
]
