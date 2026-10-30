"""Premade material configurations for illustrative tissue models."""

from __future__ import annotations
import MomentRT.materials.liquids as liquids  # noqa: F401 -- public submodule
import MomentRT.materials.skin as skin  # noqa: F401 -- public submodule

__all__ = ["liquids", "skin"]
