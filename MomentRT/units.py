"""TypedUnit registry and validation API shared by MomentRT objects.

``validate_units`` is a compatibility alias for TypedUnit's direct quantity
validator with the shared registry already selected. It preserves supplied
quantities and their units; conversions remain explicit at numerical boundaries.
"""

from __future__ import annotations
from functools import partial
import TypedUnit
from pint import Quantity, UnitRegistry

ureg: UnitRegistry = TypedUnit.ureg
validate_quantity = TypedUnit.validate_quantity
validate_annotation = TypedUnit.validate_annotation
unit_validator = TypedUnit.unit_validator
validate_units = partial(TypedUnit.validate_quantity, registry=ureg)

__all__ = [
    "ureg",
    "Quantity",
    "unit_validator",
    "validate_annotation",
    "validate_quantity",
    "validate_units",
]
