"""Physics constants and isentropic constitutive relations.

All physics scalars live here. Pure numba functions for ρ(q²), M(q²), etc.
are unit-tested against hand calculations.

Reference: design.md §2 (Governing Equations)
"""

from . import isentropic

__all__ = ["isentropic"]
