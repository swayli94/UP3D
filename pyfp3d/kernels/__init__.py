"""Element-wise residual and Jacobian assembly kernels.

Numba-jitted kernels for:
  - Laplace residual (P1)
  - Isentropic density upwinding (P3)
  - Newton linearization (P4)

All kernels use SoA arrays, colored element loops with @prange.

Reference: design.md §6 (Assembly Strategy)
"""

__all__ = []
