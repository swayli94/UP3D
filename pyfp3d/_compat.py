"""Numpy version shims (phase two GS0.1).

The library used `np.trapezoid` directly, which only exists on numpy >= 2.0,
while `pyproject.toml` declares `numpy>=1.23` -- so on a numpy 1.x install
`pyfp3d.post.surface.cl_kj_3d` raised AttributeError and four tests failed at
import/call time (audit 2026-07-28 §6.2). numpy 1.x spells the same function
`np.trapz`; numpy 2.0 renamed it and dropped the old name.

Import from here instead of reaching for either name directly:

    from pyfp3d._compat import trapezoid
"""

import numpy as np

if hasattr(np, "trapezoid"):        # numpy >= 2.0
    trapezoid = np.trapezoid
else:                              # numpy 1.x
    trapezoid = np.trapz

__all__ = ["trapezoid"]
