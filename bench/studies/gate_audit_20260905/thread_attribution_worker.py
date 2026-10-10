"""One D05 M0.80/a1.25 medium leg, under whatever thread caps the parent set.
Env vars must be set BEFORE numba/BLAS import, hence a separate process."""
import json, os, sys
sys.path.insert(0, "/home/lrz/codes/UP3D")
import numpy as np
from tests.D.test_D05_euler_naca0012 import _one, _with_cp_rms
d = _with_cp_rms(_one("medium", 0.80, 1.25), "n0012_m0800_a1.25")
print("RESULT " + json.dumps({
    "cl": d["cl"], "cd": d["cd"], "cp_rms_upper": d["cp_rms_upper"],
    "x_shock": d["x_shock"], "residual": d["residual"],
    "converged": bool(d["converged"]),
    "numba": os.environ.get("NUMBA_NUM_THREADS"),
    "blas": os.environ.get("OPENBLAS_NUM_THREADS"),
}))
