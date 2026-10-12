"""Same leg, but with n_newton_max as a second axis, plus a failure
classification -- is it budget-limited or a limit cycle?"""
import json, os, sys
sys.path.insert(0, "/home/lrz/codes/UP3D")
import numpy as np
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
from tests.D.test_D05_euler_naca0012 import RECIPE
from tests.conftest import REPO_ROOT

nmax = int(sys.argv[1])
mc, wc = cut_wake(read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh"))
kw = dict(RECIPE); kw["n_newton_max"] = nmax
r = solve_newton_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25, **kw)
phi = np.asarray(r["phi"]); dz = float(np.ptp(mc.nodes[:, 2]))
f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                            phi, alpha_deg=1.25, u_inf=1.0, s_ref=dz, m_inf=0.80)
h = np.asarray(r["residual_history"], float)
tail = h[-10:] if len(h) >= 10 else h
print("RESULT " + json.dumps({
    "n_newton_max": nmax, "blas": os.environ.get("OPENBLAS_NUM_THREADS"),
    "cl": float(f["cl"]), "converged": bool(r.get("converged")),
    "residual": float(h[-1]), "n_steps": int(len(h)),
    "descent10": float(tail[0] / tail[-1]) if tail[-1] > 0 else float("inf"),
    "tail_monotone": bool(np.all(np.diff(tail) <= 0)),
    "n_limited": int(r.get("n_limited", 0)), "n_floored": int(r.get("n_floored", 0)),
    "accept_reason": str(r.get("accept_reason")),
}))
