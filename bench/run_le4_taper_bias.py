"""LE-4: cut the tip taper's systematic error without losing the robustness it exists for.

Registered in this docstring before the first run.

CONTEXT the user supplied: the taper was introduced as a SIMPLIFIED MODEL, to stop the solver
diverging when tip circulation falls and the tip vortex cannot be represented by a flat
terminating sheet. It buys robustness and pays a systematic error. So its acceptance criterion
is the divergence it was built to prevent -- NOT a smaller bias. A config with less bias that
loses the high-Mach convergence is a regression, however good its cl looks.

What the implementation already says (pyfp3d/constraints/wake.py::tip_taper_factors): the edge
singularity is driven by the TRAILING vorticity gamma = -dGamma/dz, not bound Gamma. With
Gamma ~ C*sqrt(u), u = z_tip - z, a taper F ~ u^s gives gamma_eff ~ u^(s-1/2), so the edge is
regularised IFF s > 1/2. The forms bracket that exponent:

    vanish_sqrt    s = 1/2   borderline, gamma finite but nonzero -- not enough
    vanish_linear  s = 1     the docstring's own "minimal-bias form that should regularize"
    vanish_smooth  s = 2     "regularizes harder, unloads the tip MORE"   <-- PRODUCTION

So production runs the HARDEST-unloading of the regularising forms, while the code's own
documentation names s = 1 as the minimal-bias one that still works. And the bias scales as
O((r_c/b)^{3/2}), so r_c is a second, independent lever: halving it should cut the bias ~2.8x.

Two levers, one constraint. Configs (wing-body conforming, where the taper is production):

    none                    the UNBIASED reference -- what the Kutta rows want with no taper
    vanish_smooth  0.05     production baseline (B32; measured cost about -1.3 % cl_p)
    vanish_linear  0.05     s = 1 at the same width
    vanish_smooth  0.025    s = 2 at half width
    vanish_linear  0.025    s = 1 at half width -- least bias of the four

Measured per config:
  (a) BIAS         cl_p at M0.50 against `none`. This is the systematic error being cut.
  (b) ROBUSTNESS   the M0.84 Mach ramp: converged? how many clamps? THE binding criterion,
                   because it is what the taper was introduced for (B31 GB31.2a measured the
                   taper causal for curing the 0.83 dying level).
  (c) CONVERGENCE  tip-band R via the LE-2 instrument -- recorded, not decisive, since LE-3
                   already measured the taper contributes ~nothing there (1.270 -> 1.264).

VERDICT RULE, fixed now: adopt the smallest-|bias| config that still passes (b) at the
production level. If only vanish_smooth 0.05 passes (b), the bias is the PRICE OF ROBUSTNESS and
that is the answer -- a negative result, recorded, not tuned around.

Coarse first because it is cheap; the winner is confirmed at medium before any adoption.

Outputs (TRACKED): bench/gate_results/le4_taper_bias.csv
"""

import csv
import math
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.post.unified import wall_forces                         # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)

CSV = os.path.join(HERE, "gate_results", "le4_taper_bias.csv")
MDIR = "onera_m6_wingbody_conforming"
M_SUB, M_TRANS, ALPHA = 0.50, 0.84, 3.06

CONFIGS = [("none", 0.0), ("vanish_smooth", 0.05), ("vanish_linear", 0.05),
           ("vanish_smooth", 0.025), ("vanish_linear", 0.025)]
KEYS = ["level", "form", "r_c", "leg", "m_inf", "converged", "n_limited",
        "n_floored", "res_final", "m_attained", "m_max", "cl_p", "bias_pct",
        "wall_s", "note"]


def taper_for(wc, form, r_c):
    if form == "none":
        return None
    return tip_taper_factors(wc.station_z, B_SEMI, form, r_c * B_SEMI)


def sub_leg(mesh_path, form, r_c):
    """(a) the bias: a single subsonic solve, same recipe as production apart from the taper."""
    mc, wc = cut_wake(read_mesh(mesh_path))
    kw = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure")
    t = taper_for(wc, form, r_c)
    if t is not None:
        kw["tip_taper"] = t
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_SUB, alpha_deg=ALPHA, **kw)
    return mc, wc, r, time.perf_counter() - t0


def trans_leg(mesh_path, form, r_c):
    """(b) the robustness: the production M0.84 ramp, seeded exactly as conf_wingbody does."""
    mc, wc = cut_wake(read_mesh(mesh_path))
    t = taper_for(wc, form, r_c)
    t0 = time.perf_counter()
    seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART, alpha_deg=ALPHA,
                                **cap.CONF_SEED_KW)
    nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure")
    if t is not None:
        nk["tip_taper"] = t
    r = solve_newton_transonic(mc, wc, m_inf=M_TRANS, alpha_deg=ALPHA,
                               m_start=cap.WB_MSTART, dm=cap.DM, dm_min=0.01,
                               freeze_tol=1e-5, intermediate_tol=1e-4,
                               phi_init=np.asarray(seed["phi"]),
                               gamma_init=np.asarray(seed["gamma"]),
                               newton_kw=nk)
    return mc, wc, r, time.perf_counter() - t0


def append(row):
    head = not os.path.exists(CSV)
    with open(CSV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        if head:
            w.writeheader()
        w.writerow(row)


def main():
    levels = os.environ.get("PYFP3D_LE4_LEVELS", "coarse").split(",")
    for level in levels:
        mp = os.path.join(REPO, "cases", "meshes", MDIR, f"{level}.msh")
        if not os.path.exists(mp):
            print(f"{level}: mesh missing -- {mp}"); continue
        print(f"\n########## {level} ##########", flush=True)
        base = None
        for form, r_c in CONFIGS:
            tag = f"{form}/{r_c}"
            # ---- (a) bias --------------------------------------------------
            try:
                mc, wc, r, wall = sub_leg(mp, form, r_c)
                phi = np.asarray(r["phi"])
                sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
                cl = wall_forces(mc, phi=phi, alpha_deg=ALPHA, s_ref=sref,
                                 m_inf=M_SUB)["cl"]
                if form == "none":
                    base = cl
                bias = (100.0 * (cl - base) / base) if base else float("nan")
                print(f"  (a) {tag:22s} M{M_SUB}  cl_p {cl:.8f}  "
                      f"bias {bias:+.3f} %  conv={bool(r['converged'])} "
                      f"({wall:.0f}s)", flush=True)
                append(dict(level=level, form=form, r_c=r_c, leg="bias",
                            m_inf=M_SUB, converged=bool(r["converged"]),
                            n_limited=r.get("n_limited"),
                            n_floored=r.get("n_floored"),
                            res_final=float(r["residual_history"][-1]),
                            cl_p=cl, bias_pct=bias, wall_s=round(wall, 1),
                            note=""))
            except Exception as exc:                               # noqa: BLE001
                print(f"  (a) {tag:22s} ERROR {type(exc).__name__}: "
                      f"{str(exc)[:70]}", flush=True)
                append(dict(level=level, form=form, r_c=r_c, leg="bias",
                            note=f"{type(exc).__name__}: {exc}"))
            # ---- (b) robustness -- the binding criterion --------------------
            try:
                mc, wc, r, wall = trans_leg(mp, form, r_c)
                phi = np.asarray(r["phi"])
                sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
                cl = wall_forces(mc, phi=phi, alpha_deg=ALPHA, s_ref=sref,
                                 m_inf=M_TRANS)["cl"]
                m_att = float(r.get("m_last_converged", r.get("m_final", M_TRANS)))
                conv = bool(r["converged"]) and abs(m_att - M_TRANS) < 1e-9
                print(f"  (b) {tag:22s} M{M_TRANS}  conv={conv} "
                      f"m_att={m_att} lim/flr={r.get('n_limited')}/"
                      f"{r.get('n_floored')} |R|={float(r['residual_history'][-1]):.2e} "
                      f"cl_p {cl:.8f} ({wall:.0f}s)", flush=True)
                append(dict(level=level, form=form, r_c=r_c, leg="robust",
                            m_inf=M_TRANS, converged=conv,
                            n_limited=r.get("n_limited"),
                            n_floored=r.get("n_floored"),
                            res_final=float(r["residual_history"][-1]),
                            m_attained=m_att, cl_p=cl, wall_s=round(wall, 1),
                            note=""))
            except Exception as exc:                               # noqa: BLE001
                print(f"  (b) {tag:22s} M{M_TRANS} FAILED "
                      f"{type(exc).__name__}: {str(exc)[:70]}", flush=True)
                append(dict(level=level, form=form, r_c=r_c, leg="robust",
                            m_inf=M_TRANS, converged=False,
                            note=f"{type(exc).__name__}: {exc}"))
    print(f"\nwrote {CSV}")
    print("\n=== verdict rule (fixed before running) ===")
    print("  adopt the smallest-|bias| config that still CONVERGES at M0.84.")
    print("  if only vanish_smooth/0.05 converges, the bias is the price of")
    print("  robustness -- a recorded negative, not something to tune around.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
