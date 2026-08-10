"""LE-7: the decisive adoption leg -- does a NARROW taper hold the M0.88 envelope edge?

LE-5 flipped the verdict. At medium M0.88 / alpha 3.06 on the conforming wing-body:

    none/0.0             conv=False   |R| 8.00e-05   6 lim / 6 flr clamps   1357 s
    vanish_smooth/0.05   conv=True    |R| 1.78e-14   0 / 0                   467 s

So the taper still EARNS its bias at the envelope edge, and the pre-registered adoption rule --
`none` must converge wherever the taper does -- refuses blanket removal. The picture is scoped,
not uniform: the entropy correction covers the taper's robustness job at M <= 0.84 (LE-4's
single-knob A/B established that causally), but not at M0.88.

That reframes adoption from "remove it" to "make it cheaper", which is the same goal the user
set -- less systematic error -- reached without giving up envelope. LE-4 measured, at medium:

    vanish_smooth 0.05   bias -1.514 %   (production)
    vanish_smooth 0.025  bias -0.515 %   2.9x less, and it reached M0.84 at both levels
    vanish_linear 0.025  bias -0.575 %   nearly tied

The open question is whether either narrow form ALSO holds M0.88 medium. If yes, it delivers the
production envelope for a third of the systematic error, and THAT is the adoption candidate.
Both narrow forms are run because LE-4 left them nearly tied on bias, so the envelope decides.

Committed as a bench script rather than left in the scratchpad because it is decisive evidence
(discipline #3), and because the first attempt at it was lost: I ran `pkill -f
run_le5_taper_coverage` in the same compound command that wrote the script, and the pattern
matched my own shell's command line, so pkill killed the shell before the heredoc finished. Same
self-match family as using pgrep with a pattern that appears in the waiting command -- but
destructive rather than merely misleading.

Outputs (TRACKED): bench/gate_results/le7_narrow_taper.csv
"""

import csv
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

CSV = os.path.join(HERE, "gate_results", "le7_narrow_taper.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
M_TARGET, ALPHA = 0.88, 3.06
CANDIDATES = [("vanish_smooth", 0.025), ("vanish_linear", 0.025),
              ("vanish_smooth", 0.0375)]
KEYS = ["form", "r_c", "m_inf", "alpha", "converged", "m_attained", "n_limited",
        "n_floored", "res_final", "cl_p", "wall_s", "note"]


def main():
    print(f"M{M_TARGET} / alpha {ALPHA} / wing-body medium -- the envelope-edge leg")
    print("reference (LE-5): none = FAIL 6/6 clamps |R| 8.00e-05 ; "
          "smooth/0.05 = OK |R| 1.78e-14\n")
    rows = []
    for form, rc in CANDIDATES:
        mc, wc = cut_wake(read_mesh(MP))
        t = tip_taper_factors(wc.station_z, B_SEMI, form, rc * B_SEMI)
        t0 = time.perf_counter()
        try:
            seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART,
                                        alpha_deg=ALPHA, **cap.CONF_SEED_KW)
            nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure", tip_taper=t,
                      phi_init=seed["phi"], gamma_init=seed["gamma"],
                      n_picard_seed=0)
            r = solve_newton_transonic(mc, wc, m_inf=M_TARGET, alpha_deg=ALPHA,
                                       m_start=cap.WB_MSTART, dm=cap.DM,
                                       dm_min=0.01, freeze_tol=1e-5,
                                       intermediate_tol=1e-4, newton_kw=nk)
            sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
            cl = wall_forces(mc, phi=np.asarray(r["phi"]), alpha_deg=ALPHA,
                             s_ref=sref, m_inf=M_TARGET)["cl"]
            m_att = float(r.get("m_last_converged", r.get("m_final", M_TARGET)))
            conv = bool(r["converged"]) and abs(m_att - M_TARGET) < 1e-9
            wall = time.perf_counter() - t0
            print(f"  {form}/{rc:<7} conv={conv} m_att={m_att} "
                  f"lim/flr={r.get('n_limited')}/{r.get('n_floored')} "
                  f"|R|={float(r['residual_history'][-1]):.3e} cl_p {cl:.7f} "
                  f"({wall:.0f}s)", flush=True)
            rows.append(dict(form=form, r_c=rc, m_inf=M_TARGET, alpha=ALPHA,
                             converged=conv, m_attained=m_att,
                             n_limited=r.get("n_limited"),
                             n_floored=r.get("n_floored"),
                             res_final=float(r["residual_history"][-1]),
                             cl_p=cl, wall_s=round(wall, 1), note=""))
        except Exception as exc:                                   # noqa: BLE001
            wall = time.perf_counter() - t0
            print(f"  {form}/{rc:<7} DIED {type(exc).__name__}: "
                  f"{str(exc)[:70]} ({wall:.0f}s)", flush=True)
            rows.append(dict(form=form, r_c=rc, m_inf=M_TARGET, alpha=ALPHA,
                             converged=False, wall_s=round(wall, 1),
                             note=f"{type(exc).__name__}: {exc}"))
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    ok = [r for r in rows if r["converged"]]
    print("\n=== reading ===")
    if ok:
        best = min(ok, key=lambda r: r["r_c"])
        print(f"  narrowest form holding the M0.88 edge: {best['form']}/"
              f"{best['r_c']}  -> the adoption candidate, since LE-4 measured its")
        print("  bias at roughly a third of production's -1.514 %.")
    else:
        print("  NO narrow form holds M0.88 -- production vanish_smooth/0.05 is the")
        print("  minimum that buys the envelope, and its -1.514 % is the price.")
        print("  Recorded as a negative: the taper cannot be made cheaper without")
        print("  giving up the edge.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
