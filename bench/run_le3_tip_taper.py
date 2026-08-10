"""LE-3: the wing-alone recipe is MISSING the tip cure the wing-body already carries.

LE-2 localised the LE-band convergence deficit to the TIP bin (eta > 0.8) on all four 3-D
combinations -- root and inboard bins are healthy (p 1.4-1.7), the tip is the site. That is the
known P13/B31 tip sheet-termination class.

Then reading the recipes turned up a RECIPE GAP, not a code bug: run_capability_matrix's
conf_wingbody passes tip_taper (B31/B32's cure, vanish_smooth at 0.05*b_semi) while conf_wing
does not. And the measured tip bins line up with that -- wing-body (with the cure) R 1.644,
wing-alone (without) R 1.983, which on its 2x ladder is above the falsification ceiling.

★ That also means LE-1's "wing-body healthy / wing-alone deficient" comparison was CONFOUNDED:
it differed by recipe as well as by geometry. Recorded, and this round separates the two.

B31 records that tip_taper must run WITH kutta_estimator="pressure" -- that pair is what
activates the FD-verified Gamma-pin row blend, whose pin slope is frozen with its own sign
(measured diag D > 0 on the conforming meshes, so an unsigned weld would AMPLIFY mid-taper
loading instead of unloading it). So the taper cannot be A/B'd on its own against the shipped
probe recipe: three legs are needed to isolate it.

  L1  probe (the shipped wing-alone recipe)                      = the baseline as measured
  L2  pressure, no taper                                         = the estimator's own effect
  L3  pressure + tip_taper                                       = the cure

Isolating the taper is L3 vs L2, NOT L3 vs L1.

Read per eta bin, since a tip cure must show up in the TIP bin specifically -- a change in the
pooled number could come from anywhere, and B32 measured the taper also costs about -1.3 % cl_p,
so the total moving is expected and is not evidence about convergence.

Outputs (TRACKED): bench/gate_results/le3_tip_taper.csv
"""

import csv
import math
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors      # noqa: E402
from run_le2_spanwise import ETA_EDGES, le_bins                     # noqa: E402

CSV = os.path.join(HERE, "gate_results", "le3_tip_taper.csv")
M_INF, ALPHA = 0.50, 3.06
LEVELS = ("xcoarse_ss", "coarse_ss", "medium")
HS = (0.060, 0.030, 0.015)
R_FIRST = (HS[2] - HS[1]) / (HS[1] - HS[0])
R_MAX = math.log(HS[2] / HS[1]) / math.log(HS[1] / HS[0])

LEGS = [("L1_probe", {}),
        ("L2_pressure", {"kutta_estimator": "pressure"}),
        ("L3_pressure_taper", {"kutta_estimator": "pressure", "_taper": True})]


def solve(mesh_path, extra):
    mc, wc = cut_wake(read_mesh(mesh_path))
    kw = dict(cap.NEWTON_M6_RECIPE["newton_kw"])
    e = dict(extra)
    if e.pop("_taper", False):
        kw["tip_taper"] = tip_taper_factors(wc.station_z, B_SEMI,
                                            cap.CONF_TAPER[0],
                                            cap.CONF_TAPER[1] * B_SEMI)
    kw.update(e)
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)
    return mc, wc, r


def main():
    rows = []
    print(f"m6wing conforming, M{M_INF} alpha {ALPHA}, ladder {HS} "
          f"(p=1 -> R {R_FIRST:.3f}, falsified >= {R_MAX:.3f})")
    per_leg = {}
    for tag, extra in LEGS:
        print(f"\n=== {tag} ===", flush=True)
        vals, ntri, cls = {}, {}, {}
        for lv in LEVELS:
            mp = os.path.join(REPO, "cases", "meshes", "onera_m6", f"{lv}.msh")
            if not os.path.exists(mp):
                print(f"  {lv}: mesh missing"); continue
            try:
                mc, wc, r = solve(mp, extra)
            except Exception as exc:                               # noqa: BLE001
                print(f"  {lv}: ERROR {type(exc).__name__}: {str(exc)[:80]}",
                      flush=True)
                rows.append(dict(leg=tag, level=lv, note=f"{type(exc).__name__}"))
                continue
            phi = np.asarray(r["phi"])
            sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
            b, extra_cl, le_tot, n_le = le_bins(mc, phi=phi, alpha_deg=ALPHA,
                                                m_inf=M_INF, s_ref=sref)
            g = abs(sum(v for v, _n in b.values()) + extra_cl - le_tot)
            conv = bool(r["converged"])
            print(f"  {lv:11s} conv={conv} n_LE={n_le:6d} LE_tot {le_tot:+.7f} "
                  f"guard {g:.1e}", flush=True)
            vals[lv] = b
            ntri[lv] = {k: n for k, (v, n) in b.items()}
            cls[lv] = le_tot
        if len(vals) < 3:
            continue
        lx, lc, lm = LEVELS
        per_leg[tag] = {}
        for k in vals[lx]:
            d1 = vals[lc][k][0] - vals[lx][k][0]
            d2 = vals[lm][k][0] - vals[lc][k][0]
            R = d2 / d1 if d1 else float("nan")
            per_leg[tag][k] = R
            print(f"    eta [{k[0]:.1f},{k[1]:.1f})  d1 {d1:+.7f}  d2 {d2:+.7f}  "
                  f"R {R:+.4f}", flush=True)
            rows.append(dict(leg=tag, level="R", eta_lo=k[0], eta_hi=k[1],
                             d1=round(d1, 10), d2=round(d2, 10),
                             R=(None if R != R else round(R, 5)),
                             n_tri_min=min(ntri[lv].get(k, 0) for lv in LEVELS),
                             R_first_order=round(R_FIRST, 4),
                             R_falsify_ceiling=round(R_MAX, 4), note=""))
        d1 = cls[lc] - cls[lx]; d2 = cls[lm] - cls[lc]
        print(f"    LE band total  d1 {d1:+.7f}  d2 {d2:+.7f}  "
              f"R {d2 / d1 if d1 else float('nan'):+.4f}", flush=True)
        rows.append(dict(leg=tag, level="R", eta_lo=-1, eta_hi=-1,
                         d1=round(d1, 10), d2=round(d2, 10),
                         R=round(d2 / d1, 5) if d1 else None,
                         note="LE band pooled"))

    print("\n=== the isolation that matters: L3 vs L2 (taper alone) ===")
    if "L2_pressure" in per_leg and "L3_pressure_taper" in per_leg:
        for k in sorted(per_leg["L2_pressure"]):
            a, b = per_leg["L2_pressure"][k], per_leg["L3_pressure_taper"][k]
            tip = " ← TIP" if k[0] >= 0.8 else ""
            print(f"  eta [{k[0]:.1f},{k[1]:.1f})  R {a:+.4f} -> {b:+.4f}"
                  f"   {'BETTER' if abs(b - R_FIRST) < abs(a - R_FIRST) else 'worse'}"
                  f"{tip}")
    else:
        print("  (a leg failed -- see rows above)")
    with open(CSV, "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
