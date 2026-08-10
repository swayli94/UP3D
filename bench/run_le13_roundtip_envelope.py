"""LE-13: is the round-tip LS envelope regression a RECIPE calibration or a geometry cost?

Switching the wing-alone families to the round tip fixed a premise -- the LE band's convergence
order went from p = 0.37 (flat, and flat DIVERGES under refinement so an order on it is
meaningless) to p = 0.87, matching the 2.5-D NACA's 0.85. But the gated LS anchors moved, and one
of them is a capability regression rather than a number shift:

    test_m6_medium_ramp_anchor   flat: reaches M0.84    round: stops at M0.675
    test_m6_coarse_ramp_anchor   reaches M0.84 both ways, other values moved
    test_b7 coarse gate [M4]     clamps (0,0) -> (0,2)

A capability regression has to be reported as a result, not absorbed by re-anchoring. But before
it is charged to the geometry, the cheap alternative has to be excluded: this project has twice
found a stalled ramp to be a freeze/step CALIBRATION rather than a capability limit -- B17
("the conforming wing-body medium ramp needs freeze_tol raised to the wing-body churn floor
1e-6 -> 1e-5") and B32. The round tip changes the tip loading distribution, so the churn floor it
needs has no reason to be the one calibrated on the flat cap.

So: sweep the two knobs those precedents name, on the round-tip LS medium ramp.

  freeze_tol   1e-5 (the committed value), 1e-4, 1e-6
  dm           the committed value and half of it -- a finer Mach step is the other documented
               way past a level that a coarse step cannot cross

  any leg reaches M0.84  => CALIBRATION. The regression is not a geometry cost, and the anchor
      is re-measured with the corrected recipe rather than recorded as a capability loss.
  none does               => GEOMETRY COST, recorded honestly as M0.84 -> 0.675, with restoring
      the transonic envelope on the round tip becoming the named next target.

Outputs (TRACKED): bench/gate_results/le13_roundtip_envelope.csv
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
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te                      # noqa: E402
from pyfp3d.solve.newton_ls import (B_NEWTON_M6_DEFAULTS,           # noqa: E402
                                    solve_multivalued_newton_transonic)

CSV = os.path.join(HERE, "gate_results", "le13_roundtip_envelope.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wakefree", "medium.msh")
M_TARGET, ALPHA = 0.84, 3.06
#: (freeze_tol, dm)
LEGS = [(1e-5, None), (1e-4, None), (1e-6, None), (1e-5, "half")]
KEYS = ["freeze_tol", "dm", "m_final", "reached", "gamma", "m_max",
        "n_limited", "n_floored", "res_final", "wall_s", "note"]


def main():
    if not os.path.exists(MP):
        raise SystemExit(f"missing {MP}")
    print(f"round-tip LS M6 medium, target M{M_TARGET} / alpha {ALPHA}")
    print("reference: flat reached M0.84; round stops at M0.675 with the "
          "committed recipe\n")
    rows = []
    for ftol, dmmode in LEGS:
        #: ★ the ramp is called the way ls_wing calls it, READ from
        #: bench/run_capability_matrix.py:199-210 rather than through a wrapper I invented
        #: (my first draft called cap.ls_wing_ramp, which does not exist). ls_wing takes no
        #: knob arguments, so the two swept knobs are overridden on top of the committed
        #: B_NEWTON_M6_DEFAULTS -- freeze_tol 1e-3 and dm 0.05 there.
        base_dm = float(B_NEWTON_M6_DEFAULTS["dm"])
        dm = base_dm / 2.0 if dmmode == "half" else base_dm
        t0 = time.perf_counter()
        try:
            mesh = read_mesh(MP)
            te = np.array([[x_te(0.0), 0.0, 0.0],
                           [x_te(B_SEMI), 0.0, B_SEMI]])
            mvop = cap._ls_op(mesh, te, ALPHA)
            ramp = dict(B_NEWTON_M6_DEFAULTS, freeze_tol=ftol, dm=dm)
            r = solve_multivalued_newton_transonic(
                mvop=mvop, mesh=mesh, m_target=M_TARGET, alpha_deg=ALPHA,
                **cap.LS_WING_KW, **ramp)
        except Exception as exc:                                   # noqa: BLE001
            wall = time.perf_counter() - t0
            print(f"  freeze_tol {ftol:.0e} dm {dm:g}  DIED "
                  f"{type(exc).__name__}: {str(exc)[:60]} ({wall:.0f}s)",
                  flush=True)
            rows.append(dict(freeze_tol=ftol, dm=dm, reached=False,
                             wall_s=round(wall, 1),
                             note=f"{type(exc).__name__}: {exc}"))
            continue
        wall = time.perf_counter() - t0
        mf = float(r.get("m_final", float("nan")))
        lv = r["levels"][-1]
        reached = abs(mf - M_TARGET) < 1e-9
        print(f"  freeze_tol {ftol:.0e} dm {dm:g}  m_final={mf:.4f} "
              f"reached={reached} gamma={r.get('gamma_final')} "
              f"lim/flr={lv.get('n_limited')}/{lv.get('n_floored')} "
              f"|R|={float(lv['residual_norm']):.2e} ({wall:.0f}s)", flush=True)
        rows.append(dict(freeze_tol=ftol, dm=dm, m_final=mf, reached=reached,
                         gamma=r.get("gamma_final"),
                         m_max=r.get("mach2_max"),
                         n_limited=lv.get("n_limited"),
                         n_floored=lv.get("n_floored"),
                         res_final=float(lv["residual_norm"]),
                         wall_s=round(wall, 1), note=""))
    if rows:
        with open(CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        print(f"\nwrote {CSV}")
    print("\n=== reading ===")
    if any(r.get("reached") for r in rows):
        good = [r for r in rows if r.get("reached")]
        got = ", ".join(f"(ftol {r['freeze_tol']:.0e}, dm {r['dm']:g})"
                        for r in good)
        print(f"  CALIBRATION: M0.84 reached at {got}")
        print("  The regression is not a geometry cost -- the round tip needs its own")
        print("  churn floor, exactly as B17/B32 found for the wing-body. Re-anchor")
        print("  with the corrected recipe.")
    else:
        print("  GEOMETRY COST: no knob setting reaches M0.84 on the round tip.")
        print("  Recorded honestly as an envelope regression M0.84 -> 0.675, and")
        print("  restoring the round-tip transonic envelope becomes the named target.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
