"""LE-15: classify the round-tip LS envelope failure -- geometry cost, or something cheaper?

LE-14's failure-mode classifier changed a verdict that had already been written down, so it is
applied to the one failure that had not been classified yet.

What LE-14 established, and why it matters here:

  r_c 0      FAIL -> mode CLAMPING
  r_c 0.030  FAIL -> mode CLAMPING
  r_c 0.0375 FAIL -> mode BUDGET_LIMITED (descent10 = 2.002, still descending; converging legs
             run 1e9), i.e. NOT a failure at all -- it ran out of iterations

So the "dead band [0.030, 0.0375]" was TWO DIFFERENT DISEASES filed under one label, and my
argument that production r_c = 0.05 "sits on the safe side of a dead band" rests on a band that
may not exist. And none of the four candidate quantities separated the pattern, consistent with
the classifier: a single number cannot separate a mixture of diseases.

The round-tip LS medium regression (M0.84 -> M0.70) was recorded as a GEOMETRY COST on the same
kind of unclassified evidence -- four knob settings all "not reaching M0.84". If it is
budget_limited, that verdict has to be withdrawn the same way. If it is clamping, it is a genuine
geometry cost and stands.

One leg, the committed recipe, with the full diagnostic history captured and the LE-14 classifier
applied unchanged. Per-level histories are read too, because the LS ramp reports per level and a
ramp can stall at one level while earlier levels were healthy -- the level where it turns is part
of the diagnosis.

Outputs (TRACKED): bench/gate_results/le15_roundtip_mode.csv
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
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te                      # noqa: E402
from pyfp3d.solve.newton_ls import (B_NEWTON_M6_DEFAULTS,           # noqa: E402
                                    solve_multivalued_newton_transonic)
from run_le14_common_root import classify_failure                   # noqa: E402

CSV = os.path.join(_GATE, "le15_roundtip_mode.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wakefree", "medium.msh")
M_TARGET, ALPHA = 0.84, 3.06
KEYS = ["level_m", "converged_level", "failure_mode", "mode_evidence",
        "n_limited", "n_floored", "res_final", "descent10", "res_revisits",
        "n_newton", "note"]


def main():
    if not os.path.exists(MP):
        raise SystemExit(f"missing {MP} (regenerate the round-tip wake-free family)")
    print(f"round-tip LS M6 medium -> M{M_TARGET} / alpha {ALPHA}, committed recipe")
    print("LE-13 recorded this as a GEOMETRY COST on unclassified evidence; "
          "classifying it now\n")
    t0 = time.perf_counter()
    mesh = read_mesh(MP)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    mvop = cap._ls_op(mesh, te, ALPHA)
    r = solve_multivalued_newton_transonic(
        mvop=mvop, mesh=mesh, m_target=M_TARGET, alpha_deg=ALPHA,
        **cap.LS_WING_KW, **B_NEWTON_M6_DEFAULTS)
    wall = time.perf_counter() - t0
    mf = float(r.get("m_final", float("nan")))
    print(f"  m_final = {mf:.4f}  (target {M_TARGET})   {wall:.0f}s\n")

    rows = []
    #: per LEVEL -- a ramp can stall at one level with every earlier level healthy, so the
    #: level where it turns is part of the diagnosis, not a detail
    for lv in r["levels"]:
        m = float(lv.get("m_inf", float("nan")))
        hist = np.asarray(lv.get("residual_history", []), dtype=float)
        nlim = int(lv.get("n_limited") or 0)
        nflr = int(lv.get("n_floored") or 0)
        mode, ev, d10, rev = classify_failure(
            hist, np.array([]), np.array([]), 0, str(lv.get("accept_reason")),
            nlim, nflr)
        conv = bool(lv.get("converged", False))
        print(f"  level M{m:<7.4f} conv={conv!s:5s} MODE={mode:16s} "
              f"lim/flr={nlim}/{nflr} n_newton={len(hist)}", flush=True)
        if ev:
            print(f"                {ev}", flush=True)
        rows.append(dict(level_m=m, converged_level=conv, failure_mode=mode,
                         mode_evidence=ev, n_limited=nlim, n_floored=nflr,
                         res_final=(float(hist[-1]) if len(hist) else None),
                         descent10=(round(d10, 4) if d10 == d10 else None),
                         res_revisits=rev, n_newton=len(hist), note=""))
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    print("\n=== verdict on LE-13's 'geometry cost' ===")
    #: the diagnosis is the LAST level -- the one the ramp could not leave
    last = rows[-1] if rows else None
    if last is None:
        print("  no levels reported"); return 0
    print(f"  the ramp turned at M{last['level_m']:.4f}, mode = "
          f"{last['failure_mode']}")
    if last["failure_mode"] == "budget_limited":
        print("  => WITHDRAW the geometry-cost verdict. The round tip is not costing")
        print("  envelope; the ramp is running out of iterations at that level, and")
        print("  raising the budget is the test that decides the real envelope.")
    elif last["failure_mode"] == "clamping":
        print("  => GEOMETRY COST STANDS, and it is now specific: the round tip drives")
        print("  cells into the m_cap / rho_floor limiters at that level, which is a")
        print("  physical limit, not an iteration budget.")
    else:
        print(f"  => mode is {last['failure_mode']}, which is neither of the two")
        print("  readings LE-13 could have had. Recorded; the geometry-cost verdict")
        print("  is not supported as stated and needs this mode's own follow-up.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
