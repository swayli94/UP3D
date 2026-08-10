"""LE-16: does the r_c "dead band" exist at all, once the budget is raised?

LE-14's failure-mode classifier retracted the dead band as I had described it:

    r_c 0      FAIL -> CLAMPING
    r_c 0.025  OK
    r_c 0.030  FAIL -> CLAMPING
    r_c 0.0375 FAIL -> BUDGET_LIMITED   descent10 = 2.002, still descending
    r_c 0.045  OK
    r_c 0.05   OK

So "[0.030, 0.0375] is a contiguous failure band" was TWO DIFFERENT DISEASES filed under one
label. 0.030 hits the limiters; 0.0375 merely ran out of Newton iterations while still
descending (the converging legs descend ~1e9 over their last ten steps). And my argument that
production r_c = 0.05 "sits on the safe side of a dead band" rested on the band being real.

This settles it the way LE-4's budget question was settled -- by MEASURING with a raised cap
rather than adjudicating. One knob, n_newton_max 60 -> 200, everything else the committed recipe.

  0.0375 converges  => the band does NOT exist as described. It is one clamping failure at 0.030
      sitting between successes, not a contiguous region, and the "safe side" argument stays
      withdrawn. Whether 0.030 is then an isolated clamping point or the edge of something
      narrower is a separate question, and 0.0325/0.035 are run to bound it.
  0.0375 still fails => the band is real after all and the classifier's budget_limited reading
      was a symptom of it, not the cause. Recorded either way.

0.030 is re-run at the raised cap too. It classified as CLAMPING, which a budget raise should NOT
fix -- so it doubles as the control: if a bigger budget "fixes" a clamping leg, the classifier's
mode assignment is what needs re-examining, not the band.

Outputs (TRACKED): bench/gate_results/le16_deadband_budget.csv
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
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.post.unified import wall_forces                         # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from run_le14_common_root import classify_failure                   # noqa: E402

#: a subset run writes its own file so the committed 4-leg baseline is not overwritten
CSV = os.path.join(_GATE,
                   "le16_deadband_budget.csv" if not os.environ.get(
                       "PYFP3D_LE16_ONLY") else "le16_deadband_sigmafix.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
M_TARGET, ALPHA = 0.88, 3.06
N_NEWTON_RAISED = 200          # the committed CONF_RAMP_NK value is 60
#: (r_c, what LE-14 classified it as)
LEGS = [(0.0375, "budget_limited"), (0.030, "clamping"),
        (0.0325, "untested"), (0.035, "untested")]
#: ★ 2026-08-05: set PYFP3D_LE16_ONLY to a comma-separated r_c list to re-run a subset.
#: Added for the sigma-transport re-measurement: the fix
#: (docs/dev_phase_two/20260805-0200-sigma-transport-root-cause.md) can only flip the
#: verdict of a leg that was refused ON SIGMA, and 0.0325 is the only one -- the other
#: three carry clamping or limit-cycle evidence, so for them the fix is a ~1e-11
#: reassociation. Re-running all four costs 3.4 h; this leg alone is ~27 min.
_only = os.environ.get("PYFP3D_LE16_ONLY", "")
if _only:
    _want = {float(x) for x in _only.split(",")}
    LEGS = [lg for lg in LEGS if lg[0] in _want]
KEYS = ["r_c", "le14_mode", "n_newton_max", "converged", "m_attained",
        "failure_mode", "mode_evidence", "n_limited", "n_floored", "res_final",
        "descent10", "n_newton", "cl_p", "wall_s", "note"]


def main():
    print(f"wing-body medium, M{M_TARGET} / alpha {ALPHA}, "
          f"n_newton_max {N_NEWTON_RAISED} (committed 60)")
    print("question: is the r_c dead band real, or was 0.0375 just out of "
          "iterations?\n")
    rows = []
    for rc, le14 in LEGS:
        mc, wc = cut_wake(read_mesh(MP))
        t = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", rc * B_SEMI)
        t0 = time.perf_counter()
        try:
            seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART,
                                        alpha_deg=ALPHA, **cap.CONF_SEED_KW)
            nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
                      tip_taper=t, phi_init=seed["phi"],
                      gamma_init=seed["gamma"], n_picard_seed=0,
                      n_newton_max=N_NEWTON_RAISED)
            r = solve_newton_transonic(mc, wc, m_inf=M_TARGET, alpha_deg=ALPHA,
                                       m_start=cap.WB_MSTART, dm=cap.DM,
                                       dm_min=0.01, freeze_tol=1e-5,
                                       intermediate_tol=1e-4, newton_kw=nk)
        except Exception as exc:                                   # noqa: BLE001
            wall = time.perf_counter() - t0
            print(f"  r_c {rc:<7} DIED {type(exc).__name__}: {str(exc)[:60]} "
                  f"({wall:.0f}s)", flush=True)
            rows.append(dict(r_c=rc, le14_mode=le14,
                             n_newton_max=N_NEWTON_RAISED, converged=False,
                             wall_s=round(wall, 1),
                             note=f"{type(exc).__name__}: {exc}"))
            continue
        wall = time.perf_counter() - t0
        hist = np.asarray(r.get("residual_history", []), dtype=float)
        nlim = int(r.get("n_limited") or 0)
        nflr = int(r.get("n_floored") or 0)
        mode, ev, d10, _rev = classify_failure(
            hist, np.asarray(r.get("clamp_history", []), dtype=float),
            np.asarray(r.get("F_history", []), dtype=float),
            int(r.get("n_gmres_stalled") or 0), str(r.get("accept_reason")),
            nlim, nflr)
        m_att = float(r.get("m_last_converged", r.get("m_final", M_TARGET)))
        conv = bool(r["converged"]) and abs(m_att - M_TARGET) < 1e-9
        #: ★ a CONVERGED leg has no failure mode. Without this the row reads
        #: "converged=True, failure_mode=limit_cycle": the classifier is a
        #: failure classifier and its tail heuristics say something arbitrary
        #: about a resolved trajectory. Measured on r_c = 0.0325 after the
        #: sigma-transport fix (docs/dev_phase_two/
        #: 20260805-0200-sigma-transport-root-cause.md), which flipped that leg
        #: to converged and left the label behind.
        if conv:
            mode, ev = "", f"converged; accept_reason={r.get('accept_reason')!r}"
        sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        cl = wall_forces(mc, phi=np.asarray(r["phi"]), alpha_deg=ALPHA,
                        s_ref=sref, m_inf=M_TARGET)["cl"]
        print(f"  r_c {rc:<7} LE-14={le14:15s} conv={conv!s:5s} "
              f"n_newton={len(hist):3d} lim/flr={nlim}/{nflr} "
              f"|R|={hist[-1] if len(hist) else float('nan'):.3e} "
              f"MODE={mode:16s} cl_p {cl:.6f} ({wall:.0f}s)", flush=True)
        if ev:
            print(f"           {ev}", flush=True)
        rows.append(dict(r_c=rc, le14_mode=le14, n_newton_max=N_NEWTON_RAISED,
                         converged=conv, m_attained=m_att, failure_mode=mode,
                         mode_evidence=ev, n_limited=nlim, n_floored=nflr,
                         res_final=(float(hist[-1]) if len(hist) else None),
                         descent10=(round(d10, 4) if d10 == d10 else None),
                         n_newton=len(hist), cl_p=cl, wall_s=round(wall, 1),
                         note=""))
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    print("\n=== verdict on the dead band ===")
    d = {r["r_c"]: r for r in rows}
    a, b = d.get(0.0375), d.get(0.030)
    if a is None:
        print("  0.0375 did not run"); return 0
    if b is not None and b["converged"]:
        print("  ★ the CONTROL moved: 0.030 was classified CLAMPING and a budget raise")
        print("  should not fix that. The classifier's mode assignment is what needs")
        print("  re-examining, not the band. Nothing below is safe to read yet.")
        return 0
    if a["converged"]:
        mids = [d[k] for k in (0.0325, 0.035) if k in d]
        got = ", ".join(f"{m['r_c']}:{'OK' if m['converged'] else m['failure_mode']}"
                        for m in mids)
        print("  the band does NOT exist as described: 0.0375 converges once the")
        print(f"  iteration cap is raised, leaving 0.030 as a CLAMPING point between")
        print(f"  successes rather than a contiguous region.  interior: {got}")
        print("  The 'production sits on the safe side of a dead band' argument stays")
        print("  withdrawn -- there was no band to sit beside.")
    else:
        print(f"  the band IS real: 0.0375 still fails at n_newton_max "
              f"{N_NEWTON_RAISED} (mode {a['failure_mode']}), so LE-14's")
        print("  budget_limited reading was a symptom rather than the cause.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
