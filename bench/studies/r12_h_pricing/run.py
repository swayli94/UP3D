"""R12 -- price the h route for M1, on the committed NACA0012 2.5-D ladder.

Binding text: docs/dev_phase_five/20260823-0300-r12-h-pricing-prereg.md (committed
before this file existed).

The h the mechanism argument asks for (h_shock <~ 0.005-0.01 c) already exists as the
fine level. So this measures (P-H) whether h_wall is a fair proxy for the spacing in
the shock band, and (P-CONV) whether the production recipe converges there at all --
all four levels in ONE session at ONE thread count, because the two committed readings
of medium disagree and were taken in different sessions.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/r12_h_pricing/run.py
"""

import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "bench"))
RESULTS = os.path.join(HERE, "results")

#: G-RECIPE -- imported, never re-typed
from run_gs40h_strength_and_2p5d import NACA_KW                         # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.post.shock import shock_report                              # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve                      # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients                 # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                    # noqa: E402

MESHES = ROOT + "/cases/meshes/naca0012_2.5d"
LEVELS = ("xcoarse", "coarse", "medium", "fine")
H_WALL = {"xcoarse": 0.040, "coarse": 0.020, "medium": 0.010, "fine": 0.005}
#: ★ s_ref = dz (chord 1 x span dz), the project's 2.5-D convention
#: (tests/test_p4_transonic.py:77), and dz is COMPUTED from the mesh the way that test
#: computes it rather than hardcoded. A fixed s_ref would be wrong per level -- dz spans
#: 8x across this ladder, so it would make the two-level cl difference meaningless,
#: which is precisely what P-B compares. Caught by the G-ARITY dry-check.
DZ_STAT = {"xcoarse": 0.08, "coarse": 0.04, "medium": 0.02, "fine": 0.01}  # G-DZ only
M_INF, ALPHA = 0.80, 1.25
SHOCK_BAND = (0.55, 0.65)          # P-H: the band M1's shock sits in
TOL = 0.0055                        # M1's (b)(c) requirement, chord
PER_LEVEL_S = 15 * 60               # kill criterion 3
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}", flush=True)


def h_in_shock_band(mc):
    """P-H: the actual wall spacing in the shock band, from the mesh alone.

    Upper-surface wall nodes with x/c in SHOCK_BAND, sorted by x; the spacing is the
    median gap. No solve, no field -- this is a property of the mesh.
    """
    wf = mc.boundary_faces["wall"]
    wn = np.unique(wf.reshape(-1))
    p = mc.nodes[wn]
    up = p[(p[:, 1] > 0.0) & (p[:, 0] >= SHOCK_BAND[0]) & (p[:, 0] <= SHOCK_BAND[1])]
    if len(up) < 3:
        return np.nan, 0
    #: one spanwise station only -- the mesh is an extrusion, so collapse z
    x = np.unique(np.round(up[:, 0], 9))
    return (float(np.median(np.diff(np.sort(x)))) if len(x) >= 3 else np.nan), len(x)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        print(f"  G-THREADS  {v} = {os.environ.get(v, '<unset>')}")
    print(f"  G-THREADS  load average {os.getloadavg()}")
    print(f"  G-RECIPE   NACA_KW = {NACA_KW}")
    print(f"  G-FROZEN-LIB  pyfp3d/ untouched; M {M_INF} alpha {ALPHA}", flush=True)

    rows, hrow = [], []
    for lv in LEVELS:
        p = f"{MESHES}/{lv}.msh"
        if not os.path.exists(p):
            print(f"  [{lv}] mesh absent"); continue
        mc, wc = cut_wake(read_mesh(p))
        dz = float(np.ptp(mc.nodes[:, 2]))
        #: G-DZ -- the computed span must match the committed mesh statistic
        assert abs(dz - DZ_STAT[lv]) < 1e-9, f"G-DZ {lv}: {dz} vs {DZ_STAT[lv]}"
        h_s, n_st = h_in_shock_band(mc)
        hrow.append({"level": lv, "h_wall": H_WALL[lv], "h_shock": h_s,
                     "n_band_nodes": n_st, "n_nodes": len(mc.nodes),
                     "n_elements": len(mc.elements)})
        print(f"  [{lv:8}] nodes {len(mc.nodes):6d}  h_wall {H_WALL[lv]:.3f}  "
              f"h_shock {h_s:.5f} ({n_st} band nodes)", flush=True)

        t1 = time.perf_counter()
        try:
            r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **NACA_KW)
        except Exception as e:                                          # noqa: BLE001
            rows.append({"level": lv, "status": f"RAISED: {type(e).__name__}",
                         "wall_s": round(time.perf_counter() - t1, 2)})
            print(f"           ★ raised {type(e).__name__}: {e}", flush=True); continue
        w = time.perf_counter() - t1
        phi = np.asarray(r["phi"], float)
        #: G-CACHE -- phi and the diagnostics land on disk BEFORE anything is reported
        np.savez_compressed(
            os.path.join(RESULTS, f"{lv}.npz"), phi=phi,
            conv=bool(r["converged"]), nlim=int(r["n_limited"]),
            nflr=int(r["n_floored"]), gamma=np.asarray(r.get("gamma", []), float),
            res_hist=np.asarray(r.get("residual_history", []), float),
            clamp_hist=np.asarray(r.get("clamp_history", []), float))
        row = {"level": lv, "h_wall": H_WALL[lv], "n_nodes": len(mc.nodes),
               "wall_s": round(w, 2), "conv": bool(r["converged"]),
               "n_limited": int(r["n_limited"]), "n_floored": int(r["n_floored"]),
               #: ★ NOT "residual" -- that key does not exist. G-ARITY caught the
               #: .get() default that would have silently reported NaN.
               #: ★★ And `residual_unfrozen` is None on a leg that never froze, which
               #: the dry-check could not see: it exercised ONE state, and the TYPE is
               #: state-dependent. float(None) then raised in the reporting layer AFTER
               #: the 11-minute fine solve -- the logged float(None) hazard, again.
               #: The npz cache written above is what saved it.
               "res": (float(r["residual_unfrozen"])
                       if r["residual_unfrozen"] is not None else None),
               "res_last": float(np.asarray(r["residual_history"], float)[-1]),
               "n_newton": int(r.get("n_newton", -1)),
               "accept_reason": str(r.get("accept_reason", "")),
               "h_shock": h_s}
        #: G-CLAMP -- a clamped state is NOT converged (GS1.4)
        row["usable"] = bool(r["converged"]) and not r["n_limited"] and not r["n_floored"]
        try:
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], phi,
                                        alpha_deg=ALPHA, m_inf=M_INF,
                                        s_ref=dz)
            row["cl_p"] = float(f["cl"])
            #: the 2.5-D path is wall_cp_curve at mid-span, exactly as
            #: tests/test_p4_transonic.py:73 does it -- section_cp_curve raises here
            #: ("section too sparse"), since the mesh spans z in [0, dz]
            cur = wall_cp_curve(mc, phi, z=0.5 * dz, m_inf=M_INF)
            xs = shock_report(cur, M_INF)["upper"].get("x_shock")
            row["x_shock"] = float(xs) if xs is not None and np.isfinite(xs) else np.nan
        except Exception as e:                                          # noqa: BLE001
            row["post_error"] = f"{type(e).__name__}: {e}"
        rows.append(row)
        print(f"           {w:7.1f}s  conv={row['conv']} clamps="
              f"{row['n_limited']}/{row['n_floored']}  usable={row['usable']}  "
              f"res={row['res']:.3e}  cl_p={row.get('cl_p', float('nan')):.5f}  "
              f"x_shock={row.get('x_shock', float('nan')):.5f}", flush=True)
        if w > PER_LEVEL_S:
            print(f"           ★ kill criterion 3: {lv} exceeded {PER_LEVEL_S}s; "
                  "the timeout IS the pricing reading. Stopping the ladder.", flush=True)
            break

    for name, rr in (("levels.csv", rows), ("h_band.csv", hrow)):
        if rr:
            keys = sorted({k for d in rr for k in d})
            with open(os.path.join(RESULTS, name), "w", newline="") as fh:
                w_ = csv.DictWriter(fh, fieldnames=keys)
                w_.writeheader(); w_.writerows(rr)

    # ---- P-H --------------------------------------------------------------
    dev = [abs(h["h_shock"] / h["h_wall"] - 1.0) for h in hrow
           if np.isfinite(h["h_shock"])]
    inband = [h["level"] for h in hrow
              if np.isfinite(h["h_shock"]) and h["h_shock"] <= 0.01]
    _record("P-H", "h_shock vs h_wall in the shock band (mesh only, no solve)",
            "<=25% => h_wall is a fair proxy",
            "  ".join(f"{h['level']} {h['h_shock']:.5f}/{h['h_wall']:.3f}"
                      f"={h['h_shock']/h['h_wall']:.2f}x" for h in hrow)
            + f";  levels with h_shock <= 0.01 c: {inband}",
            "P-H: h_wall is a fair proxy" if dev and max(dev) <= 0.25 else
            "★ P-H: h_wall is NOT a fair proxy -- read the h_shock ladder directly")

    # ---- P-CONV -----------------------------------------------------------
    got = {r["level"]: r for r in rows if "usable" in r}
    fine = got.get("fine")
    _record("P-CONV", "does the production recipe converge unclamped at each level",
            "fine usable => the h route is OPEN;  fine not usable => blocker confirmed "
            "on current code",
            "  ".join(f"{lv}: {'USABLE' if got[lv]['usable'] else 'no'}"
                      f"({got[lv]['n_limited']}/{got[lv]['n_floored']}, "
                      f"{got[lv]['wall_s']}s)" for lv in LEVELS if lv in got),
            ("P-CONV: fine USABLE -- the h route is open" if fine and fine["usable"]
             else "P-CONV: fine NOT usable -- blocker confirmed on current code"
             if fine else "P-CONV: fine did not run"))

    # ---- P-B (RECORDED ONLY -- fold zone, kill criterion 4) ---------------
    if got.get("medium", {}).get("usable") and fine and fine["usable"]:
        a, b = got["medium"]["cl_p"], fine["cl_p"]
        _record("P-B", "two-level cl difference (medium, fine) vs M1's (b) < 3%",
                "RECORDED ONLY -- M0.80 NACA medium is in the FOLD ZONE "
                "(dcl/dM ~ 6-10); discipline #4 forbids a grid-convergence claim here",
                f"medium {a:.5f} -> fine {b:.5f} = {100*(b-a)/abs(a):+.2f}%",
                "RECORDED (no gate verdict)")
    else:
        _record("P-B", "two-level cl difference (medium, fine)", "conditional on both "
                "levels being usable", "not both usable", "NOT APPLICABLE")

    _record("P-PRIOR", "prior blockers on the h route, with magnitudes", "RECORDED",
            f"x_shock thread scatter 0.027 c = {0.027/TOL:.1f}x the +-{TOL} c "
            f"requirement; mixed-family seed change 0.1683 c = {0.1683/TOL:.0f}x",
            "RECORDED -- the answer is not yet reproducible to the tolerance it "
            "would be judged against")

    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as fh:
        w_ = csv.writer(fh); w_.writerow(["tag", "metric", "band", "measured", "verdict"])
        w_.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s total", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
