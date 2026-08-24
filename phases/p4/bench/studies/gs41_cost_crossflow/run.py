"""GS4.1 round 2 -- (3) cost measure at S2/S4 station counts, (4) crossflow record.

Binding text: phases/p4/docs/dev_phase_four/20260819-0100-gs41-cost-and-crossflow-prereg.md
(committed before this script existed). Codes C-SUB1 / C-PARTIAL / C-MISS and
the guards are quoted from it and are NOT re-specified here.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_cost_crossflow/run.py

★ REGISTERED LIMITATION (pre-registration section 1.2), carried verbatim into
every reading below: this measures the STRIP SOLVER at the right problem size.
It is NOT the whole VII pipeline on a real M6 field -- the cached phi died with
a previous session's scratchpad, and the committed wing state has no gamma.
Do not quote any number here as "whole-wing measured".
"""

import csv
import os
import statistics
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

REF_CSV = os.path.join(ROOT, "bench", "gate_results", "task3_s4_premise_check.csv")
#: pre-registration 2.1: registered BEFORE the measurement, never tuned to the
#: seconds that came out. n_sub = 2 is binding; 1 and 4 are reported.
N_SUB = (1, 2, 4)
N_SUB_BINDING = 2
RHO, MU = 1.0, 1.0e-5

SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------

def guard_frozen():
    paths = ["pyfp3d/viscous/ibl3.py", "pyfp3d/viscous/closures.py"]
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--"] + paths,
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  ibl3.py + closures.py unchanged vs HEAD: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed -- kill criterion 2/3 fires")


def guard_no_solve():
    import pyfp3d.solve.newton as N
    import pyfp3d.solve.picard as PC

    def _forbidden(name):
        def _f(*a, **k):
            raise AssertionError(f"G-NOSOLVE: {name} was CALLED -- round is "
                                 "registered zero-solve")
        return _f

    n = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            if not hasattr(mod, nm):
                raise SystemExit(f"G-NOSOLVE: {nm} missing -- library drift")
            setattr(mod, nm, _forbidden(nm))
            n += 1
    print(f"  G-NOSOLVE {n} solver entry points replaced by raising stubs -- a "
          "solve is IMPOSSIBLE, not merely detected  PASS")


def guard_ref():
    """G-REF: read the reference from the committed CSV, never from memory --
    and print the phase-3 prose value beside it so the 2.4 % fork stays visible.
    """
    rows = {r["level"]: r for r in csv.DictReader(open(REF_CSV))}
    ref = {}
    for lvl in ("S2", "S4"):
        r = rows[lvl]
        ref[lvl] = {
            "n_surf": int(r["n_node_surf"]),
            "bl_unknowns": int(r["bl_unknowns"]),
            "n_tet": int(r["n_tet"]),
            "ibl_s": float(r["ibl_s"]),
            "inviscid_s": float(r["inviscid_s"]),
            "p1": float(r["p1_ratio"]),
        }
        # cross-check the pre-registration's arithmetic on committed numbers
        assert ref[lvl]["bl_unknowns"] == 6 * ref[lvl]["n_surf"], lvl
        print(f"  G-REF     {lvl}: {ref[lvl]['n_surf']} surface nodes, "
              f"{ref[lvl]['bl_unknowns']} BL unknowns (= 6x, checked), "
              f"IBL {ref[lvl]['ibl_s']:.2f} s, inviscid "
              f"{ref[lvl]['inviscid_s']:.1f} s, P1 {ref[lvl]['p1']:.3f}")
    rep = rows.get("S2r")
    if rep:
        spread = abs(float(rep["ibl_s"]) - ref["S2"]["ibl_s"]) / ref["S2"]["ibl_s"]
        print(f"  G-REF     S2 repeat leg {float(rep['ibl_s']):.2f} s -> "
              f"reproducibility band {100*spread:.2f} %")
        ref["repro"] = spread
    print("  G-REF     ★ phase-3 VERDICT PROSE says 133.5 s / 9.96x -- that is "
          "2.4 % from its own CSV above. The CSV is this round's referent; the "
          "two must not be mixed (pre-registration 1.1).")
    return ref


# ---------------------------------------------------------------------------
# (3) cost
# ---------------------------------------------------------------------------

def _ue_families():
    """Three edge-velocity distributions for G-INDEP.

    Flat plate, favourable and adverse gradients -- chosen so the marched states
    differ substantially while the step COUNT is identical by construction.
    """
    from pyfp3d.viscous import strip2d as S
    return (("flat", S.flat_plate_ue(1.0), 2.591100),
            ("favourable", S.falkner_skan_ue(1.0 / 23.0), 2.495664),
            ("adverse", S.falkner_skan_ue(-0.02), 2.75))


def _one_strip(n_station, n_sub, ue_fn, target_H, x0=0.05, x1=1.0):
    """One chordwise strip with `n_station` stations, marched.

    Cost model of the round: total RK4 work is n_station * n_sub steps, four
    closure evaluations each. Returns (wall_seconds, ok).
    """
    from pyfp3d.viscous import strip2d as S
    ue0 = ue_fn(x0)[0]
    theta0 = 0.664 * x0 / np.sqrt(RHO * ue0 * x0 / MU)
    try:
        y0 = S.similar_seed(theta0, target_H, ue=ue0, rho=RHO, mu=MU)
    except ValueError:
        return float("nan"), False
    stations = np.geomspace(x0 * 1.2, x1, n_station)
    t0 = time.perf_counter()
    st = S.march(stations, y0, x0, ue_fn, rho=RHO, mu=MU,
                 n_substep=n_station * n_sub)
    wall = time.perf_counter() - t0
    ok = bool(np.all(np.isfinite(st.H)) and np.all(st.H > 1.0)
              and np.all(st.H < 10.0))
    return wall, ok


def gate_cost(ref):
    print("== (3) cost: strip core at the committed S2/S4 station counts ==")

    # --- G-JIT: warm up so the timer measures the solver, not the compiler ---
    t0 = time.perf_counter()
    _one_strip(40, 2, _ue_families()[0][1], 2.591100)
    cold = time.perf_counter() - t0
    t0 = time.perf_counter()
    _one_strip(40, 2, _ue_families()[0][1], 2.591100)
    warm = time.perf_counter() - t0
    print(f"  G-JIT     cold {cold:.4f} s vs warm {warm:.4f} s on an identical "
          f"40-station strip -> compile share {100*(1-warm/cold):.1f} %; every "
          "number below is WARM  PASS")

    # --- G-INDEP: does the cost depend on u_e? -------------------------------
    probe_n = 2000
    indep = []
    for name, ue_fn, tH in _ue_families():
        reps = [_one_strip(probe_n, N_SUB_BINDING, ue_fn, tH) for _ in range(3)]
        walls = [w for w, ok in reps if ok]
        if not walls:
            print(f"  G-INDEP   {name}: no converged leg")
            continue
        indep.append((name, statistics.median(walls), len(walls)))
        print(f"  G-INDEP   {name}: median {statistics.median(walls):.4f} s "
              f"over {len(walls)} reps ({probe_n} stations, n_sub="
              f"{N_SUB_BINDING})")
    if len(indep) < 2:
        _record("C-COST", "converged u_e families", ">= 2", f"{len(indep)}",
                "UNDEFINED")
        return None
    med = [m for _, m, _ in indep]
    spread = (max(med) - min(med)) / statistics.mean(med)
    ok_indep = spread <= 0.05
    _record("G-INDEP", "relative spread of cost across u_e families", "<= 5 %",
            f"{100*spread:.2f} % ({len(indep)} families)",
            "PASS" if ok_indep else "FAIL -> kill 1")
    if not ok_indep:
        _record("C-COST", "cost projection validity", "G-INDEP must pass",
                "u_e-dependent", "UNDEFINED")
        return None

    # --- the measurement at the committed station counts ---------------------
    rows = []
    ue_fn, tH = _ue_families()[0][1], _ue_families()[0][2]
    for lvl in ("S2", "S4"):
        n_st = ref[lvl]["n_surf"]
        for n_sub in N_SUB:
            reps = [_one_strip(n_st, n_sub, ue_fn, tH) for _ in range(3)]
            walls = [w for w, ok in reps if ok]
            if not walls:
                print(f"  {lvl} n_sub={n_sub}: all legs failed")
                continue
            t = statistics.median(walls)
            sp = (max(walls) - min(walls)) / t if t > 0 else 0.0
            rows.append({
                "level": lvl, "n_station": n_st, "n_sub": n_sub,
                "n_rk4_steps": n_st * n_sub, "t_strip_s": t,
                "rep_spread": sp, "n_rep": len(walls),
                "t_ibl3_s": ref[lvl]["ibl_s"],
                "t_inviscid_s": ref[lvl]["inviscid_s"],
                "p1_strip": t / ref[lvl]["inviscid_s"],
                "p1_ibl3": ref[lvl]["p1"],
                "speedup_vs_ibl3": ref[lvl]["ibl_s"] / t,
            })
            print(f"  G-RES     {lvl}: {n_st} stations x n_sub={n_sub} "
                  f"({n_st*n_sub} RK4 steps) -> {t:.4f} s "
                  f"(median of {len(walls)}, spread {100*sp:.1f} %)")

    with open(os.path.join(RESULTS, "cost.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # --- the verdict, on the binding n_sub at S4 -----------------------------
    b = [r for r in rows if r["level"] == "S4" and r["n_sub"] == N_SUB_BINDING]
    if not b:
        _record("C-COST", "binding S4 leg", "must exist", "missing", "UNDEFINED")
        return rows
    t = b[0]["t_strip_s"]
    code = ("C-SUB1" if t < 1.0 else
            "C-PARTIAL" if t <= 100.0 else "C-MISS")
    _record("C-COST", f"t_strip(S4) at n_sub={N_SUB_BINDING} "
            f"({b[0]['n_station']} stations)",
            "<1 s = C-SUB1; 1-100 s = C-PARTIAL; >100 s = C-MISS",
            f"{t:.4f} s", code)
    _record("C-RATIO", "P1_strip = t_strip / t_inviscid at S4 (vs committed "
            "P1_IBL3 = 6.073)", "RECORDED, not a criterion",
            f"{b[0]['p1_strip']:.5f}  ({b[0]['speedup_vs_ibl3']:.0f}x faster "
            f"than the 3-D IBL's {b[0]['t_ibl3_s']:.0f} s)", "RECORDED")

    # the gate's own number, across the registered resolution ladder
    ladder = {r["n_sub"]: r["t_strip_s"] for r in rows if r["level"] == "S4"}
    print("  ★ the gate's <1 s against the RESOLUTION knob: " +
          ", ".join(f"n_sub={k}: {v:.3f} s" for k, v in sorted(ladder.items())))
    return rows


def profile_closure_share():
    """Prediction 4: is the closure quadrature the cost? This decides whether
    (a2)'s cost argument stands at all."""
    print("== profile: where the strip core's time goes (RECORDED) ==")
    from pyfp3d.viscous import closures as C

    n = 20000
    st = np.array([2.0e-2, 7.9, 0.0, 0.0, C.CTAU_LAM, 0.0])
    t0 = time.perf_counter()
    for _ in range(n):
        C.closure_scalar(st, q=1.0, rho=RHO, mu=MU, turbulent=False)
    t_clo = (time.perf_counter() - t0) / n

    M = np.array([[1.0, 0.3], [0.2, 1.1]])
    rhs = np.array([1.0, 2.0])
    t0 = time.perf_counter()
    for _ in range(n):
        np.linalg.solve(M, rhs)
    t_sol = (time.perf_counter() - t0) / n

    share = 4 * t_clo / (4 * t_clo + t_sol)
    print(f"  one closure_scalar   {1e6*t_clo:8.2f} us  (x4 per RK4 step)")
    print(f"  one 2x2 linalg.solve {1e6*t_sol:8.2f} us  (x1 per RK4 step)")
    _record("profile", "closure quadrature share of one RK4 step",
            "RECORDED, no gate", f"{100*share:.1f} %", "RECORDED")
    return {"t_closure_us": 1e6 * t_clo, "t_solve_us": 1e6 * t_sol,
            "closure_share": share}


# ---------------------------------------------------------------------------
# (4) crossflow record -- committed evidence, zero recompute
# ---------------------------------------------------------------------------

def record_crossflow(ref):
    print("== (4) crossflow: what the 2-D strip DISCARDS (RECORDED, no gate) ==")
    rows = {r["level"]: r for r in csv.DictReader(open(REF_CSV))}
    out = []
    for lvl, xc_lo in (("S2", 0.982), ("S4", 0.9712)):
        r = rows[lvl]
        p2 = float(r["p2_ratio"])
        pmax = float(r["p2_max"])
        # consequence, exact arithmetic: |u_e|^2 = |A|^2 + |B|^2
        d_mag = np.sqrt(1.0 + pmax ** 2) - 1.0
        d_ang = np.degrees(np.arctan(pmax))
        out.append({
            "level": lvl, "max_B_over_max_A": p2, "pointwise_max": pmax,
            "median": float(r["p2_median"]), "p90": float(r["p2_p90"]),
            "top1pct_xc_from": xc_lo, "top1pct_xc_to": 1.0,
            "n_live": int(r["n_live"]), "n_tip_masked": int(r["n_tip_masked"]),
            "ue_magnitude_error_pct": 100 * d_mag,
            "flow_angle_deg_at_worst_node": d_ang,
        })
        print(f"  {lvl}: max|B|/max|A| = {p2:.6f}, pointwise max {pmax:.6f}, "
              f"top-1% |B| in x/c [{xc_lo}, 1.000]")
        print(f"       consequence at the WORST node: |u_e| low by "
              f"{100*d_mag:.4f} %, flow angle off by {d_ang:.2f} deg")
    with open(os.path.join(RESULTS, "crossflow_record.csv"), "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    s4 = out[-1]
    _record("crossflow", "S4 discarded crossflow, magnitude consequence",
            "RECORDED, no gate",
            f"{s4['ue_magnitude_error_pct']:.4f} % in |u_e|", "RECORDED")
    _record("crossflow", "S4 discarded crossflow, DIRECTION consequence",
            "RECORDED, no gate",
            f"{s4['flow_angle_deg_at_worst_node']:.2f} deg at the worst node",
            "RECORDED")
    print("  ★★ REFUSED CONCLUSION (pre-registration 3): these do NOT license "
          "'crossflow is negligible'. J = 0.600 vs null 0.060 = 10x makes it a "
          "REAL trailing-edge structure; small amplitude is not small influence "
          "for a quantity driving secondary flow and transition. Settling that "
          "needs a coupled measurement this round does not have.")
    return out


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t_all = time.perf_counter()
    print("== guards ==")
    guard_frozen()
    guard_no_solve()
    ref = guard_ref()

    rows = gate_cost(ref)
    prof = profile_closure_share()
    cf = record_crossflow(ref)

    total = time.perf_counter() - t_all
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  total round compute {total:.2f} s "
          "(kill criterion 4 budget: 20 min)")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:12s} [{tag}] {metric} = {measured}")
    print("\n★ REGISTERED LIMITATION: strip SOLVER at the right problem size, "
          "NOT the whole VII pipeline on a real M6 field (pre-registration 1.2).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
