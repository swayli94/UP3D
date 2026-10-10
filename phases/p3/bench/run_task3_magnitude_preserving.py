"""Magnitude-preserving softening: separate "the jump was the cause" from "less correction".

Pre-registered in phases/p3/docs/dev_phase_three/20260812-1700-magnitude-preserving-prereg.md, committed
before this code -- a NEW registration and not an addendum, because it CHANGES a criterion (F1), and
amending a criterion after its result has been read is post-hoc loosening even with a right motive.

The confound this exists to break: softening did TWO things at once -- it removed the jump AND
weakened the correction (sigma_min 0.743/0.786 -> 0.964-0.980, the softened seeds landing on
0.52-0.54, right beside the isentropic answer 0.5605). So the measured 32.11 % -> 4.28 % spread
collapse could not be attributed. The exponent q separates them: s = 1 - w^q (1 - sigma_RH), where
q = 1 IS the previous round's ramp and q < 1 pushes partial weights back toward 1, restoring the
MAGNITUDE while keeping the ramp continuous. Both endpoints stay exact, so default-inertness holds.

Registered criteria:
  F1' ★ RELATIVE to each family's own legacy baseline: spread <= baseline/3 AND <= 5 %. The old
      absolute <= 5 % scored an 8.8x DEGRADATION as a pass on the family whose baseline was 0.21 %.
  J1  some q restores sigma_min to within 0.03 of legacy while the spread still meets F1'
      -> removing the JUMP is enough; the weakening was incidental
  J2  as q falls and the magnitude returns, the spread climbs back to >= baseline/2
      -> the previous pass was BOUGHT by weakening; kill clause 3 then kills the route
  F2 / F3 / F5 unchanged
  ★ q is a CALIBRATION: the whole curve is the reading, and the position that matters is fixed by
    LEGACY's sigma_min, not chosen after the fact.

Outputs (TRACKED): bench/gate_results/task3_magnitude_preserving.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                        # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                       # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve               # noqa: E402
from pyfp3d.post.shock import shock_report                      # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients         # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting            # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_magnitude_preserving.csv")
M_INF, ALPHA, C = 0.80, 1.25, 1.5
SEEDS = (0, 5, 12)
#: eps is FIXED at the one cell that passed F1+F2 last round; q is this round's axis.
SOFT_EPS = 0.05
QS = (1.0, 0.5, 0.25, 0.1)
#: (eps, q); (0.0, 1.0) is the legacy control that F3 and every baseline is read from
LEGS = ((0.0, 1.0),) + tuple((SOFT_EPS, q) for q in QS)
RECIPE = dict(upwind_c=C, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80)
CASES = (("hybrid_R0", None), ("unstr_coarse", "cases/meshes/naca0012_2.5d/coarse.msh"))
SIGMA_MATCH_TOL, J2_FRAC, F1_ABS_PCT, F1_REL = 0.03, 0.5, 5.0, 1.0 / 3.0
SHOCK_REF, SHOCK_TOL = 0.61, 0.02
LEG_GATE_S = 600.0


def _mesh(path):
    if path is None:
        from run_task3_refinement_paradox import build
        return build(160, 0.004, None, 1.0)[0]
    return cut_wake(read_mesh(os.path.join(REPO, path)))


def run(tag, mc, wc, eps, q, seed):
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, n_picard_seed=seed,
                            sigma_soft_eps=eps, sigma_soft_q=q, **RECIPE)
    wall = time.perf_counter() - t0
    h = list(r.get("residual_history") or [])
    z = float(np.unique(mc.nodes[:, 2]).mean())
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=z, m_inf=M_INF), M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=float(np.ptp(mc.nodes[:, 2])), m_inf=M_INF)
    xs = rep["upper"].get("x_shock")
    fr = r.get("sigma_freeze_report") or {}
    return dict(case=tag, eps=eps, q=q, seed=seed, converged=bool(r["converged"]),
                res_final=(h[-1] if h else None), n_newton=len(h),
                n_limited=int(r.get("n_limited") or 0), n_floored=int(r.get("n_floored") or 0),
                cl_p=float(f["cl"]), x_shock=None if xs is None else float(xs),
                sigma_min=r.get("sigma_min"), m1_max=r.get("m1_max"),
                n_shock_cells=r.get("n_shock_cells"),
                selection_churn=fr.get("selection_churn"),
                churn_period=fr.get("churn_period"), wall_s=round(wall, 1))


def stats(rows):
    """Spread over CONVERGED legs only; UNDEFINED below two of them (the fixed defect)."""
    good = [r for r in rows if r["converged"]]
    cl = [r["cl_p"] for r in good]
    xs = [r["x_shock"] for r in good if r["x_shock"] is not None]
    sm = [r["sigma_min"] for r in good if r["sigma_min"] is not None]
    if len(cl) < 2:
        return dict(n=len(good), rel=None, dxs=None, sigma_mean=(np.mean(sm) if sm else None),
                    inband=sum(1 for v in xs if abs(v - SHOCK_REF) <= SHOCK_TOL),
                    churn=any(bool(r["selection_churn"]) for r in rows))
    rel = 100.0 * (max(cl) - min(cl)) / max(abs(np.mean([abs(v) for v in cl])), 1e-12)
    return dict(n=len(good), rel=rel,
                dxs=(max(xs) - min(xs)) if len(xs) >= 2 else None,
                sigma_mean=(float(np.mean(sm)) if sm else None),
                inband=sum(1 for v in xs if abs(v - SHOCK_REF) <= SHOCK_TOL),
                churn=any(bool(r["selection_churn"]) for r in rows))


def main():
    print(f"magnitude-preserving softening   M{M_INF}/alpha {ALPHA}/C {C}   eps {SOFT_EPS}   "
          f"q {QS}   seeds {SEEDS}   threads {os.environ['NUMBA_NUM_THREADS']}\n")
    rows = []
    for tag, path in CASES:
        mc, wc = _mesh(path)
        for eps, q in LEGS:
            for seed in SEEDS:
                row = run(tag, mc, wc, eps, q, seed)
                rows.append(row)
                print(f"  {tag:13} eps={eps:.2f} q={q:.2f} seed {seed:>2}  "
                      f"conv={str(row['converged']):5} |R|={row['res_final']:.2e} "
                      f"cl_p={row['cl_p']:>9.6f} "
                      f"x_shock={'-' if row['x_shock'] is None else format(row['x_shock'], '.4f')} "
                      f"sig_min={format(row['sigma_min'], '.5f') if row['sigma_min'] else '-'} "
                      f"({row['wall_s']:.0f}s)", flush=True)
                if row["wall_s"] > LEG_GATE_S:
                    print("    ★ cost gate exceeded -- stop"); _write(rows); return 1
    _write(rows)
    return _read(rows)


def _read(rows):
    #: --- F3 -----------------------------------------------------------------------------------
    print("\n=== F3: eps = 0 must bit-reproduce the legacy answers ===")
    ref = {("hybrid_R0", 0): (0.261367, 0.7141), ("hybrid_R0", 5): (0.258223, 0.6051),
           ("hybrid_R0", 12): (0.351467, 0.5458),
           ("unstr_coarse", 0): (0.407978, 0.6195), ("unstr_coarse", 5): (0.407108, 0.6073)}
    f3 = True
    for (tag, seed), (cl_a, xs_a) in ref.items():
        m = [r for r in rows if r["case"] == tag and r["eps"] == 0.0 and r["seed"] == seed]
        if not m or not m[0]["converged"]:
            print(f"  {tag:13} seed {seed:>2}: not converged -- cannot compare"); continue
        r = m[0]
        hit = abs(r["cl_p"] - cl_a) < 5e-6 and abs(r["x_shock"] - xs_a) < 5e-4
        f3 &= hit
        print(f"  {tag:13} seed {seed:>2}: {r['cl_p']:.6f}/{r['x_shock']:.4f} vs {cl_a}/{xs_a}"
              f"  -> {'MATCH' if hit else '★ MISMATCH'}")
    print(f"  -> F3 {'PASS' if f3 else 'FAIL'}")
    if not f3:
        print("  ★ kill clause 1: the flag is not default-inert, so every reading below is VOID.")
        return 1

    #: --- the q curve ---------------------------------------------------------------------------
    print("\n=== the q curve: does the spread stay collapsed once the MAGNITUDE is restored? ===")
    verdict = {}
    for tag, _ in CASES:
        base = stats([r for r in rows if r["case"] == tag and r["eps"] == 0.0])
        b_rel = "UNDEF" if base["rel"] is None else format(base["rel"], ".2f")
        print(f"\n  {tag}    legacy: spread {b_rel} % | sigma_min {base['sigma_mean']:.4f} "
              f"| in band {base['inband']}/{base['n']}")
        print(f"    F1' for this family: spread <= {(base['rel'] or 0) * F1_REL:.3f} % "
              f"(= baseline/3) AND <= {F1_ABS_PCT} %")
        print(f"    {'q':>6}{'conv':>6}{'spread %':>11}{'sigma_min':>11}{'d(sigma) vs legacy':>20}"
              f"{'F1prime':>9}{'churn':>7}{'in band':>9}")
        rowsq = []
        for q in QS:
            st = stats([r for r in rows if r["case"] == tag and r["eps"] == SOFT_EPS
                        and r["q"] == q])
            dsig = (None if st["sigma_mean"] is None or base["sigma_mean"] is None
                    else st["sigma_mean"] - base["sigma_mean"])
            f1p = (st["rel"] is not None and base["rel"] is not None
                   and st["rel"] <= base["rel"] * F1_REL and st["rel"] <= F1_ABS_PCT)
            rowsq.append((q, st, dsig, f1p))
            band = "{}/{}".format(st["inband"], st["n"])
            c_rel = "UNDEF" if st["rel"] is None else format(st["rel"], ".2f")
            c_sig = "-" if st["sigma_mean"] is None else format(st["sigma_mean"], ".4f")
            c_d = "-" if dsig is None else format(dsig, "+.4f")
            print(f"    {q:>6.2f}{st['n']:>4}/3{c_rel:>11}{c_sig:>11}{c_d:>20}"
                  f"{str(f1p):>9}{str(st['churn']):>7}{band:>9}")
        verdict[tag] = (base, rowsq)

    #: --- J1 / J2 / J3 --------------------------------------------------------------------------
    print("\n=== J1 / J2 / J3 (binding = hybrid_R0, the family where the disease is worst) ===")
    for tag, _ in CASES:
        base, rowsq = verdict[tag]
        undef = sum(1 for _, st, _, _ in rowsq if st["rel"] is None)
        matched = [(q, st, d, f1) for q, st, d, f1 in rowsq
                   if d is not None and abs(d) <= SIGMA_MATCH_TOL]
        if undef >= 2:
            print(f"  {tag:13} -> J3   {undef}/4 q values UNDEFINED -- curve unreadable")
            continue
        if matched and any(f1 for _, _, _, f1 in matched):
            hits = [q for q, _, _, f1 in matched if f1]
            print(f"  {tag:13} -> ★ J1   at q = {hits}: sigma_min back within "
                  f"{SIGMA_MATCH_TOL} of legacy AND F1' still met ⇒ the JUMP was the cause")
        elif matched:
            worst = max(st["rel"] for _, st, _, _ in matched if st["rel"] is not None)
            back = worst >= (base["rel"] or 0) * J2_FRAC
            lab = "★ J2" if back else "J3"
            why = ("this is >= baseline/2 ⇒ the earlier pass was BOUGHT by weakening" if back
                   else "this falls between the bands")
            print(f"  {tag:13} -> {lab}   at matched magnitude the spread is {worst:.2f} % "
                  f"vs baseline {base['rel']:.2f} % ({why})")
        else:
            sig = [f"{d:+.4f}" for _, _, d, _ in rowsq if d is not None]
            print(f"  {tag:13} -> J3   no q restores sigma_min to within {SIGMA_MATCH_TOL} "
                  f"(offsets {sig}) ⇒ the exponent cannot recover the magnitude at this eps")
    print("\n  ★ J2 is a USEFUL NEGATIVE and simultaneously kills the softening route (clause 3).")
    print("  ★ q is a calibration: read the curve, not the best cell.")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
