"""Is sigma's STRENGTH the dose behind the non-uniqueness? A theta sweep, not an ON/OFF pair.

★★ NOT RUNNABLE AS SHIPPED -- and that is the intended end state. The `sigma_scale` instrument
this script drives was TEMPORARY by registration and was DELETED at the round's close-out
(phases/p3/docs/dev_phase_three/20260812-0500-sigma-strength-verdict.md section 6, the B20 precedent: a knob
built solely to make an A/B measurable is removed on adoption). The evidence is the committed
CSV; this file records exactly how it was produced. To re-run, restore the instrument by
cherry-picking the library half of commit a07e9b3 (37 lines, two files, default bit-identical).

Pre-registered in phases/p3/docs/dev_phase_three/20260812-0300-sigma-strength-prereg.md, committed before
the instrument was implemented.

Why a sweep and not the ON/OFF result already in hand: ON/OFF is TWO POINTS, and two points cannot
separate "sigma causes the non-uniqueness" from "these are simply two different models, each with
its own solution structure". A dose-response can. There was no strength knob to sweep -- sigma is
physics (the Rankine-Hugoniot total-pressure ratio) and `entropy_correction` is a bool -- so the
registration declares a TEMPORARY instrument, `sigma_scale` (theta), with its removal criterion
fixed in advance per the B20 precedent.

★ Guards carried in from the defect this project fixed one round ago: every spread is computed over
CONVERGED states ONLY, and a spread with fewer than two converged seeds is UNDEFINED -- never
reported as a small spread. That single missing filter turned clamped garbage into "another
solution" and put two wrong numbers into a published verdict.

Outputs (TRACKED): bench/gate_results/task3_sigma_strength.csv
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
from failure_modes import classify_failure               # noqa: E402
from run_task3_refinement_paradox import build                  # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_sigma_strength.csv")
M_INF, ALPHA, C = 0.80, 1.25, 1.5
SEEDS = (0, 5, 12)
#: declared in advance (registration section 2) -- so is the seed set
THETAS = (0.0, 0.25, 0.5, 0.75, 1.0)
RECIPE = dict(upwind_c=C, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80)
#: unstructured MEDIUM is excluded WITH its reason on the record: measured last round at 0/3
#: converged for both sigma ON and OFF at this condition on 8 threads, so every theta would be
#: UNDEFINED. Registered so that dropping the hardest case is not an invisible act.
CASES = (("hybrid_R0", None), ("unstr_coarse", "cases/meshes/naca0012_2.5d/coarse.msh"))
#: registration P3 -- the committed 8-thread M1 anchors, theta = 1 must reproduce them
ANCHORS = {("unstr_coarse", 0): 0.6194677102, ("unstr_coarse", 5): 0.6072503529}
LEG_GATE_S = 600.0
MONO_TOL = 0.05      # a neighbouring theta may fall back by <= 5 % relative and still count
DOSE_MIN = 3.0       # spread(theta=1) / spread(theta=0.25) required by band T1


def _mesh(path):
    if path is None:
        (mc, wc), _ = build(160, 0.004, None, 1.0)
        return mc, wc
    return cut_wake(read_mesh(os.path.join(REPO, path)))


def run(tag, mc, wc, theta, seed):
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, n_picard_seed=seed,
                             sigma_scale=theta, **RECIPE)
    wall = time.perf_counter() - t0
    h = np.asarray(r.get("residual_history", []), dtype=float)
    nlim, nflr = int(r.get("n_limited") or 0), int(r.get("n_floored") or 0)
    mode = ""
    if not r["converged"]:
        mode = classify_failure(h, np.asarray(r.get("clamp_history", []), dtype=float),
                                np.asarray(r.get("F_history", []), dtype=float),
                                int(r.get("n_gmres_stalled") or 0),
                                str(r.get("accept_reason")), nlim, nflr)[0]
    z = float(np.unique(mc.nodes[:, 2]).mean())
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=z, m_inf=M_INF), M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=float(np.ptp(mc.nodes[:, 2])), m_inf=M_INF)
    xs = rep["upper"].get("x_shock")
    return dict(case=tag, theta=theta, seed=seed, converged=bool(r["converged"]),
                res_final=float(h[-1]) if len(h) else None, n_limited=nlim, n_floored=nflr,
                failure_mode=mode, m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                x_shock=None if xs is None else float(xs), cl_p=float(f["cl"]),
                sigma_min=(None if r.get("sigma_min") is None else float(r["sigma_min"])),
                n_shock_cells=r.get("n_shock_cells"), n_newton=len(h),
                threads=os.environ["NUMBA_NUM_THREADS"], wall_s=round(wall, 1))


def spread(rows, key):
    """CONVERGED states only, and UNDEFINED below two of them (registration P1)."""
    v = [r[key] for r in rows if r["converged"] and r[key] is not None]
    if len(v) < 2:
        return None, len(v)
    return max(v) - min(v), len(v)


def main():
    print(f"sigma STRENGTH sweep   M{M_INF}/alpha {ALPHA}/C {C}   thetas {THETAS}   "
          f"seeds {SEEDS}   threads {os.environ['NUMBA_NUM_THREADS']}\n")
    rows = []
    for tag, path in CASES:
        mc, wc = _mesh(path)
        for theta in THETAS:
            for seed in SEEDS:
                row = run(tag, mc, wc, theta, seed)
                rows.append(row)
                print(f"  {tag:13} th={theta:.2f} seed {seed:>2}  "
                      f"conv={str(row['converged']):5} |R|={row['res_final']:.2e} "
                      f"lim/flr={row['n_limited']}/{row['n_floored']} "
                      f"x_shock={'-' if row['x_shock'] is None else format(row['x_shock'], '.4f')} "
                      f"cl_p={row['cl_p']:>9.6f} "
                      f"sig_min={'-' if row['sigma_min'] is None else format(row['sigma_min'], '.5f')}"
                      f" ({row['wall_s']:.0f}s) {row['failure_mode']}", flush=True)
                if row["wall_s"] > LEG_GATE_S:
                    print(f"    ★ cost gate exceeded -- stop, no budget added")
                    _write(rows); return 1
    _write(rows)

    #: --- P3: theta = 1 must reproduce the committed 8-thread M1 anchors -----------------------
    print("\n=== P3 instrument check: theta = 1 vs the committed 8-thread M1 anchors ===")
    ok = True
    for (tag, seed), xa in ANCHORS.items():
        m = [r for r in rows if r["case"] == tag and r["theta"] == 1.0 and r["seed"] == seed]
        hit = bool(m) and m[0]["x_shock"] is not None and abs(m[0]["x_shock"] - xa) < 5e-10
        ok &= hit
        got = "MISSING" if not m else f"{m[0]['x_shock']:.10f}"
        print(f"  {tag} seed {seed}: {got} vs {xa}  -> {'MATCH' if hit else '★ MISMATCH'}")
    print(f"  -> {'PASS' if ok else 'FAIL'}   ★ kill clause 2: a FAIL here VOIDS every physical "
          f"reading below, because the instrument would have changed what it measures")

    #: --- T1 / T2 / T3 -------------------------------------------------------------------------
    print("\n=== dose-response: cross-seed spread vs theta (CONVERGED states only) ===")
    table = {}
    for tag, _ in CASES:
        print(f"\n  {tag}")
        print(f"    {'theta':>6}{'conv':>7}{'cl_p spread':>14}{'rel %':>9}{'x_shock spread':>16}"
              f"{'sigma_min':>11}")
        for theta in THETAS:
            sub = [r for r in rows if r["case"] == tag and r["theta"] == theta]
            dcl, n = spread(sub, "cl_p")
            dxs, _ = spread(sub, "x_shock")
            good = [r for r in sub if r["converged"]]
            mean = np.mean([abs(r["cl_p"]) for r in good]) if good else np.nan
            rel = None if dcl is None else 100.0 * dcl / max(abs(mean), 1e-12)
            smin = [r["sigma_min"] for r in sub if r["sigma_min"] is not None]
            table[(tag, theta)] = rel
            print(f"    {theta:>6.2f}{n:>4}/{len(sub):<2}"
                  f"{('UNDEFINED' if dcl is None else format(dcl, '.6f')):>14}"
                  f"{('-' if rel is None else format(rel, '.2f')):>9}"
                  f"{('UNDEFINED' if dxs is None else format(dxs, '.4f')):>16}"
                  f"{(format(min(smin), '.5f') if smin else '-'):>11}")

    print("\n  bands (pre-registered section 3):")
    for tag, _ in CASES:
        vals = [table[(tag, t)] for t in THETAS]
        undef = sum(v is None for v in vals)
        if undef >= 2:
            print(f"    {tag:13} -> T3   {undef}/5 thetas UNDEFINED -- RECORDED, no direction")
            continue
        seq = [(t, v) for t, v in zip(THETAS, vals) if v is not None]
        drops = [(a[0], b[0]) for a, b in zip(seq, seq[1:])
                 if b[1] < a[1] * (1.0 - MONO_TOL)]
        base = table[(tag, 0.25)]
        top = table[(tag, 1.0)]
        ratio = (None if base in (None, 0.0) or top is None else top / base)
        mono = not drops
        if mono and ratio is not None and ratio >= DOSE_MIN:
            print(f"    {tag:13} -> T1   monotone, spread(1)/spread(0.25) = {ratio:.2f} "
                  f">= {DOSE_MIN} ⇒ DOSED by sigma's strength")
        else:
            why = []
            if drops:
                why.append(f"non-monotone at {drops}")
            if ratio is not None and ratio < DOSE_MIN:
                why.append(f"ratio {ratio:.2f} < {DOSE_MIN}")
            if ratio is None:
                why.append("ratio undefined (spread at theta=0.25 is zero or missing)")
            print(f"    {tag:13} -> T2   {'; '.join(why)} ⇒ THRESHOLD-like, not dose-like: "
                  f"sigma's STRENGTH is not the driver")
    print("\n  T1 = strength IS the dose -> the direction is 'the correction's magnitude'")
    print("  T2 = threshold-like -> the candidate becomes sigma's INTERACTION (e.g. clamping)")
    print("  T3 = RECORDED, no direction")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
