"""Is the entropy correction a SELECTION operator? sigma ON vs OFF, cross-seed spread.

Pre-registered in phases/p3/docs/dev_phase_three/20260811-2300-sigma-selection-prereg.md, committed
before this file was written.

The upstream finding: at M0.80/alpha1.25 the discrete problem has SEVERAL solutions, and the
Picard seed decides which one you land on (cl spread 23-32 % at a bit-identical mesh, each state a
genuine root at |R| <= 3e-11 with 0/0 clamps). Transonic full-potential non-uniqueness under lift
is a documented property of the ISENTROPIC model, whose mechanism is the missing entropy condition
-- and this project's entropy correction exists for exactly that. So: does it actually narrow the
solution set?

Two things the registration pins so the result cannot be misread:
  * sigma is ALREADY ON in the baseline that fails, so this cannot deliver a cure -- only an
    attribution.
  * the comparison is SPREAD vs SPREAD, never solution vs solution. Turning sigma off changes the
    MODEL, so the two legs' solutions are not supposed to agree; the within-configuration
    cross-seed spread is the comparable statistic.

★ Writes to its OWN filename. `gate_results/m1_gate.csv` is the committed isentropic record whose
numbers the capability boundary cites by value (+101.4 % / 36.1 %), and `run_m1_gate.py --density
isentropic` would overwrite it. Nothing here touches it.

Outputs (TRACKED): bench/gate_results/task3_sigma_selection.csv
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

CSV = os.path.join(HERE, "gate_results", "task3_sigma_selection.csv")
M_INF, ALPHA, C = 0.80, 1.25, 1.5
#: the seed set is DECLARED IN ADVANCE (registration section 1) -- that declaration is the only
#: thing that makes "seed 12 was not cherry-picked after seeing the answer" a checkable claim
SEEDS = (0, 5, 12)
#: the `run_m1_gate.py` recipe, read from the source rather than recalled
RECIPE = dict(upwind_c=C, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
              precond="direct", direct_refactor_every=4, n_newton_max=80)
#: (tag, loader) -- unstructured = M1's own family; hybrid = this phase's generator
CASES = (("unstr_coarse", "cases/meshes/naca0012_2.5d/coarse.msh"),
         ("unstr_medium", "cases/meshes/naca0012_2.5d/medium.msh"),
         ("hybrid_R0", None))
LEG_GATE_S = 600.0
#: registration P3: the instrument check against committed 8-thread M1 values
ANCHOR = {(0,): (0.8798247106, 4388, 563), (5,): (0.5737947433, 0, 0)}


def _fmt(x, nd=4):
    """A named helper rather than a nested conditional inside an f-string: python 3.11 rejects
    same-type quotes nested in an f-string, and that syntax error cost a launch."""
    return "-" if x is None else f"{x:.{nd}f}"


def _mesh(path):
    if path is None:
        (mc, wc), _ = build(160, 0.004, None, 1.0)
        return mc, wc
    return cut_wake(read_mesh(os.path.join(REPO, path)))


def run(tag, mc, wc, sigma_on, seed):
    kw = dict(RECIPE)
    #: ★ sigma OFF uses the EXISTING library keyword -- no numerics change in this round
    if not sigma_on:
        kw["entropy_correction"] = False
    t0 = time.perf_counter()
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, n_picard_seed=seed, **kw)
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
    return dict(case=tag, sigma="ON" if sigma_on else "OFF", seed=seed, C=C,
                n_tets=len(mc.elements), converged=bool(r["converged"]),
                res_final=float(h[-1]) if len(h) else None, n_limited=nlim, n_floored=nflr,
                failure_mode=mode, m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                x_shock=None if xs is None else round(float(xs), 7),
                cl_p=round(float(f["cl"]), 6),
                sigma_min=(None if r.get("sigma_min") is None
                           else round(float(r["sigma_min"]), 6)),
                n_newton=len(h), threads=os.environ["NUMBA_NUM_THREADS"], wall_s=round(wall, 1))


def spread(rows, key):
    """Registration P2: a spread needs TWO converged seeds or it is UNDEFINED -- it must never be
    reported as a small spread (the mirror image of last round's `d2h == 0` vacuous PASS)."""
    v = [r[key] for r in rows if r["converged"] and r[key] is not None]
    if len(v) < 2:
        return None, len(v)
    return max(v) - min(v), len(v)


def main():
    print(f"sigma selection probe   M{M_INF} / alpha {ALPHA} / C {C}   seeds {SEEDS}   "
          f"threads {os.environ['NUMBA_NUM_THREADS']}\n")
    rows = []
    for tag, path in CASES:
        mc, wc = _mesh(path)
        for sigma_on in (True, False):
            for seed in SEEDS:
                row = run(tag, mc, wc, sigma_on, seed)
                rows.append(row)
                print(f"  {tag:13} sigma {row['sigma']:3} seed {seed:>2}  "
                      f"conv={str(row['converged']):5} |R|={row['res_final']:.2e} "
                      f"lim/flr={row['n_limited']}/{row['n_floored']} "
                      f"M_max={row['m_max']:.4f} "
                      f"x_shock={_fmt(row['x_shock'])} "
                      f"cl_p={row['cl_p']:>9.6f} sig_min={row['sigma_min']} "
                      f"({row['wall_s']:.0f}s) {row['failure_mode'] or ''}", flush=True)
                if row["wall_s"] > LEG_GATE_S:
                    print(f"    ★ cost gate {row['wall_s']:.0f}s > {LEG_GATE_S:.0f}s -- stop, "
                          f"no budget added")
                    _write(rows); return 1

    _write(rows)

    #: --- P3 instrument check: bit-identity against the committed 8-thread M1 values -----------
    print("\n=== P3 instrument check vs the committed 8-thread M1 numbers ===")
    ok_instr = True
    for seed, (xs_a, lim_a, flr_a) in ((0, ANCHOR[(0,)]), (5, ANCHOR[(5,)])):
        m = [r for r in rows if r["case"] == "unstr_medium" and r["sigma"] == "ON"
             and r["seed"] == seed]
        if not m:
            print(f"  seed {seed}: MISSING"); ok_instr = False; continue
        r = m[0]
        hit = (abs(r["x_shock"] - xs_a) < 5e-7 and r["n_limited"] == lim_a
               and r["n_floored"] == flr_a)
        ok_instr &= hit
        print(f"  seed {seed}: x_shock {r['x_shock']:.7f} vs {xs_a} , "
              f"clamps {r['n_limited']}/{r['n_floored']} vs {lim_a}/{flr_a}  "
              f"-> {'MATCH' if hit else '★ MISMATCH'}")
    print(f"  -> {'PASS' if ok_instr else 'FAIL'}  (a harness that does not reproduce the "
          f"committed numbers cannot be used to move a conclusion)")

    #: --- Q1 reading -------------------------------------------------------------------------
    print("\n=== Q1: does sigma narrow the solution set? (SPREAD vs SPREAD) ===")
    print(f"  {'case':14}{'sigma':>6}{'conv':>6}{'cl_p spread':>14}{'rel %':>9}"
          f"{'x_shock spread':>16}")
    stats = {}
    for tag, _ in CASES:
        for sig in ("ON", "OFF"):
            sub = [r for r in rows if r["case"] == tag and r["sigma"] == sig]
            dcl, ncl = spread(sub, "cl_p")
            dxs, _ = spread(sub, "x_shock")
            mean = np.mean([abs(r["cl_p"]) for r in sub if r["converged"]]) if ncl else np.nan
            rel = None if dcl is None else 100.0 * dcl / max(abs(mean), 1e-12)
            stats[(tag, sig)] = dict(dcl=dcl, rel=rel, dxs=dxs, n=ncl)
            print(f"  {tag:14}{sig:>6}{ncl:>4}/3"
                  f"{('UNDEFINED' if dcl is None else f'{dcl:.6f}'):>14}"
                  f"{('-' if rel is None else f'{rel:.2f}'):>9}"
                  f"{('UNDEFINED' if dxs is None else f'{dxs:.4f}'):>16}")

    print("\n  bands (pre-registered section 2; medium is binding):")
    verdict_rows = []
    for tag, _ in CASES:
        on, off = stats[(tag, "ON")], stats[(tag, "OFF")]
        #: ★ vacuity guard: the control has no signal unless sigma ON's own spread is non-trivial
        if on["rel"] is None or off["rel"] is None:
            band, why = "S3", (f"spread UNDEFINED (converged seeds ON {on['n']}/3, "
                               f"OFF {off['n']}/3) -- P2 forbids calling this a small spread")
        elif on["rel"] <= 5.0:
            band, why = "S3", (f"VACUOUS: sigma-ON's own spread is only {on['rel']:.2f} % "
                               f"(<= 5 %), so the control carries no signal here")
        else:
            ratio = off["rel"] / on["rel"]
            fewer = off["n"] < on["n"]
            if ratio >= 2.0 or fewer:
                band, why = "S1", (f"OFF/ON spread ratio {ratio:.2f}"
                                   + (f" and OFF converges on fewer seeds ({off['n']} < {on['n']})"
                                      if fewer else ""))
            elif 1 / 1.3 <= ratio <= 1.3 and off["n"] == on["n"]:
                band, why = "S2", f"OFF/ON spread ratio {ratio:.2f}, same converged-seed count"
            else:
                band, why = "S3", f"OFF/ON spread ratio {ratio:.2f} falls between the bands"
        verdict_rows.append((tag, band, why))
        print(f"    {tag:14} -> {band}   {why}")

    binding = dict((t, b) for t, b, _ in verdict_rows).get("unstr_medium")
    print(f"\n  BINDING (unstr_medium): {binding}")
    print("    S1 = sigma IS doing selection work but not enough -> 'strengthen the selection'")
    print("    S2 = sigma is NOT a selection operator -> that route dies; go to the limit_cycle")
    print("         / frozen-selection lead instead")
    print("    S3 = RECORDED, no direction claimed")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
