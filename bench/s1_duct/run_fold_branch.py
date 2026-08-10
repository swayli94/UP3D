"""GS1b.1: is there a fold? Follow the solution branch in M_inf and measure
dGamma/dM_inf and the existence boundary M_inf_max(h).

★★ RETRACTED BY GS1b.2 Q4 (2026-07-29) -- READ THIS BEFORE QUOTING ANY
`M_inf_max` FROM THIS SCRIPT OR ITS CSVs. The step refinement below is ONE-SHOT
(the `refined` flag), so this walker reports "the branch ends" at the first
failure it cannot refine past. `run_fold_anatomy.py` re-probed all three
terminations with smaller steps, seeded from the very same converged last-good
states, and EVERY ONE of them walks past:

    medium C=1.5   "ends" 0.7925 -> reaches 0.8015 (dGamma/dM peaks ~42 then
                                    FALLS to 11.4: Gamma turning over smoothly)
    fine   C=1.5   "ends" 0.7700 -> reaches 0.7790
    fine   C=3.0   "ends" 0.7500 -> reaches 0.7725

So `M_inf_max` as measured here is a CONTINUATION STEP-SIZE limit, not a limit
point, and GS1b.1's fold verdict is withdrawn. What survives is the slope
measurement on converged states (dGamma/dM_inf at fixed M_inf grows with
refinement), which is unaffected by the walker. The sound acceptance quantity is
the three-level cl spread at a fixed condition -- `run_h_convergence_vs_C.py`.

Use `run_h_convergence_vs_C.py::continue_to` for new continuation work: it halves
the step repeatedly (up to 4 times) instead of once.

Why this quantity and not a singular value (round file 1.1): the Schur
complement of the Gamma row is rendering-dependent -- the probe's |1-b| vanishes
even at SUBSONIC M0.5 where everything converges (0.0670/0.0454/0.0309) while
the pressure rendering's D_raw*S_code GROWS (5.65/7.42/9.73) on the very same
states. Two renderings, opposite trends, identical solutions. A fold indicator
must therefore live on the solution branch itself:

    a limit point means dGamma/dM_inf -> infinity and the branch TERMINATES,
    i.e. there is an M_inf_max beyond which no solution exists.

Method: fix the mesh, walk M_inf upward in steps, seeding each solve from the
previous Mach's converged state (branch following). Target-Mach direct Newton
(GS1.2b), no Mach ramp inside a step. Near the failure point the step is
refined. Clamped states are discarded (GS1.4 contract).

Outputs: results/gs1b_1_fold_branch.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting           # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
#: GS1b.2 Q2: the dissipation constant is a shock-WEAKENING knob, so sweeping it
#: at fixed h answers "does weakening the shock move the fold?" without any
#: entropy-correction code. Overridable from the command line.
UPWIND_C = 1.5
M_START = 0.72
DM_COARSE = 0.01
DM_REFINE = 0.0025
M_HARD_STOP = 0.86
#: fine is ~60-110 s per point, so it walks on the coarse step only
LEVELS = (("coarse", DM_COARSE, DM_REFINE),
          ("medium", DM_COARSE, DM_REFINE),
          ("fine", 0.02, 0.005))


def solve_at(mc, wc, m_inf, phi_init=None, gamma_init=None):
    kw = dict(m_inf=m_inf, alpha_deg=ALPHA, upwind_c=UPWIND_C, m_crit=0.95,
              freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
              direct_refactor_every=4, n_newton_max=80)
    if phi_init is not None:
        kw.update(phi_init=phi_init, gamma_init=gamma_init, n_picard_seed=0)
    return solve_newton_lifting(mc, wc, **kw)


def walk(level, dm, dm_fine, rows):
    path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
    if not path.exists():
        print(f"skip {level}: mesh missing")
        return None
    mc, wc = cut_wake(read_mesh(path))
    dz = float(np.ptp(mc.nodes[:, 2]))
    print(f"\n=== {level} ({len(mc.nodes)} nodes), dm {dm} -> {dm_fine} ===",
          flush=True)
    phi, gam = None, None
    m = M_START
    step = dm
    refined = False
    last_good = None
    while m <= M_HARD_STOP + 1e-12:
        t0 = time.perf_counter()
        try:
            r = solve_at(mc, wc, m, phi, gam)
            failed = False
        except Exception as exc:                               # noqa: BLE001
            print(f"   M={m:.4f}  EXC {type(exc).__name__}", flush=True)
            r, failed = None, True
        ok = False
        if r is not None:
            ok = bool(r["converged"]) and not r.get("clamped", False)
            rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                            m_inf=m), m)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], r["phi"],
                                        alpha_deg=ALPHA, s_ref=dz, m_inf=m)
            rows.append(dict(
                level=level, n_dof=len(mc.nodes), m_inf=round(m, 5),
                converged=r["converged"], clamped=r.get("clamped"),
                usable=ok, gamma=round(float(r["gamma"][0]), 8),
                cl_p=round(f["cl"], 6), x_shock=rep["upper"].get("x_shock"),
                m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                res_final=r["residual_history"][-1],
                n_newton=r["n_newton"],
                wall_s=round(time.perf_counter() - t0, 1)))
            print(f"   M={m:.4f} {'OK ' if ok else 'BAD'} "
                  f"conv={str(r['converged']):5s} clamp={r.get('clamped')} "
                  f"gamma={rows[-1]['gamma']:.6f} cl={rows[-1]['cl_p']:.5f} "
                  f"x_sh={rows[-1]['x_shock']} M_max={rows[-1]['m_max']} "
                  f"({rows[-1]['wall_s']}s)", flush=True)
        if ok:
            phi, gam = r["phi"], r["gamma"]
            last_good = m
            m += step
        else:
            if not refined:
                # step back to the last good Mach and refine the step
                refined = True
                step = dm_fine
                m = (last_good + step) if last_good is not None else m + step
                print(f"   -> refining step to {step} from "
                      f"M={last_good}", flush=True)
                if last_good is None:
                    break
            else:
                print(f"   -> branch ends: M_inf_max({level}) = {last_good}",
                      flush=True)
                break
    return last_good


def main():
    global UPWIND_C, LEVELS
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--upwind-c", type=float, default=UPWIND_C)
    ap.add_argument("--levels", nargs="+", default=None,
                    help="subset of coarse/medium/fine")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    UPWIND_C = a.upwind_c
    if a.levels:
        LEVELS = tuple(l for l in LEVELS if l[0] in a.levels)
    print(f"upwind_c = {UPWIND_C}, levels = {[l[0] for l in LEVELS]}",
          flush=True)
    rows = []
    mmax = {}
    for level, dm, dmf in LEVELS:
        mmax[level] = walk(level, dm, dmf, rows)

    name = a.out or (f"gs1b_1_fold_branch_C{UPWIND_C:g}.csv"
                     if UPWIND_C != 1.5 else "gs1b_1_fold_branch.csv")
    with open(OUT / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / name)

    print("\n=== branch existence boundary ===")
    for level, _, _ in LEVELS:
        print(f"  M_inf_max({level:7s}) = {mmax.get(level)}")

    print("\n=== dGamma/dM_inf along the branch ===")
    for level, _, _ in LEVELS:
        sub = sorted([r for r in rows if r["level"] == level and r["usable"]],
                     key=lambda r: r["m_inf"])
        if len(sub) < 2:
            continue
        print(f"  {level}:")
        for a, b in zip(sub[:-1], sub[1:]):
            dg = (b["gamma"] - a["gamma"]) / (b["m_inf"] - a["m_inf"])
            print(f"    M {a['m_inf']:.4f} -> {b['m_inf']:.4f}: "
                  f"dGamma/dM = {dg:8.3f}   (M_max {a['m_max']:.3f} -> "
                  f"{b['m_max']:.3f})")


if __name__ == "__main__":
    main()
