"""GS1b.2 Q5: does WEAKENING THE SHOCK restore h-convergence of the lift?

This replaces the acceptance quantity of GS1b.1. That round proposed
"M_inf_max rises at fixed h", and Q4 then measured that every M_inf_max in it was
a CONTINUATION STEP limit, not a limit point (a smaller step walks past all
three) -- so M_inf_max is not a sound quantity and the fold verdict is retracted.

What survives S1 and Q2 as measured, un-contaminated by continuation control:

    lift at fixed (M_inf, alpha) grows monotonically as the shock weakening is
    REDUCED, along either axis --
      h -> 0 at fixed C   : cl 0.3725 / 0.5234 / 0.7150  (S1, GS1.3b)
      C -> 0 at fixed h   : cl 0.4360 / 0.4637 / 0.5166 / 0.5714 for
                            C = 4 / 3 / 2 / 1.5 (Q2, medium M0.79)

Both knobs act through the same physical quantity, so route B's premise is:
a PHYSICAL (h-independent) shock weakening -- the entropy correction, whose
travel Q1 measured at 5.1-8.9 % of downstream density at M1 = 1.30-1.40 -- should
make cl h-convergent where the isentropic law does not.

That premise has a direct, cheap test which needs NO new library code: the
artificial dissipation is ALSO a shock weakening, so if weakening the shock is
capable of restoring h-convergence at all, then the three-level cl spread must
be SMALLER at C = 3.0 than at C = 1.5. If cl still diverges at C = 3.0, shock
weakening alone cannot buy h-convergence and B is in trouble before it is built.

Method: branch-continue each (level, C) to the common condition M0.7875 /
alpha 1.25 with an ADAPTIVE step (halve on failure, up to 4 halvings -- the fix
for the one-shot refinement flaw Q4 exposed in run_fold_branch.py). Clamped or
non-converged endpoints are reported as such, never quoted as solutions.

Outputs: results/gs1b_2_q5_h_convergence.csv
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
M_TARGET = 0.7875           # the S1 divergence condition (GS1.3b)
M_START = 0.72
DM0 = {"coarse": 0.01, "medium": 0.01, "fine": 0.02}
MAX_HALVINGS = 4
LEVELS = ("coarse", "medium", "fine")
C_LIST = (1.5, 3.0)
#: GS1b.3 criteria F/G/H: --entropy runs the SAME protocol with the
#: entropy-corrected density, so the three-level spread is compared like for
#: like (same C, same condition, same adaptive continuation).
ENTROPY = False


def solve_at(mc, wc, m, C, phi=None, gam=None):
    kw = dict(m_inf=m, alpha_deg=ALPHA, upwind_c=C, m_crit=0.95,
              freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
              direct_refactor_every=4, n_newton_max=80,
              entropy_correction=ENTROPY)
    if phi is not None:
        kw.update(phi_init=phi, gamma_init=gam, n_picard_seed=0)
    return solve_newton_lifting(mc, wc, **kw)


def usable(r):
    return bool(r["converged"]) and not r.get("clamped", False)


def continue_to(mc, wc, C, level, m_target):
    """Adaptive branch continuation to m_target. Returns (result, m_reached,
    n_halvings, n_solves)."""
    dm0 = DM0[level]
    m, phi, gam = M_START, None, None
    last = None
    halv = 0
    n_solve = 0
    while True:
        step = dm0 / (2 ** halv)
        m_next = min(m + step, m_target)
        r = solve_at(mc, wc, m_next, C, phi, gam)
        n_solve += 1
        if usable(r):
            phi, gam, m, last = r["phi"], r["gamma"], m_next, r
            if abs(m - m_target) < 1e-12:
                return last, m, halv, n_solve
        else:
            halv += 1
            if halv > MAX_HALVINGS:
                return last, m, halv - 1, n_solve
            print(f"      halving to {dm0 / 2 ** halv:g} at M={m:.5f}",
                  flush=True)


def main():
    global ENTROPY, C_LIST, LEVELS
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--entropy", action="store_true")
    ap.add_argument("--c-list", nargs="+", type=float, default=None)
    ap.add_argument("--levels", nargs="+", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    ENTROPY = bool(a.entropy)
    if a.c_list:
        C_LIST = tuple(a.c_list)
    if a.levels:
        LEVELS = tuple(a.levels)
    print(f"entropy_correction = {ENTROPY}, C = {list(C_LIST)}, "
          f"levels = {list(LEVELS)}")
    rows = []
    for C in C_LIST:
        print(f"\n########## C = {C} ##########", flush=True)
        for level in LEVELS:
            path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
            if not path.exists():
                print(f"  skip {level}: mesh missing")
                continue
            mc, wc = cut_wake(read_mesh(path))
            dz = float(np.ptp(mc.nodes[:, 2]))
            t0 = time.perf_counter()
            r, m_reached, halv, n_solve = continue_to(mc, wc, C, level,
                                                      M_TARGET)
            wall = time.perf_counter() - t0
            if r is None:
                print(f"  {level:7s} C={C}: no usable state at all")
                continue
            rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                            m_inf=m_reached), m_reached)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], r["phi"],
                                        alpha_deg=ALPHA, s_ref=dz,
                                        m_inf=m_reached)
            rows.append(dict(
                level=level, n_dof=len(mc.nodes), upwind_c=C,
                m_reached=round(m_reached, 6),
                at_target=abs(m_reached - M_TARGET) < 1e-12,
                usable=usable(r), gamma=round(float(r["gamma"][0]), 8),
                cl_p=round(f["cl"], 6), x_shock=rep["upper"].get("x_shock"),
                m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                res_final=r["residual_history"][-1], n_halvings=halv,
                n_solves=n_solve, entropy=ENTROPY,
                # GS1b.3 criterion H: the correction must NOT decay with h
                sigma_min=r.get("sigma_min"),
                n_shock_cells=r.get("n_shock_cells"),
                m1_detected=r.get("m1_max"),
                n_sigma_refresh=r.get("n_sigma_refresh"),
                wall_s=round(wall, 1)))
            print(f"  {level:7s} C={C}: M_reached={m_reached:.5f} "
                  f"{'(TARGET)' if rows[-1]['at_target'] else '(SHORT)'} "
                  f"gamma={rows[-1]['gamma']:.6f} cl_p={rows[-1]['cl_p']:.5f} "
                  f"x_sh={rows[-1]['x_shock']} M_max={rows[-1]['m_max']} "
                  f"halvings={halv} solves={n_solve} ({wall:.0f}s)", flush=True)

    name = a.out or ("gs1b_3_h_convergence_entropy.csv" if ENTROPY
                     else "gs1b_2_q5_h_convergence.csv")
    with open(OUT / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / name)

    print(f"\n=== three-level cl_p at M {M_TARGET} / alpha {ALPHA} ===")
    for C in C_LIST:
        sub = {r["level"]: r for r in rows
               if r["upwind_c"] == C and r["at_target"]}
        got = [sub[lv]["cl_p"] for lv in LEVELS if lv in sub]
        names = [lv for lv in LEVELS if lv in sub]
        if len(got) < 2:
            print(f"  C={C}: only {len(got)} level(s) reached the target "
                  f"({names}) -- cannot form a spread")
            continue
        line = "  ".join(f"{lv}={v:.5f}" for lv, v in zip(names, got))
        print(f"  C={C}: {line}")
        for a, b, na, nb in zip(got[:-1], got[1:], names[:-1], names[1:]):
            print(f"      {na} -> {nb}: {100.0 * (b - a) / a:+.2f} %")
        print(f"      total spread (max-min)/min = "
              f"{100.0 * (max(got) - min(got)) / min(got):+.2f} %")


if __name__ == "__main__":
    main()
