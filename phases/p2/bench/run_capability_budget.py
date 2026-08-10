"""Is ls_naca_medium's M0.78 failure a real limit, or just the iteration cap?

The matrix recorded ls_naca_medium M0.78 as NOT_CONVERGED at |R| = 1.5802e-10 against a
1e-10 tolerance -- 1.58x, while all four conforming NOT_CONVERGED rows sit at 940x-44500x.
The descent10 diagnostic says why: it descends 807x over its last ten steps and stopped at
n_newton = 80 = n_newton_max. It ran out of budget while converging cleanly; at the measured
rate (~0.476 per step) one or two more iterations should clear the tolerance.

So the honest fix is to raise the budget and MEASURE, not to adjudicate the row into a pass.

★ This is a RECIPE DEVIATION and is recorded as one. The matrix row stays exactly as it was
measured under the pre-registered recipe; this result lands in its own CSV with the
deviation named in the row, so nobody can later read a raised-budget number as if the
pre-registered recipe had produced it. If it converges, the ls_naca_medium envelope reads
"M0.75 under the pre-registered recipe, M0.78 with n_newton_max raised to 160" -- both
halves stated.

Outputs (TRACKED): bench/gate_results/capability_budget.csv
"""

import csv
import os
import sys

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

CSV = os.path.join(_GATE, "capability_budget.csv")
N_NEWTON_RAISED = 160          # from the pre-registered 80


def main():
    #: patch only n_newton_max, by wrapping the cell's own solver -- so the mesh, the
    #: far field, freeze_tol, n_seed and the ramp all stay verbatim and the deviation is
    #: exactly one number.
    orig = cap.ls_naca

    def ls_naca_bigger(mesh_path, m, alpha):
        import pyfp3d.solve.newton_ls as nls
        real_single, real_trans = (nls.solve_multivalued_newton,
                                   nls.solve_multivalued_newton_transonic)

        def patch(fn):
            def inner(*a, **kw):
                kw["n_newton_max"] = N_NEWTON_RAISED
                return fn(*a, **kw)
            return inner
        cap.solve_multivalued_newton = patch(real_single)
        cap.solve_multivalued_newton_transonic = patch(real_trans)
        try:
            return orig(mesh_path, m, alpha)
        finally:
            cap.solve_multivalued_newton = real_single
            cap.solve_multivalued_newton_transonic = real_trans

    cell = "ls_naca_medium"
    meta = {c[0]: c for c in cap.CELLS}[cell]
    _, path, geom, mdir, level, alpha, _ = meta
    print(f"=== {cell} M0.78 alpha {alpha}: n_newton_max 80 -> "
          f"{N_NEWTON_RAISED} ===", flush=True)
    row, payload = cap.measure(cell, path, geom, mdir, level, alpha,
                               ls_naca_bigger, 0.78)
    row["note"] = f"RECIPE DEVIATION: n_newton_max {N_NEWTON_RAISED} (matrix used 80)"
    head = not os.path.exists(CSV)
    with open(CSV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), extrasaction="ignore")
        if head:
            w.writeheader()
        w.writerow(row)
    print(f"  {row['status']}  conv={row.get('converged')} "
          f"n_newton={row.get('n_newton')} |R|={row.get('res_final')} "
          f"m_attained={row.get('m_attained')} M_max={row.get('m_max')} "
          f"cl_p={row.get('cl_p')} descent10={row.get('descent10')} "
          f"({row.get('wall_s')}s)", flush=True)
    if payload is not None and row["status"].startswith("CLEAN"):
        cap.save_cp(f"{cell}_nn{N_NEWTON_RAISED}", 0.78, geom, payload)
    if row["status"] == "CLEAN":
        print("  => budget-limited CONFIRMED: envelope M0.75 under the pre-registered "
              "recipe, M0.78 with the raised budget. Both halves get reported.")
    elif row["status"] == "MACH_NOT_ATTAINED":
        print(f"  => still short of M0.78 (attained {row['m_attained']}) -- the cap was "
              "NOT the binding limit; the ramp itself stops earlier.")
    else:
        print(f"  => still {row['status']} -- the iteration cap was not the whole story.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
