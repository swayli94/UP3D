"""Phase-two benchmark: the cheap, repeatable metric set every round A/Bs against.

GS0.3. This is deliberately SMALL and FAST (~3-4 min on 16 threads, coarse
meshes only) so that it can be run every round without thinking about cost.
It is not a validation suite -- it is a drift detector plus the live reading of
the product metrics that can be had cheaply.

What it measures (one row per case):

  linalg   -- cost of ONE elliptic solve with AMG+CG, AMG setup, splu cost and
              fill, residual/matrix assembly time (audit exp5). This is the
              algorithmic floor the solver should be judged against.
  airfoil  -- NACA0012 M0.80 / alpha 1.25 coarse transonic Newton at three
              artificial-dissipation constants (audit exp1/exp3/exp6): shock
              position, cl, M_max. Product metric M1's dissipation-sensitivity
              and (with the medium leg, --with-medium) its mesh-consistency.
  wing     -- ONERA M6 coarse M0.8395 / alpha 3.06 transonic Newton with the
              two preconditioners (audit exp2): wall time and the solution, to
              keep GS3.1 honest (same answer, less time).

Usage:
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
        python bench/run_bench.py                      # write results/bench_<stamp>.csv
    python bench/run_bench.py --compare bench/baseline_2026-07-28.csv
    python bench/run_bench.py --with-medium            # + the medium legs (~5 min more)

`--compare` prints a table of every metric that moved by more than its
tolerance (see TOL below) -- that table is what a round's report quotes.
"""

import argparse
import os
import csv
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from pyfp3d.kernels.jacobian import PicardOperator            # noqa: E402
from pyfp3d.mesh.reader import read_mesh                      # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                     # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve             # noqa: E402
from pyfp3d.post.shock import shock_report                    # noqa: E402
from pyfp3d.post.surface import (planform_area,               # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.linear import (apply_dirichlet,             # noqa: E402
                                 build_amg_preconditioner)
from pyfp3d.solve.newton import solve_newton_transonic        # noqa: E402

OUT = HERE / "results"

#: per-metric relative tolerance for --compare. Physics values are held tight
#: (a real change is much larger); wall times are noisy on a shared box.
TOL = {
    "cl_p": 1e-6, "cl_kj": 1e-6, "x_shock": 1e-6, "m_max": 1e-6,
    "gamma_mean": 1e-6, "cg_iters": 0.0, "n_newton": 0.0,
    "wall_s": 0.40, "t_amg_setup_s": 0.40, "t_cg_amg_s": 0.40,
    "t_splu_s": 0.40, "t_splu_solve_ms": 0.40, "t_residual_ms": 0.40,
    "t_matrix_ms": 0.40, "n_gmres": 0.25,
}
DEFAULT_TOL = 0.10


# ---------------------------------------------------------------------------
# case 1: the linear-algebra floor
# ---------------------------------------------------------------------------

def bench_linalg(mesh_rel, tag, rows):
    path = REPO / mesh_rel
    if not path.exists():
        print(f"  skip linalg/{tag}: {mesh_rel} missing")
        return
    mesh = read_mesh(path)
    mc, _ = cut_wake(mesh)
    op = PicardOperator(mc.nodes, mc.elements)
    phi = mc.nodes[:, 0].copy()
    rho = np.ones(len(mc.elements))
    op.assemble_residual(phi, rho)                     # warm the JIT

    t0 = time.perf_counter()
    for _ in range(5):
        op.assemble_residual(phi, rho)
    t_res = (time.perf_counter() - t0) / 5
    t0 = time.perf_counter()
    for _ in range(3):
        A = op.assemble_matrix(rho)
    t_mat = (time.perf_counter() - t0) / 3
    A = A.tocsr()

    ff = np.unique(mc.boundary_faces["farfield"])
    Af, b, _, _ = apply_dirichlet(A, np.zeros(A.shape[0]), ff,
                                  mc.nodes[ff, 0].copy())
    Af = Af.tocsr()
    t0 = time.perf_counter()
    _, M = build_amg_preconditioner(Af)
    t_amg = time.perf_counter() - t0
    nit = [0]
    t0 = time.perf_counter()
    _, info = spla.cg(Af, b, M=M, rtol=1e-8, maxiter=2000,
                      callback=lambda _x: nit.__setitem__(0, nit[0] + 1))
    t_cg = time.perf_counter() - t0
    t0 = time.perf_counter()
    lu = spla.splu(Af.tocsc())
    t_lu = time.perf_counter() - t0
    t0 = time.perf_counter()
    for _ in range(3):
        lu.solve(b)
    t_lu_solve = (time.perf_counter() - t0) / 3

    rows.append(dict(
        case=f"linalg/{tag}", n_dof=Af.shape[0], nnz=A.nnz,
        t_residual_ms=round(1e3 * t_res, 3), t_matrix_ms=round(1e3 * t_mat, 3),
        t_amg_setup_s=round(t_amg, 4), t_cg_amg_s=round(t_cg, 4),
        cg_iters=nit[0], cg_info=info,
        t_splu_s=round(t_lu, 4),
        splu_fill_ratio=round((lu.L.nnz + lu.U.nnz) / A.nnz, 2),
        t_splu_solve_ms=round(1e3 * t_lu_solve, 3)))
    print("  ", rows[-1], flush=True)


# ---------------------------------------------------------------------------
# case 2: the 2-D transonic dissipation sensitivity (product metric M1)
# ---------------------------------------------------------------------------

AIRFOIL_RECIPE = dict(m_start=0.70, dm=0.025, dm_min=0.003, freeze_tol=1e-6,
                      newton_kw=dict(freeze_refresh_max=8, precond="direct",
                                     n_newton_max=60))


def bench_airfoil(level, rows, m_inf=0.80, alpha=1.25, cs=(1.0, 1.5, 3.0)):
    path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
    if not path.exists():
        print(f"  skip airfoil/{level}: mesh missing")
        return
    mc, wc = cut_wake(read_mesh(path))
    dz = float(np.ptp(mc.nodes[:, 2]))
    for C in cs:
        t0 = time.perf_counter()
        try:
            r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                                       upwind_c=C, **AIRFOIL_RECIPE)
            err = ""
        except Exception as exc:                            # noqa: BLE001
            rows.append(dict(case=f"airfoil/{level}/C{C}", converged=False,
                             error=type(exc).__name__))
            print("  ", rows[-1], flush=True)
            continue
        wall = time.perf_counter() - t0
        rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                        m_inf=m_inf), m_inf)
        f = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], r["phi"],
                                    alpha_deg=alpha, s_ref=dz, m_inf=m_inf)
        rows.append(dict(
            case=f"airfoil/{level}/C{C}", n_dof=len(mc.nodes),
            converged=r["converged"], wall_s=round(wall, 2),
            n_newton=sum(len(lr["residual_history"]) - 1
                         for lr in r["level_results"]),
            cl_p=round(f["cl"], 8), cl_kj=round(2 * float(r["gamma"][0]), 8),
            x_shock=rep["upper"].get("x_shock"),
            m_max=round(float(np.sqrt(r["mach2_max"])), 6), error=err))
        print("  ", rows[-1], flush=True)


# ---------------------------------------------------------------------------
# case 3: the 3-D wing, both preconditioners (GS3.1 honesty check)
# ---------------------------------------------------------------------------

def wing_recipe(precond):
    """tests/test_p8_newton.py::NEWTON_M6_RECIPE with `precond` swapped."""
    return dict(dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
                newton_kw=dict(freeze_refresh_max=8, precond=precond,
                               direct_refactor_every=1000, n_newton_max=60,
                               farfield_spanwise_gamma=True))


def bench_wing(level, rows, m_inf=0.8395, alpha=3.06,
               preconds=("direct", "amg")):
    path = REPO / f"cases/meshes/onera_m6/{level}.msh"
    if not path.exists():
        print(f"  skip wing/{level}: mesh missing")
        return
    mc, wc = cut_wake(read_mesh(path))
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    for precond in preconds:
        t0 = time.perf_counter()
        try:
            r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                                       **wing_recipe(precond))
            err = ""
        except Exception as exc:                            # noqa: BLE001
            rows.append(dict(case=f"wing/{level}/{precond}", converged=False,
                             error=type(exc).__name__))
            print("  ", rows[-1], flush=True)
            continue
        wall = time.perf_counter() - t0
        f = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], r["phi"],
                                    alpha_deg=alpha, s_ref=s_ref, m_inf=m_inf)
        tm = r["timings_total"]
        rows.append(dict(
            case=f"wing/{level}/{precond}", n_dof=len(mc.nodes),
            converged=r["converged"], wall_s=round(wall, 2),
            n_newton=sum(len(lr["residual_history"]) - 1
                         for lr in r["level_results"]),
            n_gmres=r["n_gmres_total"],
            t_precond=round(tm.get("precond", 0.0), 2),
            t_linsolve=round(tm.get("linsolve", 0.0), 2),
            cl_p=round(f["cl"], 8),
            gamma_mean=round(float(np.mean(r["gamma"])), 8),
            m_max=round(float(np.sqrt(r["mach2_max"])), 6), error=err))
        print("  ", rows[-1], flush=True)


# ---------------------------------------------------------------------------
# driver + compare
# ---------------------------------------------------------------------------

FIELDS = ["case", "n_dof", "nnz", "converged", "wall_s", "n_newton", "n_gmres",
          "cl_p", "cl_kj", "gamma_mean", "x_shock", "m_max",
          "t_residual_ms", "t_matrix_ms", "t_amg_setup_s", "t_cg_amg_s",
          "cg_iters", "cg_info", "t_splu_s", "splu_fill_ratio",
          "t_splu_solve_ms", "t_precond", "t_linsolve", "n_threads", "error"]


#: ★ GS1b.3 (2026-07-29, measured): the timing metrics are wall clock, so a
#: baseline taken at a different THREAD COUNT reports them as moved with no code
#: change at all (8t vs the 16t runner default moved t_residual_ms by 85 %,
#: t_precond / t_linsolve by 13-27 %). The count now travels in the CSV and
#: `compare` says so up front -- the same fix bitcheck.py got the same day.
#: Separately, and recorded because it was mis-attributable: the `wing/coarse`
#: n_newton 17 -> 18 against the 2026-07-28 baseline is NOT an entropy-correction
#: effect. Running this script at commit a2cb9c3 (before any GS1b.3 library
#: change) also gives 18, with cl_p, gamma_mean, m_max and n_gmres identical. It
#: traces to GS1.4, which made solve_subsonic_lifting refuse `converged` while
#: clamped and so changed the Picard SEED the wing Newton starts from -- baseline
#: staleness from an earlier round, not this one.
def compare(new_rows, baseline_csv):
    with open(baseline_csv, newline="") as fh:
        base = {r["case"]: r for r in csv.DictReader(fh)}
    print(f"\n=== compare vs {baseline_csv} ===")
    here = os.environ.get("NUMBA_NUM_THREADS", "unset")
    there = next((r.get("n_threads") for r in base.values()
                  if r.get("n_threads")), None)
    if str(there) != str(here):
        print(f"  !! threads: baseline {there or 'unrecorded'}, this run {here}"
              f" -- TIMING metrics are not comparable across thread counts;"
              f" read cl_p / gamma_mean / m_max / n_newton instead.")
    moved = 0
    for row in new_rows:
        b = base.get(row["case"])
        if b is None:
            print(f"  NEW   {row['case']}")
            continue
        for k, v in row.items():
            if k in ("case", "error") or v is None or b.get(k) in (None, ""):
                continue
            try:
                vn, vb = float(v), float(b[k])
            except (TypeError, ValueError):
                if str(v) != str(b[k]):
                    print(f"  MOVED {row['case']:28s} {k}: "
                          f"{b[k]} -> {v}")
                    moved += 1
                continue
            tol = TOL.get(k, DEFAULT_TOL)
            scale = max(abs(vb), abs(vn), 1e-30)
            rel = abs(vn - vb) / scale
            if rel > tol:
                print(f"  MOVED {row['case']:28s} {k}: "
                      f"{vb:.6g} -> {vn:.6g}  ({100 * rel:+.1f} %, "
                      f"tol {100 * tol:.0f} %)")
                moved += 1
    print(f"  {moved} metric(s) outside tolerance")
    return moved


# ---------------------------------------------------------------------------
# case 4: M1a -- the in-envelope mesh-convergence point (GS1.5a / decision D3)
# ---------------------------------------------------------------------------

def bench_m1a(rows, with_fine=False):
    """NACA0012 M0.72/alpha1.25 (M_max ~ 1.17) -- INSIDE the measured
    h-convergent envelope (M_max <~ 1.2, GS1.3b). Target-Mach direct Newton,
    no Mach ramp (GS1.2b). The `airfoil/*` group above is M0.80, which is
    OUTSIDE the envelope and diverges under refinement -- these rows are the
    ones that are supposed to stay put."""
    from pyfp3d.solve.newton import solve_newton_lifting
    levels = ["coarse", "medium"] + (["fine"] if with_fine else [])
    for level in levels:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"  skip m1a/{level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        dz = float(np.ptp(mc.nodes[:, 2]))
        t0 = time.perf_counter()
        r = solve_newton_lifting(mc, wc, m_inf=0.72, alpha_deg=1.25,
                                 upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6,
                                 freeze_refresh_max=8, precond="direct",
                                 direct_refactor_every=4, n_newton_max=80)
        wall = time.perf_counter() - t0
        rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                        m_inf=0.72), 0.72)
        f = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], r["phi"],
                                    alpha_deg=1.25, s_ref=dz, m_inf=0.72)
        rows.append(dict(
            case=f"m1a/{level}", n_dof=len(mc.nodes),
            converged=r["converged"], wall_s=round(wall, 2),
            n_newton=r["n_newton"], cl_p=round(f["cl"], 8),
            cl_kj=round(2 * float(r["gamma"][0]), 8),
            x_shock=rep["upper"].get("x_shock"),
            m_max=round(float(np.sqrt(r["mach2_max"])), 6),
            res_final=r["residual_history"][-1]))
        print("  ", rows[-1], flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compare", metavar="CSV",
                    default=str(HERE / "baseline_2026-07-28.csv"),
                    help="baseline CSV to diff against ('' to skip)")
    ap.add_argument("--with-medium", action="store_true",
                    help="also run the medium legs (~5 min more)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    OUT.mkdir(exist_ok=True)
    rows = []
    print("=== linalg ===", flush=True)
    bench_linalg("cases/meshes/naca0012_2.5d/coarse.msh", "naca_coarse", rows)
    bench_linalg("cases/meshes/onera_m6/coarse.msh", "m6_coarse", rows)
    if args.with_medium:
        bench_linalg("cases/meshes/onera_m6/medium.msh", "m6_medium", rows)
    print("=== airfoil (product metric M1) ===", flush=True)
    bench_airfoil("coarse", rows)
    if args.with_medium:
        bench_airfoil("medium", rows, cs=(1.5, 3.0))
    print("=== m1a (in-envelope convergence point) ===", flush=True)
    bench_m1a(rows, with_fine=args.with_medium)
    print("=== wing ===", flush=True)
    bench_wing("coarse", rows)
    if args.with_medium:
        bench_wing("medium", rows)

    stamp = time.strftime("%Y%m%d-%H%M")
    out = Path(args.out) if args.out else OUT / f"bench_{stamp}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        nth = os.environ.get("NUMBA_NUM_THREADS", "unset")
        for r in rows:
            r.setdefault("n_threads", nth)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")
    print(f"host: {platform.node()} / {platform.platform()}")
    if args.compare:
        if Path(args.compare).exists():
            compare(rows, args.compare)
        else:
            print(f"(no baseline at {args.compare}; this run can become one)")


if __name__ == "__main__":
    main()
