"""
Track B / B12 headline -- lagged-LU makes the M6-medium LS Newton affordable.

B11 measured that the iterative escapes (ILU/AMG) do NOT scale on the fused
level-set matrix -- ILU diverges at 2.5D medium lifting and `factor_failed`s at
M6 medium, AMG stalls throughout -- so at medium/M6 sizes sparse-direct is the
only converging tool and the cost driver becomes the NUMBER of factorizations
(spsolve at 67,426 dofs = 17.5 s per factorization; roadmap G11.4). LS Newton
with `precond=None` factorizes on EVERY step, so on M6 medium it pays that
17.5 s once per Newton iteration.

This gated demo runs ONE M6-medium subsonic (M0.5) LS Newton A/B from a shared
Picard seed: `direct_refactor_every=1` (spsolve every step) vs a large k
(lagged-LU -- refactor once, then reuse the stale exact LU as a near-perfect
GMRES preconditioner). It commits the timing/agreement CSV that is the "M6 LS
Newton runs" evidence (G12.3).

★ Honest boundary (roadmap G12.3): at 67k dofs a single splu fits in memory, so
this is a REAL runnable medium-scale win -- not an extrapolation. But lagged-LU
still needs at least one splu that fits in memory, so it does NOT break the FINE
mesh memory wall (P9's 26 GB / 4h39m is per-factorization, not per-count); that
regime remains the Nunez symmetric-row-assignment -> AMG route (design_track_b
§5.3, not prebuilt).

Gated: set PYFP3D_TRANSONIC_GATES=1 (or pass --gated). Cap threads at 16
INCLUDING BLAS/OMP (the 16C/32T box; oversubscription costs ~33%):
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \\
    PYFP3D_TRANSONIC_GATES=1 python cases/demo/b12_lagged_lu/run_b12_m6_newton.py

Runtime ~8-15 min (one Picard seed + two M6-medium Newton solves). Requires the
M1 mesh: python cases/meshes/onera_m6/generate_onera_m6.py --level medium

Outputs: results/m6_newton_ab.csv, results/checks_m6_newton.csv.
Exit 0 iff the lagged solve converges, matches spsolve |Δgamma| < 1e-6, and
refactors fewer times than it takes Newton steps.
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import CheckList, MESH_DIR, write_csv  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te  # noqa: E402
from pyfp3d.solve.newton_ls import solve_multivalued_newton  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import (  # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

OUT = Path(__file__).resolve().parent / "results"
M1_DIR = MESH_DIR / "onera_m6"
LEVEL = "medium"
ALPHA = 3.06
M_INF = 0.5
LAGGED_K = 1000        # refactor once per solve; reuse the stale LU thereafter


def build(mesh):
    a = np.radians(ALPHA)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _newton(mesh, phi_seed, gamma_seed, k):
    """One M6-medium LS Newton solve from the shared seed. Returns
    (result, wall_s, n_dofs)."""
    mvop = build(mesh)
    t0 = time.time()
    r = solve_multivalued_newton(
        mvop, mesh, M_INF, alpha_deg=ALPHA, farfield="neumann",
        phi_init=phi_seed, gamma_init=gamma_seed, n_seed=2,
        n_newton_max=20, direct_refactor_every=k)
    return r, time.time() - t0, mvop.n_total


def main():
    gated = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1" \
        or "--gated" in sys.argv
    if not gated:
        print("SKIP: gated. Set PYFP3D_TRANSONIC_GATES=1 (see module header).")
        return 0
    path = M1_DIR / f"{LEVEL}.msh"
    if not path.exists():
        print(f"SKIP: {path} not generated (gitignored; see module header).")
        return 0

    checks = CheckList("Track B / B12 -- M6 medium LS Newton: lagged-LU escape")
    mesh = read_mesh(path)

    # shared Picard seed (spsolve; the cost is paid once and excluded from the
    # Newton A/B, which is about the factorization COUNT in the Newton phase)
    print("seeding (Picard-LS, M0.5, neumann) ...")
    t0 = time.time()
    seed = solve_multivalued_lifting(build(mesh), mesh, M_INF, alpha_deg=ALPHA,
                                     farfield="neumann", n_outer_max=15,
                                     tol_residual=1e-6)
    print(f"  seed: gamma={seed['gamma']:+.6f} wall={time.time()-t0:.1f}s")

    rows = []
    results = {}
    for k in (1, LAGGED_K):
        r, wall, n_dofs = _newton(mesh, seed["phi_ext"], seed["gamma"], k)
        results[k] = (r, wall)
        tag = "spsolve" if k == 1 else f"lagged(k={k})"
        print(f"  {tag:14s}: gamma={r['gamma']:+.8f} n_dofs={n_dofs} "
              f"newton={r['n_newton']} refactor={r['n_refactor']} "
              f"gmres={r['n_gmres_total']} stall={r['n_gmres_stalled']} "
              f"lim/flr={r['n_limited']}/{r['n_floored']} "
              f"wall={wall:.1f}s conv={r['converged']}")
        rows.append((tag, n_dofs, f"{r['gamma']:.8f}", r["n_newton"],
                     r["n_refactor"], r["n_gmres_total"], r["n_gmres_stalled"],
                     r["n_limited"], r["n_floored"], f"{wall:.2f}",
                     int(r["converged"])))

    write_csv(OUT, "m6_newton_ab.csv",
              "solver,n_dofs,gamma,n_newton,n_refactor,n_gmres,n_stalled,"
              "n_limited,n_floored,wall_s,converged", rows)

    ref, ref_wall = results[1]
    lag, lag_wall = results[LAGGED_K]
    dgamma = abs(lag["gamma"] - ref["gamma"])
    print(f"  |Δgamma| lagged vs spsolve = {dgamma:.2e}; "
          f"wall {ref_wall:.1f}s -> {lag_wall:.1f}s "
          f"(refactor {ref['n_refactor']}->{lag['n_refactor']} over "
          f"{lag['n_newton']} Newton steps)")

    checks.add("G12.3", "M6 medium spsolve Newton converged (baseline)",
               "converged" if ref["converged"] else "no", "converged",
               ref["converged"])
    checks.add("G12.3", "M6 medium lagged-LU Newton converged",
               "converged" if lag["converged"] else "no", "converged",
               lag["converged"])
    checks.add("G12.3", "lagged 0 limited / 0 floored",
               f"{lag['n_limited']}/{lag['n_floored']}", "0/0",
               lag["n_limited"] == 0 and lag["n_floored"] == 0)
    checks.add("G12.3", "lagged == spsolve", f"{dgamma:.2e}", "< 1e-6",
               dgamma < 1e-6)
    checks.add("G12.3", "lagged 0 GMRES stalls (stale LU is near-exact)",
               lag["n_gmres_stalled"], "== 0", lag["n_gmres_stalled"] == 0)
    checks.add("G12.3", "fewer factorizations than Newton steps",
               f"{lag['n_refactor']} < {lag['n_newton']}", "n_refactor < n_newton",
               lag["n_refactor"] < lag["n_newton"])
    return checks.report(OUT, "checks_m6_newton.csv")


if __name__ == "__main__":
    sys.exit(main())
