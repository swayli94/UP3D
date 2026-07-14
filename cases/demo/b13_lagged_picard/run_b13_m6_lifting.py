"""
Track B / B13 headline -- lagged-LU makes the M6-medium LS Picard affordable.

After B12 (lagged-LU on the LS Newton loop), the M6-medium cost driver is the
PICARD outer loop: one 17.5 s spsolve per outer (B11 lifting headline 454.8 s /
26 outers; the B12 demo's Newton seed 263 s / 15 outers). B11 measured the
iterative escapes to be coarse-only on the fused matrix (ILU diverges at 2.5D
medium, factor_faileds at M6 medium; AMG stalls), so at these sizes
sparse-direct is the only converging tool and the cost is the NUMBER of
factorizations.

This gated demo runs the B11-headline M6-medium subsonic lifting solve A/B --
`direct_refactor_every=1` (spsolve every outer, the B11 baseline) vs a large k
(lagged-LU: refactor once, reuse the stale exact LU as a near-perfect GMRES
preconditioner, warm-started) -- and commits the timing/count CSV (GB13.3).
It then adds the END-TO-END row: Picard seed + LS Newton (the B12 demo
pipeline) with both mechanisms on, against the ~330 s post-B12 baseline.

★ Honest boundary: lagged-LU amortizes the factorization COUNT; it still needs
>= 1 in-memory splu, so it does NOT break the FINE-mesh memory wall (that
regime stays the recorded B14 Schur design / Nunez route, design_track_b §5.3).

Gated: set PYFP3D_TRANSONIC_GATES=1 (or pass --gated). Cap threads at 16
INCLUDING BLAS/OMP (the 16C/32T box; oversubscription costs ~33%):
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \\
    PYFP3D_TRANSONIC_GATES=1 python cases/demo/b13_lagged_picard/run_b13_m6_lifting.py

Runtime ~12-20 min (two M6-medium lifting solves + one seed+Newton pipeline).
Requires the M1 mesh: python cases/meshes/onera_m6/generate_onera_m6.py --level medium

Outputs: results/m6_lifting_ab.csv, results/m6_end_to_end.csv,
results/checks_m6_lifting.csv. Exit 0 iff the lagged solve converges, matches
spsolve |Δgamma| < 1e-6, and refactors far fewer times than it runs outers.
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
LAGGED_K = 1000        # refactor once, reuse the stale LU for the whole solve


def build(mesh):
    a = np.radians(ALPHA)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _lifting(mesh, k):
    """One M6-medium lifting solve (the B11 headline condition). Returns
    (result, wall_s, n_dofs)."""
    mvop = build(mesh)
    t0 = time.time()
    r = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA,
                                  farfield="neumann", n_outer_max=60,
                                  tol_residual=1e-7, direct_refactor_every=k)
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

    checks = CheckList("Track B / B13 -- M6 medium LS Picard: lagged-LU escape")
    mesh = read_mesh(path)

    # ---- GB13.3 part 1: the B11-headline lifting solve, A/B ---------------
    rows = []
    results = {}
    for k in (1, LAGGED_K):
        r, wall, n_dofs = _lifting(mesh, k)
        results[k] = (r, wall)
        tag = "spsolve" if k == 1 else f"lagged(k={k})"
        n_outer = max(r["n_outer"], 1)
        print(f"  {tag:14s}: gamma={r['gamma']:+.8f} n_dofs={n_dofs} "
              f"outers={r['n_outer']} refactor={r['n_refactor']} "
              f"gmres={r['n_gmres_total']} stall={r['n_gmres_stalled']} "
              f"wall={wall:.1f}s per_outer={wall/n_outer:.2f}s "
              f"conv={r['converged']}")
        rows.append((tag, n_dofs, f"{r['gamma']:.8f}", r["n_outer"],
                     r["n_refactor"], r["n_gmres_total"], r["n_gmres_stalled"],
                     f"{wall:.2f}", f"{wall/n_outer:.3f}",
                     int(r["converged"])))
    write_csv(OUT, "m6_lifting_ab.csv",
              "solver,n_dofs,gamma,n_outer,n_refactor,n_gmres,n_stalled,"
              "wall_s,per_outer_s,converged", rows)

    ref, ref_wall = results[1]
    lag, lag_wall = results[LAGGED_K]
    dgamma = abs(lag["gamma"] - ref["gamma"])
    print(f"  |Δgamma| lagged vs spsolve = {dgamma:.2e}; "
          f"wall {ref_wall:.1f}s -> {lag_wall:.1f}s "
          f"(refactor {ref['n_refactor']}->{lag['n_refactor']} over "
          f"{lag['n_outer']} outers)")

    # ---- GB13.3 part 2: end-to-end seed + Newton (the B12 pipeline) -------
    print("end-to-end (Picard seed + LS Newton, both lagged) ...")
    t0 = time.time()
    seed = solve_multivalued_lifting(build(mesh), mesh, M_INF, alpha_deg=ALPHA,
                                     farfield="neumann", n_outer_max=15,
                                     tol_residual=1e-6,
                                     direct_refactor_every=LAGGED_K)
    seed_wall = time.time() - t0
    t0 = time.time()
    newt = solve_multivalued_newton(
        build(mesh), mesh, M_INF, alpha_deg=ALPHA, farfield="neumann",
        phi_init=seed["phi_ext"], gamma_init=seed["gamma"], n_seed=2,
        n_newton_max=20, direct_refactor_every=LAGGED_K)
    newton_wall = time.time() - t0
    total = seed_wall + newton_wall
    print(f"  seed  : {seed_wall:.1f}s (refactor {seed['n_refactor']}/"
          f"{seed['n_outer']} outers)  newton: {newton_wall:.1f}s "
          f"(refactor {newt['n_refactor']}/{newt['n_newton']} steps)  "
          f"total {total:.1f}s  gamma={newt['gamma']:+.8f} "
          f"conv={newt['converged']}")
    write_csv(OUT, "m6_end_to_end.csv",
              "stage,wall_s,n_iters,n_refactor,gamma,converged",
              [("seed_lagged", f"{seed_wall:.2f}", seed["n_outer"],
                seed["n_refactor"], f"{seed['gamma']:.8f}",
                int(seed["converged"])),
               ("newton_lagged", f"{newton_wall:.2f}", newt["n_newton"],
                newt["n_refactor"], f"{newt['gamma']:.8f}",
                int(newt["converged"])),
               ("total", f"{total:.2f}", "", "", "", "")])

    # ---- checks ------------------------------------------------------------
    checks.add("GB13.3", "M6 medium spsolve lifting converged (baseline)",
               "converged" if ref["converged"] else "no", "converged",
               ref["converged"])
    checks.add("GB13.3", "M6 medium lagged-LU lifting converged",
               "converged" if lag["converged"] else "no", "converged",
               lag["converged"])
    checks.add("GB13.3", "lagged == spsolve", f"{dgamma:.2e}", "< 1e-6",
               dgamma < 1e-6)
    checks.add("GB13.3", "far fewer factorizations than outers",
               f"{lag['n_refactor']} vs {lag['n_outer']}",
               "n_refactor <= n_outer/4",
               lag["n_refactor"] <= lag["n_outer"] / 4)
    checks.add("GB13.3", "end-to-end (seed+Newton) converged",
               "converged" if newt["converged"] else "no", "converged",
               newt["converged"])
    checks.add("GB13.3", "Newton gamma in the B12 lock band",
               f"{newt['gamma']:.8f}", "|Δ| < 1e-4 of 0.06685284",
               abs(newt["gamma"] - 0.06685284) < 1e-4)
    return checks.report(OUT, "checks_m6_lifting.csv")


if __name__ == "__main__":
    sys.exit(main())
