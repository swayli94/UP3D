"""AUDIT 2026-07-28 / experiment 5 -- the linear-algebra floor.

How expensive is ONE elliptic solve on the production wing mesh, using the
machinery already in the repo? This sets the floor a well-designed
Newton-Krylov full-potential solver would run at (a transonic solve is
typically 20-40 equivalent Poisson solves), so the measured end-to-end cost
can be compared against its own algorithmic ideal rather than against a
different language.

Also times the ingredients the hot loop repeats: residual assembly, Jacobian
assembly, AMG setup, and -- for contrast -- a sparse direct LU of the same
matrix, which is what the production recipes use as a preconditioner
(`precond="direct"`).

Outputs: results/exp5_linalg.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from pyfp3d.kernels.jacobian import PicardOperator                # noqa: E402
from pyfp3d.mesh.reader import read_mesh                          # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                         # noqa: E402
from pyfp3d.solve.linear import (build_amg_preconditioner,        # noqa: E402
                                 solve_cg_amg)

MESHES = [
    ("naca_medium", "cases/meshes/naca0012_2.5d/medium.msh"),
    ("m6_coarse", "cases/meshes/onera_m6/coarse.msh"),
    ("m6_medium", "cases/meshes/onera_m6/medium.msh"),
    ("wingbody_medium", "cases/meshes/onera_m6_wingbody_conforming/medium.msh"),
]
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)


def timed(fn, n=1):
    t0 = time.perf_counter()
    for _ in range(n):
        out = fn()
    return (time.perf_counter() - t0) / n, out


def main():
    rows = []
    for tag, rel in MESHES:
        path = REPO / rel
        if not path.exists():
            print(f"skip {tag}: mesh missing")
            continue
        mesh = read_mesh(path)
        t_cut, (mc, wc) = timed(lambda: cut_wake(mesh))
        n = len(mc.nodes)
        print(f"\n=== {tag}: {n} nodes / {len(mc.elements)} tets ===",
              flush=True)
        op = PicardOperator(mc.nodes, mc.elements)
        phi = mc.nodes[:, 0].copy()
        rho = np.ones(len(mc.elements))

        # warm up numba
        op.assemble_residual(phi, rho)
        t_res, _ = timed(lambda: op.assemble_residual(phi, rho), 5)
        t_mat, A = timed(lambda: op.assemble_matrix(rho), 3)
        A = A.tocsr()
        nnz = A.nnz

        # Real Dirichlet elimination on the far-field boundary, exactly as
        # the solver's own Laplace path does -- the resulting SPD block is
        # the system whose solve cost we are pricing.
        from pyfp3d.solve.linear import apply_dirichlet
        ff = np.unique(mc.boundary_faces["farfield"])
        vals = mc.nodes[ff, 0].copy()
        rhs = np.zeros(n)
        Ash, b, free, _ = apply_dirichlet(A, rhs, ff, vals)
        Ash = Ash.tocsr()
        n = Ash.shape[0]

        t_amg_setup, (ml, M) = timed(lambda: build_amg_preconditioner(Ash))

        def _cg():
            nit = [0]
            x, info = spla.cg(Ash, b, M=M, rtol=1e-8, maxiter=2000,
                              callback=lambda _x: nit.__setitem__(0,
                                                                  nit[0] + 1))
            return x, nit[0], info

        t_cg, (x, its, info) = timed(_cg)
        print(f"   CG+AMG to rtol 1e-8: {its} iters, info={info}, "
              f"{t_cg:.2f} s", flush=True)
        try:
            t_lu, lu = timed(lambda: spla.splu(Ash.tocsc()))
            lu_nnz = lu.L.nnz + lu.U.nnz
            t_lu_solve, _ = timed(lambda: lu.solve(b), 3)
        except Exception as exc:                             # noqa: BLE001
            t_lu, lu_nnz, t_lu_solve = float("nan"), -1, float("nan")
            print("  splu failed:", exc)

        row = dict(case=tag, n_dof_free=n, n_nodes=len(mc.nodes), n_tets=len(mc.elements), nnz=nnz,
                   t_cut_wake=round(t_cut, 2),
                   t_residual_ms=round(1e3 * t_res, 2),
                   t_matrix_ms=round(1e3 * t_mat, 2),
                   t_amg_setup_s=round(t_amg_setup, 3),
                   t_cg_amg_s=round(t_cg, 3), cg_iters=its, cg_info=info,
                   t_splu_s=round(t_lu, 3), splu_fill_nnz=lu_nnz,
                   splu_fill_ratio=round(lu_nnz / nnz, 1) if lu_nnz > 0 else None,
                   t_splu_solve_ms=round(1e3 * t_lu_solve, 2))
        print("  ", row, flush=True)
        rows.append(row)

    keys = sorted({k for r in rows for k in r})
    with open(OUT / "exp5_linalg.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case"] + [k for k in keys
                                                     if k != "case"])
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT / "exp5_linalg.csv")


if __name__ == "__main__":
    main()
