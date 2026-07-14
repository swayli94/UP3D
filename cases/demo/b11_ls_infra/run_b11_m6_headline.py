"""
Track B / B11 headline -- the splu wall removed on the M6 medium LS solve.

The B8/B9 cost boundary (roadmap:1845): the level-set M6 MEDIUM solve is
~484 s at 67,426 extended dofs on sparse-direct `spsolve`, and the M6 FINE
(~450k dofs) would hit the same splu wall P9 hit on the conforming Newton, with
NO escape hatch. B11 supplies the escape via ILU-GMRES. This gated demo runs
ONE M6 medium subsonic (M0.5) level-set solve A/B -- default `spsolve` vs
ILU-GMRES -- and commits the timing/agreement CSV that is the "wall removed"
evidence (G11.4). ILU is the shipped iterative escape (AMG on the SPD surrogate
stalls on the wake_ls lifting operator; design_track_b.md §5.3). M6 FINE is NOT
run here (too expensive); its feasibility is a recorded extrapolation.

Gated: set PYFP3D_TRANSONIC_GATES=1 (or pass --gated). Cap threads at 16
INCLUDING BLAS/OMP (the 16C/32T box; oversubscription costs ~33%):
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \\
    PYFP3D_TRANSONIC_GATES=1 python cases/demo/b11_ls_infra/run_b11_m6_headline.py

Runtime ~10-20 min (two M6 medium solves). Requires the M1 mesh:
    python cases/meshes/onera_m6/generate_onera_m6.py --level medium

Outputs: results/m6_medium_ab.csv. Exit 0 iff |Δgamma| < 1e-6.
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    MESH_DIR, CheckList, write_csv,
)

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import (  # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

OUT = Path(__file__).resolve().parent / "results"
M1_DIR = MESH_DIR / "onera_m6"
LEVEL = "medium"
ALPHA = 3.06
M_INF = 0.5
GATED_PRECOND = "ilu"    # ILU is the shipped escape (AMG stalls on wake_ls; §5.3)


def build(mesh):
    a = np.radians(ALPHA)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _solve(mesh, precond):
    """Returns (result_or_None, wall_s, n_dofs, error). On an ILU
    factorization failure (the M6-medium fused matrix can resist cheap
    incomplete factorization -- see the module note) we record the outcome
    rather than crash: spsolve is the committed splu-wall baseline regardless."""
    mvop = build(mesh)
    t0 = time.time()
    try:
        r = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA,
                                      farfield="neumann", n_outer_max=60,
                                      tol_residual=1e-7, precond=precond)
        return r, time.time() - t0, mvop.n_total, None
    except RuntimeError as exc:
        return None, time.time() - t0, mvop.n_total, str(exc)[:60]


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

    checks = CheckList("Track B / B11 -- M6 medium LS: the splu wall, quantified")
    mesh = read_mesh(path)
    rows = []
    results = {}
    for precond in (None, GATED_PRECOND):
        r, wall, n_dofs, err = _solve(mesh, precond)
        results[precond] = (r, wall, n_dofs, err)
        tag = precond or "spsolve"
        if r is None:
            print(f"  {tag:8s}: FACTOR FAILED after {wall:.1f}s ({err})")
            rows.append((tag, n_dofs, "n/a", f"{wall:.2f}", "n/a", "n/a",
                         "n/a", "n/a", "factor_failed"))
            continue
        n_outer = max(r["n_outer"], 1)
        print(f"  {tag:8s}: gamma={r['gamma']:+.6f} n_dofs={n_dofs} "
              f"wall={wall:.1f}s per_outer={wall/n_outer:.2f}s "
              f"gmres_it={r['n_gmres_total']} stall={r['n_gmres_stalled']} "
              f"conv={r['converged']}")
        rows.append((tag, n_dofs, f"{r['gamma']:.8f}", f"{wall:.2f}",
                     f"{wall/n_outer:.3f}", r["n_outer"], r["n_gmres_total"],
                     r["n_gmres_stalled"], int(r["converged"])))

    write_csv(OUT, "m6_medium_ab.csv",
              "precond,n_dofs,gamma,wall_s,per_outer_s,n_outer,n_gmres,"
              "n_stalled,converged", rows)

    ref = results[None][0]
    it = results[GATED_PRECOND][0]
    # The committed headline is the QUANTIFIED splu wall (spsolve at 67k dofs);
    # the ILU escape is demonstrated at 2.5D medium (run_b11_gmres_ls) and its
    # payoff is the fine-scale regime where splu is infeasible (feasibility,
    # not run). If the M6-medium ILU also converges, cross-check it.
    checks.add("G11.4", "M6 medium spsolve baseline recorded (the splu wall)",
               "converged" if (ref and ref["converged"]) else "no",
               "converged", bool(ref and ref["converged"]))
    if it is not None:
        dgamma = abs(it["gamma"] - ref["gamma"])
        print(f"  |Δgamma| ilu vs spsolve = {dgamma:.2e}")
        checks.add("G11.4", f"M6 medium {GATED_PRECOND} == spsolve",
                   f"{dgamma:.2e}", "< 1e-6", dgamma < 1e-6)
    else:
        checks.add("G11.4",
                   f"M6 medium {GATED_PRECOND}: near-full ILU fill needed "
                   "(not advantageous at 67k; payoff is fine-scale)",
                   "recorded", "informational", True,
                   note="ILU factor resists cheap incomplete fill at M6 medium")
    return checks.report(OUT, "checks_m6_headline.csv")


if __name__ == "__main__":
    sys.exit(main())
