"""
Track B / B11 demo -- GMRES+AMG on the level-set path reproduces spsolve.

The level-set drivers were hardcoded to sparse-direct `spsolve` (the "splu
wall"; roadmap B8 cost boundary). B11 adds `precond=None|"ilu"|"amg"`
(design_track_b.md §5.3): preconditioned GMRES (solve/linear.solve_gmres) on the
fused nonsymmetric matrix. This demo runs a subsonic M0.5 NACA0012 lifting solve
on the coarse and medium 2.5D meshes with spsolve vs ILU-GMRES and records, per
(mesh, precond): the extracted circulation, |Δgamma| vs the spsolve baseline,
total wall time, mean linear-solve time per outer, total GMRES iterations, and
stalls.

★ ILU is the shipped iterative escape on the LIFTING path: AMG built on the SPD
surrogate converges only on the SPD continuity/Laplace system, not on the
`wake_ls` lifting operator (its g1+g2 wake-LS + nonlinear TE-Kutta rows are
convection-like, not SPD) -- measured, GMRES stalls at the restart cap (coarse
M0.5, gamma 0.0033 vs 0.139, 455 s, all outers stalled). So this demo compares
spsolve vs ILU; AMG is exercised on the Laplace case in tests/test_b11_linear_ls.

Self-check: ILU reproduces spsolve to |Δgamma| < 1e-8 with 0 stalls.

Outputs: results/solver_ab.csv, results/gmres_ls_ab.png.
Exit 0 iff the identity holds. Standalone; runtime ~2-4 min (medium is ~1 min).
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, INK, MESH_DIR, S1_BLUE, S2_AQUA, S3_YELLOW,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import (  # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

OUT = Path(__file__).resolve().parent / "results"
NACA_DIR = MESH_DIR / "naca0012_2.5d"
LEVELS = ("coarse", "medium")
PRECONDS = (None, "ilu")   # AMG stalls on the wake_ls lifting operator (header)
ALPHA = 2.0
M_INF = 0.5
COLOR = {None: INK, "ilu": S1_BLUE, "amg": S2_AQUA}


def build(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
                       direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def main():
    apply_style()
    checks = CheckList("Track B / B11 -- GMRES+AMG vs spsolve on the LS path")
    rows = []
    data = {}   # (level, precond) -> dict
    for level in LEVELS:
        path = NACA_DIR / f"{level}.msh"
        if not path.exists():
            print(f"  skip {level}: {path} not generated (gitignored)")
            continue
        mesh = read_mesh(path)
        n_dofs = None
        ref_gamma = None
        for precond in PRECONDS:
            mvop = build(mesh)
            n_dofs = mvop.n_total
            t0 = time.time()
            r = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA,
                                          precond=precond)
            wall = time.time() - t0
            if precond is None:
                ref_gamma = r["gamma"]
            dgamma = abs(r["gamma"] - ref_gamma)
            n_outer = max(r["n_outer"], 1)
            ng = r["n_gmres_total"]
            data[(level, precond)] = {
                "n_dofs": n_dofs, "gamma": r["gamma"], "dgamma": dgamma,
                "wall": wall, "per_outer": wall / n_outer,
                "n_gmres": ng, "iters_outer": ng / n_outer,
                "stalled": r["n_gmres_stalled"], "converged": r["converged"],
            }
            tag = precond or "spsolve"
            print(f"  {level:6s} {tag:8s}: gamma={r['gamma']:+.6f} "
                  f"dg={dgamma:.1e} {wall:5.1f}s gmres_it={ng} "
                  f"stall={r['n_gmres_stalled']}")
            rows.append((level, tag, n_dofs, f"{r['gamma']:.8f}",
                         f"{dgamma:.2e}", f"{wall:.2f}", f"{wall/n_outer:.3f}",
                         ng, r["n_gmres_stalled"], int(r["converged"])))
            if precond is not None:
                checks.add("G11.3", f"{level}/{tag} == spsolve", f"{dgamma:.1e}",
                           "< 1e-8", dgamma < 1e-8)
                checks.add("G11.3", f"{level}/{tag} 0 stalls", r["n_gmres_stalled"],
                           "== 0", r["n_gmres_stalled"] == 0)

    write_csv(OUT, "solver_ab.csv",
              "level,precond,n_dofs,gamma,dgamma,wall_s,per_outer_s,"
              "n_gmres,n_stalled,converged", rows)

    # --- plot: per-outer linear cost + GMRES iters/outer -------------------
    if data:
        levels = [lv for lv in LEVELS if (lv, None) in data]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
        x = np.arange(len(levels))
        w = 0.25
        for i, precond in enumerate(PRECONDS):
            tag = precond or "spsolve"
            po = [data[(lv, precond)]["per_outer"] for lv in levels]
            ax1.bar(x + (i - 1) * w, po, w, color=COLOR[precond], label=tag)
        ax1.set_xticks(x)
        ax1.set_xticklabels(levels)
        ax1.set_ylabel("linear solve time / outer  [s]")
        ax1.set_title("cost per outer iteration")
        ax1.legend(frameon=False)
        for precond in [p for p in PRECONDS if p is not None]:
            io = [data[(lv, precond)]["iters_outer"] for lv in levels]
            ax2.plot(x, io, marker="o", color=COLOR[precond],
                     label=precond)
        ax2.set_xticks(x)
        ax2.set_xticklabels(levels)
        ax2.set_ylabel("GMRES iterations / outer")
        ax2.set_title("iterative convergence (warm-started)")
        ax2.legend(frameon=False)
        finish(fig, OUT, "gmres_ls_ab.png")

    return checks.report(OUT, "checks_gmres_ls.csv")


if __name__ == "__main__":
    sys.exit(main())
