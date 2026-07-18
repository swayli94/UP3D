"""
GB20.7 — is the M6-medium-M0.84 stall a RECIPE mismatch or a real capability loss?

B20 made the level-set density read the main field, which removed a spurious
supersonic contamination. Every 3-D case improved except one: the M6 MEDIUM
M0.84 ramp, which used to reach the target and now stalls at M0.6625 (2/5
levels) — with **0/0 clamped cells at every level** and |R| 9.2e-14 where it
does converge. It is neither diverging nor being clamped.

The hypothesis this script tests is written into B15's own recipe comment:

    freeze_tol=1e-3 sits ABOVE the churn floor (which RISES with Mach:
    <1e-6 / 8.6e-6 / 2.7e-4 at M0.60 / 0.65 / 0.70)

**That churn floor was measured on the CONTAMINATED field.** With the field now
clean, the floor should have dropped — and a freeze_tol far above it arms the
freeze while the selection is still far from settled, locking a bad assignment.
The observed symptoms fit exactly: freeze armed at every level, zero reverts,
levels accepted via `assignment_cycle` (the intrinsic-floor escape) rather than
by reaching tol, then no ability to climb.

So: sweep freeze_tol DOWNWARD and see whether the target comes back.

  reaches M0.84 at some freeze_tol  ->  RECIPE MISMATCH (the recipe was
                                        calibrated against the contamination)
  fails at every freeze_tol         ->  a REAL capability loss, which the user
                                        must weigh against B20's gains

Heavy: each variant is a full M6-medium Mach ramp (~22 min at the stalling
value; a converging one was ~11 min pre-B20). Run:
  python cases/analysis/c1_ls_jacobian_fd/run_gb207_recipe.py
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import read_mesh                        # noqa: E402
from pyfp3d.meshgen.wing3d import x_te                          # noqa: E402
from pyfp3d.solve.newton_ls import (                            # noqa: E402
    B_NEWTON_M6_DEFAULTS, solve_multivalued_newton_transonic,
)
from pyfp3d.wake import (                                       # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

RESULTS = Path(__file__).resolve().parent / "results"
MESH = REPO_ROOT / "cases" / "meshes" / "onera_m6_wakefree" / "medium.msh"
ALPHA, B_SEMI = 3.06, 1.1963
# B15's committed M6-medium call (run_demo.py Part 3), minus freeze_tol
BASE = dict(m_target=0.84, alpha_deg=ALPHA, farfield="neumann",
            n_seed=40, n_newton_max=80, tol_residual=1e-10)
SWEEP = [1e-3, 1e-4, 1e-5, 1e-6]     # 1e-3 = the committed value (the stall)


def build():
    mesh = read_mesh(str(MESH))
    a = np.radians(ALPHA)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def run(freeze_tol, mesh):
    # fresh operator per variant: the upwind/side caches are state-dependent
    _, mvop = build()
    kw = dict(B_NEWTON_M6_DEFAULTS)
    kw["freeze_tol"] = freeze_tol
    t0 = time.perf_counter()
    r = solve_multivalued_newton_transonic(mvop=mvop, mesh=mesh, **BASE, **kw)
    wall = time.perf_counter() - t0
    lv = r["levels"]
    n_conv = sum(1 for l in lv if l.get("converged"))
    mmax = float(np.sqrt(np.max(mvop.element_mach2(r["phi_ext"],
                                                   r["m_final"], 1.4, 1.0))))
    return dict(
        freeze_tol=freeze_tol, reached=bool(r["target_reached"]),
        m_final=r["m_final"], levels_converged=f"{n_conv}/{len(lv)}",
        gamma=float(r.get("gamma", np.nan)), m_max=mmax,
        res=float(r["residual_history"][-1]),
        n_lim=int(lv[-1]["n_limited"]), n_flr=int(lv[-1]["n_floored"]),
        accept=";".join(sorted({(l.get("accept_reason") or "not-converged")
                                for l in lv})),
        wall_s=round(wall, 1))


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not MESH.exists():
        raise SystemExit(f"{MESH} missing (gitignored); regenerate the M6 "
                         "wake-free family first")
    mesh, _ = build()
    rows = []
    for ft in SWEEP:
        r = run(ft, mesh)
        rows.append(r)
        tag = "COMMITTED" if ft == 1e-3 else ""
        print(f"freeze_tol={ft:.0e} {tag:9s} reached={r['reached']} "
              f"m_final={r['m_final']:.4f} levels={r['levels_converged']} "
              f"gamma={r['gamma']:.6f} M_max={r['m_max']:.4f} "
              f"res={r['res']:.1e} clamps={r['n_lim']}/{r['n_flr']} "
              f"({r['wall_s']:.0f}s) accept=[{r['accept']}]", flush=True)

    reached = [r for r in rows if r["reached"]]
    print()
    if reached:
        b = min(reached, key=lambda r: r["wall_s"])
        print(f"VERDICT: RECIPE MISMATCH — freeze_tol={b['freeze_tol']:.0e} "
              f"reaches M0.84 again (gamma {b['gamma']:.6f}, M_max "
              f"{b['m_max']:.4f}, {b['wall_s']:.0f}s). B15's committed "
              f"freeze_tol=1e-3 was calibrated against the contaminated churn "
              f"floor and is now too high.")
    else:
        print("VERDICT: NO recipe in this sweep reaches M0.84 — the stall "
              "survives freeze_tol 1e-3..1e-6. On this evidence it is a REAL "
              "capability loss on M6 medium, to be weighed against B20's "
              "gains elsewhere (NOT proof: only freeze_tol was varied).")

    out = RESULTS / "gb207_recipe_sweep.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
