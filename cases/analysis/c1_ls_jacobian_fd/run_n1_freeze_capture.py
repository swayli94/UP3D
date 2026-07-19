"""
N1 discriminator — does fixing the freeze CAPTURE restore the M6 medium M0.84 ramp?

GB20.7 (commit 37a0799) answered the freeze-tol axis: sweeping freeze_tol
1e-3..1e-6 moves the post-B20 ceiling only 0.6625 -> 0.675, so the stall is not
a recipe mismatch on that axis. The Kimi 2026-07-19 inspection (N1) found the
OTHER axis that sweep structurally cannot see: `freeze_side_state` captured the
(upstream, branch) selection on the UNPATCHED side field, while
`newton_side_data` applies `_apply_main_density` first — so on 3-D meshes the
frozen Newton finish iterated on a selection the live system would not make
(measured: 83 upstream + 9 branch differences on M6 coarse at a seeded M0.70
state, all on aux-touching mixed-plain elements;
docs/inspection/20260719-n1-freeze-probe.py).

That capture is now patched (multivalued.py `freeze_side_state`, B20/N1). This
script re-runs the SAME committed B15 M6-medium call as run_gb207_recipe.py,
post-fix, at the committed freeze_tol=1e-3 and the sweep's best 1e-5:

  reaches M0.84            ->  N1 WAS the missing mechanism; GB20.7's "real
                               capability loss" verdict is overturned.
  still stalls ~M0.66-0.68 ->  the loss is real AND the last known code-side
                               suspect is cleared; GB20.7 formally closes.

Heavy: each variant is a full M6-medium Mach ramp (~22-31 min). Run:
  python cases/analysis/c1_ls_jacobian_fd/run_n1_freeze_capture.py
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
SWEEP = [1e-3, 1e-5]     # committed value + the GB20.7 sweep's best


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
        tag = "COMMITTED" if ft == 1e-3 else "GB207BEST"
        print(f"freeze_tol={ft:.0e} {tag:9s} reached={r['reached']} "
              f"m_final={r['m_final']:.4f} levels={r['levels_converged']} "
              f"gamma={r['gamma']:.6f} M_max={r['m_max']:.4f} "
              f"res={r['res']:.1e} clamps={r['n_lim']}/{r['n_flr']} "
              f"({r['wall_s']:.0f}s) accept=[{r['accept']}]", flush=True)

    reached = [r for r in rows if r["reached"]]
    print()
    if reached:
        b = reached[0]
        print(f"VERDICT: N1 WAS THE MECHANISM — with the freeze capture "
              f"aligned to the B20 density, freeze_tol={b['freeze_tol']:.0e} "
              f"reaches M0.84 again (gamma {b['gamma']:.6f}, M_max "
              f"{b['m_max']:.4f}). GB20.7's 'real capability loss' verdict is "
              f"overturned; GB15.4 needs re-measurement, not a re-spec.")
    else:
        best = max(rows, key=lambda r: r["m_final"])
        print(f"VERDICT: STILL STALLS (best m_final {best['m_final']:.4f}) — "
              "the capability loss survives the N1 capture fix. The last "
              "known code-side suspect is cleared; GB20.7's 'real capability "
              "loss' verdict formally closes (bounds: freeze_tol swept by "
              "GB20.7, capture content fixed here; dm/n_newton_max/m_start "
              "still unswept).")

    out = RESULTS / "n1_freeze_fix_sweep.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
