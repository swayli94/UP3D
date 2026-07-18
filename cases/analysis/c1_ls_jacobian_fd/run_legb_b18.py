"""
B20 GB20.5 — THE HYPOTHESIS TEST (pre-registered in B19 Leg B).

B18 recorded a level-set wing-body spurious supersonic pocket that WORSENS with
refinement: the medium case dies at the first transonic level with an
element_mach2 (main-read) Mmax artifact of 3.964, nlim 43 / nflr 40. B19 Leg B
found that the SOLVER's density on mixed-side plain elements reads a spurious
supersonic side field (q^2 3.22 vs 1.34). Named test: rerun that exact case
with plain_density="main" and see whether the pocket moves.

This does NOT change any committed result -- it runs the B18 LS medium ramp
twice (side, main) on a fresh operator and reports Mmax / m_final / clamps.
A drop turns the hypothesis into a mechanism; no change bounds it out.

Heavy (medium wing-body LS Newton ramp). Run:
  python cases/analysis/c1_ls_jacobian_fd/run_legb_b18.py
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.meshgen.fuselage import FuselageParams               # noqa: E402
from pyfp3d.meshgen.wingbody import te_polyline                  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                         # noqa: E402
from pyfp3d.post.surface import planform_area                    # noqa: E402
from pyfp3d.post.surface_ls import cl_pressure_3d_levelset       # noqa: E402
from pyfp3d.solve.newton_ls import solve_multivalued_newton_transonic  # noqa: E402
from pyfp3d.wake import (                                        # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

RESULTS = Path(__file__).resolve().parent / "results"
LS_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6_wingbody"
ALPHA = 3.06
# B18's exact LS ramp recipe (run_demo.py LS_RAMP_KW) + case (ls_medium_06)
LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma", freeze_tol=1e-4,
                  freeze_max_clamped=8, intermediate_tol=1e-3, n_seed=30,
                  direct_refactor_every=1000, n_newton_max=80)
M_TARGET, M_START = 0.6, 0.5


def op_for(mode, mesh):
    a = np.radians(ALPHA)
    wls = WakeLevelSet(te_polyline(FuselageParams()),
                       direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls,
                               plain_density=mode)


def run(mode, mesh):
    mvop = op_for(mode, mesh)
    n_mp = int(mvop._mixed_plain_mask().sum())
    t0 = time.perf_counter()
    r = solve_multivalued_newton_transonic(mvop, mesh, M_TARGET, alpha_deg=ALPHA,
                                           m_start=M_START, dm=0.05, **LS_RAMP_KW)
    wall = time.perf_counter() - t0
    mf = r["m_final"]
    # Mmax read the SAME way B18 does: element_mach2 default (mixed_plain="main")
    mmax = float(np.sqrt(np.max(mvop.element_mach2(r["phi_ext"], mf, 1.4, 1.0))))
    lastlv = r["levels"][-1]
    return dict(mode=mode, n_mixed_plain=n_mp, m_final=mf,
                m_last=r["m_last_converged"], reached=bool(r["target_reached"]),
                mmax=mmax, nlim=int(lastlv["n_limited"]),
                nflr=int(lastlv["n_floored"]),
                res=float(r["residual_history"][-1]), wall_s=wall)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    p = LS_DIR / "medium.msh"
    if not p.exists():
        raise SystemExit(f"{p} missing (gitignored); regenerate the wing-body mesh")
    mesh = read_mesh(str(p))
    rows = []
    for mode in ("side", "main"):
        res = run(mode, mesh)
        rows.append(res)
        print(f"[{mode}] B18 medium LS ->M{M_TARGET}: m_final={res['m_final']} "
              f"m_last={res['m_last']} reached={res['reached']} "
              f"Mmax={res['mmax']:.3f} nlim={res['nlim']} nflr={res['nflr']} "
              f"res={res['res']:.1e} ({res['wall_s']:.0f}s) "
              f"mixed_plain={res['n_mixed_plain']}", flush=True)

    s, m = rows[0], rows[1]
    print(f"\nB18 committed baseline (side): Mmax 3.964, m_final 0.5, nlim 43/nflr 40")
    print(f"HYPOTHESIS TEST: Mmax {s['mmax']:.3f} (side) -> {m['mmax']:.3f} (main) "
          f"= {(m['mmax']-s['mmax'])/s['mmax']*100:+.1f}%; "
          f"m_final {s['m_final']} -> {m['m_final']}")
    out = RESULTS / "legb_b18_hypothesis.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
