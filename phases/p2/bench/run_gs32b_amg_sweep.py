"""GS3.2 route (b): can pyamg parameters alone fix the extruded-mesh anisotropy?

Pre-registered in phases/p2/docs/dev_phase_two/20260802-1000-gs32b-prereg.md, committed before this
file. Same protocol as run_gs32_premise.py -- same meshes, same Laplacian, same Dirichlet
set, same rtol, same seeded right-hand side -- so the numbers are directly comparable to
that round's 2.5-D 278 / 2-D 14.

Criterion: under 50 iterations on 2.5-D medium meets GS3.2's own AMG clause. And TOTAL
time (setup + CG) is recorded alongside, because halving the iterations while doubling
the setup is a fake win.

Outputs (TRACKED): bench/gate_results/gs32b_amg_sweep.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402
import pyamg                                                        # noqa: E402
import scipy.sparse.linalg as spla                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from run_gs32_premise import (assemble_tet_laplacian,               # noqa: E402
                              assemble_tri_laplacian, tri_from_symmetry)

OUT = os.path.join(_GATE)
os.makedirs(OUT, exist_ok=True)
RTOL, MAXIT = 1e-8, 2000

#: (tag, kwargs for smoothed_aggregation_solver). "default" reproduces today's
#: build_amg_preconditioner exactly, so it is the control.
CONFIGS = [
    ("default", {}),
    ("theta0.0", {"strength": ("symmetric", {"theta": 0.0})}),
    ("theta0.25", {"strength": ("symmetric", {"theta": 0.25})}),
    ("theta0.5", {"strength": ("symmetric", {"theta": 0.5})}),
    ("evolution", {"strength": "evolution"}),
    ("affinity", {"strength": "affinity"}),
    ("energy", {"smooth": "energy"}),
    ("energy+theta0.25", {"smooth": "energy",
                          "strength": ("symmetric", {"theta": 0.25})}),
    ("gs2sweep", {"presmoother": ("gauss_seidel", {"sweep": "symmetric",
                                                   "iterations": 2}),
                  "postsmoother": ("gauss_seidel", {"sweep": "symmetric",
                                                    "iterations": 2})}),
    ("blockgs", {"presmoother": ("block_gauss_seidel",
                                 {"sweep": "symmetric"}),
                 "postsmoother": ("block_gauss_seidel",
                                  {"sweep": "symmetric"})}),
    ("lloyd", {"aggregate": "lloyd"}),
    ("evolution+energy", {"strength": "evolution", "smooth": "energy"}),
]


def run_one(A, tag, kwargs):
    state = np.random.get_state()
    np.random.seed(0)
    try:
        t0 = time.perf_counter()
        ml = pyamg.smoothed_aggregation_solver(A, **kwargs)
        t_setup = time.perf_counter() - t0
    except Exception as exc:                                      # noqa: BLE001
        np.random.set_state(state)
        return dict(config=tag, note=f"setup failed: {type(exc).__name__}: {exc}")
    finally:
        np.random.set_state(state)
    M = ml.aspreconditioner()
    rng = np.random.default_rng(0)
    b = A @ rng.standard_normal(A.shape[0])
    it = [0]
    t0 = time.perf_counter()
    try:
        x, info = spla.cg(A, b, M=M, rtol=RTOL, maxiter=MAXIT,
                          callback=lambda _: it.__setitem__(0, it[0] + 1))
    except Exception as exc:                                      # noqa: BLE001
        return dict(config=tag, amg_setup_s=round(t_setup, 4),
                    note=f"cg failed: {type(exc).__name__}: {exc}")
    t_cg = time.perf_counter() - t0
    res = float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))
    return dict(config=tag, dofs=A.shape[0], amg_setup_s=round(t_setup, 4),
                cg_iters=it[0], cg_s=round(t_cg, 4),
                total_s=round(t_setup + t_cg, 4), info=int(info),
                rel_res=res, levels=len(ml.levels),
                complexity=round(float(ml.operator_complexity()), 3), note="")


def main():
    rows = []
    #: ★ the 3-D wing meshes are swept too. The premise round framed this as a
    #: 2.5-D disease, but the root cause turned out to be pyamg's DEFAULT
    #: strength-of-connection theta = 0.0, which sees no anisotropy anywhere -- so
    #: the 3-D meshes (13-16 iterations on the default) may improve as well, and
    #: that bears directly on M3b's under-60 s target.
    CASES = [("naca0012_2.5d", "medium"), ("naca0012_2.5d", "coarse"),
             ("onera_m6", "coarse"), ("onera_m6", "medium")]
    for family, level in CASES:
        path = os.path.join(REPO, "cases", "meshes", family, f"{level}.msh")
        if not os.path.exists(path):
            print(f"skip {family}/{level}: mesh missing")
            continue
        level = f"{family}/{level}"
        mesh = read_mesh(path)
        K3, _ = assemble_tet_laplacian(mesh.nodes, mesh.elements)
        far = np.unique(mesh.boundary_faces["farfield"].reshape(-1))
        free = np.setdiff1d(np.arange(K3.shape[0]), far)
        A = K3[free][:, free].tocsr()
        # ★ select on the FAMILY, not on the presence of a "symmetry" group: M6 has
        # one too (its root plane) but was never extruded, so lifting a "2-D
        # reference" off it compares against a mesh the 3-D problem was not built
        # from. The first version keyed on the group name and printed exactly that
        # meaningless reference for M6 -- caught by the mislabelled header.
        extruded = family.endswith("_2.5d")
        if not extruded:
            print(f"\n=== {level}: {A.shape[0]} dofs (3-D, no 2-D reference) ===",
                  flush=True)
            best = None
            for tag, kw in CONFIGS:
                r = run_one(A, tag, kw)
                rows.append(dict(level=level, dim="3D", **r))
                if r.get("note"):
                    print(f"  {tag:18s} {r['note']}", flush=True); continue
                if best is None or r["total_s"] < best["total_s"]:
                    best = r
                print(f"  {tag:18s} {r['cg_iters']:5d} it  setup "
                      f"{r['amg_setup_s']:6.3f}s  cg {r['cg_s']:6.3f}s  TOTAL "
                      f"{r['total_s']:6.3f}s  cplx {r['complexity']}", flush=True)
            if best:
                d = next(r for r in rows if r["level"] == level
                         and r.get("config") == "default")
                print(f"  ★ best by TOTAL: {best['config']} {best['total_s']:.3f}s "
                      f"vs default {d['total_s']:.3f}s = "
                      f"{d['total_s']/best['total_s']:.2f}x; iterations "
                      f"{d['cg_iters']} -> {best['cg_iters']}", flush=True)
            continue

        # the 2-D reference on the same triangulation, for context
        pts, tri, nodes2d = tri_from_symmetry(mesh)
        K2, _ = assemble_tri_laplacian(pts, tri)
        s = set(far.tolist())
        d2 = np.array([i for i, n in enumerate(nodes2d) if n in s])
        f2 = np.setdiff1d(np.arange(K2.shape[0]), d2)
        A2 = K2[f2][:, f2].tocsr()
        ref2 = run_one(A2, "2-D reference (default)", {})
        print(f"\n=== {level} (extruded): 2.5-D A is {A.shape[0]} dofs; "
              f"2-D reference {ref2.get('cg_iters')} it / "
              f"{ref2.get('total_s')}s ===", flush=True)
        rows.append(dict(level=level, dim="2D", **ref2))
        best = None
        for tag, kw in CONFIGS:
            r = run_one(A, tag, kw)
            rows.append(dict(level=level, dim="2.5D", **r))
            if r.get("note"):
                print(f"  {tag:18s} {r['note']}", flush=True)
                continue
            mark = ""
            if r["cg_iters"] < 50:
                mark = "  ★ < 50 (meets GS3.2's AMG clause)"
            if best is None or r["total_s"] < best["total_s"]:
                best = r
            print(f"  {tag:18s} {r['cg_iters']:5d} it  setup {r['amg_setup_s']:6.3f}s"
                  f"  cg {r['cg_s']:6.3f}s  TOTAL {r['total_s']:6.3f}s"
                  f"  lvls {r['levels']} cplx {r['complexity']}{mark}", flush=True)
        if best:
            d = next(r for r in rows if r["level"] == level
                     and r.get("config") == "default")
            print(f"  ★ best by TOTAL time: {best['config']} "
                  f"{best['total_s']:.3f}s vs default {d['total_s']:.3f}s "
                  f"= {d['total_s']/best['total_s']:.2f}x; iterations "
                  f"{d['cg_iters']} -> {best['cg_iters']}", flush=True)

    with open(os.path.join(OUT, "gs32b_amg_sweep.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", os.path.join(OUT, "gs32b_amg_sweep.csv"))

    print("\n=== the registered reading (2.5-D medium) ===")
    med = [r for r in rows if r["level"] == "medium" and r["dim"] == "2.5D"
           and r.get("cg_iters")]
    if not med:
        print("  no usable configuration")
        return 0
    b = min(med, key=lambda r: r["cg_iters"])
    print(f"  fewest iterations: {b['config']} at {b['cg_iters']} "
          f"(default 278-ish; 2-D reference 14)")
    if b["cg_iters"] < 50:
        print("  => (b) DELIVERS GS3.2's AMG clause by parameters alone")
    elif b["cg_iters"] < 150:
        print("  => PARTIAL: recorded; (a)'s trigger 1 becomes more likely")
    else:
        print("  => (b) FAILED: not reachable by pyamg parameters. (a) or a "
              "hand-written semi-coarsening is the only route, and (a)'s "
              "trigger 1 is met.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
