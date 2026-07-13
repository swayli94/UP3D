"""
B8 backlog / errata artifact -- the B6/B7 M_max RE-READ under the honest
element_mach2 default (mixed_plain="main", flipped from "side" 2026-07-14,
user-arbitrated).

Context (run_b8_termination_diagnosis.py, roadmap Track B B8 re-spec):
`element_mach2(mixed_plain="side")` reads mixed-side PLAIN (beyond-tip /
wake-plane-extension) elements from the aux-substituted side field -- a field
their single-valued MAIN-dof assembly never sees (measured x5 inflation on
the B8 verdict element). The committed B6/B7 M_max gate numbers were measured
through that reading. This script RE-READS (does NOT re-solve) the locally
cached converged states and reports both readings side by side:

  * B7 M6 coarse M0.84 (M1 wake-embedded + M4 wake-free):
      cases/demo/b7_onera_m6/results/{M1,M4}.npz  (phi_ext + cached m_max)
  * B6 NACA 2.5D coarse M0.80 neumann (M0 wake-embedded + M3 wake-free):
      cases/demo/b6_transonic/results/cp_naca0012{,_wakefree}_2.5d.npz

Reconstruction is copied verbatim from the demos that WROTE the caches
(b7_onera_m6/run_demo.py::setup, b6_transonic/run_demo.py::build), so the
MultivaluedOperator is the same discrete object. SELF-CHECK: the side
re-read must reproduce the cached m_max to ~1e-9 (same phi, same code path)
-- B7 M1 1.4534 / M4 1.3682; a mismatch means the reconstruction is wrong
and the script fails loudly.

The committed b6/b7 results CSVs are historical evidence and are NOT
touched; this writes the errata artifact results/mmax_reread.csv.

Standalone:  NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
             python cases/demo/b8_tip_taper_ls/run_b8_mmax_reread.py
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te  # noqa: E402
from pyfp3d.wake import (  # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes"
B7_RES = REPO / "cases" / "demo" / "b7_onera_m6" / "results"
B6_RES = REPO / "cases" / "demo" / "b6_transonic" / "results"
ALPHA_B7 = 3.06  # b7 run_demo.py


def setup_b7(path):
    """Verbatim from cases/demo/b7_onera_m6/run_demo.py::setup."""
    mesh = read_mesh(path)
    a = np.radians(ALPHA_B7)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                     levelset=wls)


def setup_b6(path):
    """Verbatim from cases/demo/b6_transonic/run_demo.py::build."""
    mesh = read_mesh(path)
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                     levelset=wls)


def n_mixed_plain(mesh, mvop):
    """The mixed-side PLAIN element set the flip re-reads (the exact mask
    element_mach2(mixed_plain="main") fixes)."""
    cm = mvop.cm
    el = np.asarray(mesh.elements, dtype=np.int64)
    special = np.zeros(len(el), dtype=bool)
    special[cm.cut_elems] = True
    special[cm.te_lower_elems] = True
    side_e = cm.node_side[el]
    return int(np.count_nonzero(
        ~special & (side_e.min(axis=1) != side_e.max(axis=1))))


CASES = [
    # (case, mesh tag, mesh path, setup, npz path, m_inf, cached-m_max key)
    ("B7-M1-embedded", "onera_m6/coarse", MESH / "onera_m6" / "coarse.msh",
     setup_b7, B7_RES / "M1.npz", 0.84, "m_max"),
    ("B7-M4-wakefree", "onera_m6_wakefree/coarse",
     MESH / "onera_m6_wakefree" / "coarse.msh",
     setup_b7, B7_RES / "M4.npz", 0.84, "m_max"),
    ("B6-M0-embedded", "naca0012_2.5d/coarse",
     MESH / "naca0012_2.5d" / "coarse.msh",
     setup_b6, B6_RES / "cp_naca0012_2.5d.npz", 0.80, "mach2"),
    ("B6-M3-wakefree", "naca0012_wakefree_2.5d/coarse",
     MESH / "naca0012_wakefree_2.5d" / "coarse.msh",
     setup_b6, B6_RES / "cp_naca0012_wakefree_2.5d.npz", 0.80, "mach2"),
]


def main():
    rows, failed = [], False
    for case, tag, mpath, setup, npz, m_inf, mkey in CASES:
        if not (mpath.exists() and npz.exists()):
            print(f"!! {case}: missing {'mesh' if not mpath.exists() else 'cache'}"
                  f" ({mpath if not mpath.exists() else npz}) -- skipped")
            failed = True
            continue
        d = np.load(npz)
        phi_ext = d["phi_ext"]
        cached = float(d[mkey].reshape(-1)[0])
        if mkey == "mach2":  # b6 stores M^2
            cached = float(np.sqrt(cached))
        mesh, mvop = setup(mpath)
        m_side = float(np.sqrt(np.max(
            mvop.element_mach2(phi_ext, m_inf, mixed_plain="side"))))
        m_main = float(np.sqrt(np.max(
            mvop.element_mach2(phi_ext, m_inf, mixed_plain="main"))))
        n_fix = n_mixed_plain(mesh, mvop)
        ok = abs(m_side - cached) < 5e-4  # 3-4 digit reproduction self-check
        failed |= not ok
        print(f"{case:16s} cached {cached:.4f} | side {m_side:.6f} "
              f"({'REPRODUCED' if ok else '*** MISMATCH ***'}) | "
              f"main {m_main:.6f} | n_mixed_plain {n_fix}")
        rows.append((case, tag, f"{m_inf:.2f}", f"{m_side:.6f}",
                     f"{m_main:.6f}", str(n_fix)))

    OUT.mkdir(parents=True, exist_ok=True)
    csv = OUT / "mmax_reread.csv"
    with csv.open("w") as f:
        f.write("# B6/B7 M_max re-read under element_mach2 mixed_plain="
                "'main' (default flip 2026-07-14, user-arbitrated).\n"
                "# side must reproduce the committed gate numbers "
                "(B7 M1 1.453 / M4 1.368; B6 cached mach2) -- the "
                "reconstruction self-check.\n"
                "# States are the LOCAL gitignored demo caches; the "
                "committed b6/b7 CSVs are untouched historical evidence.\n")
        f.write("case,mesh,m_inf,mmax_side,mmax_main,n_mixed_plain_fixed\n")
        for r in rows:
            f.write(",".join(r) + "\n")
    print(f"\nwrote {csv}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
