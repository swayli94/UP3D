"""GS3.2 premise check: is the 480-iteration AMG count caused by the EXTRUSION?

Pre-registered in docs/dev_phase_two/20260802-0600-gs32-prereg.md, committed before this
file was written. The audit charged the 2.5-D mesh's bad AMG behaviour (480 CG iterations
against 13-16 on the 3-D wing meshes) to the one-prism-layer extrusion. But a 2.5-D
mesh's anisotropy has two sources -- the IN-PLANE grading, which a true 2-D triangle mesh
shares identically, and the Z layer, which only the 2.5-D mesh has -- and the audit did
not separate them. If the in-plane grading is the cause, GS3.2 cannot deliver its promised
AMG improvement.

The discriminator needs no 2-D solver: assemble a P1 TRIANGLE stiffness matrix on the
SAME 2-D triangulation, taken straight off the 2.5-D mesh's z = 0 symmetry plane, and run
the same AMG + CG. Pure Laplacian with far-field Dirichlet on both sides, so the only
difference is dimension.

Outputs (TRACKED): bench/gate_results/gs32_premise.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402
import scipy.sparse as sp                                           # noqa: E402
import scipy.sparse.linalg as spla                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)

from pyfp3d.solve.linear import build_amg_preconditioner          # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry         # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402

OUT = os.path.join(_GATE)
os.makedirs(OUT, exist_ok=True)
RTOL = 1e-8


def tri_from_symmetry(mesh):
    """The 2-D triangulation the 2.5-D mesh was extruded FROM: the symmetry faces on
    the lower z plane, re-indexed to their own node set."""
    z = mesh.nodes[:, 2]
    z0 = z.min()
    f = mesh.boundary_faces["symmetry"]
    keep = np.all(np.isclose(z[f], z0, atol=1e-12), axis=1)
    tri = f[keep]
    nodes2d = np.unique(tri.reshape(-1))
    remap = -np.ones(len(mesh.nodes), dtype=np.int64)
    remap[nodes2d] = np.arange(len(nodes2d))
    return mesh.nodes[nodes2d][:, :2], remap[tri], nodes2d


def assemble_tri_laplacian(pts, tri):
    """Standard P1 triangle stiffness. K_ij = A * (grad_i . grad_j)."""
    p0, p1, p2 = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    e0, e1 = p1 - p0, p2 - p0
    det = e0[:, 0] * e1[:, 1] - e0[:, 1] * e1[:, 0]
    area = 0.5 * np.abs(det)
    # barycentric gradients: b_i = rot90(opposite edge) / det
    g = np.empty((len(tri), 3, 2))
    g[:, 0, 0] = (p1[:, 1] - p2[:, 1]) / det
    g[:, 0, 1] = (p2[:, 0] - p1[:, 0]) / det
    g[:, 1, 0] = (p2[:, 1] - p0[:, 1]) / det
    g[:, 1, 1] = (p0[:, 0] - p2[:, 0]) / det
    g[:, 2, 0] = (p0[:, 1] - p1[:, 1]) / det
    g[:, 2, 1] = (p1[:, 0] - p0[:, 0]) / det
    ke = area[:, None, None] * np.einsum("eik,ejk->eij", g, g)
    r = np.repeat(tri, 3, axis=1).reshape(-1)
    c = np.tile(tri, (1, 3)).reshape(-1)
    return sp.coo_matrix((ke.reshape(-1), (r, c)),
                         shape=(len(pts), len(pts))).tocsr(), area


def assemble_tet_laplacian(nodes, elements):
    B, V = precompute_element_geometry(nodes, elements)
    ke = V[:, None, None] * np.einsum("eik,ejk->eij", B, B)
    r = np.repeat(elements, 4, axis=1).reshape(-1)
    c = np.tile(elements, (1, 4)).reshape(-1)
    return sp.coo_matrix((ke.reshape(-1), (r, c)),
                         shape=(len(nodes), len(nodes))).tocsr(), V


def solve_amg_cg(K, dirichlet, tag):
    free = np.setdiff1d(np.arange(K.shape[0]), dirichlet)
    A = K[free][:, free].tocsr()
    rng = np.random.default_rng(0)
    b = A @ rng.standard_normal(A.shape[0])      # a consistent, deterministic rhs
    t0 = time.perf_counter()
    _, M = build_amg_preconditioner(A)
    t_setup = time.perf_counter() - t0
    it = [0]

    def count(_):
        it[0] += 1
    t0 = time.perf_counter()
    x, info = spla.cg(A, b, M=M, rtol=RTOL, maxiter=2000, callback=count)
    t_cg = time.perf_counter() - t0
    res = float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))
    print(f"  {tag:28s} dofs {A.shape[0]:6d}  AMG setup {t_setup:.3f}s  "
          f"CG {it[0]:4d} it / {t_cg:.3f}s  info={info}  rel_res {res:.2e}",
          flush=True)
    return dict(case=tag, dofs=A.shape[0], amg_setup_s=round(t_setup, 4),
                cg_iters=it[0], cg_s=round(t_cg, 4), info=int(info),
                rel_res=res)


def edge_stats(pts, cells, label):
    """Edge-length spread -- the audit's '300:1'. Split in-plane vs z where 3-D."""
    P = pts[cells]
    n = cells.shape[1]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    d = np.stack([P[:, i] - P[:, j] for i, j in pairs], axis=1)
    if P.shape[2] == 3:
        L_ip = np.linalg.norm(d[:, :, :2], axis=2)
        L_z = np.abs(d[:, :, 2])
        pos = L_ip[L_ip > 0]
        print(f"  {label}: IN-PLANE edge {pos.min():.5f}..{pos.max():.5f} "
              f"(ratio {pos.max()/pos.min():.0f}:1)   "
              f"Z edge {L_z[L_z>0].min():.5f}..{L_z.max():.5f}")
        return dict(inplane_min=float(pos.min()), inplane_max=float(pos.max()),
                    inplane_ratio=float(pos.max() / pos.min()),
                    z_min=float(L_z[L_z > 0].min()), z_max=float(L_z.max()))
    L = np.linalg.norm(d, axis=2)
    pos = L[L > 0]
    print(f"  {label}: edge {pos.min():.5f}..{pos.max():.5f} "
          f"(ratio {pos.max()/pos.min():.0f}:1)")
    return dict(inplane_min=float(pos.min()), inplane_max=float(pos.max()),
                inplane_ratio=float(pos.max() / pos.min()))


def main():
    rows = []
    for level in ("medium", "coarse"):
        path = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d",
                            f"{level}.msh")
        if not os.path.exists(path):
            print(f"skip {level}: mesh missing")
            continue
        mesh = read_mesh(path)
        print(f"\n=== NACA 2.5-D {level}: {len(mesh.nodes)} nodes / "
              f"{len(mesh.elements)} tets ===", flush=True)

        # --- 3-D (the audit's case) --------------------------------------
        K3, _ = assemble_tet_laplacian(mesh.nodes, mesh.elements)
        d3 = np.unique(mesh.boundary_faces["farfield"].reshape(-1))
        g3 = edge_stats(mesh.nodes, mesh.elements, "2.5-D tets ")
        r3 = solve_amg_cg(K3, d3, f"2.5-D tets ({level})")

        # --- 2-D on the SAME triangulation -------------------------------
        pts, tri, nodes2d = tri_from_symmetry(mesh)
        K2, _ = assemble_tri_laplacian(pts, tri)
        far3 = set(np.unique(mesh.boundary_faces["farfield"].reshape(-1)))
        d2 = np.array([i for i, n in enumerate(nodes2d) if n in far3])
        g2 = edge_stats(pts, tri, "2-D tris  ")
        print(f"  2-D triangulation: {len(pts)} nodes / {len(tri)} tris, "
              f"{len(d2)} far-field Dirichlet")
        r2 = solve_amg_cg(K2, d2, f"2-D tris ({level})")

        rows += [dict(level=level, dim="2.5D", **g3, **r3),
                 dict(level=level, dim="2D", **g2, **r2)]
        print(f"  ★ {level}: 2.5-D {r3['cg_iters']} it -> 2-D {r2['cg_iters']} it "
              f"({r3['cg_iters']/max(r2['cg_iters'],1):.1f}x fewer)")

    with open(os.path.join(OUT, "gs32_premise.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", os.path.join(OUT, "gs32_premise.csv"))

    print("\n=== the registered reading ===")
    for level in ("medium", "coarse"):
        a = next((r for r in rows if r["level"] == level and r["dim"] == "2.5D"),
                 None)
        b = next((r for r in rows if r["level"] == level and r["dim"] == "2D"),
                 None)
        if not a or not b:
            continue
        n = b["cg_iters"]
        print(f"  {level}: 2.5-D {a['cg_iters']} it, 2-D {n} it")
        if n < 50:
            print("    => attribution HOLDS: the extrusion is the cause and "
                  "GS3.2's premise is true")
        elif n >= 200:
            print("    => attribution FALSIFIED: the in-plane grading carries it, "
                  "and 2-D shares that grading. GS3.2 cannot deliver the AMG "
                  "improvement it promises; the gate needs rewriting.")
        else:
            print("    => RECORDED (50-200): both sources contribute; judge by "
                  "the ratio whether the build is worth it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
