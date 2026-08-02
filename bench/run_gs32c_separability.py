"""GS3.2 (c): do any of S1-S4 separate "wants evolution" from "wants the default"?"""
import os, sys
os.environ.setdefault("NUMBA_NUM_THREADS","16")
sys.path.insert(0, "/home/lrz/codes/UP3D"); sys.path.insert(0, "/home/lrz/codes/UP3D/bench")
import numpy as np, scipy.sparse as sp
from pyfp3d.mesh.reader import read_mesh
from run_gs32_premise import assemble_tet_laplacian

WANT = {"naca0012_2.5d/medium": "evolution", "naca0012_2.5d/coarse": "evolution",
        "onera_m6/coarse": "default", "onera_m6/medium": "default"}

def stats(A):
    A = A.tocsr()
    n = A.shape[0]
    d = np.abs(A.diagonal())
    out = {}
    s1, s2, s3, s4 = [], [], [], []
    ip, dat = A.indptr, A.data
    idx = A.indices
    for i in range(n):
        sl = slice(ip[i], ip[i+1])
        j, v = idx[sl], np.abs(dat[sl])
        off = j != i
        v = v[off]
        v = v[v > 0]
        if v.size < 2:
            continue
        mx, mn, tot = v.max(), v.min(), v.sum()
        s1.append(mx / mn)
        s2.append(mx / tot)
        s3.append(float(np.count_nonzero(v >= 0.25 * mx)) / v.size)
        s4.append(d[i] / tot if tot > 0 else np.inf)
    q = lambda a: (float(np.median(a)), float(np.quantile(a, 0.9)))
    out["S1_max_over_min"] = q(s1)
    out["S2_max_over_sum"] = q(s2)
    out["S3_frac_strong_t025"] = q(s3)
    out["S4_diag_dominance"] = q(s4)
    return out

res = {}
for key, want in WANT.items():
    fam, lvl = key.split("/")
    mesh = read_mesh(f"/home/lrz/codes/UP3D/cases/meshes/{fam}/{lvl}.msh")
    K, _ = assemble_tet_laplacian(mesh.nodes, mesh.elements)
    far = np.unique(mesh.boundary_faces["farfield"].reshape(-1))
    free = np.setdiff1d(np.arange(K.shape[0]), far)
    res[key] = stats(K[free][:, free])
    print(f"{key:26s} want={want:9s} " +
          "  ".join(f"{k.split('_')[0]}={v[0]:.3f}" for k, v in res[key].items()),
          flush=True)

print("\n=== separation (median), evolution-group vs default-group ===")
ev = [k for k, w in WANT.items() if w == "evolution"]
df = [k for k, w in WANT.items() if w == "default"]
for stat in res[ev[0]]:
    a = [res[k][stat][0] for k in ev]
    b = [res[k][stat][0] for k in df]
    lo_a, hi_a, lo_b, hi_b = min(a), max(a), min(b), max(b)
    overlap = not (hi_a < lo_b or hi_b < lo_a)
    gap = (lo_b / hi_a) if hi_a < lo_b else ((lo_a / hi_b) if hi_b < lo_a else 1.0)
    v = ("OVERLAP" if overlap else
         ("SEPARATES >=3x" if gap >= 3 else
          ("too narrow (<1.5x)" if gap < 1.5 else "RECORDED 1.5-3x")))
    print(f"  {stat:22s} evo [{lo_a:.3f},{hi_a:.3f}]  def [{lo_b:.3f},{hi_b:.3f}]"
          f"   gap {gap:.2f}x   {v}")
