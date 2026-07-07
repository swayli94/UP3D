"""Where are the M>2 cells in the cached P5 medium solution? (evidence for
the 2026-07-08 re-diagnosis of the OPEN medium gate; docs/roadmap.md P5).

This reads ONLY the committed cache results/<level>_solution.npz + the mesh
-- no solve, seconds to run -- and reproduces the cell-by-cell localization
that corrected the earlier "all far-field, not on the wing" reading:

  * 18 of the 26 M>2 cells sit ON the wing at the outboard trailing edge
    (dist-to-wall ~ 0, x/c ~ 1.0-1.1, z/b ~ 0.80); only 8 are at the
    far-field sphere. The single M_max=5.20 cell IS far-field.
  * the outboard-TE cells are well-shaped (shape quality ~0.65, aspect
    ratio <= 3.1) -- NOT slivers -- and sit at the steepest spanwise-Gamma
    roll-off (medium |dGamma/dz| peaks at z/b=0.80). => a trailing-edge /
    Kutta discretization singularity fed by per-station Kutta spanwise noise,
    NOT a far-field-BC-only artifact and NOT global under-convergence.

Usage:  python cases/demo/p5_onera_m6/diagnose_medium.py [coarse|medium]
Requires the .msh (regenerate via cases/meshes/onera_m6/generate_onera_m6.py)
and the committed *_solution.npz. Cap threads: NUMBA_NUM_THREADS=16.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, chord_at, x_le  # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field  # noqa: E402
from pyfp3d.post.surface import _element_gradients_and_centroids  # noqa: E402

RES = Path(__file__).parent / "results"
M_INF = 0.84


def _tet_shape_quality(nodes, elements):
    """Normalized shape quality 6*sqrt(2)*V / l_rms^3 (=1 regular tet, ->0
    sliver) and a crude edge aspect ratio, per element."""
    p = nodes[elements]
    v = np.abs(np.einsum("ni,ni->n",
                         np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                         p[:, 3] - p[:, 0])) / 6.0
    e = np.stack([p[:, i] - p[:, j] for i, j in
                  [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]], axis=1)
    el = np.linalg.norm(e, axis=2)
    q = 6.0 * np.sqrt(2.0) * v / (np.sqrt((el ** 2).mean(1)) ** 3 + 1e-300)
    return q, el.max(1) / np.maximum(el.min(1), 1e-30)


def diagnose(level: str):
    mesh = read_mesh(ROOT / f"cases/meshes/onera_m6/{level}.msh")
    mc, wc = cut_wake(mesh)
    d = np.load(RES / f"{level}_solution.npz")
    phi = d["phi"]

    grad, cen = _element_gradients_and_centroids(mc.nodes, mc.elements, phi)
    mach = np.sqrt(np.maximum(mach_squared_field(np.sum(grad * grad, 1), M_INF), 0.0))
    q, ar = _tet_shape_quality(mc.nodes, mc.elements)

    wall = np.unique(np.asarray(mc.boundary_faces["wall"], dtype=np.int64))
    ff = np.unique(np.asarray(mc.boundary_faces["farfield"], dtype=np.int64))
    d_wall = cKDTree(mc.nodes[wall]).query(cen)[0]
    d_ff = cKDTree(mc.nodes[ff]).query(cen)[0]

    print(f"\n===== {level}: {len(mc.elements)} tets, M_max={mach.max():.3f} "
          f"(stored {np.sqrt(float(d['mach2_max'])):.3f}) =====")
    idx = np.where(mach > 2.0)[0]
    if len(idx) == 0:
        print(f"  no M>2 cells (mesh median shape quality {np.median(q):.2f})")
    else:
        c = cen[idx]
        z = c[:, 2]
        xc = (c[:, 0] - np.array([x_le(v) for v in z])) / np.array([chord_at(v) for v in z])
        on_wing = d_wall[idx] < 0.3 * MAC
        at_ff = d_ff[idx] < 0.5 * MAC
        print(f"  {len(idx)} cells M>2:  {on_wing.sum()} ON wing (dwall<0.3MAC), "
              f"{at_ff.sum()} at far-field sphere (dff<0.5MAC)")
        print(f"    on-wing:  x/c {xc[on_wing].min():.2f}-{xc[on_wing].max():.2f} (TE if ~1), "
              f"z/b {(z[on_wing] / B_SEMI).min():.2f}-{(z[on_wing] / B_SEMI).max():.2f}")
        print(f"    shape quality of M>2 cells: median {np.median(q[idx]):.2f} "
              f"min {q[idx].min():.2f} (mesh median {np.median(q):.2f}); "
              f"aspect ratio max {ar[idx].max():.1f}  => not slivers")
        e_hot = idx[np.argmax(mach[idx])]
        print(f"    M_max cell: M={mach[e_hot]:.2f} at "
              f"(x={cen[e_hot, 0]:.1f}, z/b={cen[e_hot, 2] / B_SEMI:.2f}), "
              f"dff={d_ff[e_hot]:.2f} => {'FAR-FIELD' if d_ff[e_hot] < 0.5 * MAC else 'wing'}")

    # spanwise Gamma roll-off (steepest |dGamma/dz| co-locates with the cluster)
    o = np.argsort(wc.station_z)
    zz, gg = wc.station_z[o], d["gamma"][o]
    dg = np.gradient(gg, zz)
    k = np.argmax(-dg)
    print(f"  steepest spanwise Gamma roll-off at z/b={zz[k] / B_SEMI:.2f} "
          f"(dGamma/dz={dg[k]:.2f}); Gamma root {gg[0]:.3f} -> tip {gg[-1]:.4f}")


if __name__ == "__main__":
    levels = sys.argv[1:] or ["coarse", "medium"]
    for lvl in levels:
        diagnose(lvl)
