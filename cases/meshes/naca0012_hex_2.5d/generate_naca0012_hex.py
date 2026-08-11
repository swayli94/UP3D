"""HYBRID body-fitted NACA0012, 2.5-D (phase 3 task 2, route A -- v2).

Structured wall layers of FINITE height + unstructured far field: the architecture the user
specified and the production one. It replaced a v1 O-grid that stretched one grid all the way to
r_far; the quality history, measured, is what justifies each change:

    v1  O-grid to r_far, cosine surface        AR max 1305.4   min area 1.12e-07
    v2  finite layer block + gmsh far field     AR max   47.0   min area 7.62e-07
    v2b + LE-clustered / TE-COARSENED surface   AR max    5.9   min area 5.08e-06

★ v2b is the user's second point: `naca0012_coordinates` cosine-clusters BOTH ends, so the TE was
as finely spaced as the LE (ratio 0.81). On a sharp TE that buys nothing -- the region is nearly
straight, the Kutta condition rides on the wake -- while it drove the aspect ratio and, worse,
made every densely packed last-5%-of-chord station rotate its marching direction at once under
the TE +/-y lock. Redistributing to TE/LE = 7.05 at the SAME station count also made the LE
FINER (0.005249 -> 0.003251): the points moved from where they did nothing to where they matter.

★★ Sharp-TE constraint (user): the TE is a single point at (1, 0) with zero thickness, where the
averaged normal points along +x -- marching it would send the TE column straight down the wake
line. So the growth direction blends to +y on the upper side and -y on the lower over the last
5 % of chord, exactly +/-y at the TE, and the TE node is split for marching then WELDED for
topology. Verified: LE direction (-1, 0), TE (0, +1) / (0, -1) exactly.

    python cases/meshes/naca0012_hex_2.5d/generate_naca0012_hex.py --all

Outputs: <level>.msh + <level>_stats.csv + <level>_layer.png (tracked, like the unstructured
NACA family).
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import mesh_stats, write_mesh                    # noqa: E402
from pyfp3d.meshgen.extrude import extrude_single_layer                  # noqa: E402
from pyfp3d.meshgen.planar import naca0012_coordinates                  # noqa: E402
from pyfp3d.meshgen.structured import (airfoil_hybrid_2d,               # noqa: E402
                                       airfoil_layer_block_2d,
                                       airfoil_surface_distribution,
                                       layer_block_quality)

OUT_DIR = Path(__file__).resolve().parent
#: FIXED across the ladder so the refinement axis is wall resolution alone (G0 makes that real)
R_FAR, GROWTH, TE_BLEND, LE_CLUSTER = 20.0, 1.15, 0.05, 0.8
BLOCK_THICKNESS, H_FAR = 0.08, 2.0
LEVELS = {
    "coarse": dict(n_stations=80, h_wall_normal=0.008, dz=0.20),
    "medium": dict(n_stations=160, h_wall_normal=0.004, dz=0.10),
    "fine": dict(n_stations=320, h_wall_normal=0.002, dz=0.05),
}


def surface_for(level):
    p = LEVELS[level]
    #: resample from a FINE reference polyline so the distribution is set by the spacing law,
    #: not by whatever n_half happens to give
    return airfoil_surface_distribution(naca0012_coordinates(n_half=401),
                                        p["n_stations"], le_cluster=LE_CLUSTER)


def build(level: str):
    p = LEVELS[level]
    surf = surface_for(level)
    pts, tris, edges, interior, info = airfoil_hybrid_2d(
        surf, thickness=BLOCK_THICKNESS, h_wall_normal=p["h_wall_normal"], growth=GROWTH,
        te_blend=TE_BLEND, r_far=R_FAR, h_far=H_FAR)
    mesh = extrude_single_layer(pts, tris, edges, interior_edge_groups=interior,
                                dz=p["dz"], name=f"naca0012_hex_{level}")
    #: block-only quality, measured on the block alone where the (layer, station) indexing is
    #: unambiguous -- after the hybrid merge welds a node that mapping no longer holds, which is
    #: a probe error I made once already
    bp, bt, bo, bw, binfo = airfoil_layer_block_2d(
        surf, thickness=BLOCK_THICKNESS, h_wall_normal=p["h_wall_normal"], growth=GROWTH,
        te_blend=TE_BLEND)
    q = layer_block_quality(bp, bt, bw, binfo["n_stations"], binfo["n_layers"])
    return mesh, info, q, (pts, tris, edges, interior)


def write_inspection_png(path, pts, tris, edges, interior, level):
    """Headless artifact (hard rule). Panels chosen for the two things this route claims and the
    one it must not break: the far field, the TE lock, and the LE layer."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    views = [(0.5, 0.0, 21.0, "full domain: structured block + unstructured far field"),
             (0.8, 0.0, 0.30, "TE: layers lock to +/-y, block truncated, wake = green"),
             (0.0, 0.0, 0.05, "LE zoom: orthogonal graded layer")]
    for ax, (xc, yc, lim, ttl) in zip(axes, views):
        ax.triplot(pts[:, 0], pts[:, 1], tris, linewidth=0.4, color="0.55")
        for e, col, lw in ((edges["farfield"], "tab:blue", 1.6),
                           (edges["wall"], "tab:red", 1.8),
                           (interior.get("wake", np.zeros((0, 2), int)), "tab:green", 1.6)):
            for seg in (pts[np.asarray(e)] if len(e) else []):
                ax.plot(seg[:, 0], seg[:, 1], color=col, linewidth=lw)
        ax.set_xlim(xc - lim, xc + lim); ax.set_ylim(yc - lim, yc + lim)
        ax.set_aspect("equal"); ax.set_title(ttl, fontsize=10)
    fig.suptitle(f"naca0012_hex_2.5d {level}", fontsize=11)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", action="append", choices=sorted(LEVELS))
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    levels = sorted(LEVELS) if a.all else (a.level or ["coarse", "medium"])
    for level in levels:
        mesh, info, q, planar = build(level)
        write_mesh(mesh, OUT_DIR / f"{level}.msh")
        write_inspection_png(OUT_DIR / f"{level}_layer.png", *planar, level)
        st = mesh_stats(mesh)
        row = {"level": level}
        row.update({k: info[k] for k in sorted(info)})
        row.update({f"quality_{k}": v for k, v in q.items()})
        row.update(dz=LEVELS[level]["dz"], le_cluster=LE_CLUSTER,
                   n_nodes=len(mesh.nodes), n_tets=len(mesh.elements),
                   n_tris_wall=len(mesh.boundary_faces["wall"]),
                   n_tris_wake=len(mesh.boundary_faces.get("wake", [])))
        row.update({f"stats_{k}": v for k, v in st.items() if isinstance(v, (int, float))})
        with open(OUT_DIR / f"{level}_stats.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row)); w.writeheader(); w.writerow(row)
        print(f"  {level:7} stations {info['n_stations']:>4} layers {info['n_layers']:>3} "
              f"nodes {len(mesh.nodes):>7} tets {len(mesh.elements):>8} "
              f"wall {row['n_tris_wall']:>5} wake {row['n_tris_wake']:>4}  "
              f"AR max {q['aspect_max']:6.1f}  min area {q['min_area']:.2e}  "
              f"self-int {q['ray_crossings']}  seam {info['n_seam_nodes_matched']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
