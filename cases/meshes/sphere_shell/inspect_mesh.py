"""
Headless visual inspection of a sphere-shell mesh (design.md/roadmap.md
tooling policy: every visual check needs a script-driven, artifact-based
path -- no GUI-only tools).

PyVista's off-screen rendering needs OSMesa/EGL, which isn't available on
this machine (no GPU, no libOSMesa) -- so this uses matplotlib instead:
  1. wall surface triangulation (3D)
  2. a thin z~0 slice through the volume mesh, showing the radial grading
     from wall to far field (2D, edges only)
  3. a CSV of mesh-quality stats (node/tet counts, aspect ratio, volume)

Usage:
    python inspect_mesh.py coarse.msh
    python inspect_mesh.py medium.msh
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.metrics import compute_tet_volumes, compute_aspect_ratios


def inspect(msh_path: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    mesh = read_mesh(msh_path, verbose=True)
    nodes, elements = mesh.nodes, mesh.elements

    volumes = compute_tet_volumes(nodes, elements)
    aspect = compute_aspect_ratios(nodes, elements)
    wall_faces = mesh.boundary_faces.get("wall")
    farfield_faces = mesh.boundary_faces.get("farfield")
    wall_nodes = np.unique(wall_faces) if wall_faces is not None else np.array([])
    farfield_nodes = np.unique(farfield_faces) if farfield_faces is not None else np.array([])

    stats = {
        "n_nodes": len(nodes),
        "n_tets": len(elements),
        "n_wall_nodes": len(wall_nodes),
        "n_farfield_nodes": len(farfield_nodes),
        "total_volume": volumes.sum(),
        "aspect_ratio_min": aspect.min(),
        "aspect_ratio_mean": aspect.mean(),
        "aspect_ratio_max": aspect.max(),
        "aspect_ratio_p99": np.percentile(aspect, 99),
    }
    csv_path = out_dir / f"{msh_path.stem}_stats.csv"
    with open(csv_path, "w") as f:
        f.write("metric,value\n")
        for k, v in stats.items():
            f.write(f"{k},{v}\n")
    print(f"Wrote {csv_path}")
    for k, v in stats.items():
        print(f"  {k}: {v}")

    # --- Plot 1: wall surface triangulation ---
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    tris = nodes[wall_faces]
    coll = Poly3DCollection(tris, facecolor="tab:red", edgecolor="k", linewidths=0.2, alpha=0.9)
    ax.add_collection3d(coll)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title(f"{msh_path.name}: wall surface ({len(wall_faces)} triangles)")
    fig.savefig(out_dir / f"{msh_path.stem}_wall_surface.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 2: thin z~0 slice through the volume mesh (radial grading) ---
    z_tol = 0.03 * (nodes[:, 2].max() - nodes[:, 2].min())
    edges_local = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    segments = []
    for tet in elements:
        tet_nodes = nodes[tet]
        if np.all(np.abs(tet_nodes[:, 2]) > z_tol):
            continue
        for a, b in edges_local:
            pa, pb = tet_nodes[a], tet_nodes[b]
            if abs(pa[2]) <= z_tol and abs(pb[2]) <= z_tol:
                segments.append([(pa[0], pa[1]), (pb[0], pb[1])])

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, r_lim, title in [
        (axes[0], None, f"{msh_path.name}: z~0 slice (full domain)"),
        (axes[1], 4.0, f"{msh_path.name}: z~0 slice (zoom near wall, r<4)"),
    ]:
        lc = matplotlib.collections.LineCollection(segments, linewidths=0.3, colors="steelblue")
        ax.add_collection(lc)
        circle = plt.Circle((0, 0), 1.0, fill=False, color="tab:red", linewidth=1.5)
        ax.add_patch(circle)
        lim = r_lim if r_lim is not None else nodes[:, 0].max() * 1.05
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(title)
    fig.savefig(out_dir / f"{msh_path.stem}_slice.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_dir / f'{msh_path.stem}_wall_surface.png'}")
    print(f"Wrote {out_dir / f'{msh_path.stem}_slice.png'}")


if __name__ == "__main__":
    mesh_dir = Path(__file__).parent
    names = sys.argv[1:] or ["coarse.msh", "medium.msh"]
    for name in names:
        inspect(mesh_dir / name, mesh_dir)
