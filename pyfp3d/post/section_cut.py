"""
z = const section extraction for slice plots and sectional Cp curves.

Full deliverable lands in P2 (roadmap P2 deliverables): given z = const,
extract (a) a field slice (phi, |V|, later M/rho/Cp) and (b) the
wall-surface Cp(x/c) split into upper/lower curves, all headless
(matplotlib Agg + CSV, roadmap §0.1).

This module is pre-warmed at the G1.3 stage (cylinder oracle pre-study)
with the FINAL interface signature -- the z parameter is already there --
but only the degenerate single-layer path is implemented: on an M0
quasi-2D mesh every node lies on one of the two symmetry planes, each of
which carries its own conforming 2D triangulation, so a "cut" at that z is
just a subset restriction and fields plot directly via tripcolor. P2 adds
the general z = const tet-interpolation path (and the wall-Cp curve
extraction); until then any other z raises NotImplementedError.
"""

from typing import Dict, Optional

import numpy as np


class SectionData:
    """A z = const section: 2D points, triangulation, and restricted fields."""

    def __init__(self, points2d: np.ndarray, triangles: np.ndarray,
                 fields: Dict[str, np.ndarray], z: float):
        self.points2d = points2d      # (n_pts, 2)
        self.triangles = triangles    # (n_tris, 3), indices into points2d
        self.fields = fields          # name -> (n_pts,) nodal values
        self.z = z


def section_cut(mesh, point_fields: Dict[str, np.ndarray], z: float,
                atol: Optional[float] = None) -> SectionData:
    """Extract the z = const section of nodal fields.

    Args:
        mesh: pyfp3d.mesh.reader.Mesh (needs nodes + boundary_faces)
        point_fields: name -> (n_nodes,) nodal arrays to restrict/interpolate
        z: section plane coordinate (final-interface parameter; in the
           degenerate path it must coincide with one of the mesh's node
           planes)
        atol: plane-matching tolerance (default: 1e-9 * z-extent, or 1e-12
              for a zero-extent direction)

    Returns:
        SectionData

    Raises:
        NotImplementedError: z falls between node planes (the general
            tet-interpolation path is a P2 deliverable).
    """
    nodes = mesh.nodes
    if atol is None:
        extent = float(np.ptp(nodes[:, 2]))
        atol = 1e-9 * extent if extent > 0 else 1e-12

    on_plane = np.abs(nodes[:, 2] - z) <= atol
    if not np.any(on_plane):
        raise NotImplementedError(
            f"z = {z} matches no node plane (within atol = {atol:g}); the "
            "general z = const interpolation path is a P2 deliverable"
        )

    # Degenerate path: the symmetry-plane triangulation restricted to this z.
    sym = mesh.boundary_faces.get("symmetry")
    if sym is None:
        raise NotImplementedError(
            "mesh has no 'symmetry' boundary group; only the single-layer "
            "quasi-2D degenerate path is implemented before P2"
        )
    tri_on_plane = sym[np.all(on_plane[sym], axis=1)]
    if len(tri_on_plane) == 0:
        raise NotImplementedError(
            f"no symmetry-plane triangles at z = {z}; only the single-layer "
            "quasi-2D degenerate path is implemented before P2"
        )

    plane_nodes = np.unique(tri_on_plane)
    remap = np.full(len(nodes), -1, dtype=np.int64)
    remap[plane_nodes] = np.arange(len(plane_nodes))

    return SectionData(
        points2d=nodes[plane_nodes, :2].copy(),
        triangles=remap[tri_on_plane],
        fields={name: np.asarray(f)[plane_nodes].copy()
                for name, f in point_fields.items()},
        z=z,
    )


def plot_section_field(section: SectionData, field_name: str, output_path,
                       title: str = "", cmap: str = "viridis",
                       xlim=None, ylim=None) -> None:
    """Headless tripcolor plot of one section field (PNG, matplotlib Agg)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 7))
    tpc = ax.tripcolor(section.points2d[:, 0], section.points2d[:, 1],
                       section.triangles, section.fields[field_name],
                       shading="gouraud", cmap=cmap)
    fig.colorbar(tpc, ax=ax, label=field_name)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title or f"{field_name} at z = {section.z:g}")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
