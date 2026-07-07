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

    # Degenerate path: the symmetry-plane triangulation restricted to this z
    # (exact and cheap when z coincides with a node plane of a quasi-2D mesh).
    sym = mesh.boundary_faces.get("symmetry")
    if np.any(on_plane) and sym is not None:
        tri_on_plane = sym[np.all(on_plane[sym], axis=1)]
        if len(tri_on_plane) > 0:
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

    # General path (P2): marching tetrahedra with linear interpolation.
    return _section_cut_marching(mesh, point_fields, z)


def _section_cut_marching(mesh, point_fields: Dict[str, np.ndarray],
                          z: float) -> SectionData:
    """z = const cut through the volume mesh by marching tetrahedra.

    Each tet crossed by the plane contributes one triangle (3 crossed
    edges) or two (4 crossed edges, split quad); cut-point coordinates and
    fields interpolate linearly along the crossed edges -- exact for the
    P1 solution. Vertices lying numerically on the plane are nudged to one
    side (standard epsilon trick) so the case table stays two-sided.
    """
    nodes = mesh.nodes
    elements = np.asarray(mesh.elements, dtype=np.int64)
    s = nodes[:, 2] - z
    extent = float(np.ptp(nodes[:, 2]))
    eps = 1e-12 * (extent if extent > 0 else 1.0)
    s = np.where(np.abs(s) < eps, eps, s)

    s_el = s[elements]
    pos = s_el > 0
    n_pos = pos.sum(axis=1)
    cut = (n_pos > 0) & (n_pos < 4)
    if not np.any(cut):
        raise ValueError(f"plane z = {z} does not intersect the mesh")

    # Tet edges as local index pairs; a cut edge has opposite signs.
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    names = list(point_fields)
    field_mat = np.column_stack([np.asarray(point_fields[n], dtype=np.float64)
                                 for n in names]) if names else None

    pts, tris, vals = [], [], []
    for e in np.where(cut)[0]:
        tet = elements[e]
        se = s_el[e]
        cut_pts = []
        for a, b in edges:
            if (se[a] > 0) != (se[b] > 0):
                t = se[a] / (se[a] - se[b])
                p = nodes[tet[a], :2] * (1 - t) + nodes[tet[b], :2] * t
                f = (field_mat[tet[a]] * (1 - t) + field_mat[tet[b]] * t
                     if field_mat is not None else None)
                cut_pts.append((p, f))
        base = len(pts)
        if len(cut_pts) == 3:
            local_tris = [(0, 1, 2)]
        elif len(cut_pts) == 4:
            # Order the quad by angle around its centroid, then fan-split.
            P = np.array([c[0] for c in cut_pts])
            ang = np.arctan2(*(P - P.mean(axis=0)).T[::-1])
            order = np.argsort(ang)
            cut_pts = [cut_pts[i] for i in order]
            local_tris = [(0, 1, 2), (0, 2, 3)]
        else:  # pragma: no cover - excluded by the epsilon nudge
            continue
        pts.extend(c[0] for c in cut_pts)
        vals.extend(c[1] for c in cut_pts)
        tris.extend((base + i, base + j, base + k) for i, j, k in local_tris)

    points2d = np.asarray(pts, dtype=np.float64)
    fields = {}
    if names:
        V = np.asarray(vals, dtype=np.float64)
        fields = {n: V[:, k].copy() for k, n in enumerate(names)}
    return SectionData(points2d=points2d,
                       triangles=np.asarray(tris, dtype=np.int64),
                       fields=fields, z=z)


def _wall_section_points(mesh, phi, z: float, u_inf: float,
                         upper_hint, wall_tag: str, m_inf: float):
    """Triangle-wise plane cut of the wall at z = const (design.md Sec 9).

    Each wall triangle crossed by the plane contributes one point at its
    intersection-segment midpoint, carrying the triangle's own constant Cp
    (the in-plane tangential gradient IS the wall velocity under the natural
    BC). No nodal averaging, so the sharp-TE crease needs no special-casing.
    A vertex exactly on the plane is nudged to the +side so every crossed
    triangle yields exactly two crossing points. Fully general in z: no wall
    node need lie on the plane.

    Returns absolute-x midpoints (no chord normalization), the per-point Cp,
    and a boolean upper-side mask (dot(mid, upper_hint) > 0).
    """
    from pyfp3d.physics.isentropic import (
        pressure_coefficient,
        pressure_coefficient_incompressible,
    )
    from pyfp3d.post.surface import triangle_tangential_gradients

    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    grad_tri, _, _ = triangle_tangential_gradients(mesh.nodes, wall, phi)
    q2 = np.sum(grad_tri * grad_tri, axis=1) / u_inf**2
    hint = np.asarray(upper_hint, dtype=np.float64)
    hint /= np.linalg.norm(hint)

    sz = mesh.nodes[:, 2] - z
    extent = float(np.ptp(mesh.nodes[:, 2]))
    eps = 1e-12 * (extent if extent > 0 else 1.0)
    sz = np.where(np.abs(sz) < eps, eps, sz)

    xs, cps, sides = [], [], []
    for k, tri in enumerate(wall):
        se = sz[tri]
        if np.all(se > 0) or np.all(se < 0):
            continue
        seg = []
        for a, b in ((0, 1), (0, 2), (1, 2)):
            if (se[a] > 0) != (se[b] > 0):
                t = se[a] / (se[a] - se[b])
                seg.append(mesh.nodes[tri[a]] * (1 - t) + mesh.nodes[tri[b]] * t)
        if len(seg) != 2:
            continue
        mid = 0.5 * (seg[0] + seg[1])
        xs.append(mid[0])
        if m_inf > 0.0:
            cps.append(pressure_coefficient(float(q2[k]), m_inf))
        else:
            cps.append(pressure_coefficient_incompressible(float(q2[k])))
        sides.append(float(np.dot(mid, hint)) > 0.0)

    return (np.asarray(xs), np.asarray(cps),
            np.asarray(sides, dtype=bool))


def wall_cp_curve(mesh, phi, z: float, u_inf: float = 1.0,
                  upper_hint=(0.0, 1.0, 0.0), wall_tag: str = "wall",
                  chord: float = 1.0, x_le: float = 0.0,
                  m_inf: float = 0.0) -> Dict[str, np.ndarray]:
    """Sectional wall Cp(x/c) at z = const, split into upper/lower curves.

    Triangle-wise plane cut (see `_wall_section_points`): fully general in z,
    no nodal averaging. Sides split by `upper_hint` (default +y); points
    sorted by x/c. m_inf > 0 selects the exact isentropic Cp (2.5) instead of
    the incompressible one (P3). `chord`/`x_le` are caller-supplied; for a 3D
    swept wing where they vary with span use `section_cp_curve`, which derives
    them from the cut itself.

    Returns:
        dict: x_upper, cp_upper, x_lower, cp_lower (x as x/c from x_le)
    """
    xs, cps, sides = _wall_section_points(
        mesh, phi, z, u_inf, upper_hint, wall_tag, m_inf)
    xs = (xs - x_le) / chord
    iu = np.argsort(xs[sides])
    il = np.argsort(xs[~sides])
    return {
        "x_upper": xs[sides][iu], "cp_upper": cps[sides][iu],
        "x_lower": xs[~sides][il], "cp_lower": cps[~sides][il],
    }


def section_cp_curve(mesh, phi, *, eta: Optional[float] = None,
                     z: Optional[float] = None, b_semi: Optional[float] = None,
                     u_inf: float = 1.0, m_inf: float = 0.0,
                     wall_tag: str = "wall", upper_hint=(0.0, 1.0, 0.0),
                     min_points_per_side: int = 5) -> Dict[str, np.ndarray]:
    """Sectional wall Cp(x/c) at a spanwise station of a 3D wing (roadmap P5).

    Thin wrapper over the general `wall_cp_curve` cut that (a) resolves the
    cut plane from a normalized span fraction `eta` (z = eta * b_semi) or an
    absolute `z`, and (b) auto-derives the local chord and leading-edge x from
    the cut loop's own x-extent -- the right normalizer for plotting the
    mesh's Cp against literature x/c on a swept, tapered wing without wiring in
    the analytic planform. Output keys match `wall_cp_curve` and feed
    `shock_report()` directly.

    Args:
        eta: span fraction in [0, 1) (exactly one of eta/z required).
        z: absolute spanwise coordinate (alternative to eta).
        b_semi: semi-span; required when eta is given.
        min_points_per_side: guard -- raise if either surface has fewer
            crossing points (plane missed the wing or hit the flat tip cap).

    Returns:
        dict: x_upper, cp_upper, x_lower, cp_lower, and diagnostics
            chord, x_le, z, eta.

    Raises:
        ValueError: neither/both of eta/z given, eta without b_semi, or the
            cut is too sparse on a side.
    """
    if (eta is None) == (z is None):
        raise ValueError("pass exactly one of eta or z")
    if eta is not None:
        if b_semi is None:
            raise ValueError("b_semi is required when eta is given")
        z = float(eta) * float(b_semi)

    xs, cps, sides = _wall_section_points(
        mesh, phi, float(z), u_inf, upper_hint, wall_tag, m_inf)
    n_up, n_lo = int(sides.sum()), int((~sides).sum())
    if n_up < min_points_per_side or n_lo < min_points_per_side:
        raise ValueError(
            f"section at z={z:.4f} too sparse (upper={n_up}, lower={n_lo}); "
            "the plane likely missed the wing or hit the flat tip cap")

    x_le = float(xs.min())
    chord = float(xs.max() - xs.min())
    xn = (xs - x_le) / chord
    iu = np.argsort(xn[sides])
    il = np.argsort(xn[~sides])
    return {
        "x_upper": xn[sides][iu], "cp_upper": cps[sides][iu],
        "x_lower": xn[~sides][il], "cp_lower": cps[~sides][il],
        "chord": chord, "x_le": x_le, "z": float(z),
        "eta": (float(z) / float(b_semi)) if b_semi else None,
    }


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
