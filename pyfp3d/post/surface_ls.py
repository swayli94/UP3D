"""
Wall surface post-processing on the LEVEL-SET (Track B) path.

The conforming `post/surface.py` reads one nodal potential per node. On the
multivalued path a TE node carries TWO values (main = its own side, aux = the
other side), so the wall triangles must be told WHICH copy to read -- that is
D11 (design_track_b.md §4): the lower-surface triangles at the trailing edge
must read the TE node's AUX dof. Reading `phi_main` alone on both surfaces
makes the pressure integral junk (measured cl_pressure = -3.35 on NACA0012
alpha=2, vs 0.28 -- the B3 finding).

The side selection is by the wall triangle's OUTWARD normal: n_y > 0 is the
upper surface (read the upper field), n_y < 0 the lower surface. Away from the
TE the two side fields coincide, so the mapping is a no-op everywhere except
exactly where it matters.

Nothing here is imported by the conforming solver or post paths.
"""

from typing import Dict

import numpy as np

from pyfp3d.physics.isentropic import (
    GAMMA,
    pressure_coefficient,
    pressure_coefficient_incompressible,
)
from pyfp3d.post.surface import triangle_tangential_gradients, wall_outward_normals


def wall_cp_levelset(mesh, mvop, phi_ext, m_inf: float = 0.0,
                     u_inf: float = 1.0, gamma: float = GAMMA) -> Dict[str, np.ndarray]:
    """Per-wall-triangle Cp on the level-set path, with the D11 per-side mapping.

    Args:
        mesh: the (uncut) mesh -- the level-set path never duplicates nodes
        mvop: MultivaluedOperator (supplies the two side potentials)
        phi_ext: the extended solution vector
        m_inf: free-stream Mach; 0.0 selects the incompressible (Bernoulli) Cp
        u_inf, gamma: free-stream speed / ratio of specific heats

    Returns:
        dict: x (triangle-centroid x/c), cp, upper (bool mask, n_y > 0),
        area, n_out (outward unit normals), q2 (surface speed^2 / u_inf^2)
    """
    wall = mesh.boundary_faces["wall"]
    phi_up, phi_lo = mvop.side_potentials(phi_ext)
    n_out = wall_outward_normals(mesh.nodes, mesh.elements, wall)
    g_up, area, _ = triangle_tangential_gradients(mesh.nodes, wall, phi_up)
    g_lo, _, _ = triangle_tangential_gradients(mesh.nodes, wall, phi_lo)

    upper = n_out[:, 1] > 0.0
    grad = np.where(upper[:, None], g_up, g_lo)          # D11
    q2 = np.sum(grad * grad, axis=1) / u_inf**2

    if m_inf > 0.0:
        # scalar njit -- evaluate per triangle (the wall set is O(10^2-10^3))
        cp = np.array([pressure_coefficient(v, m_inf, gamma) for v in q2])
    else:
        cp = np.array([pressure_coefficient_incompressible(v) for v in q2])

    return {
        "x": mesh.nodes[wall].mean(axis=1)[:, 0],
        "cp": np.asarray(cp, dtype=np.float64),
        "upper": upper,
        "area": area,
        "n_out": n_out,
        "q2": q2,
    }


def cl_pressure_levelset(mesh, cp, area, n_out, alpha_deg: float) -> float:
    """Sectional pressure lift coefficient from the level-set wall Cp.

    The quasi-2D meshes are one span-cell thick, so the surface integral is
    normalised by the span extent to give a SECTIONAL cl (chord = 1).
    """
    a = np.radians(alpha_deg)
    lift = np.array([-np.sin(a), np.cos(a), 0.0])
    dz = float(np.ptp(mesh.nodes[:, 2]))
    return float((-(cp * area) @ n_out / dz) @ lift)


def surface_curve_levelset(cp_data, side: str = "upper"):
    """(x, cp) of one surface, sorted by x -- the input `post/shock.py`
    (`shock_metrics`) expects for shock detection."""
    m = cp_data["upper"] if side == "upper" else ~cp_data["upper"]
    x, cp = cp_data["x"][m], cp_data["cp"][m]
    order = np.argsort(x)
    return x[order], cp[order]


def cl_pressure_3d_levelset(mesh, cp, area, n_out, alpha_deg: float,
                            s_ref: float) -> float:
    """3D half-wing pressure lift coefficient from the level-set wall Cp (B7).

    The 3D analogue of `cl_pressure_levelset`: the same surface integral, but
    normalised by the half-wing PLANFORM AREA (`post/surface.planform_area`)
    instead of the quasi-2D span extent, so it pairs with `post/surface.cl_kj_3d`
    for the V6 circulation/pressure consistency cross-check (design.md Sec 9;
    roadmap G5.2 on the conforming path).

    Args:
        mesh, cp, area, n_out: the mesh and `wall_cp_levelset` outputs
        alpha_deg: incidence (lift direction in the chord-lift plane)
        s_ref: half-wing planform area
    """
    a = np.radians(alpha_deg)
    lift = np.array([-np.sin(a), np.cos(a), 0.0])
    return float((-(cp * area) @ n_out / s_ref) @ lift)


def section_cp_curve_levelset(mesh, mvop, phi_ext, *, eta=None, z=None,
                              b_semi=None, u_inf: float = 1.0,
                              m_inf: float = 0.0, gamma: float = GAMMA,
                              wall_tag: str = "wall",
                              min_points_per_side: int = 5) -> Dict[str, np.ndarray]:
    """Sectional wall Cp(x/c) at a spanwise station on the LEVEL-SET path (B7).

    The level-set analogue of `post/section_cut.section_cp_curve`: the same
    triangle-wise plane cut of the wall at z = const, with the same auto-derived
    local chord / x_le (so a swept, tapered planform normalises correctly), but
    the per-triangle velocity is taken from the D11 per-side potential --
    `mvop.side_potentials`, selected by the wall triangle's outward-normal lift
    component (n_y > 0 = upper). Away from the TE the two side fields coincide,
    so this differs from the `main_potential`-based curve only in the lower
    surface's TE triangles -- which is exactly the pair of triangles that make
    the lower-surface curve junk if read from `main` alone (the B3 finding
    recorded in this module's header).

    Note the UPPER surface needs no correction (its TE vertices carry the main =
    own-side value), so gate metrics taken on the upper surface -- the shock
    positions -- agree with `section_cp_curve(mesh, mvop.main_potential(...))`
    to the last bit. This function exists for the lower surface and for
    both-side plots.

    Args:
        eta / z: exactly one -- span fraction (needs b_semi) or absolute z
        b_semi: semi-span, required with eta
        min_points_per_side: guard against a plane that missed the wing or hit
            the flat tip cap

    Returns:
        dict: x_upper, cp_upper, x_lower, cp_lower (x as x/c from the local
        x_le), plus chord, x_le, z, eta -- key-compatible with
        `section_cp_curve`, so it feeds `post/shock.shock_report` directly.
    """
    from pyfp3d.physics.isentropic import (
        pressure_coefficient,
        pressure_coefficient_incompressible,
    )

    if (eta is None) == (z is None):
        raise ValueError("pass exactly one of eta or z")
    if eta is not None:
        if b_semi is None:
            raise ValueError("b_semi is required when eta is given")
        z = float(eta) * float(b_semi)
    z = float(z)

    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    phi_up, phi_lo = mvop.side_potentials(phi_ext)
    n_out = wall_outward_normals(mesh.nodes, mesh.elements, wall)
    g_up, _, _ = triangle_tangential_gradients(mesh.nodes, wall, phi_up)
    g_lo, _, _ = triangle_tangential_gradients(mesh.nodes, wall, phi_lo)
    upper = n_out[:, 1] > 0.0                                    # D11
    grad = np.where(upper[:, None], g_up, g_lo)
    q2 = np.sum(grad * grad, axis=1) / u_inf**2

    # Plane cut (mirrors section_cut._wall_section_points): every crossed
    # triangle contributes one point at its intersection-segment midpoint,
    # carrying that triangle's own constant Cp. A vertex exactly on the plane
    # is nudged to the + side so a crossed triangle always yields two points.
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
            cps.append(pressure_coefficient(float(q2[k]), m_inf, gamma))
        else:
            cps.append(pressure_coefficient_incompressible(float(q2[k])))
        sides.append(bool(upper[k]))          # D11 side, not a geometric hint

    xs = np.asarray(xs, dtype=np.float64)
    cps = np.asarray(cps, dtype=np.float64)
    sides = np.asarray(sides, dtype=bool)
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
        "chord": chord, "x_le": x_le, "z": z,
        "eta": (z / float(b_semi)) if b_semi else None,
    }
