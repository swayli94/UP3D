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
