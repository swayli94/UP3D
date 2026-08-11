"""
Wall post-processing entry points (conforming).

★ COLLAPSED 2026-08-10 (phase 3 task 1, ruling D5). This module was the single
upper-level entry point over TWO wake paths -- `post/surface.py` reading one nodal
potential (conforming) and `post/surface_ls.py` reading two side potentials through a
MultivaluedOperator (level-set) -- dispatching on `phi=` versus `(mvop=, phi_ext=)`.
The level-set route is abandoned, so the dispatch collapses onto its conforming half:
one path, no dispatch, and the `mvop=`/`phi_ext=` keywords are GONE rather than kept as
inert arguments that would advertise a route that no longer exists.

What this means for callers: `phi=` stays keyword-only and every other keyword is
unchanged, so conforming call sites are untouched -- which is the point, because this
round subtracts only. Outputs remain `np.array_equal` to the legacy functions
(`wall_force_coefficients`, `section_cp_curve`), which are still available and
unchanged; both branches always called the same shared cores (`surface._cp_from_q2`,
`surface._pressure_force`, `section_cut._wall_plane_crossings / _section_curve_dict`).
"""

from typing import Dict

import numpy as np

from pyfp3d.physics.isentropic import GAMMA
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import (
    _cp_from_q2,
    _pressure_force,
    smooth_wall_tangential_gradients,
    triangle_tangential_gradients,
    wall_outward_normals,
    wall_triangle_adjacency,
)


def _conforming_wall_state(mesh, phi, wall, u_inf, smooth_passes):
    """Per-triangle q2 (with the optional G6.1 smoothing) + geometric upper mask
    (n_y > 0)."""
    grad, area, _ = triangle_tangential_gradients(mesh.nodes, wall, phi)
    n_out = wall_outward_normals(mesh.nodes, mesh.elements, wall)
    if smooth_passes > 0:
        adj = wall_triangle_adjacency(wall)
        grad = smooth_wall_tangential_gradients(
            grad, n_out, area, adj, n_passes=smooth_passes)
    q2 = np.sum(grad * grad, axis=1) / u_inf**2
    upper = n_out[:, 1] > 0.0
    return q2, upper, area, n_out


def wall_cp(mesh, *, phi: np.ndarray, m_inf: float = 0.0, u_inf: float = 1.0,
            gamma: float = GAMMA, wall_tag: str = "wall",
            smooth_passes: int = 0) -> Dict[str, np.ndarray]:
    """Per-wall-triangle Cp: {x, cp, upper, area, n_out, q2}."""
    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    q2, upper, area, n_out = _conforming_wall_state(
        mesh, phi, wall, u_inf, smooth_passes)
    return {
        "x": mesh.nodes[wall].mean(axis=1)[:, 0],
        "cp": np.asarray(_cp_from_q2(q2, m_inf, gamma), dtype=np.float64),
        "upper": upper,
        "area": area,
        "n_out": n_out,
        "q2": q2,
    }


def wall_forces(mesh, *, phi: np.ndarray, alpha_deg: float = 0.0,
                u_inf: float = 1.0, s_ref: float = 1.0, m_inf: float = 0.0,
                gamma: float = GAMMA, wall_tag: str = "wall",
                smooth_passes: int = 0) -> Dict[str, np.ndarray]:
    """Pressure-integrated force coefficients: {cl, cd_pressure, cf, cp_tri}.
    `np.array_equal` to `wall_force_coefficients`."""
    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    q2, _, area, n_out = _conforming_wall_state(
        mesh, phi, wall, u_inf, smooth_passes)
    cp_tri = _cp_from_q2(q2, m_inf, gamma)
    cf, cl, cd = _pressure_force(cp_tri, area, n_out, s_ref, alpha_deg)
    return {"cl": cl, "cd_pressure": cd, "cf": cf, "cp_tri": cp_tri}


def section_cp(mesh, *, phi: np.ndarray, eta=None, z=None, b_semi=None,
               u_inf: float = 1.0, m_inf: float = 0.0, gamma: float = GAMMA,
               wall_tag: str = "wall", upper_hint=(0.0, 1.0, 0.0),
               min_points_per_side: int = 5,
               smooth_passes: int = 0) -> Dict[str, np.ndarray]:
    """Sectional wall Cp(x/c) at a spanwise station; delegates to
    `section_cp_curve`. Output keys feed `shock_report`."""
    return section_cp_curve(
        mesh, phi, eta=eta, z=z, b_semi=b_semi, u_inf=u_inf, m_inf=m_inf,
        gamma=gamma, wall_tag=wall_tag, upper_hint=upper_hint,
        min_points_per_side=min_points_per_side, smooth_passes=smooth_passes)
