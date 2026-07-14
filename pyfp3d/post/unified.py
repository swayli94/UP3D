"""
Unified wall post-processing dispatch across the conforming and level-set paths.

`post/surface.py` reads one nodal potential per node (conforming); `post/surface_ls.py`
reads two side potentials via a MultivaluedOperator (level-set / Track B, the D11
per-side mapping). The two share private cores (`surface._cp_from_q2`,
`surface._pressure_force`, `section_cut._wall_plane_crossings /
_section_curve_dict`, `surface_ls._d11_wall_state`); this module is the single
upper-level entry point over both.

Keyword dispatch selects the path (B11):
    - conforming: pass `phi=<nodal potential>`
    - level-set:  pass `mvop=<MultivaluedOperator>, phi_ext=<extended vector>`
Exactly one of the two must be supplied. Because both branches call the same
shared cores, the outputs are `np.array_equal` to the legacy functions
(`wall_force_coefficients`, `wall_cp_levelset`, `section_cp_curve`,
`section_cp_curve_levelset`), which remain available and unchanged.
"""

from typing import Dict, Optional

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
from pyfp3d.post.surface_ls import _d11_wall_state, section_cp_curve_levelset


def _check_dispatch(phi, mvop, phi_ext):
    """Validate the phi= vs (mvop=, phi_ext=) dispatch keywords."""
    if (phi is None) == (mvop is None):
        raise ValueError("pass exactly one of phi= or (mvop=, phi_ext=)")
    if mvop is not None and phi_ext is None:
        raise ValueError("mvop= requires phi_ext=")


def _conforming_wall_state(mesh, phi, wall, u_inf, smooth_passes):
    """Single-phi analogue of `surface_ls._d11_wall_state`: per-triangle q2
    (with the optional G6.1 smoothing) + geometric upper mask (n_y > 0)."""
    grad, area, _ = triangle_tangential_gradients(mesh.nodes, wall, phi)
    n_out = wall_outward_normals(mesh.nodes, mesh.elements, wall)
    if smooth_passes > 0:
        adj = wall_triangle_adjacency(wall)
        grad = smooth_wall_tangential_gradients(
            grad, n_out, area, adj, n_passes=smooth_passes)
    q2 = np.sum(grad * grad, axis=1) / u_inf**2
    upper = n_out[:, 1] > 0.0
    return q2, upper, area, n_out


def wall_cp(mesh, *, phi: Optional[np.ndarray] = None, mvop=None,
            phi_ext: Optional[np.ndarray] = None, m_inf: float = 0.0,
            u_inf: float = 1.0, gamma: float = GAMMA, wall_tag: str = "wall",
            smooth_passes: int = 0) -> Dict[str, np.ndarray]:
    """Per-wall-triangle Cp on either path (unified; the `wall_cp_levelset`-shaped
    dict {x, cp, upper, area, n_out, q2}). See the module docstring for dispatch."""
    _check_dispatch(phi, mvop, phi_ext)
    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    if phi is not None:
        q2, upper, area, n_out = _conforming_wall_state(
            mesh, phi, wall, u_inf, smooth_passes)
    else:
        q2, upper, area, n_out = _d11_wall_state(
            mesh, mvop, phi_ext, wall, u_inf, smooth_passes)
    return {
        "x": mesh.nodes[wall].mean(axis=1)[:, 0],
        "cp": np.asarray(_cp_from_q2(q2, m_inf, gamma), dtype=np.float64),
        "upper": upper,
        "area": area,
        "n_out": n_out,
        "q2": q2,
    }


def wall_forces(mesh, *, phi: Optional[np.ndarray] = None, mvop=None,
                phi_ext: Optional[np.ndarray] = None, alpha_deg: float = 0.0,
                u_inf: float = 1.0, s_ref: float = 1.0, m_inf: float = 0.0,
                gamma: float = GAMMA, wall_tag: str = "wall",
                smooth_passes: int = 0) -> Dict[str, np.ndarray]:
    """Pressure-integrated force coefficients on either path (unified; the
    `wall_force_coefficients`-shaped dict {cl, cd_pressure, cf, cp_tri}).

    The conforming branch is `np.array_equal` to `wall_force_coefficients`;
    the level-set branch's `cl` equals `cl_pressure_3d_levelset` (s_ref =
    planform area) or `cl_pressure_levelset` (s_ref = span extent)."""
    _check_dispatch(phi, mvop, phi_ext)
    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    if phi is not None:
        q2, _, area, n_out = _conforming_wall_state(
            mesh, phi, wall, u_inf, smooth_passes)
    else:
        q2, _, area, n_out = _d11_wall_state(
            mesh, mvop, phi_ext, wall, u_inf, smooth_passes)
    cp_tri = _cp_from_q2(q2, m_inf, gamma)
    cf, cl, cd = _pressure_force(cp_tri, area, n_out, s_ref, alpha_deg)
    return {"cl": cl, "cd_pressure": cd, "cf": cf, "cp_tri": cp_tri}


def section_cp(mesh, *, phi: Optional[np.ndarray] = None, mvop=None,
               phi_ext: Optional[np.ndarray] = None, eta=None, z=None,
               b_semi=None, u_inf: float = 1.0, m_inf: float = 0.0,
               gamma: float = GAMMA, wall_tag: str = "wall",
               upper_hint=(0.0, 1.0, 0.0), min_points_per_side: int = 5,
               smooth_passes: int = 0) -> Dict[str, np.ndarray]:
    """Sectional wall Cp(x/c) at a spanwise station on either path (unified).

    Delegates to `section_cp_curve` (conforming, geometric `upper_hint` side)
    or `section_cp_curve_levelset` (level-set, D11 physical side -- `upper_hint`
    is ignored there). Output keys are identical and feed `shock_report`."""
    _check_dispatch(phi, mvop, phi_ext)
    if phi is not None:
        return section_cp_curve(
            mesh, phi, eta=eta, z=z, b_semi=b_semi, u_inf=u_inf, m_inf=m_inf,
            wall_tag=wall_tag, upper_hint=upper_hint,
            min_points_per_side=min_points_per_side, smooth_passes=smooth_passes)
    return section_cp_curve_levelset(
        mesh, mvop, phi_ext, eta=eta, z=z, b_semi=b_semi, u_inf=u_inf,
        m_inf=m_inf, gamma=gamma, wall_tag=wall_tag,
        min_points_per_side=min_points_per_side, smooth_passes=smooth_passes)
