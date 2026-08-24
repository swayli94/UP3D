"""M0.3 球扰流的共享装配 —— 被 A 类（位相同/对称性）与 C 类（Prandtl–Glauert）两半共用。

★ C 半是本项目**唯一已存在的可压缩渐近检验**（PG），此前被我的归类漏掉过 ——
它不是精确解，只能判阶次/趋势，**不能当真值**。
"""
import csv
import numpy as np
import pytest
import matplotlib
import matplotlib.pyplot as plt
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import pressure_coefficient
from pyfp3d.post.surface import wall_tangential_gradient_quadratic
from pyfp3d.solve.picard import (
    solve_laplace,
    solve_laplace_lifting,
    solve_subsonic,
    solve_subsonic_lifting,
)


def run_sphere_pair(mesh_path):
    """Incompressible + compressible solves on one sphere mesh, wall Cp
    from the quadratic tangential recovery for both."""
    mesh = read_mesh(mesh_path)
    nodes, elements = mesh.nodes, mesh.elements
    wall = mesh.boundary_faces["wall"]
    wn = np.unique(wall)
    ff = np.unique(mesh.boundary_faces["farfield"])
    phi_ff = nodes[ff, 0]  # freestream Dirichlet, same for both solves

    r_inc = solve_laplace(nodes, elements, ff, phi_ff, rtol=1e-11, maxiter=3000)
    g = wall_tangential_gradient_quadratic(nodes, wall, r_inc["phi"])
    cp_inc = 1.0 - np.sum(g[wn] ** 2, axis=1)

    r_c = solve_subsonic(
        nodes, elements, ff, phi_ff, m_inf=M_INF_SPHERE,
        phi_init=nodes[:, 0].copy(), rtol=1e-12,
    )
    g = wall_tangential_gradient_quadratic(nodes, wall, r_c["phi"])
    q2 = np.sum(g[wn] ** 2, axis=1)
    cp_c = np.array([pressure_coefficient(q, M_INF_SPHERE) for q in q2])

    r = np.linalg.norm(nodes[wn], axis=1)
    cos_theta = nodes[wn, 0] / r
    return {
        "cos_theta": cos_theta, "cp_inc": cp_inc, "cp_c": cp_c,
        "result_c": r_c,
    }

M_INF_SPHERE = 0.3
@pytest.fixture(scope="module")
def sphere_medium(request):
    from tests.conftest import REPO_ROOT
    return run_sphere_pair(
        REPO_ROOT / "cases" / "meshes" / "sphere_shell" / "medium.msh"
    )
