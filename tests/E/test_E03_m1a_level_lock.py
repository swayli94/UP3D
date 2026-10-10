"""E03 -- M1a 逐级已提交 cl 的回归锁（E 类）。

★ 比的是 LOCK[level]["cl"]，即**代码自己提交过的值** ⇒ E 类。
★★ NACA0012 M0.80 的 CFL3D Euler 参照（2D-2）到位后**升格 D 类**。
"""
import numpy as np
import pytest
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
from tests._m1a_case import ALPHA, CL_RTOL, D2_REL_MAX, LADDER, LOCK, M_INF, RATIO_MAX, _solve, mesh_dir


@pytest.mark.parametrize("level", LADDER)
def test_m1a_level_lock(mesh_dir, level):
    """Each in-envelope level reproduces its committed cl, converges to a
    genuine solution, and carries no clamped cells."""
    r, cl = _solve(mesh_dir, level)
    assert r["converged"], (
        f"{level} did not converge: |R| = {r['residual_history'][-1]:.2e}")
    assert r["n_limited"] == 0 and r["n_floored"] == 0, (
        f"{level} carries clamps ({r['n_limited']} limited / "
        f"{r['n_floored']} floored) -- a clamped state is not a solution")
    ref = LOCK[level]["cl"]
    assert abs(cl - ref) <= CL_RTOL * abs(ref), (
        f"{level} cl {cl:.6f} vs committed {ref:.6f} "
        f"({100 * (cl - ref) / ref:+.3f} %)")
    assert float(np.sqrt(r["mach2_max"])) < 1.25, (
        "M_max moved out of the in-envelope range this lock is about")
