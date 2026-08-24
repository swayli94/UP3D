"""C04 -- M1a 三级网格收敛性质（C 类，只判收敛性质那一半）。

★★ 使用者裁决 2026-08-24（边界 ①）：整机的网格收敛性质**没有解析解**，归 C 并**只判阶次/收缩**那一半
-- 单调 + 收缩比 < 0.7 + 末步 < 5 %，**不比任何绝对值**。
★ 共享装配在 tests/_m1a_case.py（_solve 无缓存，与拆分前行为相同）。
"""
import numpy as np
import pytest
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
from tests._m1a_case import ALPHA, CL_RTOL, D2_REL_MAX, LADDER, LOCK, M_INF, RATIO_MAX, _solve, mesh_dir


def test_m1a_three_level_convergence(mesh_dir):
    """★ The re-spec'd M1a criterion (2026-08-05): on (xcoarse, coarse, medium) the
    sequence must be MONOTONE, CONTRACTING, and settled to under 5 % on the last step.

    Each of the three discriminates out of envelope -- at M0.80 the same sequence gives
    d2 = -0.0666, ratio -0.6003 and 16.33 % -- while "all three converged" does NOT,
    since the cold-start seed fallback makes M0.80 converge too. So monotonicity is
    what carries the statement, and it is asserted first.

    What this criterion no longer claims, and why, is in the LOCK comment above and in
    phases/p2/docs/dev_phase_two/20260805-2330-m1a-respec.md sec 5.
    """
    cl = {}
    for level in LADDER:
        r, c = _solve(mesh_dir, level)
        assert r["converged"] and r["n_limited"] == 0 and r["n_floored"] == 0, (
            f"{level}: the sequence needs three genuine solutions")
        cl[level] = c
    d1 = cl["coarse"] - cl["xcoarse"]
    d2 = cl["medium"] - cl["coarse"]
    assert d1 > 0.0 and d2 > 0.0, (
        f"NOT monotone: d1 = {d1:+.6f}, d2 = {d2:+.6f}. This is the criterion that "
        f"separates in-envelope from out (M0.80 gives d2 = -0.0666)")
    ratio = d2 / d1
    assert 0.0 < ratio < RATIO_MAX, (
        f"convergence ratio {ratio:.4f} outside (0, {RATIO_MAX}) "
        f"(0.5065 measured in envelope; -0.6003 at M0.80)")
    rel = abs(d2 / cl["medium"])
    assert rel < D2_REL_MAX, (
        f"last-step change {100 * rel:.3f} % >= {100 * D2_REL_MAX:.0f} % "
        f"(4.22 % measured in envelope; 16.33 % at M0.80)")
