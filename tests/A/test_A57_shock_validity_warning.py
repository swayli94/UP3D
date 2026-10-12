"""A57 — 风险 3 的警告触发锁(2026-10-12,第 50 轮 T2)。

design.md §12 风险 3:强激波(局部 M > 1.3)时等熵跳变低估激波强度 ——
缓解是 `newton.py` 的 `_FP_SHOCK_VALIDITY_MACH = 1.35` 运行时警告。
**此前没有任何测试断言它会触发**(W3.2 索引实测:风险 3 = 「有运行时缓解、
无门」)。本文件补上这条机制测试。

★ 覆盖事实(登记,不扩权):警告**只在坡道驱动器** `solve_newton_transonic` 里,
单 Mach 的 `solve_newton_lifting` 不警告。
"""
import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import solve_newton_transonic
from tests.conftest import REPO_ROOT

NACA = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"


@pytest.fixture(scope="module")
def coarse():
    p = NACA / "coarse.msh"
    if not p.exists():
        pytest.skip("naca0012_2.5d/coarse.msh not generated")
    return cut_wake(read_mesh(str(p)))


def test_the_warning_fires_beyond_the_validity_mach(coarse, capsys):
    """风险 3:M 爬过 1.35 的坡道必须把「等熵跳变不适用」打到 stdout。
    红了 = 警告机制死了或被静默改阈值 —— 设计文档承诺的缓解消失了。"""
    mc, wc = coarse
    solve_newton_transonic(mc, wc, m_inf=0.84, alpha_deg=2.0, m_start=0.70,
                           dm=0.07, newton_kw=dict(n_newton_max=80))
    out = capsys.readouterr().out
    assert "validity WARNING" in out, (
        "M 爬到 0.84(局部 M 实测 ~1.6)而 validity 警告没有出现 —— "
        "风险 3 的运行时缓解消失了")
    assert "1.35" in out


def test_the_warning_does_not_fire_subcritical(coarse, capsys):
    """G-DOMAIN 双向:低速坡道**不许**出现警告(警告乱报 = 也没有用)。"""
    mc, wc = coarse
    solve_newton_transonic(mc, wc, m_inf=0.50, alpha_deg=2.0, m_start=0.30,
                           dm=0.10, newton_kw=dict(n_newton_max=80))
    out = capsys.readouterr().out
    assert "validity WARNING" not in out, (
        "M ≤ 0.50 的纯亚声速坡道打出了 validity 警告 —— 阈值/触发条件坏了")
