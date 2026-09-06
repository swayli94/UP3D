r"""F09 — **BLAS 归约顺序不得影响答案**(W2/H30,2026-09-06)。

## 为什么有这道门

审计实测(`bench/studies/gate_audit_20260905/results/thread_attribution_D05_M080.csv`,
D05 的 M0.80/α1.25 medium,9 条腿单变量):

| 分组 | 每档内不同的 cl 值 |
|---|---|
| 按 **BLAS** 线程数 | `{1:1, 8:1, 16:1}` —— 每档**唯一** |
| 按 numba 线程数 | `{1:3, 8:2, 16:2}` —— 每档**多值** |

⇒ **结果是 BLAS 线程数的函数,numba 不是自变量**(1/8/16 逐位相同)。
同配置重复跑逐位相同 ⇒ **不是竞态**,是多线程 BLAS 的**归约顺序**。

这几道门不是效率门,**结果依赖线程数本身就是错误**(使用者裁决 2026-09-06)。
`pyfp3d/__init__.py` 因此把 BLAS 固定成单线程。

## ★★★ 判据是**行为**的,不是「环境变量被设了」

后者按构造恒真,是本仓库 F06 那族「字面量对字面量、改了参考文件什么都不会红」
的形状。所以本门**外部设 `OPENBLAS_NUM_THREADS=8` 与 `=16` 各跑一次真实求解,
要求结果逐位相同** —— 子进程,因为线程数必须在 numpy 导入前生效。

★ 载体取**便宜且对 BLAS 敏感**的一条:2.5-D NACA0012 coarse 的亚声速 Newton。
用 M0.80 那条(真正暴露问题的)要 200+ 步,太贵;而只要固定生效,
**任何**算例都该位确定 —— 位确定是全局性质,不是逐算例性质。
"""
import json
import os
import subprocess
import sys

import pytest

from tests.conftest import REPO_ROOT

_PROBE = r'''
import json, os, sys
sys.path.insert(0, %r)
import pyfp3d                      # <- 必须最先,它负责固定 BLAS
import numpy as np
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting
p = os.path.join(%r, "cases", "meshes", "naca0012_2.5d", "coarse.msh")
mc, wc = cut_wake(read_mesh(p))
r = solve_newton_lifting(mc, wc, m_inf=0.50, alpha_deg=2.0, upwind_c=1.5,
                         m_crit=0.95, precond="direct", n_picard_seed=5,
                         n_newton_max=60, tol_residual=1e-10)
f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                            np.asarray(r["phi"]), alpha_deg=2.0, u_inf=1.0,
                            s_ref=float(np.ptp(mc.nodes[:, 2])), m_inf=0.50)
print("OUT " + json.dumps({
    "cl": repr(float(f["cl"])),
    "resid": repr(float(np.asarray(r["residual_history"], float)[-1])),
    "blas_seen": os.environ.get("OPENBLAS_NUM_THREADS"),
}))
'''


def _run(blas):
    """真跑一次求解,外部把 BLAS 设成 `blas`。"""
    env = dict(os.environ)
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[v] = str(blas)
    env["NUMBA_NUM_THREADS"] = "4"
    env["PYTHONNOUSERSITE"] = "1"
    out = subprocess.run(
        [sys.executable, "-c", _PROBE % (str(REPO_ROOT), str(REPO_ROOT))],
        capture_output=True, text=True, env=env, timeout=1800)
    line = [l for l in out.stdout.splitlines() if l.startswith("OUT ")]
    assert line, f"探针没有产出(blas={blas}):\n{out.stdout[-2000:]}\n{out.stderr[-2000:]}"
    return json.loads(line[-1][4:])


@pytest.fixture(scope="module")
def pair():
    p = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh"
    if not p.exists():                       # 与 W0.1 同一条约定
        pytest.skip("naca0012_2.5d/coarse.msh not generated")
    return _run(8), _run(16)


class TestBlasReductionOrderDoesNotMoveTheAnswer:
    def test_the_probe_really_saw_different_external_settings(self, pair):
        """前提断言:两条腿的外部 BLAS 设置**确实不同**。

        ★ 少了它,两条腿若因为任何原因跑在同一设置下,下面那条就**按构造恒真**。
        ★ 注意这里读的是**子进程内**的值 —— `pyfp3d` 用 `setdefault`,所以
          外部显式设的 8 / 16 会被保留,正是本门需要的对照。
        """
        a, b = pair
        assert a["blas_seen"] != b["blas_seen"], (
            f"两条腿看到的 BLAS 设置相同({a['blas_seen']}) —— 对照没建立起来")

    def test_the_answer_is_bit_identical_across_blas_settings(self, pair):
        """★★★ 本门的判据:**逐位相同**。

        用 `repr(float)` 比字符串,而不是数值容差 —— 这里要的就是**位相同**,
        任何容差都会把「归约顺序改变了答案」这件事放过去。
        ★ G-DOMAIN:相反结果 = 两者不同 ⇒ **固定没生效**(多半是 import 顺序:
          numpy 在 `pyfp3d` 之前被导入),处置是改用**运行时**限制线程池,
          **不是放宽本判据**。
        """
        a, b = pair
        assert a["cl"] == b["cl"], (
            f"BLAS 归约顺序改变了答案:cl {a['cl']} (blas={a['blas_seen']}) "
            f"vs {b['cl']} (blas={b['blas_seen']})\n"
            "★ 这不是容差问题 —— 固定 BLAS 的目的就是让它位确定。"
            "多半是 numpy 在 pyfp3d 之前被导入,导致 setdefault 太晚。")
        assert a["resid"] == b["resid"], (
            f"残差不逐位相同:{a['resid']} vs {b['resid']}")
