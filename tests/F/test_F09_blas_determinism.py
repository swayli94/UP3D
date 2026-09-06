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

## ★★★ 载体必须是**实测敏感**的,这一条是踩出来的

本门第一版用 coarse M0.50 作载体,**G-TEETH 当场失败**:停用固定它照样绿。
再测 coarse M0.80/α1.25 与 coarse M0.75/α0 —— **也都逐位相同**。
⇒ 我在第一版 docstring 里写的「位确定是全局性质,不是逐算例性质」是**推理,
而且被测量否掉了**:OPENBLAS 敏感性**不是普遍性质**,它只出现在那几条
**病态腿**上(medium 的 M0.80 / M0.803,迭代路径对 1e-15 扰动敏感)。
⇒ 载体只能用 **medium M0.80/α1.25**,并**故意用 `n_newton_max = 80`**
(它在这个预算下不收敛 —— 本门判的是**位相同**,与收敛无关;实测
OPENBLAS 1 vs 8 给 0.34009652307781807 vs 0.34011432026323)。

## 判据为什么能有牙

`pyfp3d/__init__.py` **强制覆盖**(不是 `setdefault`),所以探针即便外部设
`OPENBLAS_NUM_THREADS=8/16`,库也会压成 1 ⇒ 两条腿必须逐位相同。
**G-TEETH 用逃生舱 `PYFP3D_ALLOW_BLAS_THREADS=1`**:放行外部值 ⇒ 两条腿走
不同的 OPENBLAS ⇒ 本门必红。**于是验证「固定会不会失效」不需要改库文件。**
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
p = os.path.join(%r, "cases", "meshes", "naca0012_2.5d", "medium.msh")
mc, wc = cut_wake(read_mesh(p))
r = solve_newton_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25, upwind_c=1.5,
                         m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
                         precond="direct", direct_refactor_every=4,
                         n_picard_seed=5, n_newton_max=80)
f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                            np.asarray(r["phi"]), alpha_deg=1.25, u_inf=1.0,
                            s_ref=float(np.ptp(mc.nodes[:, 2])), m_inf=0.80)
print("OUT " + json.dumps({
    "cl": repr(float(f["cl"])),
    "resid": repr(float(np.asarray(r["residual_history"], float)[-1])),
    "blas_seen": os.environ.get("OPENBLAS_NUM_THREADS"),
    "blas_ext": os.environ.get("PYFP3D_PROBE_EXT_BLAS"),
}))
'''


def _run(blas):
    """真跑一次求解,外部把 BLAS 设成 `blas`。"""
    env = dict(os.environ)
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[v] = str(blas)
    env["PYFP3D_PROBE_EXT_BLAS"] = str(blas)   # 记录外部值,供前提断言比对
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
    p = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh"
    if not p.exists():                       # 与 W0.1 同一条约定
        pytest.skip("naca0012_2.5d/medium.msh not generated")
    return _run(8), _run(16)


class TestBlasReductionOrderDoesNotMoveTheAnswer:
    def test_the_probe_really_saw_different_external_settings(self, pair):
        """前提断言:两条腿的外部 BLAS 设置**确实不同**。

        ★ 少了它,两条腿若因为任何原因跑在同一设置下,下面那条就**按构造恒真**。
        ★ 注意这里读的是**子进程内**的值 —— `pyfp3d` 用 `setdefault`,所以
          外部显式设的 8 / 16 会被保留,正是本门需要的对照。
        """
        a, b = pair
        assert a["blas_ext"] != b["blas_ext"], (
            f"两条腿的**外部** OPENBLAS 设置相同({a['blas_ext']}) —— 对照没建立起来")
        assert a["blas_seen"] == b["blas_seen"] == "1", (
            f"库没有把 OPENBLAS 压成 1(实际 {a['blas_seen']} / {b['blas_seen']})"
            " —— 强制覆盖失效,多半是 numpy 在 pyfp3d 之前被导入")

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
