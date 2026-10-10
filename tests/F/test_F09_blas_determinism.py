r"""F09 — **BLAS 归约顺序不得影响答案**(W2/H30 2026-09-06;**重写 2026-10-10**)。

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

## ★★★ 2026-10-10 重写:旧版只测了「固定**会生效**的那种 import 顺序」

两份独立进展审计(20261009-2200 / 20261010-1038)实测:OpenBLAS **只在 numpy
首次加载时读一次**环境变量,而测试模块普遍先 `import numpy` 再 `import pyfp3d`
⇒ 套件进程里 env 固定**写给空气**(库实际跑外部值/默认 24)。
旧版 F09 的子进程第 55 行写着 `import pyfp3d  # <- 必须最先` —— 它只验收了
pyfp3d-先导入的顺序,而前提断言读的 `os.environ` 正是固定代码自己刚写进去的值
(**按构造恒真**)。于是固定失效了一个多月,门一直绿。

⇒ 修复改两层,门也改两层:

1. **库**:`pyfp3d/__init__.py` 追加**运行时固定**(numpy 已在场时,从
   `/proc/self/maps` 找已加载的 libopenblas,调 `openblas_set_num_threads(1)`)。
2. **门**:
   - 探针跑 **两种 import 顺序**(`pyfp3d`-first 与 `numpy`-first)——
     后者就是 pytest 的真实顺序,旧版从未测过;
   - **前提断言改读库的实际线程数**(ctypes 读 `openblas_get_num_threads()`,
     不读 `os.environ` —— 「环境变量不是库的状态」)。

## ★★★ 判据是**行为**的,不是「环境变量被设了」

后者按构造恒真,是本仓库 F06 那族「字面量对字面量、改了参考文件什么都不会红」
的形状。所以本门**外部设 `OPENBLAS_NUM_THREADS=8` 与 `=16` 各跑一次真实求解,
要求结果逐位相同** —— 子进程,因为线程数必须在 numpy 导入前生效。

## ★★★ 载体必须是**实测敏感**的,这一条是踩出来的

本门第一版用 coarse M0.50 作载体,**G-TEETH 当场失败**:停用固定它照样绿。
再测 coarse M0.80/α1.25 与 coarse M0.75/α0 —— **也都逐位相同**。
⇒ OPENBLAS 敏感性**不是普遍性质**,它只出现在那几条**病态腿**上(medium 的
M0.80 / M0.803,迭代路径对 1e-15 扰动敏感)。
⇒ 载体只能用 **medium M0.80/α1.25**,并**故意用 `n_newton_max = 80`**
(它在这个预算下不收敛 —— 本门判的是**位相同**,与收敛无关;实测
OPENBLAS 1 vs 8 给 0.34009652307781807 vs 0.34011432026323)。

## 判据为什么能有牙

`pyfp3d/__init__.py` **强制覆盖** env 并在 numpy 已在场时**运行时钉库**,
所以探针即便外部设 `OPENBLAS_NUM_THREADS=8/16`,库也被压成 1 ⇒ 各腿必须
逐位相同。**G-TEETH 用逃生舱 `PYFP3D_ALLOW_BLAS_THREADS=1`**:放行外部值 ⇒
两条腿走不同的 OPENBLAS ⇒ 本门必红(前提断言与判据都红)。
**于是验证「固定会不会失效」不需要改库文件。**
"""
import json
import os
import subprocess
import sys

import pytest

from tests.conftest import REPO_ROOT

#: ★ 探针:ORDER ∈ {pyfp3d_first, numpy_first}。numpy-first 就是 pytest 的
#: 真实顺序(2026-10-10 审计:D05/D06/D08/A02/A03 全部先 import numpy)。
#: blas_lib 读**库的实际线程数**(/proc/self/maps + ctypes),不读 os.environ ——
#: 读 env 的那条前提断言按构造恒真(固定代码自己写进去的值),已随 2026-10-10
#: 重写废除。
_PROBE = r'''
import json, os, sys
sys.path.insert(0, %r)
ORDER = %r
if ORDER == "numpy_first":
    import numpy as np          # <- pytest 的真实顺序:库在 pyfp3d 之前加载
import pyfp3d
import numpy as np
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting


def blas_lib_threads():
    """已加载的每个 libopenblas 的实际线程数(库的状态,不是环境变量)。"""
    import ctypes
    paths = set()
    with open("/proc/self/maps") as fh:
        for line in fh:
            parts = line.split()
            if parts and "libopenblas" in parts[-1]:
                paths.add(parts[-1])
    return {os.path.basename(p): ctypes.CDLL(p).openblas_get_num_threads()
            for p in sorted(paths)}


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
    "order": ORDER,
    "cl": repr(float(f["cl"])),
    "resid": repr(float(np.asarray(r["residual_history"], float)[-1])),
    "blas_lib": blas_lib_threads(),
    "blas_ext": os.environ.get("PYFP3D_PROBE_EXT_BLAS"),
}))
'''

_ORDERS = ("pyfp3d_first", "numpy_first")
_BLAS_SETTINGS = (8, 16)


def _run(order, blas):
    """真跑一次求解,import 顺序为 `order`,外部把 BLAS 设成 `blas`。"""
    env = dict(os.environ)
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[v] = str(blas)
    env["PYFP3D_PROBE_EXT_BLAS"] = str(blas)   # 记录外部值,供前提断言比对
    env["NUMBA_NUM_THREADS"] = "4"
    env["PYTHONNOUSERSITE"] = "1"
    out = subprocess.run(
        [sys.executable, "-c", _PROBE % (str(REPO_ROOT), order, str(REPO_ROOT))],
        capture_output=True, text=True, env=env, timeout=1800)
    line = [l for l in out.stdout.splitlines() if l.startswith("OUT ")]
    assert line, (f"探针没有产出(order={order}, blas={blas}):\n"
                  f"{out.stdout[-2000:]}\n{out.stderr[-2000:]}")
    return json.loads(line[-1][4:])


@pytest.fixture(scope="module")
def legs():
    p = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh"
    if not p.exists():                       # 与 W0.1 同一条约定
        pytest.skip("naca0012_2.5d/medium.msh not generated")
    return {(o, b): _run(o, b) for o in _ORDERS for b in _BLAS_SETTINGS}


class TestBlasReductionOrderDoesNotMoveTheAnswer:
    def test_the_probe_really_saw_different_external_settings(self, legs):
        """前提断言:每种顺序下,两条腿的外部 BLAS 设置**确实不同**。

        ★ 少了它,两条腿若因为任何原因跑在同一设置下,下面两条就**按构造恒真**。
        ★ `blas_ext` 是 `PYFP3D_PROBE_EXT_BLAS` —— 一个 pyfp3d **不碰**的变量,
          所以它能存活与固定语义无关(旧版注释把这事记成 setdefault 的功劳,是错的:
          固定是强制覆盖,见 `pyfp3d/__init__.py`)。
        """
        for o in _ORDERS:
            a, b = legs[o, 8], legs[o, 16]
            assert a["blas_ext"] != b["blas_ext"], (
                f"order={o}:两条腿的**外部** OPENBLAS 设置相同({a['blas_ext']})"
                " —— 对照没建立起来")

    def test_the_library_is_actually_pinned(self, legs):
        """★★★ 新前提断言(2026-10-10):**库的实际线程数**被压成 1。

        读的是 `openblas_get_num_threads()`(库的状态),**不是** `os.environ`
        (固定代码自己写进去的值)。少这一条,旧版失效方式(numpy 先导入 ⇒
        env 写了而库没动)在本门里**看不见**。
        """
        for (o, b), leg in legs.items():
            assert leg["blas_lib"], f"order={o} blas={b}:没找到已加载的 libopenblas"
            for libname, n in leg["blas_lib"].items():
                assert n == 1, (
                    f"order={o} blas={b}:{libname} 实际跑 {n} 线程 —— "
                    "运行时固定没生效(多半是 /proc/self/maps 没找到库)")

    def test_the_answer_is_bit_identical_across_blas_settings(self, legs):
        """★★★ 本门的判据:**逐位相同** —— 每种 import 顺序内,跨外部 BLAS 档。

        用 `repr(float)` 比字符串,而不是数值容差 —— 这里要的就是**位相同**,
        任何容差都会把「归约顺序改变了答案」这件事放过去。
        ★ G-DOMAIN:相反结果 = 两者不同 ⇒ **固定没生效**,处置是修固定机制,
          **不是放宽本判据**。
        """
        for o in _ORDERS:
            a, b = legs[o, 8], legs[o, 16]
            assert a["cl"] == b["cl"], (
                f"order={o}:BLAS 归约顺序改变了答案:cl {a['cl']} vs {b['cl']}\n"
                f"  库线程数:{a['blas_lib']} vs {b['blas_lib']}")
            assert a["resid"] == b["resid"], (
                f"order={o}:残差不逐位相同:{a['resid']} vs {b['resid']}")

    def test_the_two_import_orders_agree_with_each_other(self, legs):
        """两种顺序都钉到 BLAS=1 ⇒ 同一外部档下,两顺序的答案也必须逐位相同。

        ★ 这一条锁的是「修复没有把 pyfp3d-first 的老路径改坏」:
        env 路径与运行时路径必须收敛到**同一个** BLAS=1,而不是两个各自"1"的变体。
        """
        for b in _BLAS_SETTINGS:
            a, c = legs["pyfp3d_first", b], legs["numpy_first", b]
            assert a["cl"] == c["cl"] and a["resid"] == c["resid"], (
                f"blas={b}:import 顺序改变了答案:cl {a['cl']} vs {c['cl']}, "
                f"resid {a['resid']} vs {c['resid']}")
