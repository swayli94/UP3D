"""
pyFP3D: 3D Unstructured-mesh Full-Potential Transonic Flow Solver

A Python + Numba implementation of the full-potential (FP) method for steady 
external flows over wings at transonic Mach numbers (0.3–0.87).

Core modules:
  - mesh: I/O, topology, metrics, coloring
  - physics: Isentropic constitutive relations
  - kernels: Residual assembly, artificial compressibility
  - solve: Linear and nonlinear solvers
  - post: VTK output, visualization

Reference: docs/design.md
Roadmap: phases/p1/docs/roadmap.md
"""

__version__ = "0.0.1"
__author__ = "swayli94"

import os

#: ★★★ W2/H30(2026-09-06,使用者裁决):**把 BLAS 固定成单线程。**
#:
#: 实测(`bench/studies/gate_audit_20260905/results/thread_attribution_D05_M080.csv`,
#: D05 的 M0.80/α1.25 medium,9 条腿单变量):
#:   * 按 **BLAS** 线程数分组,每档只有一个 cl 值   {1:1, 8:1, 16:1}
#:   * 按 numba 线程数分组,每档有多个 cl 值        {1:3, 8:2, 16:2}
#:   ⇒ **结果是 BLAS 线程数的函数;numba 不是自变量**(1/8/16 逐位相同,
#:     与 `tests/A/test_A02` 早就断言的「着色装配跨线程数位可复现」一致)。
#:   ⇒ 同配置重复跑逐位相同 ⇒ **不是竞态**,是多线程 BLAS 的**归约顺序**。
#:
#: 这几道门不是效率门,**结果依赖线程数本身就是错误**(使用者裁决)。固定成 1
#: 之后没有自由度,结果**位确定**;而 numba 的并行(真正的热点:装配、迎风扫描、
#: σ 输运)**一点不受影响** —— 实测 A/B/C/F 四层墙钟只涨 **4.5 %**。
#:
#: ★★★ **单变量分离的结果(2026-09-06 补做)**:自变量**只有
#: `OPENBLAS_NUM_THREADS`**。固定 `NUMBA_NUM_THREADS=8`、分离 OMP 与 OPENBLAS:
#:   OMP 1→8 而 OPENBLAS 固定 ⇒ cl **逐位不变**(0.34009652307781807)
#:   OPENBLAS 1→8 而 OMP 固定 ⇒ cl **变**(0.34011432026323)
#: ⇒ `OMP_NUM_THREADS` 无影响。★ 第一次归因时我把 OMP/OPENBLAS/MKL **三个一起
#: 改**,那是「把两个旋钮当一个拧」——结论侥幸没错,但当时并未被分离验证过。
#: 其余几个变量一并固定只是**无害的保险**,不是实测的自变量。
#:
#: ★★ **强制覆盖,不是 `setdefault`。** 目标恰恰是「外部环境不该决定结果」,
#: 而 `setdefault` 会把决定权留给外部环境,等于没固定。
#: ⇒ 逃生舱 `PYFP3D_ALLOW_BLAS_THREADS=1`:**只给故意做跨 OPENBLAS A/B 的探针用**
#: (归因就是这么做的)。`tests/F/test_F09` 的 G-TEETH 也用它 —— 于是验证「固定
#: 会不会失效」**不需要改库文件**。
#: ★★ 这几行只在 numpy **尚未导入**时有效,所以必须待在本文件顶部,而本文件
#: **不 import numpy**。真正的验收是**行为的**,见 `tests/F/test_F09`;
#: 「环境变量被设了」不是判据(那是 F06 那族字面量对字面量的错误)。
if os.environ.get("PYFP3D_ALLOW_BLAS_THREADS") != "1":
    for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[_v] = "1"

# Development mode flag: PYFP3D_NOJIT=1 disables Numba JIT for debugging
NOJIT = os.environ.get("PYFP3D_NOJIT", "0") == "1"
