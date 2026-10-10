# `bench/` 留下的 12 个 —— **按角色命名，不再按阶段**

★★ **2026-08-24**：`bench/` 从 **99 个脚本 / 48 个 studies** 收敛到 **12 + 7**。
判据是三条（缺一不可，见 [phase-6 普查](../docs/dev_phase_six/20260824-0600-bench-audit.md)）：
① `tests/` 的传递闭包 ② 有 cadence 的 runner / 工具 ③ 产品指标门。

| 文件 | 角色 | 备注 |
|---|---|---|
| `recipes.py` | **配方** | ★ 2026-08-24 新建（使用者裁决 a）：`NEWTON_M6_RECIPE` 原先住在 `tests/test_p8_newton.py`，而 `bench/` 与 `cases/demo/` 从那里 import —— **方向反了**。谁都从这里读 |
| `failure_modes.py` | **共享库** | ★★★ 原名 `run_le14_common_root.py` —— **被 12 个脚本 import 的失败分类器，名字却是 phase-2 的轮次号「LE-14」**。`classify_failure` 在这里 |
| `capability_matrix.py` | **共享库** | 原名 `run_capability_matrix.py`；★ 去掉 `run_`，它是库不是 runner |
| `usability.py` | **共享库** | 答案锚点（R23）；锁在 `tests/test_r23_usability.py` |
| `bitcheck.py` | **工具** | 位相同比较 |
| `run_capability_locks.py` | ★ **runner** | 快层入口，**每次收口跑**。`run_` 前缀在这里是对的 |
| `run_m1_gate.py` | **产品指标门** | M1；被 `tests/test_m1_respec.py` 加载 |
| `run_m3_budget.py` | **产品指标门** | M3。★ 曾计划拆成门+库两半，**复查后撤回** —— 19 个 importer 里 18 个本身要归档 |
| `run_g82_anchor_check.py` | **产品指标门** | G8.2 锚点复核 |
| `stamp_mesh_axes.py` | **工具** | 网格清单轴戳 |
| `run_bench.py` | **工具** | 通用基准 |
| `__init__.py` | — | ★ 使 `bench` 成为 package，`from bench.recipes import ...` 才成立 |

★ `gate_results/` **留在这里**：能力边界与快层按路径引用其中的 CSV；归档脚本**横向**去取它。
