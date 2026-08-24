# `bench/studies/` —— **只剩被活测试实际加载的 7 个**

★★ **2026-08-24 归档**：原本 48 个目录，其中 **41 个是一次性研究**，已随各自阶段移入
`phases/p*/bench/studies/`。留在这里的 **7 个**，每一个都**被 `tests/` 里的某个文件在运行时
`spec_from_file_location` 加载**（不是提及，是执行）：

| 目录 | 被谁加载 |
|---|---|
| `v5_1c_above_band_window/` | `tests/test_v5_above_band_seed.py` |
| `v5_1d_near_band_window/` | `tests/test_v5_near_band_seed.py` |
| `v5_2_rae2822/` | `tests/test_meshgen_rae2822.py` |
| `v5_tight_coupling/` | `tests/test_v5_tight_scaled.py` |
| `v6_1_wake_sheet/` | `tests/test_v6_wake_sheet.py` |
| ★ `v3_loose_coupling/` | **传递依赖** —— `v5_2_rae2822/run.py` `_load` 它 |
| ★ `v5_1b_scaled_newton/` | **传递依赖** —— `v5_1c/v5_1d` 用它的 `_p_series` |

## ★★★ 这两条传递依赖是本轮的教训，也是上一轮的

2026-08-10 那次整理写着「**保留规则原来只算了一跳**」，因此 15 个目录搬回原位。
本轮**同一个坑又出现了一次**：种子 5 个，传递闭包 **7 个** ——
`v3_loose_coupling` 与 `v5_1b_scaled_newton` **一跳规则会漏掉**。

★★ 而漏掉它们的**不是规则，是仪器**：这些测试的路径是
`os.path.join(os.path.dirname(__file__), "..", "bench", "studies", "<name>", "run.py")`
—— **逐段拼出来的**，一个找 `"bench/studies/..."` 字面量的正则**在原理上看不见它**。
第一次归档时收集直接报 2 个 `FileNotFoundError`，才暴露出来。

⇒ **口径**：判断「谁读谁」时，**路径可以是拼出来的**；只认字面量的判据不完整。
本轮最终用的是「目录名以任何形式出现在**去 docstring 后的代码**里」，再取传递闭包。

## 找已归档的研究

| 阶段 | 位置 | 内容 |
|---|---|---|
| phase 2 / Track A·V | [`phases/p2/bench/`](../../phases/p2/bench/) | `a1..a4`、`b9`、`b26`、`c1`、`p14`、`v1`、`v2`、`v5_3/5_4/5_5`、`v5_ibl_floor`、`v5_m6_bridge`、`v6_2` 等 |
| phase 3 | [`phases/p3/bench/`](../../phases/p3/bench/) | `run_task3_*`（26）· `run_hex_*`（3）· `run_gs40*`（10） |
| phase 4（GS4.1） | [`phases/p4/bench/`](../../phases/p4/bench/) | `gs41_*`（14 个 studies） |
| phase 5（M1） | [`phases/p5/bench/`](../../phases/p5/bench/) | `r11..r22`、`m1a/m1b/m1e`、`n1/n2`（17 个 studies） |

★ **证据仍被 git 跟踪** —— 归档是目录重组，不是取消证据。
★★ **但归档脚本不保证能跑**：见 [`phases/README.md`](../../phases/README.md) 的可运行性口径与
**恢复点 commit**。
