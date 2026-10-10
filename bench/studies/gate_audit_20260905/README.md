# `gate_audit_20260905` — 2026-09-05 独立门审计的两个探针

判定 note：[`docs/dev_phase_six/20260905-0200-independent-gate-audit.md`](../../../docs/dev_phase_six/20260905-0200-independent-gate-audit.md)

    python bench/studies/gate_audit_20260905/run_probes.py            # 两个都跑
    python bench/studies/gate_audit_20260905/run_probes.py --probe A  # ~30 s
    python bench/studies/gate_audit_20260905/run_probes.py --probe B  # ~8 min @8t

★ `pyfp3d/` **未改动**；两个探针都从 `tests/` 里 import 被审的门自己的常数与配方
（C06 的几何常数、D05 的 `RECIPE`），**不重打一遍** —— 免得探针与它测量的门口径分叉。

★ **G-CADENCE**：这是 `bench/` 下的一次性研究，**不在任何周期上**。它的结论若要
持续有效，必须按 note §10 建议 3 变成 `tests/` 里的门。

## `results/farfield_realization_c06.csv` — 探针 A

C06 带环量圆柱有闭式解 ⇒ `cl_rel` 是**真误差**。固定壁面间距、扫
`r_far ∈ [20, 55]`（该范围内远场早已收敛 ⇒ `r_far` 是**惰性旋钮**，唯一改变的是
gmsh 返回哪张三角剖分）。列：`h_wall, r_far, n_nodes, n_wall, cl, cl_exact,
cl_rel, cp_rms, cp_max`。

读数：h 0.02 极差 **0.10 pp**（远场确实收敛了）· h 0.04 极差 **2.09 pp**
⇒ **coarse 档的误差由抽到哪张网格决定，不是由 h 决定。**

★ 探针原本是去证实"`r_far = 15` 截断限制了 C06 的 ~1 % 残差"的 —— **它否掉了那个
假设**，而在否掉的过程中量到了上面那件更要紧的事。

## `results/nonuniqueness_naca0012_alpha0.csv` — 探针 B

NACA0012、**α = 0**、`RECIPE` 逐字取自 D05。四条 `arm`：

| `arm` | 变的是 |
|---|---|
| `cold` | 冷启动，M 扫描 |
| `ramp` | 向上暖启动链（前一个 M 的收敛解作初值） |
| `cold_no_entropy` | 同 `cold`，`entropy_correction=False` |
| **`seed_sweep`** | ★ **固定 M = 0.86，只改 `n_picard_seed`（只改初值）** |

**只读 `converged == 1` 的行。** 未收敛的行 `n_limited = n_floored = 0`
⇒ 不是钳制模式；但本探针**不读 `accept_reason`**，所以不给它们命名模式。

读数（coarse，`seed_sweep`，四条全部收敛到 |R| ≈ 3e-13、零钳制）：
**cl ∈ [−0.0059, +0.4360]** —— 同一个离散问题、只改初值的四个不同解。
参照物是 D05 自己的 α=0 基线 |cl| ≈ 1.2e-3（列 `cl_over_alpha0_baseline`）。

文献：Steinhoff & Jameson AIAA-81-1019；Salas, Jameson & Melnik NASA TP-2385。
