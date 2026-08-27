# `cases/gates/` —— C / D 类门的**被跟踪**证据

★★★ **打开这一个文件就能把全部 C/D 门的图肉眼过一遍。**

| | |
|---|---|
| 谁写它 | `tests/C/*` 与 `tests/D/*`，经 `gate_evidence_dir` fixture |
| 何时写 | ★ 只在 **`PYFP3D_GATE_FIGURES=1`** 时；平时跑**不写盘**，断言对着已提交的 `summary.csv` |
| 为什么这样 | 图与断言来自**同一次计算**（构造保证）；平时不写 ⇒ 一次普通 `pytest` **不脏工作树** |
| 刷新 | `PYFP3D_GATE_FIGURES=1 pytest tests/C tests/D` —— 数字变了要走再基线的**勘误清单**纪律 |

★★ **A / B / E / F 不在这里，而且不该在**：它们的判据**就是断言本身**，没有「肉眼检查」这个环节。
2026-08-24 之前那些图写进 gitignored 的 `artifacts/`，于是
**规则要求「每个视觉门留下产物」而 `.gitignore` 保证它永不进 HEAD** —— 按纪律 3 从来就不是证据，
却有 11 处文档把它当证据引用。`artifacts/` 与它的 fixture **已一并删除**。

★ 目录里出现的 `G1.3` / `G2.4` 这类内层名字是**历史门号**，刻意保留 —— 它把新门号接回旧证据链。

---

## C01_laplace_mms

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C01_laplace_mms/summary.csv`](C01_laplace_mms/summary.csv) — 4 行，列 `n, h, l2_error, cg_iterations`

![mms_convergence_loglog.png](C01_laplace_mms/mms_convergence_loglog.png)

*C01_laplace_mms/mms_convergence_loglog.png*


## C02_cylinder_wall_correction

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C02_cylinder_wall_correction/G1.3/summary.csv`](C02_cylinder_wall_correction/G1.3/summary.csv) — 5 行，列 `level, h_wall, max_cp_err_uncorrected, mean_cp_err_uncorrected, max_cp_err_corrected, mean_cp_err_corrected, max_rhs, spanwise_floor_uncorrected, spanwise_floor_corrected, n_cg_iterations, residual_norm`

![cp_theta_overlay.png](C02_cylinder_wall_correction/G1.3/cp_theta_overlay.png)

*C02_cylinder_wall_correction/G1.3/cp_theta_overlay.png*

![error_vs_h.png](C02_cylinder_wall_correction/G1.3/error_vs_h.png)

*C02_cylinder_wall_correction/G1.3/error_vs_h.png*

![normal_deviation.png](C02_cylinder_wall_correction/G1.3/normal_deviation.png)

*C02_cylinder_wall_correction/G1.3/normal_deviation.png*

![section_symmetry_plane.png](C02_cylinder_wall_correction/G1.3/section_symmetry_plane.png)

*C02_cylinder_wall_correction/G1.3/section_symmetry_plane.png*


## C03_laplace_sphere

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C03_laplace_sphere/summary.csv`](C03_laplace_sphere/summary.csv) — 5 行，列 `metric, value`

![sphere_cp_meridian.png](C03_laplace_sphere/sphere_cp_meridian.png)

*C03_laplace_sphere/sphere_cp_meridian.png*


## C05_nozzle_quasi1d

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C05_nozzle_quasi1d/summary.csv`](C05_nozzle_quasi1d/summary.csv) — 4 行（3 级阶梯 + **反例腿**），列 `nx, h, n_max, converged, reason, n_newton, residual, x_shock, x_s_exact, err_x, err_cells`

★★★ **项目里唯一带激波的精确解。** 左图三级 u(x) 全部贴在精确准一维解上，激波稳定落在
精确位置 x_s = 12 的**上游**约 0.7 个单元；右图 `|err_x|` 沿一阶线走（实测阶 **0.927 / 0.834**），
而 `|err_cells|` 几乎是平的（0.651 / 0.685 / 0.768）—— **位置在收敛，捕捉宽度保持亚单元**。

★★★ **CSV 第 4 行不是阶梯的一部分，它是本门最值钱的一条**：从 x_s = 8 的解出发
（边界数据仍是 x_s = 12 的真值），nx = 200 收敛到 **`reason=tol`、|R| = 1.3e-15**（机器零），
而激波停在 **x = 7.98**，错 **−40.17 个单元**。而 `verify_uniqueness` 证明**连续解唯一**
（Δφ 19.3186 > 17.6307）⇒ 那是**离散算子的伪根**。**在一个有真值的案例上，
`converged` 被证明不是正确性证书** —— `bench/usability.py` 在翼型上只能拿共识当锚点。

★★ **h = 0.05 那条腿 `reason=cap`，本门仍然用它 —— 许可是门自己证的**：迭代上限
20 → 40 → 80 时残差只在 ~2e-8 的地板上晃（4.45e-08 → 1.89e-08），而 **x_shock 只动
1.3e-06 = 2.6e-05 个单元** ⇒ 对**位置**可用、对**残差**不可用。

![c05_nozzle_shock.png](C05_nozzle_quasi1d/c05_nozzle_shock.png)

*C05_nozzle_quasi1d/c05_nozzle_shock.png*


## C06_lifting_cylinder

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C06_lifting_cylinder/summary.csv`](C06_lifting_cylinder/summary.csv) — 3 行，列 `level, h_wall, n_nodes, n_wall, cp_max, cp_rms, cl_p, cl_exact, cl_rel, residual, n_kutta`

★★★ **本门的核心是「给进 Γ、量出 cl」** —— 圆柱没有尖尾缘，Kutta 选不出环量，所以
Γ 被**规定**（`gamma_fixed`，实测 `n_kutta_updates == 0`）；右图里 `cl rel err` 那条线
走的是**表面 Cp 恢复 + 面积分**，对的是 Kutta–Joukowski 的 `2Γ/(U c)`，两条路毫无交集。

★★ **右图标题「all scales refined together」是判据的一部分，不是描述**：实测把
`h_far` 钳成固定值（只加密壁面）会让 Cp 阶从 **1.686 塌到 0.325**、cl 误差**完全不动**
（1.124 % → 1.136 %）—— 与 P11 在球上测到的固定体网格污染地板同一形状。

★ 右图 `cl rel err` 在两个粗级之间**是平的**（3.393 % → 3.477 %，阶 −0.035）：
xcoarse 只有 598 节点、`h_far` 12 而域半径 15，仍在**前渐近区**。**如实画出、记为
RECORDED，不纳入收敛判据** —— 收缩只从 coarse → medium 判（实测 3.06×）。

![c06_lifting_cylinder.png](C06_lifting_cylinder/c06_lifting_cylinder.png)

*C06_lifting_cylinder/c06_lifting_cylinder.png*


## C07_karman_trefftz

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C07_karman_trefftz/summary.csv`](C07_karman_trefftz/summary.csv) — 3 行，列 `level, h_wall, n_nodes, n_wall, n_kutta, gamma, gamma_exact, gamma_rel, cl, cl_exact, cl_rel, cp_rms, band_LE, band_MID, band_TE, n_LE, n_MID, n_TE, match_max`

★★★ **唯一能把 Kutta 条件对精确 Γ 检验的门。** C06 的圆柱没有尖尾缘、Γ 是**规定**的；
KT 有**有限尾缘角**，Kutta **必须自己选出** Γ（`n_kutta` 列 = 2，不是 0），而精确值有解析式。
medium：**Γ 2.803 % / cl 3.045 %**，两条独立路线（Kutta 行 vs 表面压力积分）**相差 8 %**。

★★★ **右图的标题就是本轮的发现：积分量收敛，前缘点值不收敛 —— 而且在同一个解上。**
Γ 阶 2.109 / 1.493、cl 阶 1.957 / 1.468、Cp **MID** 带 2.893 / 1.761；而 **LE 带阶 0.085**、
**TE 带 0.53 / 0.65**。两个慢带都有已定位的成因：
· **LE** = 本项目 **G1.6 那一族**（P1 场在前缘吸力峰上的固有能力），**首次在有精确解的
二维翼型上复现**；排除过两个混淆项 —— 折线几何点数 ×20 只让 LE 阶 0.241→0.287，
oracle 采样加密 5×（匹配距 3.5e-05→6e-06）只让 LE RMS 0.1594→0.1574。
· **TE** = **有限楔角本身**：精确解 |q| 以 (ζ−b)^(τ/π) 趋零，τ=10° 时指数只有 0.0556，
θ=1e-5 处精确 |q| 还有 0.49 ⇒ **精确解自己就没有干净的驻点**。

★★ **紫线（LE 带）只有两个点，不是漏画**：xcoarse 上 x<0.05 内**一个壁面节点都没有**
（`n_LE` 列 = 0）。**自变量改变了样本是否存在** ⇒ LE 带只判「显著慢于 MID」（实测 20.6×），
不判它自己的阶。

![c07_karman_trefftz.png](C07_karman_trefftz/c07_karman_trefftz.png)

*C07_karman_trefftz/c07_karman_trefftz.png*


## C08_ringleb

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C08_ringleb/summary.csv`](C08_ringleb/summary.csv) — 5 行（`upwind_c` 扫描 + `m_crit` 腿），列 `leg, h, upwind_c, m_crit, converged, reason, n_newton, residual, r_exact, rms, rms_sup, rms_sub, n_sup, n_cell, m_num_max, m_ex_max`

★★★ **唯一 2-D 跨声速精确解，而且是全速势方程本身的解。** 与 C05 的分工是它存在的理由：
喷管里**有激波** ⇒ 数值解本该与等熵精确解不同，那里只能判**激波位置**；Ringleb 的超声速
口袋是**无激波、光滑**的 ⇒ **人工密度加的每一分耗散都是纯误差**。

★★★ **左图把归因画死了**：速度误差**几乎只存在于精确声速线（黑）圈出的口袋里** ——
口袋内低估、紧下游高估，其余全场近零。右图是**剂量-响应**：超声速区 RMS
**C=1.5 → 0.1001**、C=0.5 → 0.0477、C=0.25 → 0.0376（收缩 **2.66×**）。
**证明那 10 % 是人工密度而不是离散化的，是这条剂量-响应，不是论证。**

★★ **红叉是价格**：`upwind_c = 0` **不收敛**（60 步撞顶，|R| 2.5e-07）⇒ 那份耗散买的是
稳健性，现在有价格。**绿三角**：`m_crit` 从 0.95 抬到 1.20 把超声速误差降到 **0.0293**
（3.4×）且仍收敛，代价是 Newton 步数 11 → 37 —— **RECORDED 而非建议改默认**，
因为本门只测了一个**无激波**算例。

★ **收敛阶判在别处**：解的误差由模型误差主导、**不随 h 收敛**（6.4e-2 → 5.9e-2 → 7.3e-2），
所以阶判在 **`r_exact` 列**（把精确解代进离散残差，实测 **2.32 / 2.28 阶**）——
那一半**不需要求解**。拿解的误差判阶会把模型误差误报成离散化缺陷。

![c08_ringleb.png](C08_ringleb/c08_ringleb.png)

*C08_ringleb/c08_ringleb.png*


## C09_prandtl_glauert_peak

**判据读的 CSV**（平时跑就是对着它断言）：
- [`C09_prandtl_glauert_peak/G3.1/summary.csv`](C09_prandtl_glauert_peak/G3.1/summary.csv) — 1 行，列 `cp_peak_incompressible, cp_peak_pg_corrected, cp_peak_compressible, rel_diff_pct, n_picard, picard_converged, mach2_max`

![v3_1_sphere_cp.png](C09_prandtl_glauert_peak/G3.1/v3_1_sphere_cp.png)

*C09_prandtl_glauert_peak/G3.1/v3_1_sphere_cp.png*


## D01_naca0012_incompressible_panel

**判据读的 CSV**（平时跑就是对着它断言）：
- [`D01_naca0012_incompressible_panel/G2.2/summary.csv`](D01_naca0012_incompressible_panel/G2.2/summary.csv) — 496 行，列 `x, phi_minus, phi_plus, jump`
- [`D01_naca0012_incompressible_panel/G2.3/summary.csv`](D01_naca0012_incompressible_panel/G2.3/summary.csv) — 2 行，列 `level, cl_pressure, cl_gamma, cl_panel_ref, rel_err_pct, n_kutta_updates`
- [`D01_naca0012_incompressible_panel/G2.3/v2_3_cp_curves.csv`](D01_naca0012_incompressible_panel/G2.3/v2_3_cp_curves.csv) — 408 行，列 `x_c, cp, surface`
- [`D01_naca0012_incompressible_panel/G2.4/summary.csv`](D01_naca0012_incompressible_panel/G2.4/summary.csv) — 1 行，列 `cl_pressure, cl_gamma, rel_diff_pct`
- [`D01_naca0012_incompressible_panel/G2.5/summary.csv`](D01_naca0012_incompressible_panel/G2.5/summary.csv) — 4 行，列 `level, case, max_abs_w_over_uinf, p99_abs_w_over_uinf, rms_abs_w_over_uinf`

![v2_1_residual_heatmap.png](D01_naca0012_incompressible_panel/G2.1/v2_1_residual_heatmap.png)

*D01_naca0012_incompressible_panel/G2.1/v2_1_residual_heatmap.png*

![v2_2_wake_jump.png](D01_naca0012_incompressible_panel/G2.2/v2_2_wake_jump.png)

*D01_naca0012_incompressible_panel/G2.2/v2_2_wake_jump.png*

![v2_3_cp_vs_panel.png](D01_naca0012_incompressible_panel/G2.3/v2_3_cp_vs_panel.png)

*D01_naca0012_incompressible_panel/G2.3/v2_3_cp_vs_panel.png*

![v2_4_cl_crosscheck.png](D01_naca0012_incompressible_panel/G2.4/v2_4_cl_crosscheck.png)

*D01_naca0012_incompressible_panel/G2.4/v2_4_cl_crosscheck.png*

![v2_5_spanwise_w_heatmap.png](D01_naca0012_incompressible_panel/G2.5/v2_5_spanwise_w_heatmap.png)

*D01_naca0012_incompressible_panel/G2.5/v2_5_spanwise_w_heatmap.png*


## D02_naca0012_m05_panel

**判据读的 CSV**（平时跑就是对着它断言）：
- [`D02_naca0012_m05_panel/G3.2/summary.csv`](D02_naca0012_m05_panel/G3.2/summary.csv) — 2 行，列 `level, cl_pressure, cl_gamma, cl_ref_pg, cl_ref_kt, cl_ref_mid, rel_err_mid_pct, n_picard, n_solves, mach2_max`
- [`D02_naca0012_m05_panel/G3.2/v3_2_cp_curves.csv`](D02_naca0012_m05_panel/G3.2/v3_2_cp_curves.csv) — 408 行，列 `x_c, cp, surface`

![v3_2_cp_mach_vs_reference.png](D02_naca0012_m05_panel/G3.2/v3_2_cp_mach_vs_reference.png)

*D02_naca0012_m05_panel/G3.2/v3_2_cp_mach_vs_reference.png*

![v3_3_picard_residual.png](D02_naca0012_m05_panel/G3.2/v3_3_picard_residual.png)

*D02_naca0012_m05_panel/G3.2/v3_3_picard_residual.png*


## D03_naca0012_m080_shock

**判据读的 CSV**（平时跑就是对着它断言）：
- [`D03_naca0012_m080_shock/G4.1/summary_coarse.csv`](D03_naca0012_m080_shock/G4.1/summary_coarse.csv) — 27 行，列 `quantity, value`

