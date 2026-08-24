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

