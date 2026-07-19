# pyFP3D 能力盘点与批判性审查（2026-07-15）

> ⚠ **SUPERSEDED（快照声明，2026-07-19 追加）**：本文是 2026-07-15 的快照，其后
> P14/B9/B14–B21/A1–A3 大量落地，多处结论已被覆写——尤其"LS medium M0.84 已演示"
> 一项经历了 B20 重基线回退（GB20.7）→ **B21 恢复**（N1 freeze 捕获修复，γ 0.088343）
> 的两次翻转。现状以 docs/roadmap/track_*.md 台账 + docs/overview.md 为准；
> 本文仅作历史参考，引用其数字前先对台账。

**性质**：横向审查报告（非 gate 文档，不改变任何 roadmap 状态）。是
[capability_review_2026-07-14.md](capability_review_2026-07-14.md) 的滚动更新——那份撰写后一天内又落地了
B11/B12/B13（2026-07-14）与 **M6 medium level-set 工作流（2026-07-15）**，其
§0/§1/§2.3/§3/§4/§5.1 已被这些工作实锤性地覆写（见 §7）。

**本次新增重点（用户指令）**：
1. 把测试算例中的**亚声速 / 跨声速情形明确分列**（§2）；
2. 给出**各情况的计算时间**，并**对照目标计算时间**判读（§2、§8）——
   目标 = 2D < 单核 1 min，3D < 1~4 核 10 min。

**方法**：三路独立审计交叉验证——① 测试与 demo 工件逐个核对（含 committed
CSV 里的实测墙钟）、② 五文档的时效性与相容性核对、③ 直接读代码建立两路径能力
矩阵与孤儿功能排查。凡代码与文档不符以代码为准；凡断言无 committed 工件按项目
自己的规则降级。所有结论落到 `file:line` 或 CSV。

**★ 本报告发布当日的追加（B15，2026-07-15）：** 报告写完后，用户裁决并落地了
**Track B B15（LS Newton 跨声速续接 + N5 冻结选择）**，它直接改写了本文 §2.3/§4/§6.1
的三条判断，并给出**首个针对"慢"的根因结论与量化提速**：
1. **24–38 min 的根因不是内层线性求解器，而是 Picard 的激波位置残差 plateau** —— M6
   medium 的 0.80/0.84 两级根本不收敛、烧满 200 outer 预算（`tol_residual` 已被钉在
   plateau 之上的 1e-5）。这是 Picard 的**固有属性**，调参无法绕开。
2. **Newton ramp 消除 plateau**：NACA coarse M0.80/α1.25（B6 gate 工况）
   **Picard 44.0 s → |R| 1.55e-5（3/5 级不收敛，不是解）** vs
   **Newton 8.1 s → |R| 3.1e-12 严格收敛（5.2×）**；
   **加 `intermediate_tol` 6.8 s（6.5×）且最终 γ 逐位相同**（0.212445）。
   *（A3 勘误 2026-07-18：本行三个墙钟值是 CSV 前的试跑值，且 44.0/8.1 = 5.4 而非
   文中的 5.2×。已提交 CSV 的真值为 **41.9 s / 7.5 s（5.6×）/ 6.5 s（6.5×）**。
   本文件是有日期的快照，正文保持原样，仅加此勘误；权威值见
   `cases/demo/b15_ls_newton_ramp/results/checks.csv`。）*
3. **LS Newton 的 `freeze`/ramp 不再是空白**（原 §6.1 "完全无代码"条目作废）。
   ★ 但 conforming 的冻结触发器**不可平移**（`live_stalled` 会冻结"仍在移动"的指派
   ⇒ 发散/revert）；LS 侧只能用 `freeze_tol` 武装。详见 roadmap B15 / design_track_b §14。
**仍未解决**：LS Newton 的 γ 比同网格 conforming-Newton 真值低 **7.4%**（B6 记为 ~13%）
—— 这是**离散化差异**，B15 使其**首次可在两侧都严格收敛下干净测量**，但未消除。

**⚠ 计算时间的口径警告（贯穿全文，最重要的一条批判）**：本仓库**所有**已测
墙钟都是在 **16C/32T 共享工作站上以 16 线程（含 BLAS/OMP）** 跑出来的
（`NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16`，见
`cases/demo/m6_medium_ls_workflow/run_demo.py:35`，MEMORY `cap-parallelism-16-threads`）。
**用户目标口径是 2D 单核 / 3D 1~4 核**。两者差一个数量级的核数。因此下文所有
时间在对照目标时都要记住：**目标核数下的真实时间会远大于表中数字**，当前没有
任何单核 / 4 核基准存在。这不是小事——它意味着"求解器有多快"这个问题**至今
没有在目标硬件口径下被测过**。

---

## 0. 一页结论（更新版）

**仓库仍是两条刻意平行的求解路径 + 一个零代码的粘性支柱；但 2026-07-14/15 的
B11/B12/B13 + M6 工作流已经把"level-set 路径卡在 spsolve 墙、无 3D Newton、
只有 coarse"这一条旧结论显著改写。**

- **conforming 路径**：全部求解器实力（全耦合 Newton、amg/ilu/direct + lagged-LU、
  Mach 续接 + freeze、Γ(z) 涡远场、翼尖 taper、Cp 平滑）。3D 跨声速真解只有它
  给（M6 medium Newton **250 s @16 线程**，cl_KJ 0.2692）。固有病：尾迹是网格
  实体（α 扫掠要重划网格、机身相交处理不了、翼尖自由边 P13）。
- **level-set 路径（Track B）**：全部工作流实力（隐式 Kutta、wake-free 网格、
  三种远场、α 免重网格机制）。**旧的"硬编码 spsolve / 无 3D Newton / 只有
  coarse"结论已过期**：B11 给了 `precond="ilu"/"amg"`（ILU 是有效逃生口），
  B12 给了 Picard/Newton 的 lagged-LU 并**在 M6 medium（3D）上真跑过 LS
  Newton**（M0.5，7 步，γ 逐位一致），B13 把 lagged-LU 上到 Picard 外层
  （M6 medium lifting **447.6→68.3 s，6.55×**），2026-07-15 的 M6 工作流**首次
  在 medium 网格上给出 LS 跨声速解**（M0.84）。**仍缺**：LS Newton 无 Mach-ramp
  包装、`freeze` 仍是空占位、**从无 fine 级研究**。
- **Track V（边界层修正）零代码**——全库 grep `viscous|IBL|transpiration|
  boundary_layer|displacement_thickness` 无命中。用户目标的第二支柱整体空白。
- **网格收敛**只在三处成立：MMS、2D 亚声速升力（G9.2）、3D 亚声速圆帽 M6
  （G13.3，p=2.31，全库唯一 3D Richardson）。**3D 跨声速在两种翼尖几何上都失败，
  至今无 M0.84 合法外推值。**
- **计算时间对照目标：全面不达标，且从未在目标核数口径下测过**（§2/§8）。跨声速
  2.5D 是 **10~40 min**（目标 2D 单核 1 min）；3D 跨声速 medium 是
  **4 min（conforming Newton@16 线程）~ 24–38 min（LS@16 线程）**（目标 3D 1~4 核
  10 min），fine 不可达。

---

## 1. 现在能算什么（含亚/跨声速与实测时间）

按"有 committed 工件"划界。凡 Picard 跨声速结果带 **P4/P5 勘误**：Picard 工程
收敛态**不是离散方程的解**（Newton 残差 coarse 2.2e-4 / M6 ~8e-6），真解须 Newton。
**所有时间 @16 线程**（口径警告见文首）。

| 问题 | 声速区 | 路径 | 网格 | 结果 | 实测时间@16线程 | 工件 |
|---|---|---|---|---|---|---|
| 自由流保持（含尾迹 cut） | — | 两条 | 全部 | 残差 8.8e-14 | 秒级 | `p1_laplace/` |
| Laplace MMS（阶 1.96/1.94LS） | — | 两条 | 三级 | 收敛 | 秒级 | `p1_laplace/`,`test_b2` |
| NACA0012 2.5D 升力 | **亚(M0.3/0.5)** | 两条 | coarse–fine | cl 0.28437∈[PG,KT]，LS 同网格 Γ<1% | 秒级 | `p3_subsonic/`,`b3_levelset_lifting/` |
| NACA 2.5D 网格收敛 | **亚** | conf | 三级 | 误差 2.71→0.33→0.03%（干净） | — | `p9_grid_discrimination/` |
| NACA 2.5D 跨声速 | **跨(M0.74–0.82)** | conf Picard+Newton | coarse+medium | Newton 真解 coarse M0.80 shock0.658/cl0.459；medium M0.7875 0.674/0.523；**medium M0.80 无孤立解(FP折叠)** | P4 medium G4.1 **~16 min**；P8 Newton coarse 分钟级 | `p4_transonic/`,`p8_newton/` |
| NACA 2.5D 跨声速，LS 路径 | **跨(M0.70/0.80)** | LS Picard(+Newton 2.5D) | **coarse** | coarse M0.80 达标(vs 同网格 Newton 真值)；medium 折叠区 Picard 停 Γ−18.8% | B6 全部解 **~25 min 上限** | `b6_transonic/` |
| ONERA M6 M0.84 | **跨** | conf Picard | coarse+medium | cl_p 0.2419/0.2453（Picard 质量） | P5 medium 从头 **~75.6 min(4539 s)** | `p5_onera_m6/` |
| ONERA M6 M0.84 真解 | **跨** | conf Newton | medium | cl_KJ 0.2692 vs Tranair 0.288，gap0.019 | **250.8 s**(mesh7.3+solve239.8) | `p8_newton/results/g82_m6_medium.csv` |
| ONERA M6 M0.5 亚声速 | **亚** | conf Newton | coarse/fine | — | coarse **7.6 s**；fine **5294 s(13级)** | `p9.../g91_m6_levels.csv` |
| ONERA M6，LS 路径 亚 | **亚(M0.5)** | LS Picard/Newton | **coarse+medium** | medium cl 0.212；Newton γ 逐位=spsolve | 见 §2 lagged-LU 行 | `b7_`,`b11_`,`b12_`,`b13_`,`m6_medium_ls_workflow/` |
| ONERA M6，LS 路径 跨 | **跨(M0.84)** | LS Picard | **coarse + (NEW)medium** | coarse cl_KJ0.2765/0.2710；**medium cl0.2765(wake-free)/0.2789(嵌入)，有界工程收敛** | coarse 18.6–22.8 min；**medium 38.4 min(wake-free)/24.5 min(嵌入)** | `b7_onera_m6/`,`m6_medium_ls_workflow/` |
| 3D 亚声速网格收敛（圆帽 M6 M0.5） | **亚** | conf Newton | 三级自相似 | **p=2.31，cl→0.2050——全库唯一 3D Richardson** | 三级重跑（heavy） | `p13.../g133rt_richardson.csv` |
| 远场 A/B（15–120c） | **亚(M0)** | LS | dual-mesh | vortex 域鲁棒0.45%；neumann O(Γ/R) | B5 全重解 ~15 min | `b4p5_farfield/` |
| 球面 Cp（光滑曲壁） | **亚** | conf | medium | **11.6% vs 2% gate，未达标** | — | `p1_laplace/`(strict xfail) |
| 翼身组合体 | — | — | 仅网格 | **无任何求解**（B9 未开工） | — | `cases/meshes/onera_m6_wingbody/` |

**物理包络（design.md §2/§12）**：等熵无旋无粘；M∞ 0.3–0.87；激波须弱（局部法向
M≲1.3）；FP 非唯一性带 M0.82–0.85 是模型固有的（medium NACA M0.80 无孤立解已
实测兑现）。不适用：分离、强激波、钝体。

---

## 2. ★ 亚声速 / 跨声速算例总表 + 计算时间（本次重点）

### 2.1 测试文件的声速区分类（`tests/`，50 文件 / 337 原始函数 / 375+18+2 展开）

**亚声速测试（M≲0.5，无激波机器 / `upwind_c=0`）**
- conforming：`test_p2_kutta_naca0012.py`（不可压升力，NACA coarse/medium）、
  `test_p3_subsonic.py`（球面 M0.3；M0→Laplace 逐位）、`test_p3_naca0012_m05.py`
  （NACA M0.5 vs 面元）、`test_laplace_sphere.py`（球面 Cp，含 G1.6 xfail）、
  `test_m0_cylinder.py`/`test_wall_correction_cylinder.py`（不可压圆柱 Cp oracle，
  后者 1 xfail）。
- level-set：`test_b3_lifting.py`（M0.0 & M0.5）、`test_b4_te_control_volume.py`
  （M0.0）、`test_b45_farfield.py`（M0 远场 A/B）、`test_b8_tip_taper_ls.py`/
  `test_b8_span_blend.py`（M6 M0.5 翼尖）、`test_b11_post_unified.py`（M0.5 后处理
  两路径对拍）、`test_b11_linear_ls.py`（M0.5 Laplace/lifting precond）、
  `test_b13_lagged_picard.py`（**M0.5** lagged-LU Picard）。

**跨声速测试（M≥~0.7，激波/迎风/Mach ramp）**
- conforming：`test_p4_transonic.py`（M0.80 coarse smoke + medium gate + sweep，
  **全 gated**）、`test_p5_onera_m6.py`（M0.84 M6 gate，gated）、`test_p8_newton.py`
  （G8.1 M0.80/M0.7875 + **G8.2 M6 medium M0.84**，gated）、`test_p8_jacobian.py`
  （converged-pocket FD，1 gated）、`test_p13_tip_taper.py`（1 条 driver-reach
  M0.84）。
- level-set：`test_b6_transonic.py`（M0.70/0.80，5 条跨声速）、`test_b6_newton.py`
  （M0.65→0.70 seed+Newton；M0.80 gated）、`test_b7_onera_m6.py`（**M6 M0.84
  gate**，gated）、`test_b12_lagged_lu_ls.py`（**M∞=0.7** lagged-LU Newton）。

**既非亚也非跨（MMS/自由流/拓扑/单元/插桩）**：Track M 全部（`test_m0/m1/m2/m5`）、
`test_mesh_*`、`test_p6_*`、`test_p7_diff_flux.py`、`test_b1/b2`、`test_io_vtk`、
`test_post_surface` 等。

**gated 门槛**：8 个文件带 `PYFP3D_TRANSONIC_GATES=1`（`p4`,`p5`,`p8_jacobian`,
`p8_newton`,`b6_transonic`,`b6_newton`,`b7`,`b11_linear`）；几乎所有 M6/翼身测试
还会因 `.msh` gitignored 而 skip。这是 18 skipped 的主体。

### 2.2 计算时间总表（实测，全部 @16 线程；来源 = CSV 或脚本 docstring）

| 算例 | 声速区 | 网格 | 路径 | 时间@16线程 | 来源 |
|---|---|---|---|---|---|
| **全 pytest 套件** | — | — | 全 | **1068 s @8 线程**（前 882.98 s @16；G8.3 301.66 s） | CLAUDE.md 台账 |
| NACA 2.5D 升力 | **亚** | coarse | 两条 | **秒级** | `p3_subsonic/` |
| B3 LS 升力 demo | **亚** | coarse/medium | LS | ~2 min | `b3_levelset_lifting/run_demo.py:42` |
| B5 远场 A/B（4 域尺寸重解） | **亚** | NACA | LS | ~15 min | `b4p5_farfield/run_demo.py:44` |
| M6 亚声速 Newton | **亚** | coarse | conf | **7.6 s** | `g91_m6_levels.csv` |
| M6 亚声速 Newton | **亚** | fine(457k dof) | conf | **5294 s / 13 级**（M0.5，非升到 M0.84） | `g91_m6_levels.csv` |
| M6 LS Newton — spsolve | **亚** | medium(67k dof) | LS | **145.57 s**（7 步 7 refactor） | `b12_lagged_lu/results/m6_newton_ab.csv` |
| M6 LS Newton — **lagged-LU(B12)** | **亚** | medium | LS | **66.66 s = 2.18×**（1 refactor+30 GMRES，γ 逐位） | 同上 |
| M6 LS lifting Picard — spsolve | **亚** | medium | LS | **447.61 s**（26 外层，17.2 s/层） | `b13_lagged_picard/results/m6_lifting_ab.csv` |
| M6 LS lifting Picard — **lagged-LU(B13)** | **亚** | medium | LS | **68.31 s = 6.55×**（2 refactor，1 stall=安全网） | 同上 |
| M6 LS seed+Newton 流水线(B13) | **亚** | medium | LS | **111.91 s**（seed42+newton70）vs ~330 pre-B13 | `b13.../m6_end_to_end.csv` |
| **M6 工作流 亚·嵌入** | **亚(M0.5)** | medium | LS | **68.7 s**（cl0.2129） | `m6_medium_ls_workflow/summary.csv` |
| M6 工作流 亚·conforming | **亚(M0.5)** | medium | conf | 缓存(cl0.2126) | 同上 |
| **NACA 2.5D 跨声速 P4 medium G4.1** | **跨(M0.80)** | medium | conf Picard | **~16 min**（12,931 迭代） | `p4_transonic/run_demo.py:24` |
| P4 G4.3 鲁棒性扫掠 | **跨** | coarse | conf | **~22 min（10 例）** | `run_demo.py:329` |
| P4 heavy demo（4+5） | **跨** | medium+扫掠 | conf | **~40 min** | `run_demo.py:30` |
| **NACA 2.5D 跨声速 B6（全解）** | **跨(M0.70–0.80)** | coarse | LS | **~25 min 上限** | `b6_transonic/run_demo.py:55` |
| **M6 M0.84 conforming Newton(G8.2)** | **跨** | medium | conf | **250.8 s ≈ 4.2 min** | `g82_m6_medium.csv` |
| M6 M0.84 conforming Picard(P5 从头) | **跨** | medium | conf | **~75.6 min（4539 s）** | `g82_m6_medium.csv` `p5_picard_solve_s` |
| M6 M0.84 conforming fine(G13.2) | **跨** | fine(2.5M) | conf | **~40–45 min**（先前错 precond 跑到 1h16m 被杀） | `run_g132_transonic.py:41` |
| **B7 M6 coarse 嵌入(M1)** | **跨** | coarse | LS Picard | **22.8 min**（cl_KJ0.2765，0 lim/flr） | `b7_onera_m6/results/summary.csv` |
| B7 M6 coarse wake-free(M4) | **跨** | coarse | LS Picard | **18.6 min**（cl_KJ0.2710） | 同上 |
| **M6 工作流 跨·wake-free** | **跨(M0.84)** | medium | LS | **2304.7 s ≈ 38.4 min**（cl0.2765，M_max2.455，有界工程收敛） | `m6_medium_ls_workflow/summary.csv` |
| **M6 工作流 跨·嵌入** | **跨(M0.84)** | medium | LS | **1469.9 s ≈ 24.5 min**（cl0.2789，M_max2.195） | 同上 |
| M6 工作流 跨·conforming | **跨(M0.84)** | medium | conf | 缓存(P5 Picard 场，cl0.2499) | 同上 |

**M6 fine direct 陷阱（复发性成本，MEMORY 记录）**：`precond="direct"` 单次 splu
在 ~457k dof 上 **4h39m / 26 GB RSS 不返回**（`run_g132_transonic.py:90`）——这正是
B12/B13 lagged-LU 存在的理由，也是 M6 fine 必须用 `amg`+η=1e-8 的原因。

### 2.3 对照目标计算时间的即时判读

| 目标 | 目标口径 | 当前最快实测 | 判读 |
|---|---|---|---|
| **2D 亚声速** | 单核 <1 min | 秒级 @16 线程 | 表面达标，但**未在单核测过**；亚声速本就便宜，风险低 |
| **2D 跨声速** | 单核 <1 min | @16 线程 **16–40 min**（NACA medium）；coarse Newton 分钟级 | **⛔ 差 1–2 个数量级**，且这还是 16 线程；单核会更糟。这是最刺眼的差距 |
| **3D 亚声速** | 1~4 核 <10 min | @16 线程：conf coarse 7.6 s；LS medium lagged-LU **66–68 s**；conf fine 88 min | medium 若换算到 4 核约仍在 ~几分钟量级，**有希望但未在 4 核口径验证**；fine 远超 |
| **3D 跨声速** | 1~4 核 <10 min | @16 线程：conf Newton medium **4.2 min**；LS medium **24–38 min**；fine 不可达 | conf Newton medium 在 16 线程下 4.2 min——**换到 4 核几乎肯定破 10 min**；LS 破得更多；fine 无解 |

**结论**：唯一接近目标的是 **conforming Newton 的 M6 medium 跨声速（4.2 min@16
线程）**，但（a）它在目标 1~4 核口径下几乎必然超 10 min，（b）它是 conforming
路径（无任意几何入口），（c）它是"真解"但网格只到 medium。跨声速 2.5D（"2D"目标
的替身）在 16 线程下都要 10~40 min，**离单核 1 min 差两个数量级**。**当前没有任何
一个跨声速算例在目标核数口径下达标，也没有任何单核/4 核基准。** 若目标时间是硬
指标，这是最大的、被计时口径长期掩盖的差距。

---

## 3. 数值格式清单（更新标注）

- **空间**：Galerkin P1 tet（顶点中心中位对偶 FV）。**只有 tet**——"2D"实为单层
  挤出 2.5D（`meshgen/extrude.py`），无真 2D 核。
- **人工密度迎风**：默认多跳有向 walk（`kernels/upwind.py`）；streamline-Gaussian
  kernel（`mode="kernel"`）opt-in **孤儿**，从未进配方，Newton 侧直接 raise
  （结构上进不了 Newton）。
- **尾迹/Kutta 两个不同数学对象**：conforming = 节点复制 + 每站 secant 压力相等
  Kutta（斜率 b≈0.93 故必须 secant）；level-set = TE duplication + 切元辅助 DOF +
  B4 非线性 TE 压力相等 |q_u|²=|q_l|²（Γ 是解模态，无 secant，P5 st133 类失稳结构
  上不可能）。翼尖：conforming Γ(tip)=0 拓扑强制 + G13.2 紧支撑 taper（`vanish_smooth`,
  r_c=0.05b）；LS Γ→0 涌现（±3e-4）但片终止环奇异 **characterized-not-cured**（B8，
  诚实指数 +0.62/+0.37 = 与 conforming +0.52 同对象）。
- **求解器**：见 §4 矩阵。★ **相对 2026-07-14 的净变化**：LS 路径不再"硬编码
  spsolve"——B11 加 `precond=None|"ilu"|"amg"`（ILU 有效、AMG 在 wake_ls lifting
  上实测 STALL 只留给 Laplace）、B12 加 Newton lagged-LU 且在 M6 medium 3D 真跑、
  B13 加 Picard 外层 lagged-LU。
- **后处理**：conforming `post/surface.py`（`smooth_passes` 默认 0，只用于 Cp 曲线，
  用于力反使 V6 恶化）；LS `post/surface_ls.py`（D11 逐侧映射强制，否则 cl=−3.35）；
  两路径经 `post/unified.py`（B11）统一分派，逐位一致。★ `element_mach2` 默认已从
  `"side"` **翻转为 `"main"`**（`multivalued.py:624`，×5 度量伪影修复），B6/B7 M_max
  已重读（B7-M1 committed 1.453 本身是伪值单元，诚实 1.392）。

---

## 4. 两路径能力矩阵（post-B11/B12/B13，代码核对）

✅ 完整 ◐ 部分 ❌ 无。证据 = `file:line`（直接核对）。**粗体 = 相对 2026-07-14 的变化。**

| 能力 | conforming | level-set | 后果 |
|---|---|---|---|
| 升力 Kutta | ✅ 每站 secant `picard.py:367` | ✅ 隐式压力相等 `picard_ls.py:247` | 不同模型，B8 证不可互移 |
| 跨声速 Mach ramp | ✅ Picard+Newton `continuation.py:78`,`newton.py:806` | ✅ **Picard + Newton（B15 新增 `solve_multivalued_newton_transonic`）** | ★ B15 补齐（下方） |
| Newton | ✅ 全功能 | ✅ **全功能（B15：freeze-selection 落地，`freeze` 不再是空占位）** | ★ B15 补齐 |
| 预条件/规模 | ✅ amg/ilu/direct+lagged-LU `newton.py:329` | ◐ **`None`/ilu/amg（B11）+ Picard/Newton lagged-LU（B12/B13）** `picard_ls.py:274,279`,`newton_ls.py:108,114` | **medium 不再受 spsolve 墙限**；fine 仍未跑 |
| 远场选项 | ◐ 仅 Dirichlet+vortex(+Γ(z)) `picard.py:506,405` | ✅ vortex/neumann/freestream `picard_ls.py:263` | 互补空洞，方向相反 |
| Γ(z) 锥化涡远场 | ✅ | ❌（用 neumann 绕开） | |
| 翼尖 taper | ◐ **仅 Newton 入口** `newton.py:343`（Picard 无） | ✅ **Picard+Newton 均收** `picard_ls.py:273`,`newton_ls.py:107`（B8） | conforming Picard 3D 仍带未治愈翼尖边 |
| 表面 Cp 平滑 | ✅ opt-in `surface.py:318` | ✅ **有对应物** `surface_ls.py:36`+`unified.py:62`（B11） | 不再是空洞 |
| 3D 实跑 | ✅ Picard+Newton（medium） | ◐ **Picard(coarse+medium) + Newton(medium 亚声速,B12)**；跨声速仅 Picard | LS Newton 跨声速 3D 仍缺 |
| α 扫掠免重网格 | ❌（尾迹是网格实体） | ◐ `WakeLevelSet.update_direction` `levelset.py:100` 仅构造器+B1 测试调用，**未接入 solve 流程** | 卖点能力仍是孤岛 |
| wake-free 任意网格摄入 | ❌ | ✅（M3/M4/M2 家族） | LS 核心工作流优势 |
| 多尾迹/翼身求解 | ❌ | ❌ `CutElementMap`/`MultivaluedOperator` 各**只收单个 levelset** `cut_elements.py:105`,`multivalued.py:58` | B9 未开工 |
| 非升力 Newton 入口(wc=None) | ❌ `newton.py:306`（wc 必传） | ❌ | G10.1 开放 |

**默认值陷阱（同概念不同入口默认不同，易误用；代码核对）**：
- LS Picard `farfield="vortex"` (`picard_ls.py:263`) vs LS Newton `"neumann"`
  (`newton_ls.py:95`)——有 B6 依据但无警告。
- `solve_multivalued_lifting` 默认 `upwind_c=0.0`（`picard_ls.py:264`，**静默关闭
  跨声速**）vs conforming `solve_subsonic_lifting` 1.5（`picard.py:385`）vs LS Newton 1.5。
- `wing3d.py` `tip_cap="flat"`（`wing3d.py:215`，已证发散几何）vs `wingbody.py`
  `"round"`（`wingbody.py:116`）。
- conforming Newton `precond="amg"` vs LS Newton `precond=None`(=spsolve)。

---

## 5. 网格收敛：哪里成立、哪里不成立（无变化）

**成立**：① MMS（阶 1.96）；② 2D 亚声速升力（G9.2 2.71→0.33→0.03%，尖 TE 无升力
地板）；③ **3D 亚声速圆帽 M6 M0.5**（G13.3，p=2.31，cl→0.2050，全库唯一 3D
Richardson，花三个 phase 挣到：G13.2 taper + M1b h_far 钳位修复 + M5 圆帽）。

**不成立（明确失败）**：① **3D 跨声速 M0.84 两种几何双败**——平帽 fine 收敛但序列
0.2593→0.2652→0.2866 上升非渐近；圆帽 fine 的 Mach 续接死在 M=0.75，无 M0.84 fine
态 ⇒ **P9 判定带至今未触发**；② G9.1 原始 M6 序列 fine 非离散解 + 网格族钳位缺陷
双重作废；③ fold 区（NACA M≈0.78–0.80）dcl/dM≈6–10 结构上不可比；④ **LS 路径至今
无任何 fine 级收敛研究**（B13 让 medium 快了 6.55×，但 fine 仍未跑）。

---

## 6. 没实现 / 孤儿 / 阴性清单（更新）

### 6.1 完全无代码
| 项 | 设计程度 | 备注 |
|---|---|---|
| **Track V 粘性/边界层（V1–V4）** | 设计完整（IBL3 六方程、transpiration BC、GV gate） | `viscous/` 不存在，grep 0 命中。**第二支柱整体空白**。VII 预期使 CL 下降至 ~0.26–0.27，不能指望它解释 0.019 的离散精度差 |
| **B9 多尾迹/翼身 LS 求解** | gate 已定义，B9=NEXT | 仅网格；`CutElementMap`/`MultivaluedOperator` 只收单 levelset，无组合机制；无 `cases/demo/b9*`、无测试驱动翼身网格过升力解 |
| **钝尾缘** | 无 roadmap 条目 | grep `blunt` 0 命中；NACA/M6 均削尖。**用户目标中此项连排期都没有** |
| **真 2D 核** | 无 | 仅 2.5D 挤出替身；所有 kernel 只吃 tet |
| **LS Newton 的 Mach-ramp 包装 + freeze** | freeze 是空占位 `newton_ls.py:105` | LS 跨声速 3D 真解的前置 |
| **LS fine 级研究** | — | B11/B12/B13 解锁了 medium，fine 仍未跑（无工件） |
| G10.1 非升力 Newton / P11 曲壁元 / P12 backlog / B10 自由尾迹(搁置) | gate 级 | B10 搁置意味着刚性平面尾迹尖端物理**永不被模型级治愈**，只被 taper/圆帽工程压制 |

### 6.2 实现了但没人用（孤儿，caller 核对）
| 功能 | 位置 | 状态 |
|---|---|---|
| `ptc_dtau` | `newton.py:335` | 零外部调用 |
| `upwind_c_post` | `newton.py:816` | 零外部调用（M6 配方从未用） |
| `pseudo_dt`（retired） | `picard.py:400` | 无 caller 传值 |
| streamline-Gaussian kernel | `upwind.py mode="kernel"` | 仅测试+1 demo，无 solver caller，Newton 侧 raise |
| `newton_ls(freeze=)` | `newton_ls.py:105` | 空占位，函数体从不引用 |
| `WakeLevelSet.update_direction` | `levelset.py:100` | 仅构造器+B1 测试，未接 solve |
| `wall_correction.py` 整模块 | `solve/` | oracle 否决，无复用者 |
| `tanh_half`/`span_blend`/`te_weld_coo` | `wake.py:184`,`multivalued.py:58,334` | 阴性标本；`te_weld_coo` 已被 tip_taper blend **接线（非孤儿）** |
| `smooth_wall_tangential_gradients` | `surface.py:147` | opt-in，默认 0（LS 侧 B11 已有对应物） |

### 6.3 有工件的阴性 / 开放
1. **G1.6 球面 Cp 11.6%**（strict xfail）：光滑曲壁变分犯罪；h 加密饱和 ~3.6%；唯一
   活路 = 曲壁元 + gate 重定义。**含义：对"光滑曲面为主"几何（机身主导）壁面 Cp
   有 ~10% 级系统误差风险**——直接关系翼身推进；B9 gate 已加机身 Cp 护栏条款。
2. **3D 跨声速网格收敛双败**（§5）+ 圆帽 fine tip-TE 超限（已定位 20/20 最快单元在
   尖 tip TE）。
3. **B8 LS 片终止奇异** characterized-not-cured（用户裁决关闭；gate 内 mechanism
   probe 一项 NOT MET，标"kept by design"）。
4. **B6 medium 折叠区**：LS Picard 停 −18.8%；LS Newton 折叠解升力 −13% 未分摊。
5. **B7 只有 coarse**；顶部 Mach 级 |R|~4-6e-6 有界非收敛（gate 只断言 bounded）。
   **M6 工作流的 medium 跨声速同样是"有界工程收敛"**（残差 plateau 1e-5..1e-4，
   非严格收敛；跨声速方法 A/B 差 **10.65%**，靠 <15% 的松判据通过）。
6. **P4/P5 Picard 收敛态非离散解**（勘误在案，引用须带"Picard 质量"限定）。
7. FP 折叠（M0.80 medium 无孤立解）是模型极限非 bug。

---

## 7. 文档时效性与相容性核查（"检查文档是否最新且相容"的答复）

**结论：roadmap.md / agent-rules.md / demo_report.md 已跟上 B11/B12/B13；
STATUS.md 与 CLAUDE.md 陈旧；roadmap 缺 2026-07-15 工作流条目；本审查前身
`capability_review_2026-07-14.md` 有 ~8 处已被覆写。**

### (A) 各文档时效性
| 文档 | B11 | B12 | B13 | M6 工作流(07-15) | 判定 |
|---|---|---|---|---|---|
| roadmap.md | ✓ | ✓ | ✓ | **✗ 缺** | 基本最新，**唯缺 2026-07-15 工作流条目**（仅在 agent-rules+demo_report） |
| agent-rules.md | ✓ | ✓ | ✓ | ✓（line 3） | **完全最新** |
| demo_report.md | ✓ | ✓ | ✓ | ✓ | **完全最新** |
| STATUS.md | ✗ | ✗ | ✗ | ✗ | **陈旧** |
| CLAUDE.md | 部分 | ✗ | ✗ | ✗ | **陈旧** |

- **STATUS.md**（快照标 2026-07-14）：line 11 `当前工作重心 = P13`（应为 B9=NEXT）；
  line 68 同一格内 `下一步=B8` 与 `=B9` 并存、无 B11/B12/B13；line 124 基线
  `358+17+2`（陈旧）。**建议整体刷新或明确标注 superseded。**
- **CLAUDE.md**：line 11 文档地图把 Track B 记为"B1–B5+B7 closed / B6 open"
  （漏 B8/B11/B12/B13）；line 58 基线 `375+18+2` 显式只加了 B11 的 +9+8，
  **未含 B12(+4)/B13(+5)**，故基线相对 B12/B13 陈旧。**建议更新 Track-B 状态行
  + 基线纳入 B12/B13。**
- roadmap 应补 `m6_medium_ls_workflow`（2026-07-15）条目——目前只在 agent-rules
  与 demo_report 有。

### (B) 跨文档矛盾
1. 当前阶段：STATUS(P13) vs agent-rules/roadmap(B9=NEXT)。
2. Track B "下一步"：STATUS 同格 B8/B9 并存；CLAUDE.md 暗示 B6 是前沿。权威答案：
   B8 closed，B9 next。
3. 基线三方不一致：STATUS `358+17+2` vs CLAUDE `375+18+2` vs
   roadmap/demo_report/agent-rules（含 B12/B13 增量）。仅后三者互洽至 B13。
4. **`capability_review_2026-07-14.md` 内部自相矛盾**：§2.5(line 69) 说
   `element_mach2` 默认还是 `"side"`（未翻转），而 §5.2(line 142)/§7(line 179/183)
   说 2026-07-14 已翻转为 `"main"`。代码为 `"main"`（`multivalued.py:624`），故
   line 69 内部错误。
5. "0.019 gap = resolution / P11 refuted" 降级 **已核实一致**：唯一残留 `refuted`
   在 CLAUDE.md:27 是对被审计出的 prose-only 断言的**历史叙述**，非活动断言；所有
   活动语句已降为"strongly indicated, NOT earned"。**残留**：agent-rules line 17
   与 STATUS line 65 仍有加粗小标题"THE 0.019 GAP IS RESOLUTION"紧接降级正文，
   标题措辞未软化（宜一并软化以免误读）。

### (C) 文档 vs 代码抽查（全部 file:line 核对通过）
| 断言 | 结果 | 位置 |
|---|---|---|
| `solve_multivalued_{laplace,lifting,newton}` 有 `precond=None\|ilu\|amg` | ✓ | `picard_ls.py:110,274,376`,`newton_ls.py:108,156` |
| B12 `direct_refactor_every` in newton_ls | ✓（默认1=逐位spsolve，rtol1e-8） | `newton_ls.py:114,303-332` |
| B13 `direct_refactor_every` in lifting | ✓（默认1，rtol **1e-10**≠B12的1e-8，刻意） | `picard_ls.py:279,514-531` |
| `newton_ls freeze=` 仍空占位 | ✓（签名+docstring 声明，函数体从不引用） | `newton_ls.py:105,128` |
| `element_mach2` 默认 | ✓ `="main"` | `multivalued.py:624` |
| `WakeLevelSet.update_direction` 仅 B1 测试调 | ◐ **不精确**：构造器 `__init__` 也调；"未接 solve 流程"成立，"仅 B1 测试调"字面不准 | `levelset.py:94,100` |

### (D) `capability_review_2026-07-14.md` 已陈旧的断言（撰写后一天被 B11/B12/B13 覆写）
1. §0/§2.3：LS"硬编码 sparse-direct/无 AMG 逃生口/LS Newton 无预条件/无 lagged-LU/
   从未跑过 3D"——**4 点陈旧**（B11 precond、B12 lagged-LU、B12 M6 medium 3D 真跑
   LS Newton）。仅"无 Mach-ramp 包装""freeze 未实现"仍真。
2. §0：LS"3D 只在 coarse 用 Picard 验证过一次(B7)"——**陈旧**（2026-07-15 M6
   medium 亚+跨；B12 M6 medium Newton）。
3. §1 表：LS 跨声速"仅 coarse"——**陈旧**（M6 medium M0.84 已示范）。
4. §3 矩阵：LS Newton splu/无 3D/无 lagged-LU；"预条件 ❌ 硬编码 spsolve，fine 不可
   达"——**precond/3D/lagged-LU 三点陈旧**；仍真：无 ramp、freeze 未实现、无 fine。
5. §5.1：LS GMRES+AMG "B3+ 递延"——**已由 B11 交付**。
6. §5.1：LS Newton freeze/lagged-LU/ramp/3D——**lagged-LU+3D 已做**，仅剩 freeze+ramp。
7. §5.2：`newton_ls(freeze=)` 断言仍对，**行号陈旧**（:87→:105）。
8. §4：LS"受 spsolve 墙限制"——"无 fine 研究"仍真，但"spsolve 墙"框架已被 B13 的
   6.55× 覆写（medium 不再 spsolve-bound）。

（前身 §7 item 1–6 标注的已修复项——P1 护栏、B8 NOT-MET、G13.2/3 降级、LS 指数勘误、
M2 census、element_mach2 翻转——本次复核仍成立，非陈旧。）

### 本次审查未做的编辑
以上均为**发现报告**，未改动任何文档（考虑到 phase-scoped 提交纪律 + 可能有并发
会话在改同一 repo）。**可立即行动的高优先级修复**（如需，可代为执行）：① 刷新
STATUS.md（重心行/Track B 行/基线）；② CLAUDE.md line 11 Track-B 状态 + 基线纳入
B12/B13；③ roadmap 补 2026-07-15 M6 工作流条目；④ 为前身 review 加 dated addendum
指向本文；⑤ 软化两处"0.019 GAP IS RESOLUTION"加粗标题。

---

## 8. 对照最终目标的差距表（含计算时间维度）

目标：非结构网格通用全速势 **+ 边界层积分** 求解器；2D/2.5D/3D 任意几何（翼型、
机翼、翼身、多段翼型、钝尾缘）；亚/跨声速、有/无升力；**快**（2D 单核<1min，
3D 1~4 核<10min）。

| 目标要素 | 现状 | 差距 |
|---|---|---|
| 非结构全速势核 | conforming 3D Newton 成熟（真解级） | ✅ 基本到位 |
| **边界层修正** | 设计完整、**零代码** | ⛔ 整支柱空白 |
| 2D | 无真 2D，仅 2.5D 替身 | ⚠ 若 2.5D 可接受则无差距 |
| 2.5D | 两路径成熟 | ✅ |
| 3D | conf 全能力；LS 亚 medium(Newton) + 跨 medium(Picard) | ◐ LS 差跨声速 Newton + fine |
| 任意几何 | 生成器仅 NACA/M6/旋成体族；翼身求解(B9)未开工；G1.6 光滑壁未解 | ◐→⛔ |
| level-set 尾迹 | 亚 medium 成熟；跨 medium(有界)；尖端 characterized-not-cured；无 fine | ◐ |
| **网格收敛** | 2D 亚 ✓、3D 亚圆帽 ✓；**3D 跨双几何败，无 M0.84 外推**；LS 无 fine | ◐ 完成一半 |
| 尖尾缘 | ✓（2D 无地板已证；3D 尖 tip TE 是跨声速 fine 拦路虎） | ✅/⚠ |
| **钝尾缘** | 无代码无排期 | ⛔ |
| **多尾缘/多尾迹** | 无代码（单 levelset 硬约束） | ⛔ |
| **计算时间（快）** | 全部 @16 线程；跨声速 2.5D 16–40min，3D 跨 medium 4.2min(conf Newton)~38min(LS)，fine 不可达；**无任何单核/4 核基准** | ⛔ **对照目标核数口径全面不达标且未测** |

**优先级（按目标倒推）**：
1. **B9（翼身 LS 求解 + 多尾迹机制）** + **LS Newton 跨声速（ramp+freeze）** 是打通
   "任意几何 × level-set × 跨声速真解"的咽喉；两者的共同前置（LS fine 级 GMRES/AMG）
   B11 已交付基础，但 fine 从未跑过。
2. **Track V 从零起步**——第二支柱，用户目标的一半。
3. **计算时间必须在目标核数口径（单核/4 核）下重新基准化**——当前 16 线程数字对
   "快"这个目标没有说服力，且掩盖了跨声速 2.5D 差两个数量级的事实。若"快"是硬指标，
   现有 Picard/Newton 在跨声速下的迭代成本（NACA medium 上万次 Picard 迭代 / M6
   Newton splu 填充）需要专门的性能 phase，而非附带产物。
4. 钝尾缘 / 多段翼型目前在权威文档之外，若确属目标需先立项。

---

*审计方法与三路原始报告（测试+时间、文档时效、代码矩阵）结论均落到 file:line 或
committed CSV。本文为其批判性综合，是 `capability_review_2026-07-14.md` 的
2026-07-15 滚动更新，非 gate 文档。*
