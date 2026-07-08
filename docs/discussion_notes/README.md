# Discussion Notes 索引

> 本文件管理 `discussion_notes/` 下所有讨论文档的状态。
> 每次新增或更新文档时，同步更新此索引。

## 文档列表

| # | 文件名 | 日期 | 主题 | 状态 | 备注 |
|---|--------|------|------|------|------|
| 1 | `20260707_1505_levelset_wake_design.md` | 2026-07-07 15:05（更新 22:45） | Level-Set 尾流 + 多值有限元实现方案 | ✅ 有效 | 方案 B 选定；6 阶段 B1-B6；已更新 P4 修复后的接口和性能基线 |
| 2 | `20260707_2118_ibl_viscous_coupling_design.md` | 2026-07-07 21:18（更新 22:45） | IBL 粘性修正耦合方案（VII） | ✅ 有效 | 方案 VI-1 选定；4 阶段 V1-V4；已确认与 P4 修复完全兼容 |
| 3 | `20260707_2313_newton_solver_strategy.md` | 2026-07-07 23:13（更新 07-08 01:05） | Newton 求解器实现策略 | ❌ 弃用 | 策略 B + N-Γ 分裂路线基于错误前提（López 实际用策略 A）；被文档 4 替代 |
| 4 | `20260708_1646_newton_solver_strategy.md` | 2026-07-08 16:46 | Newton 求解器实现策略（修订版） | ✅ 有效（核心成立，3 处已修正） | 策略 A（完整 Jacobian）+ 全耦合 Newton + P6 硬前置：核心架构经 PDF 核对全部成立；但尖-TE 前提、load-stepping 调度表、可微性=可微选择 3 处被 PDF 修正，见 §"文档 4 的 PDF 交叉核对"，修正已落实到 design.md §3.1/§8.1 + roadmap P6/P7 |

## 状态定义

- ✅ **有效**：当前适用，未被否定或替代
- ⚠️ **部分有效**：核心思路有效，但部分内容需要更新以匹配代码变化
- ❌ **弃用**：被后续文档或代码变更否定，不再适用
- 📋 **草案**：初步构思，尚未确定方案

## 文档关系

```
20260707_1505_levelset_wake_design.md（方案 B：Level-Set 尾流）
    │
    ├── B1: 切单元识别
    ├── B2: 多值 FE 装配
    ├── B3: 罚函数 Kutta  ←──── 20260707_2118 依赖此阶段
    ├── B4: 跨声速          ←──── 与 20260708_1646 的 Newton 路径互补
    ├── B5: 多尾流
    └── B6: 曲线尾流 / 自由尾迹
            │
            └── V4: 尾流面 IBL 修正 ← 依赖 B3 + V1-V3

20260707_2118_ibl_viscous_coupling_design.md（方案 VI-1：Transpiration IBL）
    │
    ├── V1: IBL 求解器（独立，可与 B/N 并行）
    ├── V2: 壁面 Transpiration BC（依赖 V1）
    ├── V3: VII 迭代（依赖 V2）
    └── V4: 尾流面 IBL 修正（依赖 V3 + B3）

20260708_1646_newton_solver_strategy.md（Newton 求解器，对应 P6→P7）
    │
    ├── 前置：P5 TE 奇异性缓解 → P6 可微通量（**硬依赖**——策略 A 需要）
    ├── N0: 密度截断光滑化
    ├── N1: P6 可微通量（switching function + 上游选择）
    ├── N2: Newton Jacobian 装配（策略 A：项 1+2+3，stencil 宽一层）
    ├── N3: GMRES + AMG 线性求解
    ├── N4: 全耦合 Newton 驱动器 + 亚声速验证（Γ 作为未知数，无 secant）
    ├── N5: 跨声速验证 + load stepping 参数调度 ←── 与 B4 互补
    └── N6: ONERA M6 + 性能 gate

20260707_2313_newton_solver_strategy.md（❌ 已弃用，被 20260708_1646 替代）
```

## 与当前代码的适配度（2026-07-08 17:30 快照，HEAD: 0b0749e）

### 仓库变化（a3f75ef → 0b0749e）——P5 **已关闭**（medium PASS）

P5 ONERA M6 验证阶段**关闭**（2026-07-08）：
- **coarse (55.5k)**：PASS，M_max 1.398，0 floored/limited，CL 0.2419。
- **medium (350.7k)**：**现 PASS**，M_max 1.995，**0 floored/limited**，CL 0.2453，Kutta |F| 5.8e-4。
- **根因二次翻案（T1–T4，`INVESTIGATION_kutta_closure.md`）**：medium 失败**不是** TE 离散奇异性——而是**单站 Kutta 闭合失败**（st133，z/b=0.801，欠环量 32%；T3 决定性：只把该站 Γ 设为其 Kutta 目标即把 18→0 cell 团簇消掉，1/h 奇异性不可能依赖 Γ）+ 远场 2D 点涡分支射线伪影（翼尖外）。之前笔记里"TE P1 离散奇异性"的归因**已被推翻**。
- **修复（recipe 级，默认关闭，非 3D 路径 bit-identical）**：`farfield_spanwise_gamma=True`（Γ(z) 渐缩远场涡）+ `n_kutta_polish=4`（continuation 后定 Γ Kutta 闭合抛光，`omega_rho_polish=0.5`，secant-free 收缩）。
- **V6<1% 重新定义为 post-P6 目标**（O(h) 尖 TE/LE P1 + sawtooth 离散地板；G5.2 re-spec）。
- `continuation.py`：forward rtol 1e-10→1e-7，5.5× 加速，M_max 一致到 5 位（已在 a3f75ef）。

### 文档 4 的 PDF 交叉核对（2026-07-08，直接读 `references/Dissertation_Inigo_Lopez.pdf`）

在把 P6/P7 写入 design.md/roadmap.md 前，直接对博士论文正文+附录做了核对。**核心架构判断全部成立**（策略 A / 全耦合 Newton / 上游耦合 B.4 非零 / 开关函数导数纳入 / 上游选择冻结），但发现文档 4 有 **3 处需修正**（已在 design.md §3.1/§8.1、roadmap P6/P7 落实）：

| # | 文档 4 声称 | PDF 实际（页码/公式） | 修正 |
|---|-----------|---------------------|------|
| 1 | López 用**钝 TE**（Eq 4.2），故 UP3D 尖 TE 是 Newton 障碍，需 N0 光滑化+可能钝化 | 正文 p63：Eq 4.2 **把钝 TE 改尖**（NACA0012 原式 4.1 是钝的），López 用**尖 TE** 仍达二次收敛 | 尖 2D TE **不是** Newton 障碍；N0 光滑截断降为**可选**（收敛态 0 floored/limited，硬截断在收敛态非激活）；只有 3D 翼尖 TE 几何奇异性需截断（Fig 3.5, §3.4） |
| 2 | §9.2 调度表：每步随 M∞ 调 M_crit(0.99→0.90) 和 μc(↑)，"参考 Table 4.7-4.8" | Table 4.7/4.8：单个 case 内 M_crit/μc **固定**、只 ramp M∞；Table 4.13 M6：M_crit **恒 0.95**、μc 先恒 2.0 再在定 M∞ 后**降** 2.0→1.6 | 无 per-step M_crit 扫掠；μc 调度是**到达目标 M∞ 后的降耗散**，非 ramp 中增 |
| 3 | Defect B（可微性）需要"可微通量"含可微上游选择 | Appendix B：推导**冻结** u(e)，只对 ρ_up/μ/ρ 求导（B.3–B.8） | P6 只需在**固定选择下**可微，无需可微化整数选择——这也把 sawtooth（选择跳变的空间伪影）与 Newton 可微性**分离**为两个独立缺陷 |

### 文档 1：Level-Set 尾流方案（`20260707_1505_levelset_wake_design.md`）

**整体评估：✅ 有效——P5 进展不改变方案 B 的设计，但 B4 阶段需关注 TE 奇异性。**

| 检查项 | 状态 | 说明 |
|--------|------|------|
| 前置假设：当前尾流用 conforming 面 + master-slave | ✅ 仍然成立 | `wake_cut.py` + `constraints/wake.py` 仍是此架构 |
| 方案 A/B/C 对比逻辑 | ✅ 有效 | 三方案对比不受 P5 进展影响 |
| Davari 2019 + Núñez 2022 理论引用 | ✅ 有效 | 公式编号引用已对照原文核实 |
| 模块设计（新增 `wake/` 目录） | ✅ 有效 | 当前代码结构未变 |
| B1-B6 阶段划分 | ✅ 有效 | 阶段依赖关系不受 P5 进展影响 |
| B4 跨声速阶段 | ⚠️ 需关注 | P5 发现 TE P1 奇异性在 fine mesh 恶化；B4 实现时需同步处理 TE 单元（与 P6/P7 修复路线协同） |
| `constraints/wake.py` 的 h_j 批处理 | ✅ 已落地 | commit e0353d9 中完成 |

### 文档 2：IBL 粘性修正耦合方案（`20260707_2118_ibl_viscous_coupling_design.md`）

**整体评估：✅ 有效——与当前代码结构完全兼容。**

| 检查项 | 状态 | 说明 |
|--------|------|------|
| 文献澄清（Davari 2019 vs Núñez 2022） | ✅ 有效 | UP3D 仍是 body-fitted 几何 |
| 方案 VI-1（Transpiration BC）选型 | ✅ 有效 | 不依赖尾流表示方式 |
| 壁面 Transpiration 实现路径 | ✅ 有效 | `kernels/residual.py` 壁面 Neumann 项增加源项 |
| 尾流面 IBL 修正 | ✅ 有效（条件性） | 依赖方案 B 的 B3 阶段 |
| VII 迭代循环设计 | ✅ 有效 | 与势流求解器接口清晰 |
| V1-V4 阶段划分 | ✅ 有效 | V1 可独立开发，不阻塞 |
| P5 TE 奇异性影响 | ✅ 无直接影响 | IBL 的 δ\* 修正不涉及 TE 单元梯度；但 VII 验证需等 P5 gate 通过 |

### 文档 3：Newton 求解器策略（`20260707_2313_newton_solver_strategy.md`）

**整体评估：❌ 弃用——基于 López 博士论文 Appendix B 交叉验证，核心判断有误。被文档 4 替代。**

| 检查项 | 状态 | 说明 |
|--------|------|------|
| Picard vs Newton 算法对比 | ✅ 有效 | 当前 Picard 实现与文档描述一致 |
| N-Γ 分裂策略 | ❌ 弃用 | 保留 secant 外层，H3 耦合失稳风险未消除；改为全耦合 Newton |
| design.md (6.3) Newton Jacobian 公式 | ❌ 弃用 | 策略 B 为错误推荐；López 实际用策略 A（Appendix B 证明） |
| P6 前置依赖 | ❌ 弃用 | P6 是 P7 的硬前置依赖（策略 A 需要可微通量） |
| Jacobian 策略选择（A/B/C） | ❌ 弃用 | 推荐改为策略 A；策略 B 降为 fallback |
| GMRES + AMG 线性求解器选择 | ✅ 有效 | design.md §7 已规划 |
| `damping_theta` 作为全局化策略 | ✅ 有效 | P4 已落地 |
| Eisenstat-Walker inexact Newton | ✅ 有效 | 当前 `forcing=0.0`，Newton 阶段需打开 |
| 性能预估 | ❌ 弃用 | 过于乐观（未考虑 load steps）；需修正 |
| **P5 TE 奇异性修复** | ⚠️ 仍需处理 | 文档 4 新增 N0 密度截断光滑化阶段 |

### 文档 4：Newton 求解器策略修订版（`20260708_1646_newton_solver_strategy.md`）

**整体评估：✅ 有效——基于 López 2021 博士论文 Appendix B 交叉验证重写。**

| 检查项 | 状态 | 说明 |
|--------|------|------|
| 策略 A（完整 Jacobian 含上游耦合） | ✅ 正确 | López Appendix B Eq.(B.3)-(B.6) 证明上游耦合存在 |
| 全耦合 Newton（Γ 作为未知数） | ✅ 正确 | 消除 secant，根除 H3 耦合失稳 |
| P6 硬前置依赖 | ✅ 正确 | 策略 A 需要可微通量 + 上游选择 |
| TE 奇异性光滑化（N0 阶段） | ✅ 新增 | 密度截断 C¹ 光滑化，Newton 可微性前提 |
| Load stepping 参数调度 | ✅ 新增 | Mcrit/μc 随 Mach 步调整（López Table 4.7-4.8） |
| 尾流 Jacobian 对比（master-slave vs least-squares） | ✅ 新增 | UP3D 方式更简单，支持全耦合 Newton |
| 性能预估（含 load steps） | ✅ 修正 | ONERA M6: 12 steps × ~7 Newton/step ≈ 84 迭代 |
| 收敛阶预期 | ✅ 修正 | 策略 A → 严格二次；López Table 4.9 验证 |

## 更新历史

| 日期 | 操作 | 说明 |
|------|------|------|
| 2026-07-07 22:40 | 创建索引 | 初始版本，含 2 个文档，适配度分析基于 HEAD e0353d9 |
| 2026-07-07 22:45 | 更新文档 1 | P4 修复后适配：§2.4 注释、§8 B4 性能基线、§10 接口继承；状态从 ⚠️ 部分有效 → ✅ 有效 |
| 2026-07-07 22:45 | 更新文档 2 | 加附录确认与 P4 修复完全兼容，无需修改 |
| 2026-07-07 22:45 | 更新索引 | 同步状态变化 |
| 2026-07-07 23:13 | 新增文档 3 | Newton 求解器实现策略；5 阶段 N1-N5；更新依赖关系图 |
| 2026-07-07 23:18 | 更新文档 3 | 阶段重编号 P6→P7；新增 P6 可微通量前置依赖说明；gate 编号 G6→G7 |
| 2026-07-07 23:18 | 更新索引 | 同步阶段编号变化和依赖关系图 |
| 2026-07-07 23:22 | 同步 baseline | 文档 1、2 的 baseline 从 e0353d9 更新到 e3d0386；内容无变化 |
| 2026-07-07 23:38 | 补全文档 3 适配度 | README 索引补全文档 3 的适配度分析小节（13 项检查） |
| 2026-07-07 23:42 | LaTeX 格式化 | 三个文档的数学符号统一为 LaTeX `$...$` / `$$...$$` 格式；代码块内不动 |
| 2026-07-07 23:50 | 公式代码块修正 | 将 32 个误写为代码块的公式转为 `$$...$$` LaTeX 格式（文档 1: 21, 文档 2: 9, 文档 3: 2） |
| 2026-07-08 00:30 | LaTeX 格式全面修正 | 用户在独立对话中产出了格式规范的三个文档，已替换工作区文件；表格对齐、公式编号、符号一致性全面修正 |
| 2026-07-08 01:05 | 更新文档 3 | P6 前置依赖修正（非硬依赖）；新增 Jacobian 三策略分析（A/B/C）；策略 B 为推荐路线；stencil 不变；性能预估修正 |
| 2026-07-08 16:10 | 同步 baseline | 仓库更新到 a3f75ef；P5 coarse PASS / medium FAILS；TE P1 奇异性根因确认；spanwise-Γ smoothing 反驳；B4 ⚠️ TE 奇异性；P6/P7 优先级提升 |
| 2026-07-08 16:46 | 弃用文档 3 | López 博士论文 Appendix B 交叉验证：策略 B 论据错误（López 用策略 A）；N-Γ 分裂保留 secant 风险；P6 是硬前置依赖 |
| 2026-07-08 16:46 | 新增文档 4 | Newton 求解器策略修订版；策略 A + 全耦合 Newton；P6 硬前置；N0 密度截断光滑化；load stepping 参数调度；性能预估修正 |
| 2026-07-08 17:30 | 同步 baseline + PDF 核对 | 仓库更新到 0b0749e（P5 **已关闭**，medium PASS，TE-奇异性归因翻案为单站 Kutta 闭合失败）；直接读 López 博士论文 PDF 核对文档 4：核心（策略 A/全耦合/B.4 上游耦合/开关函数导数/上游选择冻结）全部成立，发现 3 处需修正（尖 TE、调度表、可微性≠可微选择）；已将 P6/P7 设计落实到 design.md §3.1–3.2/§6.3/§8.1/§13 + roadmap P6/P7 阶段与 ledger |
