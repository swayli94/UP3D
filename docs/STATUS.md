# pyFP3D 开发现状（快照）

> **快照日期：2026-07-12。** 本文件是给人读的高层概览，**不是**权威进度源：
> 阶段/gate 状态以 [roadmap.md](roadmap.md)（含进度台账）为准，当前阶段细节以
> [agent-rules.md](agent-rules.md) 的 "Current phase" 行为准，证据在
> [demo_report.md](demo_report.md)。若本文件与它们冲突，以它们为准。
> 每次 gate 关闭后按 CLAUDE.md 工作流更新 roadmap/agent-rules 时，顺手刷新本文件。

## 一句话状态

求解器主线（Track P）P0–P9 全部关闭（P10 的 G10.2/G10.3 亦关）。**P9 判别阶段给出
了否定性但决定性的结论：3D 升力缺口不是壁面单元的问题，而是刚性平面尾迹在翼尖自由
边上的涡片边缘奇点**（★2026-07-12 更正：这是**尾迹模型**缺陷，其修复是尾迹卷起/显式翼尖涡，**不是 Track B**——Track B 只换尾迹表示不换模型；`cases/demo/tip_edge_singularity/` 在 M0.5 无限幅器下实测该奇点在 conforming 与 level-set 两条路径上都随加密发散，同网格 coarse→medium 峰值 ×1.38/×2.28，机翼对照平）；P11 曲面壁元作为"3D 升力修复"不再被
支持（其 G1.6 球面 Cp 理由仍然成立，待用户仲裁）。**当前工作重心 = Track B**：
B1–B5 + **B7** 已关（level-set 尾迹 + 隐式 Kutta 已产生升力，无 Γ secant）；
**★ B7（ONERA M6 三维 gate）2026-07-12 一次通过**——双网格达标，且 B6 的**升力反转在
三维复现**：以 conforming **Newton 真解**（cl_KJ 0.2692）为基准，LS Picard 落在
**+2.7%(M1) / +0.7%(M4 无尾迹工作流网格)**，而 conforming Picard（P5, 0.24788）低了
**−8.6%**；三维专属机械（斜交标架、展向裁剪 ⇒ Γ(tip)→3e-4 离散归零）**无需新求解器
代码**，唯一缺口是后处理（`section_cp_curve_levelset` / `cl_pressure_3d_levelset`）；
**三维远场 = `neumann`，P5 的 Γ(z) taper 在 B 路径上结构性地不需要**（B 路径的涡是
展向均匀、支割线恒为 y=0 的二维涡，两个独立错法都已实测；neumann 不带涡故二者皆
不可能）。诚实缺口：顶部 Mach 级停在 Picard 残差尾（|R| 4–6e-6，有界物理、指标在带内，
gate 断言"有界"而非 `converged`）；**LS Newton 在 M6 上推迟**（朴素 splu，需 lagged-LU）。
下一步 = B8（多尾迹）。**B6（跨声速）仍进行中**——
coarse M0.80 gate 已达成（基线改为同网格 Newton 真解）、**LS Newton 已交付并 FD 验证、
在折叠区取得二次收敛真解**（工作流网格 M3 medium |R| 1.5e-12），medium 定量收尾余两项
（用户决定暂不追）。

## 阶段一览

| 阶段 | 状态 | 要点 |
|------|------|------|
| P0–P3 | ✓ | 基础设施、Laplace、尾迹割缝/Kutta、亚声速压缩（P1 仅 G1.6 未关，见下） |
| P4 | ✓（带勘误） | 人工密度跨声速 Picard；★2026-07-11 勘误：Picard 态**不是**离散方程的解（Newton 残差 2.2e-4），P4 gate 降级为 Picard 质量/鲁棒性 gate（warm-start 引擎） |
| P5 | ✓（带注记） | ONERA M6 M0.84/α3.06；同类欠收敛注记（真解 cl 高 +5.8%/+7.9%，激波位置不变） |
| P6 | ✓ | 表面 Cp 锯齿 = 壁面梯度**恢复**伪迹（非通量），法向门控平滑修复 |
| P7 | ✓ | 冻结选择下的 ∂ρ̃/∂φ 精确灵敏度，FD 验证 3e-10~6e-9 |
| P8 | ✓ | 全耦合 (φ_red, Γ) Newton：精确 Jacobian + 直接步 + 停滞自适应冻结；G8.2 M6 medium 249 s（G10.2 后 ~145 s）、G8.3 全套件 302 s |
| P9 | ✓（2026-07-11） | 判别阶段。**G9.2 干净 PASS：2D 尖 TE 升力误差 2.71%→0.33%→0.03% ⇒ 尖 TE 不构成升力地板**（"尖 TE/LE 平壁"归因的 2D 支腿没了）。**G9.1 INVALID：M6 fine 不是离散解**（无限幅 M_max 1.40→2.13→**7.93**，9 个单元触 M_cap=3 ⇒ 无三点 Richardson）。★**发散奇点位于尾迹片的自由翼尖边**（9 个封顶单元全在 z/b 0.998–1.000、TE 之后、弦平面内）= 刚性平面尾迹的经典涡片边缘奇点 ⇒ **尾迹模型缺陷（修复=卷起/翼尖涡，非 Track B——见 tip_edge_singularity demo），不是壁元缺陷**。G9.3：0.019 缺口按原提法**不可拆分**（份额如实记 n/a）。顺带测得：`precond="amg"` + 紧 EW forcing（η=1e-8）在每个规模都更快且有效 |
| P10 | ◐ | G10.2 ✓（分裂裁决：M6 +41% 提速已晋升；折叠区禁用松容差）；G10.3 ✓（裁决：**保留 Mach ramp**）；G10.1（非升力 Newton 入口）未动、无顺序约束 |
| P11 | 未开（条件性） | 曲面/等参壁元。**G9.3 裁决：作为"3D 升力修复"不被 P9 支持**（曲面壁元去不掉尾迹片边缘奇点）；仅剩 G1.6（光滑曲壁上的球面 Cp，另一机理）为有效理由 — **待用户仲裁** |
| P12 | 未开 | Backlog：离散伴随、VII transpiration、混合单元/BO 标定 |
| P13 | ◐（G13.1 ✓ 2026-07-13） | **翼尖/尾迹边缘奇点 — 定性 + 尾迹模型收尾**（新阶段，P9 直系后裔；追加不重编号）。**G13.1 CLOSED**（demo `tip_edge_singularity` 10/10，亚声速 M0.5 无跨声速机械）：★奇点是 **1/√r 平板边缘型，非"1/r"**（conforming 三点 coarse/medium/fine 峰值 Mach 对数斜率 **p=0.59∈[0.4,0.65]**，峰值 0.712→0.981→1.510；1/r 线涡应 p=1）；★驱动量 = **脱落涡强 dΓ/dz（翼尖最大 ~10× 中展），非束缚 Γ**（Γ 翼尖→0 是必要非充分）；★**模型非表示**（conforming ×1.38 与 level-set ×2.28 同网格皆发散，LS 斜率 1.34≥conforming 0.52；同盒 p95/mean 平；峰值在 TE 之后弦平面 ⇒ 非壁元、P11 修不了）；★**conforming fine M0.5 不收敛**（限幅/落地 ~1.4k NaN）= G9.1 的亚声速对照。**修复 = 尾迹模型改（卷起/显式翼尖涡，G13.2 未来，实现交 Track B 重定 B9）；G13.3（未来）= 模型到位后重做 M6 三点 Richardson**。已纠正 roadmap:1036/demo_report:1333 "1/r"→"1/√r" |
| M0/M1 | ✓ | 准 2D 挤出网格族 + M6 后掠尾迹网格族（.msh gitignored，脚本再生） |
| M3 / M4 | ✓ | 无尾迹面("O 型")网格族（2026-07-11,Track B 双网格规则用）：M3 = NACA 准 2D（走廊扇形覆盖 α 扫掠;coarse 已提交）、M4 = ONERA M6（复用 M1 走廊尺寸场,单元数与 M1 差 6–9% ⇒ B7 受控 A/B;.msh 全部 gitignored） |
| Track B | ◐ | **B1–B5 ✓ + B7 ✓（B1 2026-07-11;B2/B3/B4/B5/B7 2026-07-12）;B6 跨声速进行中（2026-07-12）;下一步 = B8**。level-set 尾迹产生升力,隐式 Kutta,无 Γ secant。**B4：尾迹 LS 对常数跳跃恒零（单位分解,1.9e-16）⇒ "g₂ 即 Kutta" 错、退休;真 Kutta = 非线性 TE 压力相等 |q_u|²=|q_l|²,q 在壁面邻接控制体恢复**（全扇形 +11~15%、壁面邻接 <1%）。Γ 与 conforming 同网格 <1%,无尾迹 M3 差 0.3%。**B5：远场 A/B — 选项 a（Dirichlet+涡）保持默认（亚声速）**;选项 b（López Neumann 出口、无涡）截断 O(Γ/R),15c −4%。**B6 已测得三个机理结论（design_track_b.md §10）**：逐侧人工密度落地（亚临界严格无操作,激波区与 conforming 同构）;P4 全场阻尼**不可平移**（Jacobi 光滑子掐死作为解模态的 Γ）⇒ 阻尼局部化到 ν>0 行;★折叠区 option a 实时涡反馈增益 >1（Γ 单调穿过真解后爆掉,lagged 外层映射无不动点）⇒ **跨声速配方 = Neumann 出口**,coarse M0.80 收敛到距 Newton 真解 −6%（比 conforming Picard 自己的 stall 态近得多）;★A/B 反转——原始 Picard-vs-Picard 差（M0.75 +10.5%）经同网格 **Newton 仲裁**证明是 **conforming Picard 的欠环量偏置**（−4~−8%,P4 勘误在弱激波区的定量化）,**LS Picard 距 Newton 真解仅 +0.25%~+1.0%** 且随加密收敛（隐式 Kutta 无可早停的 Γ 外层）。**gate 基线已改（用户裁决 2026-07-12）为同网格 Newton 真解**;**coarse M0.80 gate 达成**（M0 Γ −7.9% / M3 +0.9%,激波 0.644/0.678,demo 14/14 PASS）;**medium M0.7875 = 折叠区**——Picard 只到有界 stall（Γ −18.8%）。**★ LS Newton（`solve/newton_ls.py`，§5.5/§10.6）已交付+FD 验证 1.3e-9**：精确 Jacobian = Picard 矩阵 + 逐侧 Term 2/3 + TE Kutta 精确二次导数,尾迹 LS 行线性,无 Γ DOF/Woodbury。**在折叠区拿到机器精度、二次收敛的真离散解（0 lim/flr）**——medium M0.7875 **M3 无尾迹（工作流网格）|R| 1.5e-12** Γ 0.2292、coarse M0.80 双网格皆收敛。**两个诚实缺口**：M0 嵌入 medium 活选择 Newton 极限环于 3e-6（P8/N5 近平局 churn 的 LS 形态,需接入冻结选择）;收敛的 LS 折叠解比 conforming Newton 低 ~13%（真离散差异,待厘清 neumann 截断/切层 O(h)/人工密度网格依赖）。设计 [design_track_b.md](design_track_b.md)。B6:coarse ✓、LS Newton 交付+折叠区工作流网格达真解、medium 定量收尾余两项。**★ B7（ONERA M6 三维 gate）2026-07-12 一次通过**（design_track_b.md §11;demo 35/35 PASS）:M∞0.84/α3.06 coarse、`farfield="neumann"`、ramp 0.60→0.84@dm0.04;**M1 嵌入** cl_KJ 0.2765／激波 0.635/0.588/0.449／Γ 0.1076→−0.0003／M_max 1.453／**0 lim,flr**（22.7 min）,**M4 无尾迹** 0.2710／0.634/0.584/0.454／0.1055→+0.0003／1.368／**0 lim,flr**（18.4 min）;V6 1.77%/1.97%;双网格 A/B 2.0%。**★ B6 的升力反转在三维复现**:以 conforming **Newton 真解** 0.2692 为基准,LS Picard **+2.7%(M1)／+0.7%(M4)**,而 conforming Picard（P5, 0.24788）**−8.6%** ⇒ 按 P5 设 gate 等于因"更接近真解"而扣分,故升力锚定 Newton;**无尾迹工作流网格 M4 是两者中更准的**。**★ 三维远场 = neumann,P5 的 Γ(z) taper 结构性不需要**（B 路径的涡展向均匀、支割线恒为 y=0:(a) α 对准的尾迹面与之非共面 ⇒ 出口承载无切割支撑的跳跃,出口 M 0.958 vs neumann 0.513;(b) 即使共面,标量 Γ 也匹配不了 Γ(z)→0 且翼尖外无切割,出口 M 0.825 = P5 branch-ray 伪迹;neumann 不带涡故二者皆不可能）。**★ Γ(z) 本身展向光滑、无需任何平滑**（计划外结论）:归一化二阶差分 0.0079/0.0091 vs conforming P5 的 0.0970 = **光滑 11–12 倍**——conforming 每个 TE 站点解一个独立 secant（st133 那套机械;`INVESTIGATION_gamma_smoothing.md` 曾试图平滑其抖动而失败）,隐式 Kutta 没有逐站循环（Γ 是**一个解模态**）⇒ P5 的展向 Γ 问题在此**结构上不可能发生**,而非被修好。**★ 三维专属机械无需新求解器代码**（B1 的斜交标架+展向裁剪已生效,Γ(tip)→~3e-4 离散归零）,唯一缺口=后处理 `post/surface_ls.py::section_cp_curve_levelset`（D11 逐侧;**上表面与 main 场逐位相同** ⇒ gate 激波指标不受影响）+ `cl_pressure_3d_levelset`。成本远低于风险估计(~0.6 s/outer @ ~12k 三维 DOF ⇒ ~20 min/解)。**诚实缺口**:顶部 Mach 级停在 Picard 残差尾（|R| 4–6e-6,600-outer 上限;但每级有界物理 0 lim/flr、指标在带内 ⇒ gate 断言"有界"而非 `converged`）;**LS Newton 在 M6 推迟**（朴素 splu,需 lagged-LU);激波比 P5 靠后 0.02–0.04c（在带内）;仅 coarse。**下一步 = B8（多尾迹）** |
| M2 / Track V | 未开 | 翼身组合体宜与 Track B 同排（B8);Track V = VII 粘性耦合，已设计未动工 |

## 长期挂起项（勿反复重提）

- **G1.6 球面 Cp <2%**（strict xfail，11.6%）：已定根因 = 平坦面片壁上的自然边界条件
  （变分罪）；h 加密、恢复调参、Nitsche、边界数据修正**均已用证据排除**。
  唯一在案路线 = Option C 重定义 + P11 曲面壁元。
- **V6 <1%**（CL_p 对 CL_KJ 的 O(h) 离散地板）：挂到 P11，P9 正在检验其归因。
- 后掠 TE Kutta 探针跨站共享（P5 记录的鲁棒性隐患，未修）。

## 当前最大困难（按层级）

1. **模型极限（不是 bug）：FP 非唯一性折叠。** NACA medium 在 M≈0.79 附近进入
   全势方程的折叠区（dcl/dM ≈ 6–10）：同工况不同网格不可比、warm-start 有陷阱
   （G10.2b 负结果）、M0.80 无可达孤立解、无 ramp 直接求解 class C（G10.3）。
   对策已固化为纪律：折叠区只做单网格回归锁、Mach ramp 保留、松容差禁用。
2. **尾迹模型：刚性平面尾迹的翼尖涡片边缘奇点（P9 新定性，P13/G13.1 定量，已取代
   "平坦 P1 壁面"归因）。** 0.019 的 M6 升力缺口（cl_KJ 0.2692 vs Tranair/KRATOS 0.288）
   曾归因于尖 TE/LE 的 P1 壁面梯度；P9 推翻了它：2D 尖 TE 干净收敛（0.03%），而 3D 的
   fine 网格根本不是离散解——奇点在**尾迹片的自由翼尖边**（真实流动在此卷成翼尖涡，
   而刚性平面片直接"截断"）。**P13/G13.1 定量：驱动量是脱落涡强 dΓ/dz（翼尖最大，非
   束缚 Γ——Γ 翼尖→0），奇点是 1/√r 平板边缘型（conforming 三点 p≈0.59，非"1/r"）,
   在 conforming 与 level-set 两路径皆随加密发散（模型非表示）。** 修复 = 尾迹模型改
   （卷起/显式翼尖涡，P13/G13.2，实现交 Track B 重定 B9）。**任何 3D 网格收敛断言都必须
   等翼尖/尾迹处理落地。**
3. **Track B 的下一个未知数：跨声速（B6）。** 切单元有**两个速度状态**（上/下），
   人工密度必须逐侧求值，且上游 walk 的邻接图要**限制在同侧**（尾迹是滑移线，密度
   信息不穿过切向间断）——design_track_b.md §5.2/D10（DN1 曾错判为"upwind 不用改"）。

## 运行/成本纪律（易踩的坑）

- 16 线程封顶且**必须同时盖住 BLAS/OMP**（`NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
  OPENBLAS_NUM_THREADS=16`；漏了 BLAS 会慢 ~33% 并挂掉 G8.2 断言）。
- 贵重工件不随手重算（P4 heavy demo ~40 min、G4.1 medium ~17 min、P5 medium
  45–75 min、M6 fine 网格数分钟）；已提交的 CSV/PNG 是权威。
- fine 网格与解算 npz 一律 gitignored，本地缓存、demo 缺则重算。
- 回归基线：**276 passed + 17 skipped + 2 xfailed**（B7 起,2026-07-12,
  实测 719.29 s @16 线程;较 B6 的 270+12+2 增加 `test_b7_onera_m6.py`（6 快测 +5
  gated——gated 是 M6 三维解,每个约 20 min,故 skipped 12→17）;G8.3 的 302 s 仍是
  CI 参考;含 B1–B7 的 Track B 测试,其中若干条在无尾迹网格未本地生成时跳过;
  重 gate 走 `PYFP3D_TRANSONIC_GATES=1`）。
  内核/装配改动后先跑 `tests/test_v0_freestream.py`。
- P9 demo 等重活在跑时,新测试/求解一律降到 8 线程(NUMBA/OMP/OPENBLAS 同盖)。
