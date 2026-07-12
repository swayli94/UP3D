# pyFP3D 开发现状（快照）

> **快照日期：2026-07-12。** 本文件是给人读的高层概览，**不是**权威进度源：
> 阶段/gate 状态以 [roadmap.md](roadmap.md)（含进度台账）为准，当前阶段细节以
> [agent-rules.md](agent-rules.md) 的 "Current phase" 行为准，证据在
> [demo_report.md](demo_report.md)。若本文件与它们冲突，以它们为准。
> 每次 gate 关闭后按 CLAUDE.md 工作流更新 roadmap/agent-rules 时，顺手刷新本文件。

## 一句话状态

求解器主线（Track P）P0–P9 全部关闭（P10 的 G10.2/G10.3 亦关）。**P9 判别阶段给出
了否定性但决定性的结论：3D 升力缺口不是壁面单元的问题，而是刚性平面尾迹在翼尖自由
边上的涡片边缘奇点**——即 Track B 存在的理由；P11 曲面壁元作为"3D 升力修复"不再被
支持（其 G1.6 球面 Cp 理由仍然成立，待用户仲裁）。**当前工作重心 = Track B**：
B1–B5 已关（level-set 尾迹 + 隐式 Kutta 已产生升力，无 Γ secant），**正在开展 B6
（跨声速 + Mach 续接）**。

## 阶段一览

| 阶段 | 状态 | 要点 |
|------|------|------|
| P0–P3 | ✓ | 基础设施、Laplace、尾迹割缝/Kutta、亚声速压缩（P1 仅 G1.6 未关，见下） |
| P4 | ✓（带勘误） | 人工密度跨声速 Picard；★2026-07-11 勘误：Picard 态**不是**离散方程的解（Newton 残差 2.2e-4），P4 gate 降级为 Picard 质量/鲁棒性 gate（warm-start 引擎） |
| P5 | ✓（带注记） | ONERA M6 M0.84/α3.06；同类欠收敛注记（真解 cl 高 +5.8%/+7.9%，激波位置不变） |
| P6 | ✓ | 表面 Cp 锯齿 = 壁面梯度**恢复**伪迹（非通量），法向门控平滑修复 |
| P7 | ✓ | 冻结选择下的 ∂ρ̃/∂φ 精确灵敏度，FD 验证 3e-10~6e-9 |
| P8 | ✓ | 全耦合 (φ_red, Γ) Newton：精确 Jacobian + 直接步 + 停滞自适应冻结；G8.2 M6 medium 249 s（G10.2 后 ~145 s）、G8.3 全套件 302 s |
| P9 | ✓（2026-07-11） | 判别阶段。**G9.2 干净 PASS：2D 尖 TE 升力误差 2.71%→0.33%→0.03% ⇒ 尖 TE 不构成升力地板**（"尖 TE/LE 平壁"归因的 2D 支腿没了）。**G9.1 INVALID：M6 fine 不是离散解**（无限幅 M_max 1.40→2.13→**7.93**，9 个单元触 M_cap=3 ⇒ 无三点 Richardson）。★**发散奇点位于尾迹片的自由翼尖边**（9 个封顶单元全在 z/b 0.998–1.000、TE 之后、弦平面内）= 刚性平面尾迹的经典涡片边缘奇点 ⇒ **尾迹模型缺陷（Track B 的靶心），不是壁元缺陷**。G9.3：0.019 缺口按原提法**不可拆分**（份额如实记 n/a）。顺带测得：`precond="amg"` + 紧 EW forcing（η=1e-8）在每个规模都更快且有效 |
| P10 | ◐ | G10.2 ✓（分裂裁决：M6 +41% 提速已晋升；折叠区禁用松容差）；G10.3 ✓（裁决：**保留 Mach ramp**）；G10.1（非升力 Newton 入口）未动、无顺序约束 |
| P11 | 未开（条件性） | 曲面/等参壁元。**G9.3 裁决：作为"3D 升力修复"不被 P9 支持**（曲面壁元去不掉尾迹片边缘奇点）；仅剩 G1.6（光滑曲壁上的球面 Cp，另一机理）为有效理由 — **待用户仲裁** |
| P12 | 未开 | Backlog：离散伴随、VII transpiration、混合单元/BO 标定 |
| M0/M1 | ✓ | 准 2D 挤出网格族 + M6 后掠尾迹网格族（.msh gitignored，脚本再生） |
| M3 / M4 | ✓ | 无尾迹面("O 型")网格族（2026-07-11,Track B 双网格规则用）：M3 = NACA 准 2D（走廊扇形覆盖 α 扫掠;coarse 已提交）、M4 = ONERA M6（复用 M1 走廊尺寸场,单元数与 M1 差 6–9% ⇒ B7 受控 A/B;.msh 全部 gitignored） |
| Track B | ◐ | **B1–B5 ✓（B1 2026-07-11;B2/B3/B4/B5 2026-07-12）;B6 跨声速进行中（2026-07-12）**。level-set 尾迹产生升力,隐式 Kutta,无 Γ secant。**B4：尾迹 LS 对常数跳跃恒零（单位分解,1.9e-16）⇒ "g₂ 即 Kutta" 错、退休;真 Kutta = 非线性 TE 压力相等 |q_u|²=|q_l|²,q 在壁面邻接控制体恢复**（全扇形 +11~15%、壁面邻接 <1%）。Γ 与 conforming 同网格 <1%,无尾迹 M3 差 0.3%。**B5：远场 A/B — 选项 a（Dirichlet+涡）保持默认（亚声速）**;选项 b（López Neumann 出口、无涡）截断 O(Γ/R),15c −4%。**B6 已测得三个机理结论（design_track_b.md §10）**：逐侧人工密度落地（亚临界严格无操作,激波区与 conforming 同构）;P4 全场阻尼**不可平移**（Jacobi 光滑子掐死作为解模态的 Γ）⇒ 阻尼局部化到 ν>0 行;★折叠区 option a 实时涡反馈增益 >1（Γ 单调穿过真解后爆掉,lagged 外层映射无不动点）⇒ **跨声速配方 = Neumann 出口**,coarse M0.80 收敛到距 Newton 真解 −6%（比 conforming Picard 自己的 stall 态近得多）;★A/B 反转——原始 Picard-vs-Picard 差（M0.75 +10.5%）经同网格 **Newton 仲裁**证明是 **conforming Picard 的欠环量偏置**（−4~−8%,P4 勘误在弱激波区的定量化）,**LS Picard 距 Newton 真解仅 +0.25%~+1.0%** 且随加密收敛（隐式 Kutta 无可早停的 Γ 外层）。**gate 基线已改（用户裁决 2026-07-12）为同网格 Newton 真解**;**coarse M0.80 gate 达成**（M0 Γ −7.9% / M3 +0.9%,激波 0.644/0.678,demo 14/14 PASS）;**medium M0.7875 = 折叠区**——Picard 只到有界 stall（Γ −18.8%）。**★ LS Newton（`solve/newton_ls.py`，§5.5/§10.6）已交付+FD 验证 1.3e-9**：精确 Jacobian = Picard 矩阵 + 逐侧 Term 2/3 + TE Kutta 精确二次导数,尾迹 LS 行线性,无 Γ DOF/Woodbury。**在折叠区拿到机器精度、二次收敛的真离散解（0 lim/flr）**——medium M0.7875 **M3 无尾迹（工作流网格）|R| 1.5e-12** Γ 0.2292、coarse M0.80 双网格皆收敛。**两个诚实缺口**：M0 嵌入 medium 活选择 Newton 极限环于 3e-6（P8/N5 近平局 churn 的 LS 形态,需接入冻结选择）;收敛的 LS 折叠解比 conforming Newton 低 ~13%（真离散差异,待厘清 neumann 截断/切层 O(h)/人工密度网格依赖）。设计 [design_track_b.md](design_track_b.md)。B6:coarse ✓、LS Newton 交付+折叠区工作流网格达真解、medium 定量收尾余两项;B7 M6 三维在后 |
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
2. **尾迹模型：刚性平面尾迹的翼尖涡片边缘奇点（P9 新定性，已取代"平坦 P1 壁面"
   归因）。** 0.019 的 M6 升力缺口（cl_KJ 0.2692 vs Tranair/KRATOS 0.288）曾归因于
   尖 TE/LE 的 P1 壁面梯度；P9 推翻了它：2D 尖 TE 干净收敛（0.03%），而 3D 的
   fine 网格根本不是离散解——奇点在**尾迹片的自由翼尖边**（真实流动在此卷成翼尖涡，
   而刚性平面片直接"截断"）。**任何 3D 网格收敛断言都必须等翼尖/尾迹处理落地。**
3. **Track B 的下一个未知数：跨声速（B6）。** 切单元有**两个速度状态**（上/下），
   人工密度必须逐侧求值，且上游 walk 的邻接图要**限制在同侧**（尾迹是滑移线，密度
   信息不穿过切向间断）——design_track_b.md §5.2/D10（DN1 曾错判为"upwind 不用改"）。

## 运行/成本纪律（易踩的坑）

- 16 线程封顶且**必须同时盖住 BLAS/OMP**（`NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
  OPENBLAS_NUM_THREADS=16`；漏了 BLAS 会慢 ~33% 并挂掉 G8.2 断言）。
- 贵重工件不随手重算（P4 heavy demo ~40 min、G4.1 medium ~17 min、P5 medium
  45–75 min、M6 fine 网格数分钟）；已提交的 CSV/PNG 是权威。
- fine 网格与解算 npz 一律 gitignored，本地缓存、demo 缺则重算。
- 回归基线：**268 passed + 10 skipped + 2 xfailed**（B6 起,2026-07-12;
  +9 B6 快测 +2 gated,较 B5 的 259+8+2;实测 682.77 s @16 线程但紧接重算之后
  冷缓存,G8.3 的 302 s 仍是 CI 参考;含 B1–B6 的 Track B 测试,其中若干条在
  无尾迹网格未本地生成时跳过;重 gate 走 `PYFP3D_TRANSONIC_GATES=1`）。
  内核/装配改动后先跑 `tests/test_v0_freestream.py`。
- P9 demo 等重活在跑时,新测试/求解一律降到 8 线程(NUMBA/OMP/OPENBLAS 同盖)。
