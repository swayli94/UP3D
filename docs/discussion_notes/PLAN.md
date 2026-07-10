# UP3D 开发计划（PLAN.md）

> **本文件**是 UP3D 项目的统一开发计划，整合 `docs/roadmap.md` 的阶段定义与进度、各 design note 的理论分析与技术决策，按开发逻辑排序。
>
> **维护规则**：每次 roadmap 阶段状态变化或 design note 新增/更新时，同步更新本文件。
>
> **权威性**：`docs/roadmap.md` 是 Track P/M/B/V 阶段状态与 gate 定义的**唯一权威**
> （2026-07-10 起 Track B/V 已并入 roadmap）；本文件是跨线整合视图，冲突时以 roadmap 为准。
>
> **交叉引用约定**（各 DN 文件与本文件同目录 `docs/discussion_notes/`）：
> - `roadmap Pn` → `docs/roadmap.md` §Pn
> - `DN1 §x` → `20260707_1505_levelset_wake_design.md` §x
> - `DN2 §x` → `20260707_2118_ibl_viscous_coupling_design.md` §x
> - `DN4 §x` → `20260708_1646_newton_solver_strategy.md` §x
> - `DN5 §x` → `20260708_1744_implementation_audit.md` §x
> - `DN6 §x` → `20260709_0145_3d_vii_implementation_analysis.md` §x
>
> **命名消歧**：Track V 的 V1–V4 是**开发阶段**编号，与 design.md §10 验证阶梯的
> 用例编号 V0–V6（如 "V6 < 1%" 的 Γ 一致性指标）无关；Track V 的 gate 记作 GV<阶段>.<n>。

---

## 总体架构

UP3D 是三维非结构网格全速势方程（FPE）数值模拟工具，目标是支持亚声速—跨声速工况的气动分析，并逐步引入粘性耦合（VII）能力。

**四条开发线**：

| 开发线 | 目标 | 当前状态 |
|--------|------|---------|
| **Track M** — 网格生成 | 为求解器提供非结构网格，支持后续复杂几何 | M0–M1 已关闭；M2 待开发 |
| **Track P** — 势流求解器 | FPE 的准确、高效求解 | P0–P6 已关闭；P7–P10 待开发 |
| **Track B** — Level-Set 尾流 | 嵌入式尾流替代 conforming 网格尾流 | B1–B6 全部待开发（设计完成） |
| **Track V** — 粘性耦合 (VII) | IBL3 6 方程 + transpiration 耦合 | V1–V4 全部待开发（设计完成） |

**依赖关系**：

```
Track M:  M0 ✓ → M1 ✓ → M2 ☐
                      Ↄ (为 P/B/V 提供网格)
Track P:  P0 ✓ → P1 △ → P2 ✓ → P3 ✓ → P4 ✓ → P5 ✓ → P6 ✓ → P7 ✓ → P8 ☐ → P9 ☐ → P10 ☐
                                                                    ↑
Track B:  (B1 → B2 → B3 → B4 → B5)   [B6 搁置]          B3 ⇢ P8（可选 Kutta 接口，非依赖；
                                                         P8 Newton 先在 conforming 尾流上落地）
Track V:  (V1 → V2? → V3 → V4)                          V3 依赖 P8; V1 依赖 P6
```

---

## 进度总览

| 阶段 | 状态 | 关闭日期 | 备注 |
|------|------|---------|------|
| M0 | ✅ | 2026-07-06 | 网格生成基础设施 |
| M1 | ✅ | 2026-07-07 | ONERA M6 翼网格 |
| M2 | ☐ | — | 翼身组合体尾流-机身交接 |
| P0 | ✅ | 2026-07-06 | 脚手架 |
| P1 | △ | — | G1.6 球 Cp 待 P9 |
| P2 | ✅ | 2026-07-06 | 尾流切割 + Kutta |
| P3 | ✅ | 2026-07-07 | 亚声速可压缩 |
| P4 | ✅ | 2026-07-07 | 跨声速人工密度 |
| P5 | ✅ | 2026-07-08 | ONERA M6 3D 验证 |
| P6 | ✅ | 2026-07-08 | Cp 锯齿消除 |
| P7 | ✅ | 2026-07-10 | frozen-walk ∂ρ̃/∂φ,FD 验证 3–5e-10 |
| P8 | ☐ | — | 全耦合 Newton |
| P9 | ☐ | — | 曲面壁面单元 |
| P10 | ☐ | — | Backlog |
| B1–B5 | ☐ | — | Level-Set 尾流（设计完成） |
| B6 | ⏸ | — | 曲线尾流（搁置） |
| V1 | ☐ | — | IBL3 松耦合 |
| V2 | ☐ | — | Quasi-simultaneous（可选） |
| V3 | ☐ | — | 紧耦合 Newton |
| V4 | ☐ | — | 尾流面 IBL 修正 |

**测试基线**：165 passed + 4 skipped + 2 xfailed = 171 collected（HEAD after P7；4 skipped = `PYFP3D_TRANSONIC_GATES=1` 门控的重载 gate；M6 `.msh` 缺席时另有 13 个 M1 测试转 skip）

---

## Track M — 网格生成

### M0 — 网格生成基础设施

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-06) |
| **roadmap** | `roadmap Track M` |
| **交付物** | `meshgen/planar.py` (NACA0012 2.5D + 圆柱), `meshgen/extrude.py` (prism→tet split) |
| **网格族** | NACA0012: 16.4k/61.8k tets |

### M1 — ONERA M6 翼网格

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-07) |
| **roadmap** | `roadmap Track M` |
| **交付物** | `meshgen/wing3d.py` (ONERA M6 半翼) |
| **网格族** | ONERA M6: 55.5k/350.7k/2513k tets |

### M2 — 翼身组合体（Backlog）

| 项 | 内容 |
|---|---|
| **状态** | ☐ Backlog |
| **roadmap** | `roadmap Track M` (M2) |
| **目标** | 翼身组合体网格生成，处理尾流-机身交接 |
| **备注** | 需要 Track B 的嵌入式尾流支持（尾流面与机身相交） |

---

## Track P — 势流求解器

### P0 — 仓库脚手架 + 网格基础设施

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-06) |
| **roadmap** | `roadmap P0` |
| **交付物** | `mesh/reader.py`, `metrics.py`, `coloring.py`, `physics/isentropic.py`, `post/vtk_out.py` |
| **验证** | G0.1–G0.4 单元测试通过；87 passed + 2 xfailed |

### P1 — Laplace 求解器 (ρ ≡ 1)

| 项 | 内容 |
|---|---|
| **状态** | △ 进行中（G1.1/G1.2 关闭；G1.6 待 P9 解决） |
| **roadmap** | `roadmap P1` |
| **交付物** | `kernels/residual.py`, `solve/linear.py`, `solve/picard.py`, `post/surface.py` |
| **关键发现** | G1.6 球 Cp 误差 ~11.6% 根因是 flat-facet P1 壁面的变分犯罪（非恢复问题），需 curved elements（P9）解决；DP1 选定 "> 5%" 分支 |
| **阻塞项** | G1.6 pending Option C re-spec → 实际由 P9 curved elements 解决 |

### P2 — 尾流切割、环量、Kutta 条件（Laplace 基础）★ 关键阶段

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-06) |
| **roadmap** | `roadmap P2` |
| **交付物** | `mesh/wake_cut.py` (节点复制 + flood-fill 侧分类), `constraints/wake.py` (master-slave T^TAT 消元), `constraints/dirichlet.py` (远场 + 涡修正), `solve/picard.py::solve_laplace_lifting` (secant Kutta) |
| **验证** | G2.1–G2.5 全部通过；100 passed + 2 xfailed |
| **关键设计决策** | TE 节点复制（非单值），避免 spurious TE suction ~Γ²/h |

### P3 — 亚声速可压缩

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-07) |
| **roadmap** | `roadmap P3` |
| **交付物** | 高速装配路径 (`kernels/gradient.py`, `jacobian.py`, `PicardOperator`), `solve/picard.py::solve_subsonic_lifting` (嵌套 Picard + secant Kutta), PG-scaled 涡远场, 可压缩 Cp |
| **验证** | G3.1 cl 0.32% (< 2%); G3.2 cl −0.33%, 15 iterations; G3.3 M∞=0 bit-identity；117 passed |

### P4 — 跨声速：人工密度

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-07, 当天修复重测) |
| **roadmap** | `roadmap P4` |
| **交付物** | `solve/continuation.py::TRANSONIC_DEFAULTS` (damping_theta), `solve/picard.py` 阻尼参数化 |
| **验证** | G4.1 medium 激波 x/c 0.633, M_max 1.366, cl 0.349；G4.2 亚临界 bit-identity；G4.3 10-case sweep；136 passed |
| **关键问题** | medium 首次运行发散 → 根因: pseudo-time damping 弱化 → 改用 damping_theta=0.2 |

### P5 — 3D 验证：ONERA M6

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-08) |
| **roadmap** | `roadmap P5` |
| **交付物** | `post/surface.py::cl_kj_3d`, `solve/continuation.py` rtol 参数化, `farfield_spanwise_gamma` (Γ(z) taper), `n_kutta_polish` (fixed-Γ Kutta polish) |
| **验证** | V6 1.82% medium (3% floor, <1% deferred to P9)；CL 0.2453；140 passed |
| **关键教训** | medium 失败根因三次翻转：最终确认为单站 Kutta 闭合失败 (st133, z/b=0.801, 32% under-circulated)，而非 TE 奇异性或 secant 耦合失稳。教训：correlation is not causation；1/h 奇异性不可能依赖于 Γ |
| **遗留** | V6 < 1% deferred to P9 (sharp-TE P1 wall gradient floor)；swept-TE Kutta probe 共享 (st133/134, 记录未修) |

### P6 — 表面压力恢复（Cp 锯齿消除）★ 精度阶段

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-08) |
| **roadmap** | `roadmap P6` |
| **交付物** | `post/surface.py::smooth_wall_tangential_gradients` (Jacobi 边邻居平均, TE crease gating), `smooth_passes` 参数贯穿 `wall_cp_curve`/`section_cp_curve`/`wall_force_coefficients` |
| **关键发现** | Cp 锯齿是 wall-gradient recovery 伪影（P1 常梯度在三角形上不连续），**不是**人工密度通量的问题。修复在后处理，不动求解器 |
| **验证** | G6.1 锯齿 0.0758→0.00023 (330×); G6.2 物理不变; G6.4 bit-identity (smooth_passes=0)；157 passed |
| **遗留** | G6.3 V6 re-measure under smoothing → P9 (residual floor) |

### P7 — 可微人工密度通量（Newton 前置）

| 项 | 内容 |
|---|---|
| **状态** | ✅ 已关闭 (2026-07-10) |
| **roadmap** | `roadmap P7` |
| **design note** | DN4 §3.2 (strategy A 全 Jacobian), DN5 §4.3 (Jacobian 敏感度), DN5 §11.1 |
| **交付物** | `kernels/upwind.py::rho_tilde_sensitivities_sweep` + `UpwindOperator.rho_tilde_sensitivities`（frozen u(e) 下分支式 (s_e, s_u)，floor/self-upstream 平坦分支→0）; `physics/isentropic.py::mach_squared_derivative_wrt_q_sq` |
| **验证** | G7.3 JVP-vs-FD（对 SHIPPED `rho_tilde_sweep`，冻结选择）：cube 3–5e-10（8 测试，含 PYFP3D_NOJIT）；NACA coarse 构造场 3.5e-9；**收敛 G4.1 M0.80 场 5.7e-9**（pocket 1189 accel + 977 shock-point）。`max_ε` 未实现（范围决策）——前向通量逐字节不变，G7.1/G7.2 by construction；165 passed |
| **关键发现** | (1) design.md §6.3 符号 FD 裁决：`dμ/dM² = +M_c²/M⁴`（原 "−" 为笔误，已改）；(2) max 折点陷阱：FD 跨折点读到分支平均（~1e-5，非 bug）；可分离场在结构化/prism-split 网格上整片元素落在平局上——未来 FD 检验必须用打破简并的通用场 |
| **Gate** | G7.1/G7.2 held by construction; G7.3 ✅ 3–5e-10 (< 1e-6) |

### P8 — 全耦合 Newton + 性能 ★ 待开发

| 项 | 内容 |
|---|---|
| **状态** | ☐ 待开发（设计完成 2026-07-08） |
| **roadmap** | `roadmap P8` |
| **design note** | DN4 全文 (strategy A + 全耦合 + load stepping), DN5 §11.2 (Newton 实现细节), DN5 §11.3 (线性求解器), DN6 §12 (与 VII 的交互) |
| **目标** | 增广 (φ_red, Γ) 全耦合 Newton 求解，替代 Picard + secant |
| **关键设计** | (1) 完整 Jacobian "strategy A" — 保留 ∂μ/∂φ + 上游耦合 → 严格二次; (2) 全耦合 (φ_red, Γ) — Γ-Jacobian 近乎现成 (`wake.py::self._h`); (3) load stepping — M∞ ramp at fixed M_crit/μ_c; (4) sharp TE 不是 Newton 障碍 (López Eq.4.2 sharpened TE, 仍二次收敛) |
| **子阶段** | N0(optional density smoothing) / N1(P7 flux) / N2(Jacobian) / N3(GMRES+AMG) / N4(coupled driver) / N5(transonic+load-stepping) / N6(M6+perf) |
| **Gate** | G8.1 二次收敛 + FD-verified Jacobian; G8.2 M6 medium < 5 min; G8.3 suite < 10 min |
| **前置** | P7 |

### P9 — 曲面 / 等参壁面单元 ★ 待开发

| 项 | 内容 |
|---|---|
| **状态** | ☐ 待开发 |
| **roadmap** | `roadmap P9` |
| **目标** | 解决两个 P1 flat-facet 壁面精度瓶颈 |
| **解决问题** | (i) G1.6 球 Cp < 2% (当前 ~11.6%, 变分犯罪); (ii) V6 < 1% floor (sharp-TE/LE P1 wall gradient; M6 CL 0.245 vs 0.288) |
| **Gate** | G9.1 球 Cp < 2%; G9.2 V6 < 1% on M6 medium |

### P10 — Backlog (post-v1.0)

| 项 | 内容 |
|---|---|
| **状态** | ☐ Backlog |
| **roadmap** | `roadmap P10` |
| **项目** | (1) Discrete adjoint = P8 Newton Jacobian 转置; (2) VII hook: transpiration BC (→ Track V); (3) Mixed prism/tet + (C, M_c, ω) BO calibration |

---

## Track B — Level-Set 嵌入式尾流

> **设计文档**：DN1 全文
> **动机**：当前 conforming 网格尾流（`wake_cut.py` + master-slave）有三个结构性限制——(1) 网格必须预嵌入尾流面; (2) 改攻角需重新生成网格; (3) 钝尾缘无明确锚定点。方案 B 用 level-set 尾流 + 多值 FE (CutFEM) 解决。
> **状态**：设计完成，实现全部待开发

| 阶段 | 内容 | 状态 | 依赖 | design note |
|------|------|------|------|-------------|
| **B1** | Level-set 尾流 + 切单元识别（验证基础设施） | ☐ | — | DN1 §8 Phase B1 |
| **B2** | 多值 FE 装配（核心数值验证） | ☐ | B1 | DN1 §8 Phase B2 |
| **B3** | 罚函数 Kutta + 升力求解 | ☐ | B2 | DN1 §8 Phase B3 |
| **B4** | 跨声速 + Mach continuation | ☐ | B3 | DN1 §8 Phase B4 |
| **B5** | 多尾流验证（多段翼/翼身组合体） | ☐ | B4 | DN1 §8 Phase B5 |
| **B6** | 曲线尾流 / 自由尾迹 | ⏸ 搁置 | — | DN1 §8 Phase B6; 决策见 DN2 §4.5.6 |

**B6 搁置原因**（DN2 §4.5.6）：精度收益微小 (O(0.1%))、工程代价不成比例、与 Newton 紧耦合结构性冲突、López 先例、优先级排序。保留 `update_direction()` 接口能力以备未来。

**Track B 与 Track P 的关系**：B1–B4 完成后可以替代当前 conforming 尾流路径；B3 为 P8 Newton 提供罚函数 Kutta 的可能接口；但 Track B 不阻塞 Track P 的 P7–P10（当前 conforming 尾流继续工作）。

**时序护栏（2026-07-10，roadmap Track B 同步记录）**：P8 全耦合 Newton 的 Γ-Jacobian 块建立在 master-slave 消元上（`wake.py::self._h`），而 B3 罚函数 Kutta 会取消 Γ 自由度——两套 Newton 设计不要并行推进。P8 先在 conforming 尾流上落地；level-set 路径的 Newton 是 B4 之后的重推导。

---

## Track V — 粘性耦合 (VII)

> **设计文档**：DN2 全文（IBL 耦合方案）、DN6 全文（3D VII 实现路径分析）
> **动机**：纯全速势方程无法捕获粘性效应（摩擦阻力、激波-边界层相互作用、升力修正）。Drela IBL3 6 方程 + transpiration BC 是选定的耦合方案。
> **状态**：设计完成，实现全部待开发

### V1 — IBL3 求解器 + 松耦合（独立验证）

| 项 | 内容 |
|---|---|
| **状态** | ☐ 待开发 |
| **design note** | DN2 §8.2 (松耦合循环), DN2 §10 Phase V1, DN6 §7.2 (松耦合流程), DN6 §8 (transpiration 实现), DN6 §10 (IBL 求解器设计) |
| **目标** | 势流 → IBL3 → transpiration → 重解势流，固定点迭代 |
| **交付物** | `viscous/ibl3.py` (6 方程表面 Galerkin FE), `viscous/transpiration.py` (δ*→ṁ), `viscous/coupling.py` (VII 松耦合控制器), 修改 `kernels/residual.py` 壁面 Neumann + transpiration, 修改 `constraints/wake.py` 尾流面 RHS + δ*_wake |
| **验证** | 2D NACA 0012 亚声速 δ* vs XFOIL; 松耦合 5-10 次迭代收敛; δ*=0 时 bit-identical |
| **前置** | P6 (smooth wall gradient 作为 IBL 输入) |

### V2 — Quasi-simultaneous coupling（可选中间步骤）

| 项 | 内容 |
|---|---|
| **状态** | ☐ 可选 |
| **design note** | DN2 §10 Phase V2, DN6 §15.5 item 6 (BLWF58 Hilbert 积分参考) |
| **目标** | Hilbert 积分预估势流响应，加速松耦合收敛 |
| **交付物** | `viscous/hilbert.py` (Hilbert 积分算子, 非结构网格) |
| **验证** | 松耦合迭代数减少 30-50%; 接近分离时收敛性改善 |
| **备注** | 如果 V1 已够快（5-10 次迭代），可跳过直接做 V3 |

### V3 — 紧耦合 Newton

| 项 | 内容 |
|---|---|
| **状态** | ☐ 待开发 |
| **design note** | DN2 §8.3 (紧耦合 Newton), DN2 §10 Phase V3, DN6 §7.3 (紧耦合流程), DN6 §12 (与 P8 Newton 的交互) |
| **目标** | 增广 (φ, Γ, δ, A, B, Ψ, C_τ1, C_τ2) 同时 Newton 求解 |
| **关键** | J_φ,δ* = ∂ṁ/∂δ* (transpiration Jacobian); J_δ*,φ = ∂R_δ/∂φ (势流速度对 BL); GMRES + 块预处理 (势流 AMG + BL ILU) |
| **验证** | Newton 二次收敛; 2D NACA 0012 跨声速 VII 修正激波位置 vs 实验; 3D ONERA M6 CL 从收敛的无粘值**下降**、向实验 ≈0.26–0.27 靠近（DN6 §13.3）。注意：0.245 vs 0.288 的差距是 P9 离散精度问题（sharp-TE/LE P1 壁面梯度），**不是**粘性能弥补的——粘性只会把 CL 往下修 |
| **前置** | P8 (Newton Jacobian 框架) + V1 |

### V4 — 尾流面 IBL 修正（V1 的自然延续）

| 项 | 内容 |
|---|---|
| **状态** | ☐ 待开发（非独立阶段） |
| **design note** | DN2 §4 (δ* 如何影响尾流面), DN2 §10 Phase V4, DN6 §9 (尾流面 IBL 修正) |
| **目标** | IBL3 方程从物面延续到尾流面，输出 δ*_wake 作为势流尾流面 RHS 源项 |
| **关键认识** | 不是独立于 V1 的工作——Drela IBL3 的 6 方程在物面和尾流面上是同一套方程，区别只在闭合关系。V1 实现时就应预留尾流面 6 个未知量 |
| **TE 处理** | DN2 §4.5：局部基自适应吸收几何折转（§4.5.1-4.5.3）；直尾流 + 质量穿透松弛（§4.5.6）；不做尾流几何松弛 |

### VII 的关键理论决策汇总

| 决策 | 内容 | design note |
|------|------|-------------|
| IBL 方程选择 | Drela IBL3 6 方程 (δ, A, B, Ψ, C_τ1, C_τ2) | DN2 §2 |
| δ* 耦合方式 | Transpiration BC（方案 VI-1），不动网格 | DN2 §3.2, §7 |
| 尾流面修正 | δ*_wake 作为尾流面 RHS 质量源 | DN2 §4.2 |
| TE 拐折处理 | Drela 局部基自适应 + 闭合关系切换，不需特殊 TE 方程 | DN2 §4.5.1-4.5.3 |
| 尾流偏折 | 直尾流 + 质量穿透松弛（López 路线），不做几何松弛 | DN2 §4.5.6 |
| 耦合策略 | 松→紧渐进：V1 松耦合 → V2? quasi-simultaneous → V3 紧耦合 Newton | DN2 §8, DN6 §7 |
| Level-set 用于 IBL | 不用于物面 δ* 修正（UP3D 是 body-fitted）；仅用于尾流面（Track B） | DN2 §6 |

### VII 的适用范围与限制

- **适用**：亚声速附着流；跨声速弱激波（激波前附着）；薄翼、小到中等攻角
- **不适用**：大攻角分离流；强激波诱导分离；钝体绕流
- **design note**：DN2 §9

---

## 跨线依赖关系总览

```
P6 ✓ (smooth wall gradient) ──────→ V1 (IBL 输入)
P7 ✓ (∂ρ̃/∂φ frozen walk) ────────→ P8 (Newton Jacobian)
P8 ☐ (Newton 框架) ───────────────→ V3 (紧耦合 Newton)
B3 ☐ (罚函数 Kutta) ──────────────→ P8 (可选 Kutta 接口)
P9 ☐ (curved elements) ───────────→ P1 G1.6 (球 Cp) + P5 V6 (< 1%)
```

**关键路径**：P7 → P8 → V3（Track V 的最终目标紧耦合 Newton）

**可并行路径**：
- Track B 的 B1–B2（切单元识别 + 多值装配，纯 Laplace 验证）可与 P7–P8 并行；**B3 起涉及 Kutta/Newton 结构，受上文时序护栏约束**（P8 先落地）
- V1 可在 P6 完成后独立开始（不依赖 P7–P8）；但 V1 是 Track-P 量级的独立求解器工程（6 方程非线性表面 FE + 闭合关系），按整阶段排期，不是顺手任务
- P9 可与 P7–P8 并行

---

## 关键设计决策记录

| # | 决策 | 日期 | design note | 要点 |
|---|------|------|-------------|------|
| 1 | TE 节点复制 | 2026-07-06 | roadmap P2 | 避免单值 TE 的 spurious suction ~Γ²/h |
| 2 | damping_theta 替代 pseudo_dt | 2026-07-07 | roadmap P4 | 网格/激波无关的阻尼 |
| 3 | P5 根因：单站 Kutta 闭合失败 | 2026-07-08 | roadmap P5 | 非 TE 奇异性、非 secant 耦合失稳 |
| 4 | Cp 锯齿是 recovery 伪影 | 2026-07-08 | roadmap P6 | 修复在后处理，不动求解器 |
| 5 | Newton 策略 A（完整 Jacobian） | 2026-07-08 | DN4 | 保留 ∂μ/∂φ + 上游耦合 → 严格二次 |
| 6 | IBL3 6 方程 | 2026-07-09 | DN2 §2 | 替代 Green's 3-eq；表面 Galerkin FE |
| 7 | Transpiration BC（方案 VI-1） | 2026-07-09 | DN2 §7 | 不动网格，只改 BC |
| 8 | 直尾流 + 质量穿透，不做几何松弛 | 2026-07-10 | DN2 §4.5.6 | 5 条理由；B6 搁置 |
| 9 | Level-set 仅用于尾流面 | 2026-07-09 | DN2 §6 | UP3D 是 body-fitted，物面不用 level-set |

---

## 更新历史

| 日期 | 内容 |
|------|------|
| 2026-07-10 | 创建 PLAN.md，整合 roadmap 进度与 design notes 理论决策 |
| 2026-07-10 | Track B/V 并入 roadmap.md（权威）+ design.md §4/§5/§11 指针；本文件修正：交叉引用路径（design_notes→同目录）、V3 验证方向（VII 使 CL 下降，0.245→0.288 是 P9 离散差距）、依赖图 "P8 依赖 B3"→可选接口、测试基线补 4 skipped、新增 B3/P8 时序护栏与 V1 工作量提示、V1–V4 与验证阶梯 V0–V6 命名消歧 |
| 2026-07-10 | **P7 关闭**：frozen-walk ∂ρ̃/∂φ（`rho_tilde_sensitivities_sweep`）+ G7.3 FD 验证（cube 3–5e-10 / 构造场 3.5e-9 / 收敛 G4.1 场 5.7e-9，gate 1e-6）；`max_ε` 未实现（G7.1/G7.2 by construction）；design.md §6.3 dμ/dM² 符号 FD 裁决为 "+"；max 折点 FD 陷阱记录（可分离场平局）；基线 165+4+2 |
