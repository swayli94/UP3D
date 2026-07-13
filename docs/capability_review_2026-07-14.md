# pyFP3D 能力盘点与批判性审查（2026-07-14）

**性质**：横向审查报告（非 gate 文档，不改变任何 roadmap 状态）。
**方法**：五路独立审计交叉验证——① roadmap/agent-rules 状态台账、② design.md + design_track_b.md 数值方法清单、③ demo_report + `cases/demo/` 证据链逐工件核对、④ **直接读代码**建立两路径能力矩阵与孤儿功能排查、⑤ 讨论笔记（git 历史）与 Track V 设计审计。凡代码与文档不符，以代码为准；凡断言无 committed 工件，按项目自己的规则（"A claim without a committed artifact is not evidence"）降级处理。
**注**：`docs/discussion_notes/` 已于 commit `0e4895a`（本报告撰写当天、由并发会话）整目录删除，其内容仅存于 git 历史；本报告引用处均已注明。

---

## 0. 一页结论

**这个仓库目前不是一个统一的求解器，而是两条刻意平行、能力互补但互不打通的求解路径，外加一个完全空白的粘性支柱。**

- **conforming 路径**（节点复制尾迹 + 每站 secant Kutta）拥有全部"求解器实力"：全耦合 Newton、AMG/lagged-LU 预条件、Mach 续接、Γ(z) 涡远场、翼尖 taper、表面 Cp 平滑——3D 跨声速真解（M6 medium 249 s，cl_KJ 0.2692）只有它能给。但它的尾迹是网格实体：α 扫掠要重划网格、尾迹-机身相交处理不了、翼尖自由边规则是它的固有病（P13 一整个 phase 就是治它）。
- **level-set 路径**（Track B）拥有全部"工作流实力"：隐式 Kutta（无 secant、Γ(z) 平滑 11–12×）、wake-free 网格（真正的任意几何入口）、三种远场选项、α 免重网格的机制——但它的线性代数**硬编码 sparse-direct**（M6 fine ~45 万 dof 不可达，无 AMG 逃生口）、LS Newton 是精简版（无预条件选项、无 Mach ramp 包装、freeze 是占位参数、**从未跑过 3D**）、3D 只在 coarse 网格上用 Picard 验证过一次（B7）。
- **Track V（边界层修正）零代码**。全库 grep 无 viscous/IBL/transpiration 命中。设计（Drela IBL3 六方程 + transpiration BC，GV gate 定义齐全）停留在 roadmap/design 文本。
- **网格收敛性**目前只在三处真正成立：MMS、2D 亚声速升力（G9.2）、**3D 亚声速圆帽 M6（G13.3，p=2.31，全仓库唯一的 3D Richardson）**。3D **跨声速**网格收敛在两种翼尖几何上都失败（平帽序列非渐近；圆帽 fine 到不了 M0.84），至今没有 M0.84 的合法外推值。
- 用户目标清单（2D / 钝尾缘 / 多尾迹 / 粘性修正）中有四项**没有任何代码**，不是"没验证"，是"没实现"。

---

## 1. 现在能算什么（已验证的能力包络）

按"有 committed 工件的验证"划界。凡 Picard 跨声速结果，注意 P4/P5 勘误：**Picard 工程收敛态不是离散方程的解**（Newton 残差 2.2e-4 coarse / ~8e-6 M6），真解须 Newton。

| 问题 | 路径 | 验证网格 | 结果与误差 | 工件 |
|---|---|---|---|---|
| 自由流保持（含尾迹 cut） | 两条 | 全部 | 残差 8.8e-14 | `cases/demo/p1_laplace/` |
| Laplace MMS | 两条 | 三级 | 阶 1.96 / 1.94(LS) | 同上；`tests/test_b2_multivalued.py` |
| 不可压/亚声速升力 2.5D NACA0012 | 两条 | coarse–fine | cl 0.28437 ∈ [PG, KT]，−0.33% 距中点；LS 同网格 Γ 差 <1% | `p3_subsonic/`、`b3_levelset_lifting/` |
| 2D 亚声速网格收敛 | conforming | 三级 | 误差 2.71→0.33→0.03%（G9.2，干净） | `p9_grid_discrimination/` |
| 跨声速 2.5D NACA M0.74–0.82 | conforming Picard+Newton | coarse+medium | Newton 真解 coarse M0.80（shock 0.658/cl 0.459）、medium M0.7875（0.674/0.523）；**medium M0.80 无孤立解（FP 折叠）** | `p4_transonic/`、`p8_newton/` |
| 跨声速 2.5D，LS 路径 | LS Picard(+Newton 2.5D) | coarse ✓ / medium ◐ | coarse M0.80 达标（vs 同网格 Newton 真值）；medium 折叠区 Picard 停在 Γ −18.8%，LS Newton 可达机器精度折叠解但升力比 conforming Newton 低 ~13%（未分摊） | `b6_transonic/` |
| ONERA M6 M0.84/α3.06 | conforming Picard | coarse+medium | cl_p 0.2419/0.2453（Picard 质量） | `p5_onera_m6/` |
| ONERA M6 M0.84 Newton 真解 | conforming Newton | medium（249 s） | cl_KJ 0.2692 vs Tranair/KRATOS 0.288，gap 0.019 | `p8_newton/results/g82_m6_medium.csv` |
| ONERA M6 M0.84，LS 路径 | LS Picard | **仅 coarse** | cl_KJ 0.2765(嵌入)/0.2710(wake-free) vs Newton 0.2692（+2.7%/+0.7%）；Γ(z) 平滑 11–12× | `b7_onera_m6/` |
| 3D 亚声速网格收敛（圆帽 M6 M0.5） | conforming Newton | 三级自相似 | **p=2.31，cl_KJ(h→0)=0.2050 —— 全库唯一 3D Richardson** | `p13.../g133rt_richardson.csv` |
| 远场模型 A/B（域尺寸 15–120c） | LS | dual-mesh | vortex 域鲁棒 0.45%；neumann 截断 O(Γ/R)（15c −4%）；freestream 紧域发散 | `b4p5_farfield/` |
| 球面 Cp（光滑曲壁） | conforming | medium | **11.6% vs 2% gate，未达标**（平面壁元变分犯罪；h 加密饱和 ~3.6%） | `p1_laplace/`（strict xfail） |
| 翼身组合体 | — | 网格已交付 | **无任何求解**（B9 未开工） | `cases/meshes/onera_m6_wingbody/` |

**物理包络（design.md §2/§12）**：等熵、无旋、无粘；M∞ 0.3–0.87；激波须弱（局部法向 M≲1.3，M_max>1.35 运行时警告）；FP 非唯一性带 M0.82–0.85 是模型固有的（medium NACA M0.80 无孤立解已实测兑现，roadmap:1000-1016）。不适用：分离、强激波、钝体。

---

## 2. 数值格式清单（含默认/可选/孤儿标注）

### 2.1 空间离散
- Galerkin P1 四面体（= 顶点中心中位对偶 FV），design.md §6。**只有 tet 路径**——"2D"实为单层挤出的 2.5D（`meshgen/extrude.py`），不存在真 2D 核。
- 人工密度迎风（Hafez/Holst flux-biasing）：ρ̃ = ρ − ν·Δℓ·(∂ρ/∂ℓ)_up，ν = C·max(0, 1−M_c²/M²)。**默认 = 多跳有向 walk**（`kernels/upwind.py`）；streamline-Gaussian kernel（`mode="kernel"`）是 opt-in **孤儿**——从未进入任何配方，且其 Newton 敏感度直接 raise（`upwind.py:906`），永远进不了 Newton。
- 限制器：q² 限速 M_cap=3 + ρ̃ floor 0.05。"收敛"判定拒绝任何 limited/floored>0（Newton 侧）。

### 2.2 尾迹与 Kutta（两个不同的数学对象）
- **conforming**：尾迹为网格内部面，i⁺/i⁻ 节点复制（TE 必须复制，P2 证据），[φ]=Γ 主从消元；Kutta = 每站压力相等的 **secant/Aitken 外层**（映射斜率 b≈0.93 故必须 secant）。翼尖自由边单值 ⇒ Γ(tip)=0 拓扑强制，其 Γ_last~√h 尾迹边奇异由 **G13.2 紧支撑展向 taper**（`vanish_smooth`, r_c=0.05b）治愈（代价 cl −1.1~−1.6%）。
- **level-set**（design_track_b.md）：TE duplication + 切元逐节点辅助 DOF；尾迹 LS BC 两分量（g₁ 法向质量 + g₂ 流向对流，展向刻意留空）；Kutta = **B4 非线性 TE 压力相等** |q_u|²=|q_l|²（精确因式分解、冻结 s̄ 线性化、壁面邻接控制体强制——全扇形误差 +11~15%）。Γ 是解模态，无 secant——P5 的 st133 类失稳**结构上不可能**。翼尖 Γ→0 是涌现的（±3e-4），但**片终止环奇异未治愈**（B8 characterized-not-cured：诚实指数 +0.62/+0.37，与 conforming +0.52 同量级同对象；两条约束侧修法均实测死路，治愈需改终止处函数空间）。

### 2.3 求解器
- **Picard**（两条路径）：ρ̃ 外迭代 + `damping_theta` 局部阻尼（conforming 全场 / LS 必须 `damping_scope="supersonic"`——全场阻尼掐死 Γ 解模态，B6 实测）+ Mach ramp（`solve_transonic_lifting` / `solve_multivalued_transonic`）。
- **conforming Newton**（P8，全功能）：frozen-selection 精确 Jacobian（Term1+2 融合 CSR + Term3 active-set COO；∂μ/∂M² 符号勘误已修）、全耦合 (φ_red,Γ) 精确 δΓ 消元（Kutta |F|~1e-16）、precond amg/ilu/direct、lagged-LU（`direct_refactor_every`，真 3D splu 填充 ~100× 的解药）、EW forcing、stall-adaptive freeze/refresh、Mach 续接 + `intermediate_tol` 松弛（fold 区禁用）。**M6 fine 必须 amg+η=1e-8**（direct 是 4h39m/26GB 陷阱）。
- **LS Newton**（`newton_ls.py`，精简版）：残差-Jacobian FD 验证 1.3e-9，折叠区达机器精度；但 **splu 硬编码、无 precond 选项、无 lagged-LU、无 EW、无 Mach ramp 包装器、`freeze=True` 是无实现的保留字段**（`newton_ls.py:102`），且**从未在 3D 上运行**（调用方仅 2.5D 测试；B7 明确 DEFER）。

### 2.4 边界条件
- 固壁：自然弱式 ρ∂φ/∂n=0。**已知病**：平面壁元在光滑曲壁上 Cp ~11.6%（G1.6，边界数据修正整族被 G1.3/G1.4 oracle 排除，仅剩 P11 曲壁元/Option C 路线，未开工）。
- 远场（**两路径不对称**）：conforming 只有 Dirichlet+2D 涡（PG 缩放）+ 可选 `farfield_spanwise_gamma` Γ(z) 锥化（3D 必需，否则 branch-ray 伪影）；**没有 neumann/freestream 选项**。LS 有三选（vortex 默认/neumann/freestream），**没有 Γ(z) 锥化**（3D 用 neumann 绕开，B7 论证结构上不需要）。
- Track V transpiration BC：设计在案（RHS-only，δ*=0 位相同），**未实现**。

### 2.5 后处理
- conforming：`post/surface.py` 壁面梯度恢复 + 可选法向门控平滑（`smooth_passes`，默认 0；锯齿是恢复伪影非通量伪影——P6 推翻旧归因；平滑只用于 Cp 曲线，**用于力反而使 V6 恶化**）。
- LS：`post/surface_ls.py`（D11 逐侧映射，强制——否则 cl=−3.35 垃圾）。**无平滑对应物**（0 处 smooth）。
- ⚠ `element_mach2` 默认 `mixed_plain="side"` 是 B8 定性的 **×5 度量伪影来源**；诚实读法 `mixed_plain="main"` 已实现但默认未翻转——**B6/B7 的 M_max 锁读的仍是伪影度量**（backlog 在案未落地）。

---

## 3. 两路径能力对照矩阵（碎片化的实锤）

✅ 完整 ◐ 部分 ❌ 无。证据为代码位置（第 4 路审计直接核对）。

| 能力 | conforming | level-set | 后果 |
|---|---|---|---|
| 升力 Kutta | ✅ 每站 secant | ✅ 隐式压力相等 | 不同模型，B8 已证不可互相移植 |
| 跨声速 Mach ramp | ✅ Picard + Newton | ◐ 仅 Picard；**Newton 无 ramp** | LS medium 折叠区定量闭合被卡死 |
| Newton | ✅ 全功能 | ◐ 精简版（splu、无 freeze、无 3D） | LS 无法复现 G8.2 类 3D 真解 |
| 预条件/规模 | ✅ amg/ilu/direct+lagged-LU | ❌ 硬编码 spsolve；**fine 不可达** | LS 网格收敛研究只能到 medium |
| 远场选项 | ◐ 仅 Dirichlet+vortex(+Γ(z)) | ✅ vortex/neumann/freestream | 互补空洞，方向相反 |
| Γ(z) 锥化涡远场 | ✅ | ❌（用 neumann 绕开） | |
| 翼尖 taper | ◐ **仅 `solve_newton_lifting` 接受 `tip_taper`；Picard 入口无此参数** | ◐ 两套机械均为阴性结果标本（默认 None） | conforming Picard 3D 仍带未治愈翼尖边 |
| 表面 Cp 平滑 | ✅ opt-in | ❌ 无对应物 | LS 侧若现锯齿无工具 |
| 3D 实跑 | ✅ Picard+Newton（medium） | ◐ 仅 Picard、仅 coarse | |
| α 扫掠免重网格 | ❌（尾迹是网格实体） | ◐ `WakeLevelSet.update_direction` 存在但**仅 B1 测试调用，从未接入任何 solve 流程** | 卖点能力是孤岛 |
| wake-free 任意网格摄入 | ❌ | ✅（M3/M4/M2 家族） | LS 的核心工作流优势 |
| 非升力 Newton 入口 | ❌（`wc` 必传；G10.1 开放） | ❌ | |
| 翼身组合体求解 | ❌ | ❌（仅网格） | B9 未开工 |

**默认值陷阱（同一概念在不同入口默认不同，容易误用）**：
- LS Picard `farfield="vortex"` vs LS Newton `farfield="neumann"`（有 B6 依据但无警告/文档）；
- `solve_multivalued_lifting` 默认 `upwind_c=0.0`（**静默关闭跨声速能力**）vs conforming 默认 1.5；
- `wing3d.py` 默认 `tip_cap="flat"`（已证发散的几何）vs `wingbody.py` 默认 `"round"`。

---

## 4. 网格收敛性：哪里成立、哪里不成立

**成立（有合法三点 Richardson）**
1. MMS（G1.1，阶 1.96）。
2. 2D 亚声速升力（G9.2：2.71→0.33→0.03%）——尖 TE 对 2D 升力**无地板**。
3. **3D 亚声速、圆帽 M6、M0.5**（G13.3 subsonic：p=2.31，cl→0.2050）。这是花掉三个 phase 才挣到的：需同时修好 ① 尾迹翼尖边（G13.2 taper）、② h_far 网格族钳位（M1b，**该缺陷作废了此前所有 M6 三点外推，含 G9.1**）、③ 平帽壁棱（M5 圆帽）。

**不成立（明确失败，非"没跑"）**
1. **3D 跨声速 M0.84：两种几何都失败。** 平帽：fine 收敛（0.2866）但序列 0.2593→0.2652→0.2866 上升非渐近 ⇒ 单点值非外推；圆帽：fine 的 Mach 续接**死在 M=0.75**，不存在 M0.84 fine 态（site=设计尖锐的 tip TE，圆帽放大而非制造它）。⇒ **P9 判定带至今未触发，"0.019 gap = 分辨率" 只有方向没有外推数**。
2. G9.1 原始 M6 序列：fine 非离散解（M_max 7.93 极限环）+ 网格族本身有钳位缺陷，双重作废。
3. fold 区（NACA M≈0.78–0.80）：dcl/dM≈6–10，O(h) 表现为 O(0.01) 马赫移位，**结构上不可比**——只能单网格回归锁。
4. LS 路径：受 spsolve 墙限制，**没有任何 fine 级收敛研究**；其翼尖终止环指数（诚实 +0.62）与 conforming 同量级但未治愈。

---

## 5. 没实现 / 没验证 / 阴性的完整清单

### 5.1 完全没有代码（设计程度各异）
| 项 | 设计程度 | 备注 |
|---|---|---|
| **Track V 粘性/边界层修正（V1–V4 全部）** | 设计完整（IBL3 六方程、transpiration BC、GV gate、松→紧路线；DN2/DN6→roadmap:1942-2010） | `viscous/` 不存在，全库 0 命中。**用户目标的第二支柱整体空白**。预期方向：VII 使 CL **下降**至 ~0.26–0.27——0.245 vs 0.288 的差距是离散精度问题，不能指望粘性去解释 |
| **B9 多尾迹/翼身 LS 求解** | gate 已定义，2026-07-14 解锁为 NEXT | `CutElementMap`/`MultivaluedOperator` 只接受**单个** levelset，无组合机制；`levelset.py:40` 那句 "one instance per wake" 只是 docstring 愿望 |
| **钝尾缘** | 无设计条目 | 全库无支持；NACA/M6 均 foilmod 削尖；López 原文也是削尖的（DN4 修正 1）。**用户目标中此项连 roadmap 都没排** |
| **真 2D 求解** | 无 | 只有 2.5D 单层挤出替身；所有 kernel 只吃 tet |
| **G10.1 非升力 Newton 入口** | gate 开放 | `wc=None` 不被接受 |
| **P11 曲壁元 / G1.6 Option C** | 条件性排期 | 升力理由已被 P9/P13 削弱，仅剩球面 Cp 一案 |
| **P12 backlog**（离散伴随、VII hook、混合单元） | 列表级 | |
| **B10 曲面/自由尾迹（roll-up）** | **搁置**（straight-wake 误差 O(θ²)~0.1%） | 注意：这意味着刚性平面尾迹的尖端物理（G9.1 的根源）**永远不会被模型级治愈**，只被 taper/圆帽工程化压制 |
| LS 路径 GMRES+AMG | B3+ 递延 | fine LS 一切工作的前置 |
| LS Newton 的 freeze/lagged-LU/ramp/3D | freeze 是保留字段 | B6 medium 与 B7 quantitative 闭合的前置 |

### 5.2 实现了但没验证 / 没人用（孤儿功能）
| 功能 | 位置 | 状态 |
|---|---|---|
| `ptc_dtau`（伪瞬态续接） | `newton.py:335` | **零调用、零测试** |
| `upwind_c_post`（López 2.0→1.6 后退火） | `newton.py:816` | 零外部调用——M6 配方从未用它 |
| `pseudo_dt`（retired 全局伪时间） | `picard.py` | 零测试，文档称 fallback 实无覆盖 |
| streamline-Gaussian kernel 通量 | `upwind.py:744` | opt-in，从未进配方；Newton 侧直接 raise，结构上进不了 Newton |
| `newton_ls(freeze=)` | `newton_ls.py:87` | 无实现的占位参数 |
| `WakeLevelSet.update_direction`（α 免重网格） | `levelset.py:100` | 仅 B1 测试调用，未接入任何 solve 工作流 |
| `wall_correction.py` 整模块 | `solve/` | 路线已被 oracle 否决；"kept as reusable" 至今无复用者 |
| `tanh_half` taper、`span_blend`、`te_weld_coo` | `wake.py:184`、`multivalued.py:57/336` | 阴性结果标本（有测试锁位相同，保留合理，但须知它们**不是**能力） |
| `element_mach2(mixed_plain="main")` | `multivalued.py:623` | ~~默认未翻转~~ **2026-07-14 默认已翻转为 "main"**，B6/B7 M_max 已重读（见 §7 第 6 条） |
| `smooth_wall_tangential_gradients` | `surface.py:147` | opt-in、默认 0、无 LS 对应物 |

### 5.3 验证失败 / 开放问题（有工件的阴性）
1. **G1.6 球面 Cp 11.6%**（strict xfail）：光滑曲壁上的变分犯罪；Nitsche/惩罚/边界数据修正/h 加密全部实测排除，唯一活路 = 曲壁元 + gate 重定义（未开工）。**含义：当前求解器对"光滑曲面为主的几何"（如机身主导构型）的壁面 Cp 有 ~10% 级系统误差风险**——这与翼身组合体推进方向直接相关，但 M2/B9 计划里没有对应护栏。
2. **3D 跨声速网格收敛**（§4，两种几何双败）+ 圆帽 fine 的 tip-TE 超限（`run_g133_roundtip_transonic_locate.py` 已定位 20/20 最快单元在尖 tip TE 上）。
3. **B8 LS 片终止奇异**：characterized-not-cured（用户裁决关闭）；gate 内 mechanism probe 一项 `NOT MET`。
4. **B6 medium 折叠区**：LS Picard 停在 −18.8%；LS Newton 折叠解升力 −13% 未分摊（离散差异 vs 缺陷不明）。
5. **B7 只有 coarse**；顶部 Mach 级 |R|~4-6e-6 有界非收敛（gate 只断言 bounded）。
6. **P4/P5 的 Picard "收敛"态非离散解**（勘误在案）——引用这些 gate 数字时必须带 "Picard 质量" 限定。
7. FP 折叠（M0.80 medium 无孤立解）是**模型极限**而非 bug——接近 0.82–0.85 带的算例天然超出可交付范围。

---

## 6. 碎片化诊断（回应"功能局部交叉、没形成统一能力"）

用户的判断**在结构上成立**，且可以说得更精确：

1. **平行宇宙是刻意的隔离策略**（Track B 规则"conforming 求解器数值逐字节不动"），它保护了回归基线，但代价已经显形：§3 矩阵里**没有一条路径拥有全套能力**。要在 2026-07 的今天算"一个任意几何、fine 网格、跨声速、真解级"的 3D 算例，两条路径都做不到——conforming 缺任意几何入口，LS 缺 Newton/AMG/fine。
2. **同一物理概念有 2–3 份互不相认的实现**：翼尖处理（conforming taper / LS 行混合 / LS span-blend，前者是治愈、后两者是阴性标本）；远场（两套互补但不重叠的选项集）；Cp 恢复（有平滑/无平滑）。B8 的教训是这类"移植"不是机械工作——**两路径的尖端奇异根本是不同对象**，所以碎片不是简单合并能消的，但选项集（远场、预条件、平滑）的不对称纯属工程债。
3. **插桩、配方、开关的可达性极不均匀**：`tip_taper` 只有 Newton 入口收（Picard 3D 照样带病）；`intermediate_tol` 只经一个间接路径可达；四个死参数（§5.2 前四行）连测试都没有。
4. **能力孤岛**：α 免重网格（LS 的核心卖点之一）机制存在却从未接进求解流程；kernel 通量做完即弃。

**若要收敛为统一能力，最短路径大致是**（供讨论，非计划）：
a) LS 路径补 GMRES+AMG（解锁 fine 与 3D Newton，B6/B7/B8 三处的共同前置）；
b) LS Newton 补 Mach ramp + freeze（向 conforming Newton 的 N5 机械看齐）；
c) 决定 conforming 路径的终局角色（若 LS 是终局，conforming 的 Γ(z)/taper/平滑不必反向移植，只作真值锚）；
d) 清理死参数与默认值陷阱（`upwind_c=0`、两个 farfield 默认、flat/round 默认），翻转 `element_mach2` 默认并重读 B6/B7 M_max；
e) B9 之前给 G1.6 风险一个护栏（机身 Cp 的系统误差量级评估）。

---

## 7. 文档可信度批判（docs 里哪些不能直接信）

按项目自己的证据标准核查，以下各条引用时须打折：

1. **P1 从未真正闭合**。§0 工作规则要求 medium gate 通过才算 phase 完成，G1.6 从未通过（strict xfail），ledger 也如实记 `☐ in progress`——但 P2–P13 全部越过它标 ✓。这是长期悬置的规则-状态矛盾（经 DP1 改判合法化，但读者容易误读 P1 为已完成）。**裁决 2026-07-14：B9 gate 已新增机身表面 Cp 护栏条款**（G1.6 类光滑壁误差在翼身验证前必须量化记录，roadmap B9 gate 第三项）；P1 账面重组（Option C）未排期。
2. **B8 的 "✓ CLOSED" 是裁决性关闭**：gate 内 `[ ] ❌ mechanism probe — NOT MET` 仍在，closed 的语义是 characterized-not-cured，不是达标。~~agent-rules 段首 "B8 ◐" 与段尾/ledger "✓" 并存~~（已于同日重写消除）。**裁决 2026-07-14：NOT MET 行维持现状**，其旁已加 "kept by design, do not fix" 注记防止未来会话误改。
3. ~~G13.2 强断言与 G13.3 回撤并存~~ **已解决 2026-07-14（用户裁决）**：roadmap（G13.2 条目、G13.3 条目、P11 ledger 行、证据块引言）、agent-rules、STATUS.md 中的 "REFUTED / 0.019 gap 是分辨率" 全部降级为 **"strongly indicated, NOT earned"**（预注册带以 Richardson 外推值触发，两种几何均无合法外推）；demo_report 原本就是弱口径，未动。
4. ~~G9.1/G13.1 文本中的 LS 指数勘误未落地~~ **已解决 2026-07-14**：勘误注记已加在 roadmap（G13.1、P9 两处）与 demo_report（三处）的原断言旁，均指回 B8 re-spec 块（诚实 +0.62 ≈ conforming +0.52；"两路径都发散"的定性结论保留）。
5. ~~M2 census prose-only~~ **已解决 2026-07-14**：census 落成 `tests/test_m2_wingbody.py` 的 7 条新测试（现 20 passed）+ `cases/meshes/onera_m6_wingbody/ls_ingest_census.csv`。实测**逐项证实** f3c7989 散文（coarse 1,415/76/z∈[0.15,1.1963]/交界内侧 0/TE 集与 wall_nodes 口径无关，α=0 口径），并新增 medium 级（29,108/150）。
6. ~~B6/B7 M_max 读伪影度量~~ **已解决 2026-07-14**：`element_mach2` 默认已翻转为 `mixed_plain="main"`（"side" 保留 opt-in，demo 复现脚本显式钉住），B6/B7 M_max 从缓存态重读（`mmax_reread.csv`，side 值复现 committed 到 6 位）——**B7-M1 的 committed M_max 1.453 本身就是伪值单元，诚实 1.392**；M4 与两个 2.5D 态逐位相同，gate 带无需移动。仍挂 backlog：`element_densities` 的 junk-weight 修复。
7. design.md/design_track_b.md 内含大量已标注 superseded/retired/erratum 的段落（§4.1 旧机制、"g₂ 即 Kutta"、Option A 推荐语、§6.3 符号等——第 2 路审计列了 20+ 条），阅读时必须以标注为准，不能按章节字面采信。
8. `docs/discussion_notes/` 已删（0e4895a，用户确认）；CLAUDE.md / agent-rules / roadmap / design.md / design_track_b.md 中指向它的路径链接已于 2026-07-14 同步更新（改为"已删除 + `git show 8aa4aee:docs/discussion_notes/…` 找回"的历史注记；正文中残余的 "DN1 §…" 类字样为内容引证，非路径链接）。

---

## 8. 对照最终目标的差距表

目标（用户表述）：非结构网格通用全速势 + 边界层修正求解器；2D/2.5D/3D 任意几何；level-set 尾迹；网格收敛性；尖/钝/多尾缘。

| 目标要素 | 现状 | 差距定级 |
|---|---|---|
| 非结构全速势核 | conforming 3D Newton 成熟（真解级、249 s/M6 medium） | ✅ 基本到位 |
| **边界层修正** | 设计完整、**零代码** | ⛔ 整支柱空白（Track V V1 起步即是 IBL3 求解器 + 松耦合，前置 P6 已备） |
| 2D | 无真 2D，只有 2.5D 替身 | ⚠ 若 2.5D 替身可接受则无差距，否则需新核 |
| 2.5D | 两条路径均成熟 | ✅ |
| 3D | conforming 全能力；LS 仅 coarse Picard | ◐ LS 侧差 AMG+Newton+fine |
| 任意几何 | 生成器仅 NACA0012/M6/旋成体参数族；wake-free 摄入机制是对的入口，但翼身求解（B9）未开工；光滑曲壁 Cp 有 G1.6 未解病 | ◐→⛔ |
| level-set 尾迹 | 亚声速 medium 级成熟（Γ<1%、隐式 Kutta 结构优势实证）；跨声速 coarse 级；尖端终止 characterized-not-cured；无 fine | ◐ |
| **网格收敛性** | 2D 亚声速 ✓、3D 亚声速圆帽 ✓（p=2.31）；**3D 跨声速双几何失败，无 M0.84 外推值**；LS 侧无 fine 研究 | ◐ 核心指标只完成一半 |
| 尖尾缘 | ✓（2D 无升力地板已证；3D 尖 tip TE 是跨声速 fine 的当前拦路虎） | ✅/⚠ |
| **钝尾缘** | 无代码、无 roadmap 条目 | ⛔ 未排期 |
| **多尾缘/多尾迹** | 无代码（单 levelset 硬约束）；B9 是占位 | ⛔ |

**优先级含义（按目标倒推）**：B9（多尾迹机制 + 翼身求解）与 LS-AMG 是打通"任意几何 × level-set"的咽喉；Track V 是第二支柱的从零起步；钝尾缘目前在所有权威文档之外，若确属目标需先立项；3D 跨声速网格收敛的下一个具体拦路虎是圆帽 fine 的尖 tip-TE 超限（已定位，未治）。

---

*审计方法与各路原始报告：roadmap 台账、design 数值清单、demo 证据链核对、代码矩阵、讨论笔记审计各一份，结论均以 file:line 落地；本文为其批判性综合。发现的两个可立即行动项：① M2 LS 摄入 census 补工件；② `element_mach2` 默认翻转 + B6/B7 M_max 重读（均已在 roadmap backlog，本报告仅提升其可见性）。*
