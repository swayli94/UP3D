# Track V 粘性耦合——设计、现状与边界

> **本文件是 Track V 的纲领**：设计要点、代码地图、实测能力与边界、开放项。
> 它不是档案。完整的设计决策与逐轮执行记录在归档原件
> [phases/p1/docs/design_track_v.md](../phases/p1/docs/design_track_v.md)（891 行，
> phase 1 时期，叙事止于 GV5.6/GV6.2，2026-07-25）；`pyfp3d/viscous/` 源码注释里
> 的 § 引用全部指向该原件。现行能力声明（本文件引用的每个门号的当前读数）以
> [phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md](../phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md)
> （2026-10-11）为准，冲突时它胜出。
>
> 装订参考：M. Drela, *Three-Dimensional Integral Boundary Layer Formulation for
> General Configurations*, AIAA 2013-2437（下称 D13，方程号用论文原号；本地副本
> `reference/`，该目录被 gitignore，干净克隆需自取）。

---

## 1. 定位

Track V 交付 `pyfp3d/viscous/`：Drela IBL3 **6 方程积分边界层**，壁面三角网上
Galerkin P1 面元 FE（**不做流线积分**），经 **transpiration 边界条件**与全势流
（FP）求解器耦合——边界层的位移效应建模为壁面质量通量吹气分布 ṁ：不动网格、
只改右端项；ṁ ≡ 0 时与无粘路径**逐位一致**（GV2.1(b) 纪律）。

与无粘内核的关系：无粘 FP 内核是已定量闭合的部分（对 CFL3D Euler 同层级参照，
M1 靶子工况 cl 差 −2.16 %，能力边界 §4）；Track V 在其上加位移厚度反馈。两条
通道：**物理通道** δ* → ṁ → RHS（`transpiration.py`）；**数据通道** u_e 恢复
（A4 纪律：medium 光滑壁 ≈2.5 % 峰值相对误差、前缘/驻点带 4–7 %，
`phases/p1/cases/analysis/a4_ue_error_band/VERDICT.md`）。尾迹 δ* 修正（V6）
骑同一通道。

范围之外（能力边界 §7.6）：本库没有 RANS/LES、没有非定常；全速势是等熵+无旋，
激波熵产生是模型形式误差，不是粘性耦合能修的。

## 2. 设计概要

### 2.1 方程与闭包

每面节点 6 个主未知量（D13 §III.A）：δ（厚度尺度）、A/B（顺流/横流壁面剖面
斜率）、Ψ（剖面扭转 crossover）、Cτ1/Cτ2（外层应力尺度）。6 个积分方程（稳态
= D13 (24)(26)(28)(29)(31)：x-动量、z-动量、动能、横向曲率、两条剪应力输运），
每条残差在其节点的**局部笛卡尔基**内构造（面内旋转不变 ⇒ 后缘 kink 无需特殊
方程）。x-动量在 2-D 极限即 von Kármán 动量积分方程。剖面族：层流修正
Bernstein 多项式；湍流 Spalding 壁面律 + Coles 尾迹；密度 Crocco–Busemann、
黏度 Sutherland（可压路径保留，不可压 R≡1）。

闭包（`closures.py`）：每节点输入 6 未知量 + 边数据 (q, ρ, μ, M)，输出 30 个
量（积分厚度、Cf/耗散系数、应力输运因子、派生量）+ 对 6 未知量的**解析**
(30,6) 导数——η 积分号下解析微分，与取值共用同一求积引擎（FD 可验至求积
精度）。层流/湍流分支按节点 regime 标志切换；**强制转捩**（x ≥ x_tr 切湍流分支，
湍流侧 Cτ 以平衡值播种）。自由转捩（e^N）是记录的 follow-up，**从未实现、
从未设门**——D13 门只用 xtr005 工况即因此（见 §4）。

D13 欠定义处的实现决策（决策表当前形态，全文与逐条出处在归档原件 §3–§4）：

- **D-CT**（剪应力输运闭包）：K̄/P/D 全部由剖面族积分出来，不引入新拟合式；
  唯一自由标定 = 耗散长度 L/δ（经典 Bradshaw 外区值 0.09），标定约束 = 2-D
  约化回到 XFOIL/MSES 的 Cτ lag 形式。
- **D-HB**（数值扩散）：三角元各向同性 h̄ = h·I（h = √(2A)）+ **顺流张量**
  ε_s·max(q)·h̄·(s1·∇G)(s1·∇N)（GV1.1(e) 暴露 2h 网格模态后落地；标定膝点
  ε_s = 0.02，求解器默认 eps_diff_s = 0.02，0.0 回收旧格式）。
- **D-QUAD**（求积）：单元 3 点边中点规则；η 向闭包积分**已由 Gauss–Legendre
  换成 θ 换元** η = (1−cos θ)/2——2026-08-22 F 系列勘误：原规则的
  「多项式精确、单测定」两半都被实测为假（`(4η(1−η))^{3/2}` 端点平方根奇异，
  没有任何点数精确），换元把奇异吸收掉（同一积分 n=16：1.9e-06 → 6.7e-16；
  近壁 δ⁺=1e3、n=24：4.4e-01 → 3.5e-03，**仍不足 ⇒ F6 未关闭**）；层流表
  8 → 24 点，湍流 24 点**故意未改**。判定见
  `phases/p5/docs/dev_phase_five/`（F1b/F3/F6 轮）。
- **D-STRESS-2**：2-D 极限（B=Ψ=Cτ2=0）严格退化为单应力 lag 方程，作结构锁
  单测。
- 审计状态（phase-5 第六轮 F8）：D-CT-1/D-CT-2/D-TR/D-BC/D-PT 五个标识符
  **无定义记录**，按诚实原则留在原地并附说明——引用前先查它是否 resolve
  （`closures.py`/`ibl3.py` 文件头审计块）。

### 2.2 离散与求解

- Galerkin P1 面元：闭包量在**节点**求值后 P1 插值到 Gauss 点（D13 §III.D.2
  同款，Jacobian 因此可解析链式求导）；物理通量散度保持**强形式**；仅扩散项
  分部积分，扩散密度逐方程取其守恒量 (Mx, Mz, e, k◦, kτ1, kτ2)。入流
  Dirichlet（翼型 = 前缘带 x/c ≤ 0.02 逐节点 Blasius 种子），出流/侧缘自然。
- 装配：3 节点贪心着色 + serial-color/prange，bit-deterministic；符号 CSR
  （6×6 块/node 对）一次 + 单元散射映射；解析单元 Jacobian 18×18，FD 验证
  （项目纪律）；V1 规模 scipy spsolve 直解。
- 全局化：backward-Euler 伪时间作用于**物理守恒密度**（不用状态向量——那会
  破坏横流方程对 (B,Ψ,Cτ2) 的齐次性），几何 CFL 递增回收稳态 Newton；
  backtracking merit = 伪时间残差 F_pt（Stage 4 实测修正，归档原件 §9.1）。
- 数据布局：节点表 group-aware（"wall"/"wake" 槽位同构，V6 复用同一布局）；
  `SurfaceMesh.volume_node_of` 主映射钩子，使建于未切壁面的 IBL 面网可从
  cut-mesh 解取 u_e。

### 2.3 Transpiration 与 wake sheet 通道

- ṁ = div_Γ(ρ_e·u_e·δ*)（缺陷质量通量的面散度；2-D 规则 v_n = d(u_eδ*)/ds
  的 3-D 推广）。壁面 RHS = ṁ 的 Galerkin 载荷散射到体网格右端项；**符号约定
  钉死**：ṁ > 0 = 吹气（质量出物体），弱形式右端为其**负**载荷（GV2.1(a) 对
  解析 Fourier 圆柱解钉定）。
- u_e 恢复按 A4 分区：前缘带 linear+crease-gated smoothed，其余 quadratic
  lstsq（`transpiration.py::edge_velocity_per_zone`）。
- wake sheet（V6，`wake_sheet.py`）：TE 汇流 δ*_TE = δ*_upper + δ*_lower（θ 同
  理），随后直尾迹指数松弛 δ*_wake(s) = θ_TE + (δ*_TE − θ_TE)·exp(−s/L_rel)
  （L_rel = 钉死的模型选择，灵敏度是登记项）；**每面拷半份**散射——Tᵀ 折叠把
  两侧载荷加进主行，半份使折叠后的跳量等于源面强度本身。flag 默认 OFF =
  legacy 逐位一致。

### 2.4 三条求解路径（3-D 耦合）

1. **松耦合**（V3，`coupling.py`）：FP → u_e → IBL3 → δ* → ṁ → RHS → FP 的
   Picard 外循环，δ* 欠松弛；IBL3Solver 每外迭代重建并暖启动状态。现状见 §4
   ——亚声速翼型收敛，跨声速反馈太弱，举力机翼上跑不起来。
2. **准同时耦合**（V4）：2026-07-22 用户决策**跳过**（skip 判据按字面满足；
   闭体尾部增益 >1 证据是它的存在理由，重开触发自当日起**挂起**）。其思想
   （局部 BL 响应）后来在 GV5.6 以预条件代数形式回归并被实测判死（见 §4）。
3. **紧耦合 = 增广 Newton**（V5，`tight.py` + `tight_driver.py`）：状态
   x = (φ_free, Γ, U)，残差三块（F_φ 带活 transpiration 链、F_Γ 未消元 Kutta
   行、F_BL 稳态 IBL）；全部耦合块**精确装配**（J_φ,BL = −Tᵀ·W·S·P·L·D、
   J_BL,φ = J_e·D_ue·G、J_φφ 增补经 ρ_e·u_e 链）；splu 直解 + 可选行/列均衡
   与 LM 阻尼（三旗标默认全关 = legacy 逐位）；P8/P14 回溯 + 探针守卫沿用。
   交付时 FD 干净（全系统 2.7e-9 @ k=1 态），但撞上了 IBL 地板（见 §4）。

另有二维条带路径（GS4.1，`strip2d.py` + `closures_2d.py`）：沿流向逐站行进、
不组装全局系统，与 3-D 系统**不共享离散、不共享闭包**——它是独立的参照臂，
不是 3-D 求解路径的替代品。

## 3. 代码地图（`pyfp3d/viscous/`，11 个模块约 7500 行）

| 模块 | 一句话职责 |
|---|---|
| `closures.py` | D13 闭包：剖面族、积分厚度、Cf/耗散系数、剪应力输运；30 输出 + 解析导数。**3-D IBL 的唯一权威**（冻结研究分支）。 |
| `ibl3.py` | 3-D IBL 求解核：SurfaceMesh 上 Galerkin P1 装配、解析 Jacobian、守恒密度伪时间全局化 Newton。 |
| `surface_mesh.py` | 壁面紧凑面网：DOF 编号、P1 几何表（面积/法向/梯度/局部基/求积/着色/符号 CSR）；group-aware 节点表 + `volume_node_of` 钩子。 |
| `transpiration.py` | V2 通道：δ*→ṁ 面散度算子、壁面 RHS 装配（符号钉死、ṁ≡0 逐位为零）、A4 分区 u_e 恢复。 |
| `coupling.py` | 松耦合驱动（纯编排）：外循环、钉扎带与 ṁ 掩蔽、翼型/机翼 case 构建（`build_airfoil_case`/`build_wing_case`）。 |
| `tight.py` | 增广系统的固定/线性稀疏算子：J_φ,BL、G（φ→u_e 恢复图）、D_ue（边数据导数）。 |
| `tight_driver.py` | 增广 (φ_free, Γ, U) Newton 驱动：未消元块系统、splu + 可选均衡/阻尼、`step_solve` 注入点。 |
| `wake_sheet.py` | V6 保形尾迹面 δ* 源：TE 汇流 + 指数松弛、半份双面散射；默认 OFF。 |
| `strip2d.py` | GS4.1 二维条带核：逐站行进；闭包**逐字复用** `closures.py`（G-REUSE），自身不含任何闭包公式。 |
| `closures_2d.py` | GS4.1 二维 Drela-Giles 拟合相关式闭包（状态 (θ,H)；层流 + 湍流 + 剪应力滞后，常数带 XFOIL file:line 出处）。**条带核的唯一权威**；与 `closures.py` 互不得 import（G-AUTHORITY 有套件锁）。 |
| `__init__.py` | 空。 |

测试分布（截至 2026-08-28 侦察，`20260828-0100-ibl-xfoil-recon.md` §0）：
18 个活测试碰 `viscous/`，全在 `tests/A`，测**接线与数值机理**（闭包公式、
Jacobian FD 精确性、ṁ=0 严格为零、站点表闭环）；**答案锚在 D 层**（phase 6
建的 D08/D09/D10/D13，见 §4）——「机器接对了」不是「答案对了」（R23 教训）。

## 4. 已测能力与边界（全部实测，带门号）

### 4.1 phase 1：V1–V6（2026-07-22 → 07-25，逐轮记录见归档原件 §9–§21）

- **V1 standalone 核**：GV1.1 **9 PASS / 2 FAIL**。两条 FAIL 都在 band (a)
  （Blasius 平板）：H 带 3.76 %（带 ≤2 %）、δ* 指数 0.5288（带 [0.48, 0.52]），
  同一根因 = **闭包族自洽不动点 H\* ≈ 2.7083 ≠ Blasius 2.59**——FE 与闭包
  自身的 2-D march 吻合 2.2e-4 ⇒ 是闭包族性质，不是离散误差。门(b)(c)(d)(e)
  全过（含准 2-D 锁、加密阶 ≈1.0）。
  （`phases/p2/bench/studies/v1_ibl3_standalone/VERDICT.md`）
- **V2 transpiration 通道**：GV2.1 **23 PASS / 0 FAIL**（符号约定、ṁ≡0 逐位、
  三驱动贯通）。
- **V3 松耦合**：GV3.1/3.2 2P/4F/23R、GV3.3 0P/2F/7R。翼型侧成立——NACA0012
  M0.5/α2° medium 5 次外迭代（ω=1.0）收敛到 ‖Δδ*‖/‖δ*‖ < 1e-3；跨声速
  M0.72 记录点（Newton 驱动）4 次迭代、cl 0.3764。
  闭体侧 GV3.3 FAIL：尾锥 BL 特征线汇到尾极点、无自然出流；钉扎+ṁ 掩蔽+
  FP 护栏后数值干净，但残余不稳定是尾锥**真实**反馈（增益 >1，Veldman 类）
  ——这是 V4 的存在理由与 V5 的动机证据。
- **V4**：跳过（用户决策；重开触发挂起至今）。
- **V5 紧耦合**：增广 Newton 交付且 FD 干净，但四轮窗口实测（GV5.1b/1c/1d）
  证明：**IBL 稳态块有公式层地板**（coarse 3.15e-6 / medium 1.71e-6；残差住
  TE 带 x/c 0.96–1.00 的 (B,δ) 方程；原始条件数墙 4e10 大半是缩放 artifact，
  行/列均衡后余 1e5–1e7 的 (A,Ψ) 真刚度），地板**之上和紧邻处都没有二次
  盆地**（中程还有一层下降屏障）。GV5.5 的 TE 出流行替换**不动地板反而更差**
  （m2 = 5554×/245998× 地板；峰值原系统残差不在 TE 而在 LE 吸力区 x/c≈0.027），
  flag 保持默认 OFF。应用侧：**GV5.2**（RAE2822 跨声速松环）band (b) FAIL——
  每个算出的激波都在实验带下游 0.06–0.10 c，4 腿中 3 腿 ≤10 外迭代不收敛；
  **GV5.3**（M6 Cp）honest FAIL——松环粘性修正几乎不动上翼面 Cp（pooled RMS
  0.1288 → 0.1299 反升），Δcl −2.20 % 方向物理正确但在 A4 输入地板之下
  （input-limited）。结论（登记原文）：**跨声速 VII 的下一次读数属于
  tight/augmented 耦合，不再调松环**。成本侧：**GV5.4** 增广步 7.53× 无粘步
  （RECORDED，参考带 ≤~2×），块预条件 NOT-WORKING；「measure before Schur」
  答案 = J_BL,BL 直接消元只要 1.8 s，瓶颈在 Krylov 收敛不在分解。**GV5.6**
  Schur-aware corrected-AMG：coarse 工作、medium 灾难性劣化（rel_res 0.664）
  ⇒ 该形式在翼尺度**实测死刑**。**V5 CLOSED 2026-07-25**。
  （幸存证据：`phases/p2/bench/studies/v5_*/`、`phases/p1/cases/analysis/v5_6_schur_prec/`）
- **V6 wake sheet**：GV6.1 6 PASS（机制）；GV6.2 实测效应 Δcl +0.00015、TE 区
  max|ΔCp| 0.00250——**显著性条件未满足** ⇒ 求解版尾迹 IBL（producer ii）
  **未开**；V6 CLOSED 2026-07-25。

### 4.2 phase 4：GS4.1 二维 strip IBL（2026-08-18 → 08-21，17 轮，
记录见 `phases/p4/docs/dev_phase_four/`）

- 条带核（`strip2d.py`，复用 `closures.py`）对 Blasius **V-FAIL**：H +4.52 %
  带内、`c_f`√Re_x +6.94 % 带外；归因到闭包族而非离散（自收敛阶 4.04/3.86
  = RK4 设计阶，离散误差 1.3e-13 比模型偏置低 11 个数量级；代数不动点
  2.708292 与 GV1.1 三种独立构造五位吻合）。
- 相关式闭包路线（`closures_2d.py`）：Blasius 按构造是不动点（−0.026 % vs
  型线族 +4.52 % = 175.6×）；成本 0.5040 s（全分辨率阶梯 < 1 s；生产级 BL
  从无粘成本的 607 % 降到 0.33 %，= 3-D IBL 的 1/1815）。注意：过 Blasius 是
  **转录检查**（拟合对自己的训练数据），不是模型验证。
- 湍流闭包对着 XFOIL 源码逐行转录（常数带 file:line；五处整项缺失补回；
  `GCC/Re_θ` 修错）；**滞后方程 L-LAG PASS**（对 XFOIL 自己的 CT 松弛曲线，
  41/43 站，弛豫段中位 0.48 % vs 平衡臂 1.36 %）；转捩种子链对 XFOIL
  **0.12 %/0.15 %**（G20，XT 落站技巧，粘性套件第一条外部期望值锁）。
  转捩区差距分解定案：XFOIL 自身一步离散 ≈41 % + 层流闭包族差 ≈59 %——
  都被理解、都是刻意 ⇒ **不是缺陷**。
- 收口裁决（2026-08-21，使用者，逐字）：「**我们不引入边界层剖面的实验参照，
  以 XFOIL 为准，边界层修正的目的也是让 Cp 正确，独自追求边界层正确没有
  意义。**」⇒ 可执行规矩：**一条边界层量与 XFOIL 的差异，只有在被证明会移动
  Cp 时才算债务**（举证责任在重开方）；13 个开口处置到 4 个；零判据放宽、
  不撤销任何 FAIL。

### 4.3 phase 6：D 类门（现行能力声明；证据 `cases/gates/`）

- **D08**（FP+IBL vs CFL3D RANS，NACA0012，SST–SA 分歧为带）：**没有一条腿
  落在带内**。BLAS=1 档读数：xtr005 medium cl 0.2833（**+10.3 %**，带上沿外）、
  xtr030 medium 0.2821（**+6.1 %**）；加密使 cl 远离带。而 Cp 几乎贴带
  （medium 对 SST 的 RMS = 带宽的 0.99×/0.91×）——分歧是一个**铺开的小偏置
  积分出的 cl 差**，不是局部大偏差。记录方向；带内包含是将来的判据。
  （参考侧注意：xtr030 的 SA 边自身收敛性 unstable，但带宽 1.8–2.7 % 远大于
  任何网格差。）
- **D09**（RAE2822）：**能力边界——耦合环跑不起来**。4 腿全部死在**第 0 次
  外迭代**：`mach_squared_field` 返回负值（9/206 站位 M²<0，min −63.90；
  max |u_e| = 4.86 倍来流；全在上表面激波上下游两处）→ NaN 进闭包崩。
  FP 解本身健康（D06 同网格：收敛、激波 0.6049）⇒ 定位到
  `run_loose_coupling` 的**站位边缘速度提取**，**不归因**；「分侧错误」假设
  已被测量否掉（`coupling.py:334` 的全局 y 基准是未被触发的隐患，单独登记）。
  参考侧完好待用（|SST−SA| = cl 上 9.3 %/11.0 %）。
- **D10**（ONERA M6）：**能力边界——三维粘性耦合在举力机翼上跑不起来**。
  coarse 死在**第 1 次外迭代**（FP driver 不收敛，mdot_max = 6.2e-2），与
  GV5.0 机制一致（根部上翼面 TE 分离斑驱动 δ*↔ṁ↔u_e 失控）。参考侧完好
  （|SST−SA| = cl 3.07 %，是参考自身网格差的 9.5 倍——主导不确定度是湍流
  模型不是网格）。
- **D13**（松耦合 IBL vs XFOIL 6.99，NACA0012 M0.5/Re3e6/α2°/xtr005）：
  (a) 上游 δ*（x/c ≤ 0.30）**只门上表面**（实测比 0.917–1.285，带
  [0.80, 1.40]）；(b) cl 严格夹在无粘 0.2921 与 XFOIL 粘性 0.2691 之间
  （coarse 0.2753 = +2.3 %、medium 0.2833 = +5.3 %）。锁成 RECORDED 的结构
  缺陷：(c) **δ* 比值沿弦单调下降到尾缘 ~0.5**（下游发育不足）；(d) **加密使
  (b)(c) 同向变差**——δ* 分布本身两级几乎重合 ⇒ 不是网格分辨率问题，是
  **耦合强度随网格变化**；(e) **两面不对称度不足**（XFOIL 上/下 δ* 比
  1.26–1.49，我们 0.90–1.12，前半弦甚至倒过来）⇒ 耦合对两面压力梯度差异
  响应不足。(c)(e) 与 c_f 偏高同指一个根因：**湍流闭包的 H 族**（摩擦高 ⇒
  H 偏小 ⇒ 边界层过饱满 ⇒ δ* 长得慢；侦察
  `20260828-0100-ibl-xfoil-recon.md` §2–§3）。
- **BLAS 多吸引子与暖态纪律**（能力边界 §0/§5.2、D08 docstring 登记）：
  库运行时把 BLAS 钉成单线程（H30 裁决：结果依赖线程数本身即错误）。
  D08/D13 的松耦合腿**跨 BLAS 档存在多吸引子**（xtr030 medium：cl
  0.2570 → 0.2821 = +9.7 %，rel_to_band 从 −1.6 % 翻到 +6.1 %，外迭代
  7 → 4）；固定 BLAS=1 下两次独立运行**逐位复现**。裁决：① 这些腿的值锁钉在
  BLAS=1 吸引子上（维持现状）；② 多吸引子是**开放研究项**（可能是其他 bug
  或不足引发）。引用这两门的逐点数字**必须带 BLAS 档位**（证据 CSV 有
  `blas_threads` 列）；证据生成/验证必须在暖缓存态。另：D13 medium 与 D08
  xtr005 medium 是**同一求解**，两门两级 cl 逐位一致。

## 5. 已知问题与开放项（下一阶段候选）

以能力边界文档 §9 的登记为准（工程/研究两分），Track V 相关的在册项：

1. **松耦合跨声速反馈太弱**（GV5.2/GV5.3 联合结论，D09/D10 是其现行门形状）：
   位移厚度反馈在 M ≥ 0.725 拉不回激波；登记的数据主张 = 下一次跨声速 VII
   读数走 tight/augmented 耦合。阻碍：增广 Newton 的成本（7.53–7.87× 无粘步）
   与 IBL 地板。
2. **IBL 地板 = 公式层问题**（TE 带 (B,δ) 方程，机制上 = row 0 x-动量 +
   row 2 动能）：GV5.5 升级阶梯（upwind 边界通量变体、闭包正则化）**登记
   未开**；盆地搜寻已穷尽（GV5.1b/1c/1d 三类种子都降不进地板带）。
3. **Schur-aware 预条件**：corrected-AMG 形式翼尺度死刑（GV5.6）；剩余登记
   路线 = Ĉ 内换节点间 (A,Ψ)-structured M_BB（站位条块逆 / ILU 型近似逆）
   或放弃 corrected-AMG 形式，连同 EW-forcing 变体，登记未开。
4. **xtr030 多吸引子**（裁决②，研究项）：松耦合环跨配置多吸引子的成因
   ——可能是其他 bug 或不足引发。
5. **W3.3 ZPG 湍流平板锚**（§9 工程项）：零压梯度湍流平板的锚定。
6. **D09 站位边缘速度提取缺陷**：已定位（`run_loose_coupling` 的站位 u_e
   提取给出 M²<0），未归因；是「有粘耦合在 RAE2822 上跑不起来」的直接
   瓶颈。
7. **F6 未关闭**：湍流近壁 η 求积分辨率（θ 换元后 δ⁺=1e3 仍 3.5e-03）。
8. **e^N 自由转捩**：从 V1 起登记的 follow-up；D13 只用 xtr005 即因无 e^N
   （xtr030 上 XFOIL 的 e^N 自然转捩 x/c = 0.2668 早于强制点 0.30，不可比）。
9. **V4 重开触发挂起**（2026-07-22 起；GV5.x 各门失败均按预注册不触发）。
10. **V6 producer (ii)**（求解版尾迹 IBL）：GV6.2 显著性条件未满足 ⇒ 未开，
    重开条件在册。
11. D13 的 RECORDED 锁 (c)(d)(e)（下游漂移 / 加密变弱 / 不对称度不足）：
    根因指针在湍流闭包 H 族；按 GS4.1 收口裁决，这些差异**被证明会移动 Cp**
    才升格为债务——举证责任在重开方。

## 6. 出处

- **设计与逐轮执行全文**（V1–V6，决策表、方程推导、每个 gate 的叙事）：
  [phases/p1/docs/design_track_v.md](../phases/p1/docs/design_track_v.md)；
  gate 定义与状态：[phases/p1/docs/roadmap/track_v.md](../phases/p1/docs/roadmap/track_v.md)。
  代码注释里的 § 引用指向归档原件。
- **phase 1 门禁幸存证据**：`phases/p2/bench/studies/v1_ibl3_standalone/`、
  `v5_ibl_floor/`、`v5_m6_bridge/`、`v5_3_m6_cp/`、`v5_4_cost/`、`v5_5_te_floor/`、
  `v6_2_measured_effect/`；`phases/p1/cases/analysis/v3_fuselage_smoke/`、
  `v5_6_schur_prec/`、`v6_0_design_adjudication/`。
- **phase 4（GS4.1）**：`phases/p4/docs/dev_phase_four/`（入口 progress.md +
  每轮预注册/判定）。
- **phase 5（F 系列求积勘误）**：`phases/p5/docs/dev_phase_five/`
  （20260822 的 F1b/F3/F6 等轮）。
- **phase 6（D 类门）**：测试 docstring 即记录（`tests/D/test_D08_rans_naca0012.py`、
  `test_D09_rans_rae2822.py`、`test_D10_rans_onera_m6.py`、
  `test_D13_ibl_vs_xfoil.py`）；已提交证据 `cases/gates/D08_rans_naca0012/`、
  `cases/gates/D13_ibl_vs_xfoil/`（D09/D10 是能力边界门，证据形状 = 探针捕获
  的失败方式，记于 docstring）；参照数据 `cases/reference_data/`（ground truth，
  永不修改）；侦察 `phases/p6/docs/dev_phase_six/20260828-0100-ibl-xfoil-recon.md`；
  能力声明
  `phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md`。
