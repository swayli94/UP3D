# Track V — 粘性/无粘交互（VII）：数值设计参考

> 状态：**设计参考文档**（2026-07-22，V1 开工时建立，随实现推进补充实测记录）。
> 阶段/gate/进度以 [roadmap/track_v.md](../phases/p1/docs/roadmap/track_v.md) 为准，冲突时 roadmap 胜出。
> 绑定参考：M. Drela, *Three-Dimensional Integral Boundary Layer Formulation for
> General Configurations*, AIAA 2013-2437（本仓库
> `docs/references/Drela_2013_IBL3_general_configurations.pdf`，以下简称"D13"，方程号
> 用论文原号）。D13 欠定义处的实现决策逐条记录于 §4。

---

## 1. 定位

Track V 交付 `pyfp3d/viscous/`：Drela IBL3 6 方程积分边界层，壁面三角网上
Galerkin P1 面元 FE（**无流线积分**），经 transpiration BC 与 FP 求解器耦合
（V2/V3/V5；V1 只做 standalone 核）。本文记录方程体系、离散方案、闭包实现决策与
数据布局设计点。

## 2. 方程体系（D13 已通读核实）

### 2.1 主未知量与坐标基

每面节点 6 个主未知量（D13 §III.A）：

| 未知量 | 符号 | 作用 |
|---|---|---|
| υ1 | δ | 厚度尺度 |
| υ2 | A | 顺流剖面形状（壁面流线方向斜率 U′(0)） |
| υ3 | B | 横流剖面形状（W′(0)） |
| υ4 | Ψ | 剖面扭转（crossover） |
| υ5 | Cτ1 | 顺流外层应力尺度 |
| υ6 | Cτ2 | 横流外层应力尺度 |

每条残差在其所属节点 i 的**局部笛卡尔基** (x̂ᵢ, ŷᵢ=法向, ẑᵢ) 内构造（D13 §II.B、
§III.D.1；面内旋转不变 ⇒ TE kink 无需特殊方程）。局部基构造沿用
`post/surface.py` L651 模式：ŷ = 面积加权顶点法向；x̂ = normalize(seed−(seed·ŷ)ŷ)，
seed 取全局 X（|n_x|>0.9 时取 Y）；ẑ = ŷ×x̂。所有向量点积在全局 XYZ 基完成。

剖面基向量（D13 (38)）：ŝ1 = qᵢ/qᵢ（沿 EIF 速度），ŝ2 = (ŝ1×n̂_w)/|…|，n̂ = ŝ2×ŝ1。

### 2.2 6 个积分方程（稳态 + 伪时间全局化）

稳态形式（D13 (24)(26)(28)(29)(31)，∂/∂t=0；伪时间项见 §5.3）：

```
R1 (x-动量):   ∇̃·J_x − uᵢ ∇̃·M − τ_xw = 0
R2 (z-动量):   ∇̃·J_z − wᵢ ∇̃·M − τ_zw = 0
R3 (动能):     ∇̃·E − qᵢ² ∇̃·M − ρ Q·∇̃qᵢ² − 2D = 0
R4 (横向曲率): ∇̃·K◦ + E·∇̃ψᵢ + ½ρ Q×∇̃qᵢ²·ŷ − ρ Q◦·∇̃qᵢ² + D× − 2D◦ = 0
R5/R6 (应力):  ∇̃·K̄_c − S_τc = 0,  c = 1, 2
```

∇̃ 为面内梯度。defect 通量 M、J̄、E、K◦、Q、Q◦ 由积分厚度按 D13 (62) 组装；
壁面剪切 τ_w 与耗散 D、D×、D◦ 按 (63)。积分厚度 δ*₁…θ◦₂ 与系数 Cf1、Cf2、
C_D、C_D×、C_D◦ 按 (60)(61) 由剖面族 + η 向 Gauss 积分"on the fly"求得。

**2-D 一致性**：x-动量在 2-D 极限即 von Kármán 动量积分方程（D13 §II.I (36)(37)）。

### 2.3 剖面族

层流（D13 (42)–(46)，修正 Bernstein 多项式，4 参数）：
U(η)=A(1−0.6(A−3)η³)f₁(η)+f₀(η)，W(η)=B f₂(η)+Ψ f₃(η)，η=y/δ；
S=(1/Re_δ)(μ/μᵢ)dU/dη，T 同理。不可压时 (60)(61) 全部为 η 多项式 ⇒
解析闭式值用于单测交叉校验。

湍流（D13 (47)–(57)）：Spalding 壁面律（u⁺ 由 y⁺_S(u⁺) 反解，单调 ⇒ 标量 Newton）
+ Coles 尾迹 g₀=3η²−2η³；K、Υ 由边界匹配 (53)(54)；q_τ、U_τ、W_τ 按 (55)–(57)。
壁面斜率含义与层流一致：U′(0)=A（黏性单位），Cf1=2R(0)A/Re_δ 两种 regime 同式。

密度 Crocco–Busemann (58)（V1 门不可压 ⇒ R≡1，可压路径保留）；黏度 Sutherland (59)。

### 2.4 转捩

V1 ship **强制转捩**：x≥x_tr 切湍流分支，湍流侧 Cτ 以平衡值播种。
自由转捩（e^N：层流 TS 源 D13 (34)(35) + Cτ 幅值判据 §II.H.7）为记录的 follow-up，
**V1 不 gate**；闭包 API 结构上预留（`transition_mode` 标志）。

## 3. 闭包实现（`viscous/closures.py`）

### 3.1 节点闭包映射

每节点：输入 (δ,A,B,Ψ,Cτ1,Cτ2) + 外部参数 (qᵢ, ρᵢ, μᵢ, Mᵢ, x_tr 状态) →
输出 (60) 全部厚度、(61) 全部系数、派生量（θ₁₁=φ₁₁−δ*₁、H=δ*₁/θ₁₁ 等）、
应力积分（§3.2），及对 6 未知量的**解析导数**（η 积分号下解析微分，
与取值共用同一 Gauss 求积 ⇒ FD 可验至逼近误差）。层流/湍流分支按 §2.4 判据。

### 3.2 剪应力输运闭包（D13 欠定义处 ⇒ 实现决策 D-CT）

D13 给出应力方程结构 (30)–(33) 与剖面 (49)(50)，但 K̄_τ、P_τ、D_τ、L 的具体
闭式未给（仅"L 按 Clauser G-β 轨迹标定"的原则）。**决策 D-CT**：

1. 应力积分全部用同一 η-Gauss 引擎从剖面族**积分出来**，不引入新拟合式
   （τ′ 外层剖面取 (49)(50) 的 Cτ 驱动项，形状 w(η)=4η(1−η) 乘密度 R，D13 (41)）：
   - k_τc = ρq²δ Cτc·a_k，a_k=∫Rw dη（OUT_AK）；
   - K̄_c = ρq³δ Cτc·(ku₁ ŝ₁ + ku₂ ŝ₂)，ku_c=∫Rw·(U,W) dη（OUT_KU1/2，向量通量）；
   - P_τc = ρq³ Cτc·s_pc，s_pc=∫Rw·(U′,W′) dη（OUT_SP1/2，产生，(30) 的 |τ′|∂q/∂y）；
   - D_τc = ρq³ (1/c_l) √Cτmag·Cτc·s_d，s_d=∫(Rw)^{3/2} dη（OUT_SD，耗散 (33)），
     |τ′| 为**矢量幅值** ⇒ √Cτmag，Cτmag=√(Cτ1²+Cτ2²)；L=c_l δ ⇒ δ/L=1/c_l 严格。
2. 唯一自由标定 = 耗散长度 L/δ。标定约束 = **2-D 约化回到 XFOIL/MSES 的 Cτ lag
   形式**（D13 §II.H.5 自述其 2-D 退化目标）：先取经典 Bradshaw 外区值 L/δ≈0.09，
   以 2-D 平板湍流平衡（S_τ=0）给出的 Cτ_eq(H) 与 Cf(Re_θ) 落带情况验收，
   实测记录于 VERDICT；若单常数跨 Re_θ 不足，再按 G-β 原则引入 L̃(H) 并记录。
3. 源项 S_τc = 2a1(P_τc − D_τc)，a1=0.15（D13 (32)）。

GV1.1(b) 的对标基准 = **同一闭包的 2-D ODE marching 参考解**（gate 脚本内生成），
±5% 判的是 FE 实现与其自身 2-D 极限的一致性；与文献经验相关式（Schultz-Grunow/
White）的偏差作为 RECORDED 参考值，不作 pass/fail。

## 4. D13 欠定义处 / 实现决策汇总

| # | 位置 | 决策 |
|---|---|---|
| D-CT | 应力闭包 K̄/P/D/L | §3.2：剖面积分 + 单自由长度 L/δ，2-D 约化标定 |
| D-HB | 数值扩散 h̄（D13 (70)–(72) 对四边形定义） | 三角元取各向同性 h̄=h·I（h=√(2A)，对直角结构网格严格回收格距；量纲同 D13 的 h̄=长度），守恒形式 + 仅扩散项分部积分（D13 (74)），物理通量散度保持强形式 Gauss 点求值；扩散密度逐方程取其守恒量（Mx, Mz, e, k◦, kτ1, kτ2）；V_ε=ε·max_j q_j，ε∈[0.001,0.01] 旋钮，实测值诚实记录；**各向异性顺流张量已实现**（GV1.1(e) 暴露 2h 网格模态后落地原 follow-up）：+ε_s·max(q)·h̄·(s1·∇G)(s1·∇N)，s1=单元平均缘流方向（边数据，非状态 ⇒ Jacobian 一阶精确、FD 双路绿），标定 ε_s=0.02 膝点（(e) 双测度严格降 + H 阶≈1.0 + 阻尼裕度；0.01 为临界值，详见 §9.4 与 VERDICT 修订） |
| D-QUAD | 单元求积 | P1 三角 3 点边中点规则（`solve/wall_correction.py` 先例；P1 形函数梯度单元常数）；η 向闭包积分用 Gauss–Legendre（层流 8 点 / 湍流 24 点，以多项式精确性与 Spalding 分辨率单测定） |
★★★ **勘误（F8，2026-08-22）**：上一行 `D-QUAD` 的记述里，「**层流 8 点 / 湍流 24 点，以多项式精确性与 Spalding 分辨率单测定**」这句话**两半都已被实测为假**，**原文保留不改写**：① **不是多项式** —— `OUT_SD` 的被积函数是 `(R·w)^{3/2}` 而 `w = 4η(1−η)` 两端取零 ⇒ 端点平方根型奇异，**没有任何点数精确**（[第二轮](../phases/p5/docs/dev_phase_five/20260822-0400-f1b-verdict.md)，层流表已由 **8 → 24 点**）；② **"单测定"是假的** —— 全 `tests/` 从无任何断言触及 `ETA_TURB`，唯一的 Spalding 锁测的是**反演**不是**求积**（[第三轮](../phases/p5/docs/dev_phase_five/20260822-0600-f3-verdict.md)；湍流 24 点在每个测到的状态上留 **≥1.7e-03**，点数**故意未改**，见 F6）。★ **本条由 F8 补上** —— F1/F3/F7 的勘误覆盖了库注释，**漏了这份设计记录**，那正是纪律 #11 存在的理由。
★★★ **第二道勘误（F4，2026-08-22）**：上面那行连**规则本身**也已不成立 —— **不再是「η 向 Gauss–Legendre」**。第七轮把它换成 **θ 换元**：`η = (1 − cos θ)/2`、`w = ½ sin θ · W^GL`，因为 `OUT_SD` 的 `(4η(1−η))^{3/2}` 在两端有平方根型奇异 ⇒ **原规则只能代数收敛，没有任何点数精确**；换元把奇异**吸收**掉（同一积分：`n=16` 时 **1.9e-06 → 6.7e-16**）。★ 它同时把节点向两端聚集，对近壁层也有数量级改善（`δ⁺=1e3`、`n=24`：4.4e-01 → 3.5e-03），**但仍不足** ⇒ **F6 未关闭**。见[第七轮判定](../phases/p5/docs/dev_phase_five/20260822-1700-f2f4f6-verdict.md)。★★ **本行原文与上一道勘误都保留不改写** —— 而**这条勘误本身就是 F9 的第一次执行**：第六轮只勘误了「多项式精确性/单测定」那句**断言**，**没有勘误「Gauss–Legendre」这个规则名**，因为那时它还是对的。
| D-PSI | ψᵢ=atan2 分支 | 准 2-D 门 ψ≡const；一般情形节点展开（unwrap）纪律，V1 只在平滑区取证据 |
| D-STRESS-2 | Cτ2 在 (49)(50) 中不出现 | 交叉应力经横向曲率/横流扇区与第 6 方程闭环；2-D 极限（B=Ψ=Cτ2=0）严格退化为单应力 lag 方程，作结构锁单测 |

## 5. 离散与求解（`viscous/ibl3.py`）

### 5.1 Galerkin P1 面元

tent 加权 Wᵢ（D13 (73)）；节点闭包量（厚度/系数/τ_w/D/q1/q2/ρ）在**节点**求值后
P1 插值到 Gauss 点（D13 §III.D.2 同款：剖面求积只在节点做，Jacobian 因此可解析
链式求导）。物理通量散度保持**强形式**（P1 常数梯度 × 节点通量，D13 (74) 同款）；
仅扩散项分部积分：∫V_ε h̄∇̃(ρ_d)·∇Wᵢ，扩散密度 ρ_d 逐方程 = (Mx, Mz, e, k◦, kτ1, kτ2)。
边界线积分仅域边界非零：入流 Dirichlet、侧缘/对称零法向通量（自然）、出流自然。
残差分量基 = 全局 (x̂, ẑ)（V1 门全部平面算例；曲面局部基投影为记录的 follow-up）。

### 5.2 装配 / Jacobian / 线性解

- 3 节点贪心着色（`mesh/coloring.py` 同算法推广），serial-color + prange，
  bit-deterministic；SoA 预计算几何表；热核零分配；`_njit` shim（PYFP3D_NOJIT=1）。
- 符号 CSR（6×6 块/node 对）一次 + `elem_to_csr` 映射；解析单元 Jacobian 18×18
  （链式：闭包导数 → 插值 → 通量代数 → 残差）；FD 验证（项目纪律，B19/B31 惯例）。
- V1 规模（O(10³–10⁴) 面节点）：scipy spsolve；GMRES+ILU 备选记录。

### 5.3 初始化与全局化

Rayleigh 启动（D13 (78)）：δ=4√(νt₀)（t₀ 使 δ≪ 特征尺度），A=2.5，B=Ψ=0，
Cτ 小量（强制湍流段以平衡值播种）。backward-Euler 伪时间作用于**物理守恒密度**
（D13 (70)–(72) 的 ∂/∂t 项取稳态缘流：G=(Mx−u·m, Mz−w·m, e−q²·m, k◦, kτ1, kτ2)，
质量集中、Jacobian 节点块对角；**不**用状态向量伪时间——那会把 z-动量行耦合到
未知量 A，破坏横流方程对 (B,Ψ,Cτ2) 的齐次性，实测导致准 2-D 锁失效），Δt 几何
递增 →∞ 回收稳态 Newton（D13 §IV.B）；稳态门按稳态残差验收，Δt 序列诚实记录。
backtracking merit = 伪时间残差 F_pt = R + w(G−G_old)（步长本就是 F_pt 的 Newton
步，线性模型只对 F_pt 保证下降）；若改用纯稳态残差判接受，伪时间权重非可忽略时
稳步与判据失配——实测 FS 减速分支从近解种子出发每步被拒、CFL 塌到下限停滞
（Stage 4 诊断，2026-07-22 修正）。

## 6. 数据布局设计点（roadmap V1 要求记录）

1. **wake-sheet 未知量预留**：节点表 group-aware（"wall"/"wake" 槽位同构，同 6 方程
   块）；V1 只建 "wall" 组，V6 续用同一布局换 wake 闭包。
2. **master-map 钩子**：`SurfaceMesh` 保留 surface→volume 节点映射
   （`volume_node_of`），使建于**未切**壁面的 IBL 面网可从 cut-mesh（LS）解取 u_e。
3. **单组（"wall"）scope**：wall+fuselage 接缝 = wing-body，出 V1 范围。

## 7. 验证映射（GV1.1，standalone prescribed u_e）

| gate | 对应设计条款 |
|---|---|
| (a) 层流平板→Blasius（H ±2% of 2.59；δ*∝√x） | §2.3 层流族 + §5；inflow = Blasius Dirichlet 播种（x₀>0 避驻点） |
| (b) 湍流平板 Cf(Re_θ) ±5% vs 闭包自身 2-D 参考 | §3.2 D-CT；参考 = 同闭包 2-D ODE march |
| (c) 减速 u_e 分离指示（H 升）预注册带 | §5.3 全局化；不穿越 Goldstein 奇异，指示器在奇异前取证 |
| (d) 准 2-D 锁（B,Ψ,Cτ2≈0） | 结构锁：横流方程对 (B,Ψ,Cτ2) 齐次 ⇒ Newton 保零 |
| (e) 面加密 ×2 误差降、阶数记录 | §3.1 求积精度 + §5.1 |

## 8. 范围纪律

- V1 不触碰 `solve/` 任何现有路径（纯新增包）；backport 检查结论 = N/A。
- u_e 输入误差带（A4：medium ≈2.5% 峰值相对，LE/驻点带 4–7%）是 **V3+** 的对标
  纪律；V1 为 prescribed-u_e 解析输入，VERDICT 仅作前瞻注释引用，不与 V1 数字混合。
- 线程帽 16（含 BLAS/OMP）；证据 = committed artifact（CSV/PNG）。

## 9. 实现记录（Stage 4 诊断与修正，2026-07-22）

门执行中抓到并按证据修正的三处实现问题（1–3），外加一项门暴露的
离散稳定性缺陷（4）；均不改变 §2–§5 的方程与离散设计，只修正/细化
实现与参考码或诚实记录。详见 `bench/studies/v1_ibl3_standalone/VERDICT.md`。

1. **PTC backtracking merit（§5.3 细化）**。`solve()` 的步长接受判据初版用
   纯稳态残差 |R|∞，但 Newton 步是伪时间残差 F_pt = R + w(G−G_old) 的
   Newton 步——伪时间权重非可忽略时线性模型不保证 |R| 下降。实测 FS 减速
   分支（u_e=x^m, m=−0.05）从近解种子（|R|∞≈3.5e-7）出发每一步都被拒、
   CFL 塌到下限停滞。修正：merit 改 F_pt（与设计 §5.3 一致），FS 两支
   14 步收敛；平板各门无回归（`tests/test_v1_ibl3.py` 双路绿）。
2. **2-D 参考 march 的起点瞬移 bug（gate 参考码，非求解器）**。
   `march_2d` 初版把种子状态直接放在第一个记录站（xs[0]=x0+0.1）开始积
   分——方程 x-自治，整条参考轨线因此平移 0.1，首站记录的是未积分的入
   流种子。症状：GV1.1(e) 误差平台（h-无关、eps-无关，eps 扫描
   0.005→0 全等）+ 参考自身 δ* 指数 0.608。修正后参考与 FE 首站四位有
   效数字一致（H=2.6237 vs 2.6236），(e) 恢复判别力。教训：eps_diff 扫
   描先于一切归因——**扩散不是平台成因**（ν=0.005 对结果无可见影响，
   保留 eps_diff=0.005 预注册设置）。
3. **δ* 指数带外 0.0087 的物理归因（RECORDED）**。闭包族自洽不动点
   H*≈2.7083 ≠ Blasius 2.59 ⇒ 全测段被族调整瞬态主导，瞬态中 δ* 增长
   快于 √x；committed 诊断：march 参考同窗拟合 0.5149（自身近带边）、
   FE 下游半段（x>1.2，H 已近平衡）拟合 0.5101（带内），全窗 FE
   0.5288。这是 (a) H 带失败的同一根因，不是离散误差。
4. **(e) 出流条带 2h 网格模态 — 已修复（D-HB 张量 follow-up 落地）**。
   精确落点参考把比较噪声压到 RK4 精度后，(e) 误差分解为两支：入流区
   扩散误差 ∝ε·h（100→200 正常降，阶≈1）；出流区 station 交替变号的
   2h 模态，入流处播种、向下游放大、增长率 ∝1/h（ε=0 时 200×32 上
   −2.4e-6@x=1.0 → −4.4e-5@x=2.1；100×16 全场 ≤4e-6）。D-HB 各向同
   性扩散 ν=ε·q·h̄ 的逐格阻尼 ∝ε/h，在旋钮带 ε∈[0.001,0.01] 内输给
   ∝1/h 增长 ⇒ (e) 首跑判 FAIL。排查排除：交错对角网格无效、FD
   Jacobian 干净、2-D 约化通量恒等式（Jx−uM=ρq²θ、E−q²M=ρq³θ_s）机
   精成立。**修复（同日，用户指示）**：扩散项加顺流张量
   ε_s·max(q)·h̄·(s1·∇G)(s1·∇N)，s1=单元平均缘流方向（边数据 ⇒
   Jacobian 一阶精确，FD 双路绿）；标定扫描 ε_s∈{0.01,0.02,0.03,0.05}
   全部严格降，取膝点 **ε_s=0.02**（0.01 临界、阶 ds≈0.5；0.02 H 阶
   ≈1.0 且误差仅×1.3）。(e) 复测 **PASS**：errH 4.31e-4→2.21e-4→
   1.12e-4（阶 0.96/0.99）、errds 1.02e-5→6.81e-6→4.39e-6（阶
   0.59/0.63）。求解器默认 eps_diff_s=0.02（0.0 回收旧格式）。SUPG
   （一致、零建模误差）留作 V3+ 升级路线——需处理丢二阶导项或非精确
   Jacobian，超出 V1 补丁范围。

其他已记录实现点：闭包安全 floor（DELTA_MIN=1e-8、RE_D_MIN=1e-3、
a2b2 1e-12，floor 区导数严格 0 且 FD 一致）；N_OUT=30（KU1/KU2 应力通
量积分）；湍流播种 A = ½·cf·Re_δ（线性，非平方）；`run.py` 以
sys.path 锚定本 worktree 的 pyfp3d（site-packages editable 安装指向
姊妹 worktree，直接 import 会拿错包）。

## 10. V3 实现记录（松耦合，2026-07-22）

`viscous/coupling.py` 交付松耦合驱动（FP → u_e → IBL3 → δ* → ṁ → RHS →
FP），设计不变，以下为执行中实测抓到的实现/边界条件问题与修正；门
证据见 `bench/studies/v3_loose_coupling/VERDICT.md`（GV3.1/3.2）与
`phases/p1/cases/analysis/v3_fuselage_smoke/VERDICT.md`（GV3.3）。

1. **IBL3 局部基矢修复（`viscous/ibl3.py`）**。GV3.3 准备期发现 3-D 闭
   曲面上缘流/横流分量的局部基矢投影有系统误差，横流泄漏达
   max|B|/max|A| = 25.9、max|Cτ2|/max|Cτ1| = 0.15。修正后降至
   1.8e-4 / 1.6e-3，52/52 V1/V2 回归绿；GV3.3 中段机身横流比实测
   ~1e-6（FAIL 仅尾锥态缺陷，非基矢缺陷）。
2. **入流 Dirichlet 带状钉扎（翼型）**。单站位钉扎（min-q 站）在 α=0
   也会把 Newton 走进近分离奇异盆地（四腿二分定位，
   `docs/temp/v3_case_bisect.py`）；改为前缘带 x/c ≤ 0.02 全部站位逐节
   点 Blasius 种子、首个外迭代冻结（中途离散开关会搅动不动点）——恢
   复 V1 的 x0 线纪律。闭体类似：鼻极点+第一环为入流候选。
3. **cf 归一化坐标系（预注册补遗①）**。XFOIL DUMP 的 cf 列是**自由流
   归一**（`tools/xfoil/Xfoil/src/xoper.f:1970,2178`，
   CF=TAU/(0.5·QINF²)）；本闭包 OUT_CF1 是**当地归一**（D13 (61)）。
   裸比较在减速区产生系统性偏差；补遗规定统一转自由流参考系
   （cf_fs = cf_local·ρe·|ue|² 逐节点转换后站位平均）再进带，裸当地系
   比较降级为 RECORDED。教训：外部参考数据的归一化约定要先读源码再
   进预注册。
4. **闭体尾部边界条件三轮调试（GV3.3 核心）**。闭曲面 BL 特征线汇聚
   到尾极点，无自然出流（翼型有 TE 出流，旋成体没有）：
   (i) 首跑 IBL Newton "Matrix is exactly singular"（极点无出流约束）
   → 闭包 ZeroDivisionError 崩；(ii) 5% 尾带钉扎湍流种子消除奇异性，
   但冻结的肥厚种子 δ* 在尾锥汇聚几何下变成 ṁ 汇（k=1 −1.05 → k=3
   −3.3 → k=4 −5.7），汇加速尾部 u_e（1.29→1.80），正反馈增益 >1，
   k=4 FP 不收敛（φ~2e23）、k=5 爆（q~1e25）——增益 >1 时 Picard 类
   松迭代对任意 ω 不收敛（Veldman 教训，V4 准同时耦合的存在理由）；
   (iii) 改窄钉扎（极点+第一环）后尾锥 BL 自由演化并真实分离，
   Newton 迎头撞上 **Goldstein 分离奇异性**，崩得更早（k≈2）；
   (iv) **最终方案**：尾带（x/L > 0.95）Dirichlet 钉扎挡 Goldstein +
   **钉扎带 ṁ 掩蔽**（冻结种子 δ* 是边界数据不是解，其汇是人工的；
   `transpiration_from_delta_star` 后 `m_surf[pin]=0`，净源不平衡很小
   由 Dirichlet 远场吸收）+ **FP 不收敛护栏**（`dinfo.converged is
   False` → RuntimeError，垃圾不前传）。翼型路径零影响
   （outflow_pin_surf=None）。结果：10/10 外迭代无数值事件，中段机
   身轴对称优秀；但环不收敛——掩蔽去掉人工种子汇后，残余不稳定是
   尾锥**真实**反馈（经线汇聚把方位角 δ* 不对称放大成 ṁ，FP 以锥面
   u_e 畸变返回），ṁ_max k=5→10 ×5.7。实测记录即 V4 动机证据。
5. **松耦合收敛性记录（GV3.2，翼型侧反方证据）**。NACA0012 M0.5/α2°
   medium 5 次外迭代 ω=1.0（coarse 4 次）收敛到 ‖Δδ*‖/‖δ*‖<1e-3；
   跨声速 M0.72 记录点（Newton 驱动）4 次迭代无调参，cl 0.3764，IBL
   残差地板 3.2e-6。V4 跳过判据按字面已满足；与第 4 条的闭体尾部证
   据并列，**2026-07-22 用户决策：V4 跳过**（判据按字面满足；重开触发
   = V5 受挫或闭体粘性算例提前进范围）。
6. **已记录实现点**。钉子带站位在 profiles CSV 标 `pinned` 列、不计
   入任何误差统计（是边界数据不是解）；LE 区记录值因此从 3.52 修正
   为 0.274（预注册补遗②）。δ* 每外迭代做物理 floor（负值计数入
   history）。全套件基线见 `PROJECT_STRUCTURE.md` 页脚。

## 11. V5 实现记录（GV5.0 M6 松耦合桥，2026-07-23）

GV5.0 入口检查（roadmap GV5.0，RECORDED）经 `viscous/coupling.py::
build_wing_case` 把松耦合驱动铺到 3-D 升力翼壁面；门证据与诊断见
`bench/studies/v5_m6_bridge/VERDICT.md`。实现决策记录：

1. **3-D 翼 IBL 边界拓扑**。wake-cut 网格壁面：上下表面在 LE 共节点、
   TE 因切割复制 ⇒ 双侧 TE 线均为自然出流边界边；根部 z=0 截面是开
   边界边，自然零通量即对称条件。入流 = 局部 x/c ≤ 0.02 的 LE 带
   Dirichlet 钉扎（逐节点 Blasius 种子，k=1 冻结——翼型鼻端带纪律按
   剖面局部弦长推广）；转捩按 side_node ±1 各自的 x_tr/c（局部弦）。
2. **翼尖带掩蔽（复用 GV3.3 尾带机制）**。z > 0.95·b_semi（= 生产
   tip_taper 半径 r_c = 0.05·b_semi，B32）整带 Dirichlet 钉扎到逐节
   点 regime 种子，且钉扎带 ṁ 掩蔽——机制上复用
   `CouplingCase.outflow_pin_surf`（钉扎 + ṁ 置零），语义即 GV3.3 尾
   带：冻结种子 δ* 是边界数据不是解，不产生 transpiration。A4 u_e
   恢复区 = LE 带 + 翼尖带（尖缘奇异旁取 linear+smoothed 稳健路）。
   `run_loose_coupling` 走其 stations=None 分支，零改动复用。
3. **驱动 = 生产 M6 Newton 配方 + V2 external_rhs 通道**。
   `solve_newton_lifting`（farfield_spanwise_gamma、precond="direct"、
   pressure Kutta（P14）、tip_taper（B32））；冷启动按 P14 配方
   probe→pressure 链（n_picard_seed=0），外迭代热启动。
4. **桥答案（两档实测 regime）**。coarse：根部上翼面 TE 分离斑块
   （H 4–5.5）δ*↔ṁ↔u_e 增益 > 1，ṁ_max ×12.4 单调增长（GV3.3 尾部/
   Veldman 同类，首次在升力翼上测到）；medium：加密消除斑块（TE 区
   H>3.5 归零），失稳消失但留有界 δ* 极限环（2–12 %/k）不达
   tol_ds 1e-3。两档 FP 侧全程干净（热启动 Newton 2–4 步收敛）。
   ΔCL 双估计量下行（medium −2.4 % 低于 A4 2.5 % 输入带 ⇒ 输入受限）。
   δ*(z) CSV = GV5.3 带预注册喂料。

## 12. V5 实现记录（GV5.1 紧耦合增广 Newton，2026-07-23）

GV5.1（roadmap GV5.1，9 PASS / 1 FAIL / 36 RECORDED）交付精确增广
(φ, Γ, U) Newton；门证据与诊断见 `bench/studies/v5_tight_coupling/`
（PRE_REGISTRATION + Addenda 1–2、VERDICT、summary.csv、
`results/gv5_1_medium_seed_diagnosis.md`）。实现决策记录：

1. **增广系统架构**。状态 x = (φ_free, Γ, U)：φ_free 为 wake-cut 网格
   自由势 DOFs，Γ 为 Kutta 环量，U (n_s,6) 为 IBL 表面态。残差三块：
   F_φ = R_bare(φ, Γ) + Tᵀ·W·S·P·ṁ(δ*(U))（transpiration 装配走 V2
   通道，符号同 GV2.1）；F_Γ = 未消元 [J_ff B; K −I] Kutta 行
   （probe 估计器；列映射 T[:,dir_red]@V_red + G_jump 与 B 共享，
   G_jump = wake 主从指示阵）；F_BL = 稳态 IBL 残差（边数据冻结在
   pack 基态 = pre-registration 语义；对外输出的 δ*(U) 用当前边数据
   闭包，即松环一致输出字段）。
2. **精确雅可比块与求解器**。J_φ,BL = −(Tᵀ W S P L D)[free]
   （L = 表面散度向量算子，D = δ* 闭包行算子）；J_BL,φ = J_e·D_ue·G
   （J_e = IBL 边数据残差导数，G = 分区 u_e 恢复算子：LE 带
   linear+crease-gated smoothed、其余 quadratic lstsq）；J_φφ 增补 =
   dṁ/dφ 经 ρ_e·u_e 链（surface_divergence_vector_operator +
   rhou_jacobian）。线解 = splu 直解（2.5-D 密度下 ~28k² 可负担）
   + P8/P14 safety-only 回溯（max_backtracks = 30；探针守卫 =
   IBL3Solver 的 halving-on-nonfinite 习语——抛出（非物理探针在闭包
   求积内除零）或非有限的探针按 merit = +inf 继续折半，只影响探针
   拒绝，不影响已接受步）。GMRES + 块预条件（AMG-φ / ILU-BL）留
   GV5.4。
3. **边数据链与 closures douts_e 派生栈**。边基 7 标量
   (q, ρ, μ, M, û)（û = 归一化边方向，D13 链路）；closures 扩
   douts_e (30,2) = d/d(re_d, e_prime)，经 8 宽导数栈贯穿 eta 积分
   包。状态列 0–5 与重构前逐位一致（256 态探针实测 outs max diff
   0.0），douts 差 ≤ 3.6e-12（派生栈重排）。veps/veps_s 全局扩散
   尺度在 Newton 步内冻结（pre-registration decision 5），冻结 vs
   自然重算遗漏实测 ≤ 3.0e-8 scaled（FD 门，Recorded）。
4. **FD 门裁决出的三个装配 bug（Stage 2）**。drhom wa 权重
   （ibl3）、s1e 双因子扩散链（ibl3）、dR range(6)→range(7) re_d
   列（closures，污染湍流分支）。修后甜点：Stage-3 全系统 FD
   2.7e-9（k=1 态）；GV5.1 coarse 种子 2.246e-8 / 端点 2.244e-8、
   medium 种子+端点 5.074e-9（tol 1e-5），掩蔽行 0/1236 +
   0/2460。
5. **IBL 地板机理（band (b)/(c) 未达的根因）**。稳态 BL 块在
   k=1 类态附近有内禀条件数墙：cond(J_BL,BL) ~ 4e10，501/1236
   奇异值 < 1e-6·max，42 节点的 Λ/A/B 行原始 dU O(5e2–6e3)；
   standalone 伪时间解在同一地板（~1e-6 残差）100 步不收敛
   （converged=False）⇒ 地板是 IBL 公式/条件数内禀，**非紧耦合
   缺陷**。松环从不需要该态的 BL 块收敛（欠松弛滑过），增广
   Newton 要求三块同时收敛而暴露之。polish 自松环收敛种子：
   F_φ 第 1 步即 1.16e-7（medium），F_BL 自第 0 步钉在松环末态
   自带地板（medium 1.708e-6 / coarse 3.11e-6），lam → 0，无
   斜率-2 尾段（medium p = 0.02/0.50/16.07，coarse
   0.98/3.68/0.57）；N_polish = 10，N_total 14/13 vs 松环 4/5。
6. **松环 medium 不动点不可复现（finding）**。松环 medium 轨迹
   在 IBL 地板上混沌：三个代码/环境组合 → 三个不动点
   （committed cl 0.2719 / δ* 6.84e-3（n_outer 5）、HEAD 重生成
   cl 0.2814 / δ* 3.45e-3（3）、c2dc325 重生成 cl 0.2217 / δ*
   9.73e-3（6））；k=0 无粘基线逐位一致（Δ ≤ 1.3e-9），HEAD
   运行间逐位一致 ⇒ 环境/代码 1e-12 级微扰经 100 步截断的 IBL
   解在近零流形上放大到 O(0.3)。诊断
   `results/gv5_1_medium_seed_diagnosis.md`；用户裁决（2026-07-23）
   接受 HEAD 重生成种子（接线守卫 = converged + |dcl_k0| ≤ 1e-8，
   medium 实测 1.309e-9）。coarse 不动点条件良好（k=1 δ* 差
   0.14 %），其 cross-check 通过。

## 13. V5 实现记录（IBL 地板诊断 = GV5.1 follow-up，2026-07-24）

诊断研究 `bench/studies/v5_ibl_floor/`（预注册 53bf904 先于首次执行；
RECORDED 类、无 pass/fail 带；14 RECORDED；`run.py` 单 runner 从头
再生成全部 artifact，`--states`/`--phases` 支持分相位续跑；执行证据
`results/findings.md` + `results/summary.csv`）。三态：S1/S2 =
coarse/medium 松环收敛态（HEAD 重生成、三次独立再生成逐位一致，
接线守卫 |dcl_k0| ≤ 1e-8 通过），S3 = coarse k=1 fixture。S2 数字
一律带 GV5.1 §4 轨迹散布告诫，结论只立 S1/S2 共同特征。诊断结论：

1. **近零簇在松环收敛态持续存在，由湍流 (A,Ψ) 变量承载**。S1 谱与
   S3 几乎逐点重合（S1 500/1236 < 1e-6·σmax、cond 1.3e11 vs S3
   501/1236、4.1e10；S2 1082/2460、4.0e13）——非 k=1 特异。
   top-20 右奇异向量能量 A+Ψ 占 98–99 %（δ、Cτ ≈ 0），节点支撑
   遍布湍流区、mid-chord 偏重（质量 0.86/0.82）、TE 次峰
   （0.14/0.18）、LE 带机器零。
2. **原始 cond 主要是缩放 artifact**。一次行+列 2-范数均衡后 cond
   降到 2.1e4（S3）/ 7.4e5（S1）/ 1.1e7（S2），亚 1e-6·σmax 计数
   501/500/1082 → 0/0/2——无精确零方向。§12 第 5 条的
   "cond(J_BL,BL) ~ 4e10 内禀条件数墙" 按此重新理解：raw 数字
   无误，但墙体大部分是 pin 行（σ=1）与物理行范数混杂（行/列范数
   动态范围 1e4–1e6）的表象；均衡后余 1e5–1e7 的 (A,Ψ) 真刚度
   才是靶子。
3. **地板残差住 TE 带 (B,δ) 方程，且几乎全在 J 值域内**。F_BL
   支撑 x/c 0.96–1.00（S1 3.154e-6 / S2 1.710e-6 = 松环末态自带
   地板）；方程范数份额 B 0.83 / δ 0.48（A,Ψ ≈ 0.19，Cτ ≈ 0.01）；
   与 top-20 近零左奇异向量对齐 ≤ 7.7e-3（vs 良态向量 1e-3–1e-1）
   ⇒ 非值域亏缺：J "看得见" 残差，但消它要沿 (A,Ψ) 平方向运动，
   而非线性残差在那里不跟随线性模型。
4. **两条候选机理被否**。闭包地板活动集 S1/S2 全空（min δ 高于
   DELTA_MIN 2.2 个 decade、min re_d 高于 RE_D_MIN 3.4 个 decade）
   ⇒ "地板行零导数 ⇒ 精确零方向" 假说不成立，且 DELTA_MIN 灵敏度
   在这些态上恒为零（Q6(b) 探针延期——numba cache 烧死的模块
   常量、无接口暴露——以此空集为替代证据）；eps_diff ×0.5/×2
   （含 eps_diff_s 联动变体）地板仅移 −3 %/+6 %（联合 ×2 最多
   +23 %），近零计数变化 < 4 % ⇒ 非人工粘性截断地板。
5. **伪时间控制器触底 = 公式化地板经控制器表现**。自 S1 重解：
   稳态残差自第 0 迭代冻结于 3.154e-6（偶发的接受步沿近零流形
   移动 U 而残差不动），cfl 1.0 → 1e-3 = cfl_min 钉死，14/21
   步拒绝，n_fail > 10 退出 ⇒ 任意小伪时间步都找不到下降方向，
   单靠全局化过不了地板。

GV5.1b 设计输入（排序待用户裁决）：① equilibration 并入紧解
（J_BL,BL 乃至全增广系统的行/列均衡是廉价先决，直接把 cond 打到
1e5–1e7 量级）；② 靶向缩放后 (A,Ψ) 块刚度的阻尼/投影 Newton
（GV5.4 的 BL 块预条件同样应按 (A,Ψ) 结构组织）；③ band (b) 的
slope-2 窗口重定义为地板之前（地板已由本诊断钉死为公式性质，
不再是收敛判据的一部分）；④ 全局化改动单用无效（Q7 已证）。

## 14. V5 实现记录（GV5.1b scaled+damped Newton，2026-07-24）

门禁 `bench/studies/v5_1b_scaled_newton/`（预注册 8b7793f 先于首次
执行；裁决后 **2 PASS / 0 FAIL / 7 RECORDED**——band (a) medium
cond-aware 读 PASS，2026-07-24 用户裁决（VERDICT §3）；执行时读数
1 PASS / 1 FAIL / 7 RECORDED 保留在 commit 1c55906；`run.py` 单
runner 从头再生成
全部 artifact，协议 = GV5.1 amended 逐字——松环重生成种子 + 接线
守卫 |dcl_k0| ≤ 1e-8 两腿均过（coarse 1.56e-12 / medium 1.31e-9）；
执行证据 `results/summary.csv` + `results/compare.csv` + 三条
newton_history CSV，判决 `VERDICT.md`）。实现决策与结果：

1. **求解器内部实现，装配逐位不动**。`newton_tight` 新增三旗标
   `scaling / lm_damping / floor_stop`：每迭代对装好的全增广 J 做
   一遍行/列 2-范数均衡（零安全）得 R、C，解 (R·J·C + μI) δy =
   −R·F（splu，保稀疏、无法方程），δx = C·δy；μ 走确定性日程
   （μ₀ = 1e-6，线搜索拒绝 ×10 封顶 1e2 并重解，接受 ÷3 封底
   1e-12）；P8/P14 回溯 + 探针守卫逐字沿用；floor-reached 停止类
   = merit 相对下降 < 1e-4 连续 3 个接受步。三旗标默认全关 =
   legacy 路径逐位（对 committed k1seed 历史回归 rel ≤ 2e-6 通过；
   tight 舰队 28 passed 两次，执行前后各一；新测试
   `tests/test_v5_tight_scaled.py` 8 个）。
2. **band (a) 套件精确；medium 活体 e2 = 阈值校准问题，非代数错；
   裁决落地 PASS（2026-07-24 用户）**。良态合成系统上的机器精度
   恒等式 + μ 日程转移全绿。
   活体种子 J 上 e1（对角代数）≤ 2.6e-16 两级均过；e2（μ=0 阻尼
   步 vs 无阻尼 splu 步）medium 1.96e-10 超实现时自设的 ≤1e-10
   前向阈值——两矩阵数学相同，+0.0·I 引入的显式零元改变 SuperLU
   列主序，差 = 经 cond(J) ~ 1e10 放大的舍入（backward-error
   意义下恒等成立，离 cond·eps 上界还有 4 个 decade）。阈值非预
   注册、事后未动；**2026-07-24 用户裁决：cond-aware 读 PASS**
   （VERDICT §3）——e2 容差改为 tol = max(1e-10, 10·κ₁(J)·eps)
   （κ₁ 由一范数估计现算，~1e10 → ~1e-5 量级界，实测 1.96e-10
   以 ~4 个 decade 余量通过）；套件仍是 binding gate，活体检查
   按 cond-aware 容差重发。as-executed 的 1/1/7 读保留在 commit
   1c55906。
3. **band (b) 无窗可读 = 构造性结果，非机构失败**。amended 种子
   本身就是松环末态，F_BL 自第 0 迭代坐在 1.00× 诊断地板
   （coarse 3.154e-6/3.154e-6、medium 1.710e-6/1.712e-6），深在
   预注册 10× 地板带内侧——全程无 above-band 收缩段，median-p
   判读为空，走预注册 fallback（RECORDED）。fallback 读数：
   medium termination = **floor_reached** 第 5 迭代干净收官
   （merit 9.074e-11 ≈ GV5.1 的 9.025e-11），取代 GV5.1 的 10 步
   λ-collapse 爬行；coarse 末 merit 2.044e-10 < GV5.1 的
   2.068e-10 且仍在降（λ_last = 0.031 未塌）；k=1 standalone
   下潜显著更深（F_BL 3.268e-6 vs k1seed 4.726e-6 = −31 %，
   merit 2.25e-10 vs 5.28e-10 = 2.3× 低）——给足房间时缩放
   Newton 确实下潜，但 10 步慢降不是二次尾段。
4. **μ 惰性；缩放是活性成分**。三条 run 上 μ 拒绝重试合计 0 次，
   μ 自 1e-6 单调衰减到 ~5e-11/1e-8——Levenberg 对角臂在这些态
   上从未启用，与 §13 第 5 条"单靠全局化过不了地板"互证；真正
   起作用的是行/列均衡（与 §13 第 2 条"cond 主要是缩放
   artifact"互证）。阻尼代码路径保留（套件覆盖），但在
   above-band 种子出现前不指望它起作用；均衡件是 GV5.4 预条件
   的现成配料。
5. **窗口问题被重构而非回答 → GV5.1c 输入**。要在地板带之上读
   slope-2，需要 F_BL 高于地板带的种子（早松环迭代态或扰动 δ*
   态；k=1 态也只有 1.5× 地板）；在此之前"增广 Newton 在地板
   之前二次收敛" = 未检验而非证伪。破地板本身仍是公式层工作
   （TE 带 (B,δ) 方程，§13 第 3 条），排队待用户裁决。band (c)
   计数：coarse N_polish 10 vs 期望 ≤8 NOT met（记录）；medium
   5 vs ≤10 met（退化：band-entry iter 0）。下一步 = GV5.1c
   （above-band 种子）或 TE 带公式层工作，
   排序 = 用户裁决；V4 重开触发保持挂起。

## 15. V5 实现记录（GV5.1c above-band 种子：地板前 slope-2 窗实测，2026-07-24）

门禁 `bench/studies/v5_1c_above_band_window/`（预注册 1e90d59 先于首次
执行；**2 PASS / 1 FAIL / 7 RECORDED**；`run.py` 单 runner 从头再生成
全部 artifact，协议 = GV5.1 amended 逐字 + 标定 δ 扰动；接线守卫两腿均过
coarse 1.56e-12 / medium 1.54e-9；执行证据 `results/summary.csv` +
`results/compare.csv` + 两条 newton_history CSV + 两条 seed_calibration
CSV，判决 `VERDICT.md`）。**本 session 临时 8 线程约束**（用户定，仅本
session；runner 默认 16 不动，约束经环境变量落地并记入 artifact；壁时
不可与 16 线程账目直接比）。设计决策与结果：

1. **above-band 种子按预注册交付**。种子 = 松环收敛态 + 自由 BL 节点
   δ×(1+ε)（inflow Dirichlet 带不动——pin 行是边界数据）；ε 由确定性
   log10 对分（括号 [1e-8, 1e4]，inf 安全，≤20 次残差评估）标定进 T1
   窗 [5e-2, 5e-1]：两级均 2 次评估即在括号边上窗（ε = 1e4 → 种子
   F_BL 3.219e-1 coarse / 1.819e-1 medium ≈ 地板带的 1e4 倍）。记录在
   案：F_BL 对 δ 缩放的响应饱和（δ×10001 只把 max-范数抬到 ~0.2–0.3）
   ——野态下残差 max-范数不是 δ 扰动的线性函数，未追。
2. **窗口实测：地板之上处处无二次收缩**。干净下降段（coarse 前 3 步 /
   medium 前 2 步）全部 λ = 0.50 封顶，收缩恰 0.30 dex → p = 1.00 为
   构造值（回溯上限，非 Newton 渐近率）；随后轨迹**中程停滞**
   （F_BL ~ 3e-2 → 1.3e-2 / 2.2e-2，λ 塌到 1e-3–0.1，medium 第 6 迭代
   F_BL 反跳 3.06e-2 → 3.68e-2 而 merit 仍降），10 迭代内从未进地板带
   （cap 时仍距地板 4262×/12867×）。binding medium：7 个 above-band
   三元组 median p = 0.56 ∉ [1.5, 2.5] → honest FAIL（coarse 1.00
   记录）；回归斜率 0.75/0.62。停滞段的 p 值（0.00–340）是近零收缩
   比值的噪声放大，非物理。
3. **新发现：中程下降屏障**。紧 Newton 的障碍不止 ~3e-6 处的公式地板
   ——地板之上 3–4 个 decade 还有一层下降屏障（10 迭代只走到
   F_BL ~ 1e-2）。地板紧邻处是否存在局部二次盆 = 近带种子后续问题
   （候选 GV5.1d，待用户裁决）；本gate 不测它（种子按预注册是远种子）。
4. **band (a) PASS 两腿；cond-aware e2 容差本次预注册**（上轮裁决
   落地）：扰动种子 J 上 e1 2.6e-16/2.3e-16、e2 2.06e-9/2.40e-9
   （容差 3.9e-2/5.2e-2，κ₁ 一范数估计现算）、e3 3.99e-10/2.57e-10。
   套件 37 passed 两次（执行前后；tight 舰队 28 + 新 9，
   `tests/test_v5_above_band_seed.py` 合成映射测试）。μ 拒绝重试再次
   为 0（μ 自 1e-6 衰减到 5e-11）——远种子下 Levenberg 臂仍不启用，
   全部全局化由线搜索承担；行/列均衡仍是活性配料。
5. **medium 不动点在 8 线程下再次散布**（第 4 个不动点 cl 0.28245999，
   n_outer 3，未扰动种子 F_BL 1.824e-6 = 1.07× 诊断地板）——GV5.1 §4
   轨迹散布机制（并行归约顺序随线程数变）；接线守卫（|dcl_k0| ≤ 1e-8）
   按设计验证配方而非不动点。coarse 在 8 线程下逐位一致
   （cl 0.26791639，F_BL 3.153842e-6）。band (b) 判读不受影响：种子
   在任意轨迹上都 above-band（构造性质）。
6. **窗口问题状态**：远种子假设已实测为否（不再是"未检验"）；下一步
   = 用户在 GV5.1d（近带种子）/ GV5.5（已登记的 TE 带 (B,δ) 公式层
   独立项，未开工）/ GV5.2 / GV5.3 / GV5.4 间排序；V4 重开触发保持
   挂起（预注册：本 gate 失败不触发）。

## 16. V5 实现记录（GV5.1d 近带种子：地板紧邻处无二次盆，2026-07-24）

门禁 `bench/studies/v5_1d_near_band_window/`（预注册先于首次执行；
**2 PASS / 1 FAIL / 7 RECORDED**；`run.py` 单 runner 从头再生成全部
artifact，协议 = GV5.1c 逐字 + 新近带窗；helper 自 GV5.1c runner IMPORT
不镜像；接线守卫两腿均过 coarse 1.56e-12 / medium 1.54e-9；执行证据
`results/summary.csv` + `results/compare.csv` + 两条 newton_history CSV
+ 两条 seed_calibration CSV，判决 `VERDICT.md`）。**本 session 临时
8 线程约束**（用户定，仅本 session；runner 默认 16 不动，约束经环境
变量落地并记入 artifact；壁时不可与 16 线程账目直接比）。设计决策与
结果：

1. **近带种子按预注册交付**。窗 = T1 [1e-4, 1e-3] 主（medium 带的
   5.8–58× / coarse 带的 3.2–31.6×——带正上方 1–2 个 decade，且低于
   GV5.1c 中程停滞区 ~1e-2）/ T2 [1e-3, 1e-2] 升级（仅在 <3 above-band
   三元组时触发；两级 T1 均 ≥3 三元组 ⇒ T2 未触发）。标定：coarse
   ε = 10 → 种子 F_BL 1.711e-4（5.42× 带），medium ε = 56 → 6.02e-4
   （35× 带），对分 4/6 次评估。
2. **判读：地板紧邻处也不存在二次盆**（band (b)，medium binding
   honest FAIL）。coarse：一次 λ = 0.5 封顶折半（1.7e-4 → 8.7e-5）后
   爬行——λ 塌缩到 6e-5–8e-3，收缩因子 ≤ 0.03 dex/step，10 迭代走到
   7.59e-5 = **地板的 24.07×（带的 2.4×）仍未进带**；median p = 0.35
   （recorded），回归斜率 0.15。medium：第一个接受步把 F_BL **推离
   带**（6.0e-4 → 9.8e-4——merit 1.82e-5 → 1.09e-5 全靠块再平衡
   （F_φ 7.8e-4 → 4.9e-4）换得，BL 块反而增大），随后同样爬行到
   8.43e-4 = **地板的 492.6×**；median p = 1.17 ∉ [1.5, 2.5] → FAIL，
   回归斜率 0.88。两腿 termination 均为 cap；band-entry iter = none
   （本 gate 提升为关键 datum——GV5.1c 从上方从未进带，GV5.1d 从近带
   也未进带）。
3. **物理读法**：GV5.1c 的"中程屏障+下方或有盆地"图像被否定——平坦/
   锯齿 merit 邻域从 1e4× 地板一路延伸**向下到距地板 ~1.5 个 decade
   内**（coarse 24× 处仍是停滞区）。与诊断自洽：缩放后 (A, Ψ) 刚度
   1e5–1e7（findings Q2）意味着沿 Newton 方向的 merit 极平，线搜索
   只能走 1e-3–1e-4 相对步——"收敛"的是线搜索不是模型。近带处
   Newton 方向甚至不是 BL 下降方向（medium 第一步）。p 序列的
   0.06…8.45…128 散布判为平台上的舍入噪声，非收缩regime。
4. **band (a) PASS 两腿**（cond-aware e2 容差沿用上轮裁决）：扰动种子
   J 上 e1 2.6e-16/2.3e-16、e2 1.54e-11/1.20e-11（容差 9.2e-2/9.7e-2，
   ~12 个 decade 余量）、e3 6.7e-12/7.1e-12。套件 49 passed 两次
   （执行前后；tight 舰队 33 + 9 + 新 7 `tests/test_v5_near_band_seed.py`
   合成映射测试）。μ 拒绝重试**第三次为 0**（μ 1e-6 → 5e-11）——
   Levenberg 臂在三类种子（amended / above-band / near-band）下均不
   启用；行/列均衡仍是唯一活性配料。
5. **medium 不动点**：8 线程下仍是第 4 个不动点（cl 0.28245999，
   n_outer 3，未扰动 F_BL 1.824e-6 = 1.07× 地板）——与 GV5.1c 同一
   散布；coarse 逐位一致（cl 0.26791639）。近带判读引用 committed
   地板带，不受散布影响（窗以带为单位）。
6. **程序状态**：盆地搜寻穷尽（GV5.1b 全局化 / GV5.1c 远种子 /
   GV5.1d 近种子全部无法下降最后 ~1.5 个 decade）——地板及其平坦
   邻域是同一障碍，且是公式层的：**GV5.5（TE 带 (B,δ) 公式层独立项，
   机制上 = 诊断行命名的 row 0 x-动量 + row 2 动能）成为破地板唯一
   在册路线**；GV5.4 块预条件解决成本不解决地板。V4 重开触发保持
   挂起（预注册：本 gate 失败不触发）。下一步 = 用户在 GV5.5 /
   GV5.2 / GV5.3 / GV5.4 间排序。

## §17 GV5.5 —— TE 带 (B,δ) 公式层破地板：V1 TE 出流行替换执行记录（2026-07-24，2 PASS / 1 FAIL / 9 RECORDED）

独立项（2026-07-24 用户定序：GV5.1d → GV5.5 → GV5.2–5.4）。预注册
`bench/studies/v5_5_te_floor/PRE_REGISTRATION.md` 先于首行代码提交；
路线选择在开工时按登记落定：**路线 (a) TE 自然出流离散先行**，行级
变体 **V1 = TE 出流行替换**（δ 载体行 6i+0 `R = δ_i − δ_up`、H 载体行
6i+2 `R = H_i − H_up` 一阶外推，精确雅可比行，CSR pattern 内构造守卫，
默认 OFF flag `te_extrapolate`，`te_outflow_pairs` 由 case 层提供冻结
数据）。binding 指标 = **m2**（V1 终态上的**原系统**残差，防行替换
作弊）；m1 = 变体系统自身地板（记录）。全程 8 线程临时约束。

1. **V0 控制**：coarse 逐位接近 committed 地板（3.1537e-6 vs 3.154e-6，
   control_rel 1.06e-4 ≤ 1%）；medium 8 线程散布条款**第三次触发**——
   同一第 4 不动点（cl 0.28245999），floor_ref 改用种子自身 flag-OFF
   地板 1.8238e-6（与 GV5.1c/1d 的 1.824e-6 一致）。接线守卫两腿过
   （|dcl_k0| ≤ 1e-8）。
2. **band (a) PASS 两腿**：变体系统在 amended 种子处的雅可比作用 FD
   最大相对误差 1.79e-7（coarse）/ 1.09e-8（medium），均 < 1e-5——
   替换行解析精确（9 个新测试 `tests/test_v5_te_outflow.py` 在单元
   层面同样锁死：默认 OFF 逐位、行结构、FD、J_e 清零、越 pattern
   守卫、平板 smoke、strip 配对）。
3. **band (b) honest FAIL（binding，"变差"档）**：变体系统在种子处
   残差 9.821（coarse）/ 4.846（medium）——替换行测得 TE 自然跳量
   （δ_TE − δ_up、H_TE − H_up 在松环收敛态非零，O(1) 量级）。伪时间
   （诊断 Q7 协议，同一协议下 flag-OFF 系统 11–21 迭代即达地板）仅降
   ~2–5× 即**全步拒 stalled**（cfl 压到 1e-3 地板：coarse it=11 起冻结
   于 4.703；medium 磨 75 迭代至 0.986 后冻结）。m2 = 1.752e-2 =
   **5554× 地板**（coarse）/ 4.487e-1 = **245998× floor_ref**（medium，
   binding）——比预注册"≥0.9× 不动"档更差，属"变差"档。
4. **次要读数与 guards**：GV5.1b tight polish（次要读数）coarse 终态
   F_max 7.32e-5（committed GV5.1b 末值 3.07e-6 的 23.8×）；medium
   **发散**（F_max 3.98，floor_stop 未救回）——polish 路径同样不破
   地板。guards：平板 H 带 flag-ON 过（层流 [2.606, 2.687] ∈
   [2.55, 2.75]，湍流 [1.509, 1.872] ∈ [1.2, 2.0]）；松环 smoke
   flag-ON **coarse 红**（cl_rel 2.62% > 2.5% A4 带且撞 10-outer 帽未
   收敛）/ **medium 边际过**（cl_rel 2.49%，3 outer 收敛）——TE 处理
   会把松环不动点 cl 推动 ~2.5–2.6%，物理上可察觉，是 flag 保持默认
   OFF 的又一理由。
5. **物理读法**：m2 峰值**不在 TE**——anatomy（诊断同款列命名
   F_delta/F_A/F_B/F_Psi/F_Ct1/F_Ct2）显示最大原系统残差在
   **x_c ≈ 0.027 的 LE 吸力区**（coarse 节点 53 F_B −1.75e-2；medium
   节点 310 F_Psi −0.449、F_B −0.271，两侧对称），TE 行本身安静。
   硬外推约束（δ_TE = δ_up、H_TE = H_up）与耦合动量/闭包系统在条带
   分辨率下**不相容**——变体系统在 amended 种子附近没有解，伪时间
   找不到下降路径，破坏迁移到最刚的耦合块（LE）。破地板不能靠在
   TE 带局部改写边界行绕开 (B,δ) 公式层刚度。
6. **程序状态**：盆地搜寻（GV5.1b/1c/1d）+ 首个公式层处理（GV5.5
   V1）均未动地板分毫。flag `te_extrapolate` 保持默认 OFF（legacy
   逐位一致；tight 舰队 + 全套件 flag-OFF 绿）。**升级阶梯保持登记
   未开**：upwind 边界通量 (a)-变体、闭包正则化 (b)——开启与否 =
   用户裁决。V5 保持 OPEN；下一步 = GV5.2（RAE2822 跨声速 VII）/
   GV5.3（M6 Cp）/ GV5.4（cost RECORDED）（用户 2026-07-24 定序）。
   V4 重开触发保持挂起。

## §18 GV5.2 —— RAE2822 跨声速 VII 对 committed 实验：band (b) FAIL + 松环配方极限解剖（2026-07-24→25，band (a)/(c)/(d) RECORDED）

用户定序项（GV5.1d → GV5.5 → GV5.2–5.4）。预注册
`bench/studies/v5_2_rae2822/PRE_REGISTRATION.md` 先于首行代码提交，
addenda #1–#3 先于各自（重）执行提交。协议 = 松环 VII（GV3.1 配方逐字
ω = 1.0、≤ 10 outer、tol_ds = 1e-3）+ GV3.2 跨声速 Newton driver 协议
（NEWTON_ARGS 导入，不重造）；几何 = Cook/AGARD-AR-138 Table 6.1
坐标（下表面列 positive-DOWN，签名锁定：12.1 % @ 37.9 %c 厚度、
1.3 % @ 75.7 %c 弯度、0 轮廓自交）；新网格族
`cases/meshes/rae2822_2.5d/`（M0 式嵌入尾迹准 2D 配方镜像 NACA 族：
coarse 5560 节点 / 16236 四面体，medium 20790 / 61494；stats + 分层
PNG 入账）；两点 P1（M 0.725 / α 2.55）/ P2（M 0.73 / α 3.19），
Re 6.5e6，两侧强制转捩 x_tr/c 0.03；coarse 记录、medium binding。
全程 8 线程临时约束（壁时标记不可比）。

1. **band (a) TE 楔角预检 RECORDED（无 fallback）**：mesh-crease 楔角
   （A4 法，**未切**网格上测——切后 TE 边消失）9.46° coarse /
   9.92° medium vs 坐标拟合 12.91°（x ≥ 0.95 割线和；差值 = 反弯度
   在拟合窗内的曲率）。≈6° quadratic 恢复守卫两级过
   （`quadratic_available=True`）→ 预注册 linear+smoothed fallback
   未触发。执行前修了两个反弯度特有缺陷（addendum #1）：`cut_wake`
   Kutta 探针加 TE 楔 bisector-normal fallback（RAE2822 后段下表面
   在弦线上方——两翼瓣邻居都在 TE 节点 +y 侧，旧探针失效；fallback
   只在旧代码抛错处触发，旧网格位等价；回归
   `test_kutta_probes_cambered_te`）；runner Cp 分侧换 outward-normal
   惯例（D11 `wall_outward_normals`——centroid-y 分侧把 coarse 上
   9 个 aft 下表面三角错标为 upper，会污染 band (b) 激波窗与分侧
   RMS）。
2. **band (b) FAIL（medium binding）**：x_shock = 窗内
   （x/c ∈ [0.2, 0.9]）压缩分支 max dCp/dx，读松环终态上壁面 Cp；
   接受带 = 实验 bracket ± 0.03 c（G4.1 无黏带）：P1 [0.495, 0.580]、
   P2 [0.520, 0.605]。实测 coarse P1 **0.6122**、coarse P2 **0.6520**、
   medium P1 **0.6288**（k = 10 顶格终态上读）全部带外；medium P2
   不可读（§6 配方极限条款）。**每个算出的激波都在实验带下游
   0.06–0.10 c**（P1 偏 0.06–0.08，P2 偏 0.05–0.10），随 Mach/α
   恶化；coarse→medium 在 P1 反而再往后移（0.6122 → 0.6288）且 LE
   吸力峰加深——**不是粗网格伪影**。方向（无黏族激波在黏性实验
   下游）符合 full-potential 求解在位移厚度反馈拉不回激波时的
   预期；幅度 = 本 gate 的诚实阴性发现。
3. **(c) Cp RMS RECORDED**：rms_upper/rms_lower = coarse P1
   0.185/0.146、coarse P2 0.176/0.118、medium P1 0.265/0.129——
   比 A4 输入带（medium ~2.5 % peak-rel u_e）高一个量级，由 (i) 激波
   位移（RMS 穿过错位跳变积分）与 (ii) 过深 LE 吸力峰（medium P1
   上表面 x/c ≈ 0.009 处 diff −0.27）主导；下表面一致性普遍显著
   好于上表面。
4. **(d) 收敛与 guards RECORDED**：松环不动点在跨声速点**不收缩**
   ——顶格腿 ds_change_rel 振荡（medium P1：0.06 → 0.49 → 0.11 →
   0.32）而非衰减，mdot_max 增长（0.013 → 0.13），ds_neg_floored
   爬到 k = 10 的 36 节点，IBL 每个 outer 都撞 100 迭代帽。
   medium P2 outright 发散：k = 4 mdot_max = 1.59（GV3.3
   transpiration 失控类；loud-fail 守卫触发；§6 条款读 RECORDED）。
   coarse P2 中途 IBL `MatrixRankWarning`（J 精确奇异）继续跑，
   k = 10 IBL 残差 NaN（记录）。M_peak：coarse P1 1.271 @ 0.155；
   coarse P2 **1.365**、medium P1 **1.306** 均超 1.3 包线 →
   outside-envelope RECORDED（非 FAIL）。**FP 救援链
   （addenda #2/#3：严格 1e-10 → 库内 Mach 延拓 m_start = 0.70 →
   诚实守卫 stall-accept）每腿承力**：每级 2 次延拓冷启动；
   stall-accept 占比 2/10（P1 coarse）→ 10/22（P1 medium）次 FP
   调用，全部卡在 |R| ≲ 1e-9 平台（Kutta 约束已收敛、接受时零
   limiter/floor 活动，accept_reason 入 run.log）；严格 1e-10
   始终是每次调用的首选。
5. **物理读法**：失败是**配方极限**，按预注册条款记录，不是代码
   崩溃：4 腿中 3 腿 ≤ 10 outer 不收敛（一顶格、一顶格 + 奇异 IBL、
   一失控），唯一收敛腿（coarse P1：7 outer / 150 s / cl 0.9036）
   激波仍偏后 0.06 c。miss 方向 + 环不收缩 ⇒ 松环更新通道的位移
   厚度反馈在 M ≥ 0.725 太弱/太慢——与 IBL 地板发现
   （GV5.1b/5.5）一致；逼出救援链的激波单元 Newton 平台同时标出
   内层求解在 M ≥ 0.725（2.5-D）的鲁棒边界。
6. **程序状态**：本 gate 未开任何跟进项（登记未开 = 用户裁决）；
   数据主张下一次跨声速 VII 读数走 **tight/augmented 耦合**，不再
   调松环。V5 保持 OPEN；下一步 = GV5.3（M6 Cp）/ GV5.4（cost
   RECORDED）（用户 2026-07-24 定序）。全套件基线见三面 ledger 行
   （随 baseline commit 填入）。V4 重开触发保持挂起。

## §19 GV5.3 —— M6 机翼方向+量级检查对 committed Cp：band (b) honest FAIL + band (a) input-limited（2026-07-25，0P/1F/17R）

用户定序项（GV5.1d → GV5.5 → GV5.2–5.4）。预注册
`bench/studies/v5_3_m6_cp/PRE_REGISTRATION.md` 先于首行代码提交，
addendum #1 先于重执行提交。问题：在 committed 实验条件（TEST 2308
逐字 M 0.8395 / α 3.06，Re_MAC 11.72e6，两侧强制转捩 x_tr/c 0.05）下，
松环 VII（GV3.1 配方逐字 + GV5.0 翼 case 逐字：翼尖带 z > 0.95·b_semi
钉扎 + ṁ 遮蔽）是否把 (a) CL 从无黏基线向下推过 A4 地板、(b) 壁面 Cp
向 committed 7 站实验靠近？FP driver = P14 跨声速配方逐字（M0.70 probe
种子 → NEWTON_M6_RECIPE ramp 从 `tests/test_p8_newton.py` 导入，pressure
Kutta，n_picard_seed=0；暖 outer 解；FP 侧**无 tip_taper** 使 k=0 解可
直接对 committed P14 锚点）；coarse 记录、medium binding。全程 8 线程
临时约束（壁时标记不可比）。

1. **接线守卫（全过）+ 首次执行叙事**：W1 = k=0 无黏 cl 对 committed
   P14 锚点（coarse 0.2628/0.2688、medium 0.2776/0.2823，1% 容差吸收
   ΔM = 0.0005 标签差）两级过（k=0 cl_KJ 0.2685 / 0.2819）；W2 = 实验
   侧向映射（每区 max-Cp 点在 x/c < 0.05 = LE 驻点；k=0 pooled RMS 选定
   映射 0.95/0.63 << 翻转映射 2.65/2.68）两级过；W3（截面提取稀疏守卫）
   未触发（7 站两级各侧 37–58 点）。**首次执行 medium k=0 未过 W1**
   （cl 0.226）：driver 冷启动带出了 GV5.0 bridge M0.5 的短路（probe
   种子不收敛就提前返回）——medium @8 线程 M0.70 种子停在 |R|~1e-6，
   Mach ramp 根本没跑，k=0 成了半收敛 M0.70 态。实测诊断（同工作副本同
   线程数）：**同一个失败种子**跑 P14 ramp 逐级收敛（0.70: 3 Newton →
   6.4e-6 … 0.8395: 12 → 7.6e-15）落在锚定分支（0.27726/0.28188，与
   本地 P14 缓存种子解 9 位一致）——**不是 8 线程分支散布**。addendum
   #1 删短路（ramp 无条件跑；执行机制，不动 band/协议）。值得注意的是
   循环本身曾从带毒种子 9 outer 收敛（k=1 暖解跳回 cl 0.277）——循环
   动力学健康，问题只在基线种子。W1 按设计立功。
2. **band (a) RECORDED input-limited**（medium binding，cl_KJ 读，
   cl_p 一致性引用）：Δcl_KJ = 0.2819 → 0.2757 = **−2.20 %**（cl_p
   −2.40 %，两估计量一致无分歧注记）；coarse −1.03 %/−1.35 %。方向
   向下（物理上正确的黏性去弯度效应）但低于 A4 2.5 % 地板 → 预注册
   "较小移动" 条款：RECORDED 标 input-limited（非 PASS 也非 FAIL）。
   注意 medium 移动集中在最后两轮：k=1–8 钉在 Δcl ≈ −0.3 % 的紧极限
   环，**k=9–10 晚期分离斑事件**（ds_max 0.0021 → 0.0039，
   ds_change_rel 0.448，mdot 0.007 → 0.085 → 0.106）带来其余部分——
   移动在 10-outer 配方帽处仍在增长。
3. **band (b) honest FAIL**（medium binding）：5 个未遮蔽站仅根站
   η=0.20 改善（−0.0019），其余 +0.0001…+0.0034；pooled RMS
   **0.1288 → 0.1299 反升**（带要求 ≥4/5 且 pooled 下降）。coarse
   0/5、pooled +0.0024 同向。所有 |ΔRMS| < 0.05 按预注册 A4 Cp 尺度
   注释标 input-limited——**FAIL 是方向判**：黏性修正除根站边际外
   处处不把 Cp 推向实验，反向移动的幅度小。翼尖遮蔽站（η=0.96/0.99）
   按构造记录（Δ ≈ +0.001…+0.003，无大移异常）。叠加图显示原因：
   无黏基线本身就是无黏族失配（LE 吸力峰太浅：实验 ≈ −1.0…−1.2 vs
   计算 ≈ −0.8；激波偏后），该强度下的松环黏性修正几乎不扰动上翼面
   曲线（虚线贴在实线上），加的小扰动在 4/5 站落在测量值的反侧。
4. **(c) 收敛与 guards RECORDED**：两级均 ≤10 outer 不收敛（GV5.0
   （M0.5）/GV5.2（2.5-D M≥0.725）同款松环签名——coarse ds_change_rel
   振荡 0.09–0.67 伴 mdot 爬至 0.33；medium 紧极限环（ds_change_rel
   0.010–0.025 不达 tol_ds 1e-3）后接 k=9 有界分离事件）；IBL 每个
   outer 撞 100 迭代帽，残差地板 1.9e-6（medium）；横流
   max|B|/max|A| = 0.026/0.055（GV5.0 类量级，无 3-D 失控）。预注册
   FP 救援链两级承力：10/21（coarse）/ 10/22（medium）次 FP 调用在
   诚实守卫平台上 stall-accept（各 1 次延拓调用；严格解始终首选；
   逐次路径 + accept_reason 入 run.log）。壁时：coarse 1721 s +
   medium 12479 s @8 线程。
5. **物理读法**：band (b) FAIL 是 GV5.2 RAE2822 发现的 3-D 对应——
   跨声速条件下松环位移厚度反馈太弱，修不动无黏族失配（LE 吸力不足、
   激波偏后），其加的小扰动反而落在测量值反侧；band (a) 方向物理正确
   但幅度在输入地板之下，且随分离斑事件仍在增长。与 GV5.2 联合的
   结论：后续跨声速读数属于 **tight/augmented 耦合**，不再属于松环
   调参。本 gate 未开跟进项（用户裁决）。
6. **程序状态**：V5 保持 OPEN；下一步 = **GV5.4**（cost RECORDED，
   用户 2026-07-24 定序）。全套件基线见三面 ledger 行（随 baseline
   commit 填入；本 gate 无库/测试改动，计数沿用 GV5.2 的 642+25+2，
   壁时重测）。V4 重开触发保持挂起。

## §20 GV5.4 —— M6 medium 增广步成本实测：7.53× RECORDED + 块预条件 honest FAIL（2026-07-25，0P/1F/17R）

用户定序末项（GV5.1d → GV5.5 → GV5.2–5.4）。预注册
`bench/studies/v5_4_cost/PRE_REGISTRATION.md` 先于首行代码提交，
addendum #1–#4 各先于对应（重）执行提交。问题（登记原文）：M6 medium
上增广 Newton 步壁时能否读进 ≤ ~2× 无黏 Newton 步的参考带（块预条件
工作的前提下；数值照录不论落点）？附带回答登记的 "measure before
Schur" 问题：J_BL,BL 的直接消元到底多贵（决定 Schur 路线是否值得）。
系统 = W2 增广系统 124,216 DOFs（62,820 φ + 166 Γ + 61,230 BL），
种子链 = A1 conf_newton 逐字（addendum #3）。全程 8 线程临时约束
（壁时标记不可比）。

1. **实现（库改动 + 接线守卫）**：`pyfp3d/viscous/tight_driver.py` 的
   `scaled_damped_step` 增可选 `solve(A, b)` 回调、`newton_tight` 增
   `step_solve` 透传（两条路径；默认 None = splu 逐位一致——默认路径
   不触任何已提交数字）；+2 测试 `tests/test_v5_tight_scaled.py`（回调
   接线 + 默认逐位一致）。runner 把 splu 与两档块预条件 rung 作为
   `step_solve` 注入同一驱动，步级壁时逐步入 CSV。守卫：W1 = 种子态
   cl_p 对锚（addendum #4 重订 P14 探针 G8.2 锁 0.2646；#1 把容差
   只绑 medium——A1 锚是 medium 数，coarse 偏差与已提交压力网格效应
   −5.35 % 相符）；W2 = 系统规模/块界守卫；W3 = FD 点验（增广 Jacobian
   对 φ/Γ/BL 三块，中位 φ 8.7e-12 / Γ 7.3e-12 / BL 0）全过。W1 cl_p
   0.26429 vs 0.2646 = 0.116 %。
2. **band (a) RECORDED**：增广步 22.93 s / 同 session 无黏锚 3.05 s/步
   （末层 13 步均匀）= **7.53×**，高于 ≤ ~2× 参考带——登记明言
   "recorded either way"，故记 RECORDED 而非 FAIL；4/5 增广步含
   capped-GMRES 工作（rung 打帽后回退 splu 的步），VERDICT 如实注明。
   A1 committed @16t 4.35 s/step 仅作非约束对照（线程数不同不可比）。
   coarse 摸底（叙事不入账）：rung-2 working（RECORDED），ratio
   20.58×——coarse 无黏步 0.17 s 太小读不出数。
3. **band (b) honest FAIL（D5）——块预条件在 medium 不工作**：
   rung-1 块 Jacobi（AMG-φ + ILU-BL 块对角近似）**发散**（rel_res
   5.75e4；φ–BL 非对角耦合太强，块对角近似根本压不住；该步为
   least-bad 步）；rung-2 exact-BL Schur（Schur 算子上用 AMG-φ）
   仅 1/4 步收敛（2.66e-8 @277 迭代），其余三步停滞
   2.07e-7 / 6.13e-5 / 2.52e-6 vs rtol 1e-8 @300 帽——纯 AMG-φ
   循环表不出 Schur 修正 J_hB·J_BB⁻¹·J_Bh（φ 侧预条件看不到 BL
   消元对势流的反作用）。
4. **"measure before Schur" 答案**：splu(J_BL,BL) setup 仅 **1.8 s**
   ——BL 块直接消元便宜，**瓶颈在 Krylov 收敛而不在分解**（AMG-φ
   setup 0.2 s、ILU-BL 2.5 s）。即 Schur 路线的消元成本不是障碍，
   障碍是 Schur 算子的谱性质需要感知耦合结构的预条件。IBL 地板
   轨迹全程完整（merit 钉 9.35e-9，|F_BL| 2.888e-6，步 2–5
   λ ~ 1e-4 严格下降）——成本实验不扰动已入账的地板结论。
5. **执行叙事（四个 addendum，各先于重执行提交）**：#1 coarse W1
   触发（cl_p −4.97 % vs A1 锚）→ 容差只绑 medium；#2 medium M0.70
   probe 种子严格不收敛 → 种子链不收敛不 raise（照 A1 行为）；#3
   种子链定为 A1 conf_newton 逐字（消任何配方漂移）；#4 W1 锚从
   陈旧 A1 0.26918 重订 P14 探针 G8.2 锁 0.2646（A1 数是旧分支态，
   G8.2 探针锁才是现行锚）。
6. **判读与程序状态**：预注册要记的诚实负结果——翼尺度增广 Newton
   在块对角/纯 AMG-φ 类预条件下读不进 ≤ ~2× 带；要读进带需要更强的
   （Schur-aware、(A,Ψ) 结构化）约化空间预条件。跟进登记未开：
   Schur-aware 预条件阶梯 + EW-forcing 变体（均用户裁决）。**V5 五个
   gate 全部执行完 → CLOSED 2026-07-25**（close-out 待用户裁决）。
   全套件基线 642 → 644（+2 = 新 W4 回调测试）见三面 ledger 行。
   V4 重开触发保持挂起。

## §21 GV5.6 —— Schur-aware 约化空间预条件器：corrected-AMG 路线实测死刑（2026-07-25，0 PASS / 1 FAIL / 17 RECORDED）

GV5.4 登记的 follow-up，2026-07-25 用户裁决开立。预注册
`phases/p1/cases/analysis/v5_6_schur_prec/PRE_REGISTRATION.md` 先于首行代码提交
（`091f9fe`）；系统/种子/协议 = GV5.4 逐字（124,216-DOF W2 系统、A1
conf_newton 种子链、rowcol 均衡、mu ≡ 0、N = 5 实测步、同 W1/W2/W3
守卫 + D5 二元裁决）；**零库改动**（阶梯全部在 case runner 内，复用
已提交的 `step_solve` 注入点；W4 = pyfp3d/tests diff 为空，套件基线
不变 652+25+2）；runner 默认 16 线程（GV5.4 的 8t 数字仅作非约束
对照）。

1. **设计（rung 3 + escalation rung 4）**：rung 3 = GV5.4 rung-2 的
   exact-BL Schur 算子逐字，但约化空间预条件换成
   bdiag(AMG(Ŝ_φφ), M_Γ)，Ŝ_φφ = J_φφ − Ĉ_φφ 显式装入稀疏化 Schur
   修正 Ĉ = J_hB·D_BB⁻¹·J_Bh——D_BB = J_BL,BL 的每节点 6×6 块对角
   （(6,6)-BSR 视图提取；rcond 守卫 1e12 → 零安全对角 fallback；
   物理上 = quasi-simultaneous 局部 BL 响应，即被跳过的 V4 路线以
   预条件代数形式回归；(A,Ψ) 结构化指导经每节点整块兑现）。rung 4
   = 块上三角全系统预条件（y_B = lu.solve(r_B)；y_h = P_hh(r_h −
   J_hB·y_B)），两个耦合方向 + Schur-aware φ 循环都进预条件。
2. **medium（binding）band (b) honest FAIL**：预注册的可证伪假说
   以其朴素形式被**证伪**——rung 3 唯一一步 GMRES info=5 @242
   迭代、rel_res **0.664**（对照 GV5.4 rung-2 同种子：首步收敛
   2.66e-8 @277，其后停滞 2e-7..6e-5）——修正 AMG 矩阵在 medium
   灾难性劣于纯 AMG(J_φφ)；rung 4 四步全部打帽 0.68 → 1.06（单调
   恶化）⇒ NOT-WORKING。(a) RECORDED 23.88 s / 3.03 s = **7.87×**
   （≤ ~2× 带上方，照录；≈ GV5.4 的 7.53×；全部步含 capped-GMRES
   工作）。
3. **coarse/medium 分裂 = 解剖**：coarse 上 rung 3 **工作**（5/5 步
   info=0，120–209 迭代，rel_res ≤ 3.5e-9，无 escalation，
   RECORDED）——D_BB-局部修正在 coarse 够用；medium 上 Ĉ_φφ 把 φ
   块密度约翻倍（nnz(Ĉ) 597,154、nnz(Ŝ) 1,560,192），要么被丢弃的
   节点间 BL 耦合在 medium 变为本质，要么 Ŝ 失掉让纯 AMG(J_φφ) 近乎
   可用的代数性质。setup 开销无罪：t_corr 0.2 s（D_BB 求逆 + 三重
   积）、t_amg 0.2 s、t_lu 1.8 s（GV5.4 的 "measure before Schur"
   数字复现）、n_fallback = 0（10,205 个 6×6 块全部精确可逆）。
4. **守卫与轨迹**：W1/W2/W3 一次过（cl_p 0.26429 vs P14 探针锁
   0.2646 = 0.116 %；FD 中位 φ 9.6e-12 / Γ 7.3e-12 / BL 0）；IBL
   地板轨迹完整（merit 钉 9.385e-9，λ ~ 1e-4，步全接受）；无
   addendum（一次干净执行）。
5. **判读与程序状态**："AMG on a sparsified Schur" 方向在翼尺度
   实测死刑（两个块方向、D_BB-局部物理都测过）；剩余登记路线 =
   Ĉ 内换节点间 (A,Ψ)-structured M_BB（站位条块逆 / ILU 型近似逆）
   或放弃 corrected-AMG 形式，连同 EW-forcing 变体保持
   registered-not-opened（用户裁决）。**GV5.6 ✓ CLOSED
   2026-07-25**。
