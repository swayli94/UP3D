# Track V — 粘性/无粘交互（VII）：数值设计参考

> 状态：**设计参考文档**（2026-07-22，V1 开工时建立，随实现推进补充实测记录）。
> 阶段/gate/进度以 [roadmap/track_v.md](roadmap/track_v.md) 为准，冲突时 roadmap 胜出。
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
实现与参考码或诚实记录。详见 `cases/analysis/v1_ibl3_standalone/VERDICT.md`。

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
证据见 `cases/analysis/v3_loose_coupling/VERDICT.md`（GV3.1/3.2）与
`cases/analysis/v3_fuselage_smoke/VERDICT.md`（GV3.3）。

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
`cases/analysis/v5_m6_bridge/VERDICT.md`。实现决策记录：

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
(φ, Γ, U) Newton；门证据与诊断见 `cases/analysis/v5_tight_coupling/`
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

诊断研究 `cases/analysis/v5_ibl_floor/`（预注册 53bf904 先于首次执行；
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
