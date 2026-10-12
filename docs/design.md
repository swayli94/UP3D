# pyFP3D — 设计参考
## 三维非结构网格全势流跨声速外流求解器（Python + Numba）

本文是**现行设计参考**。求解器已实现，并经历了六个开发阶段（phase 1–6，
2026-07-06 → 2026-10-12；各阶段归档索引见 `phases/README.md`）。本文描述现行
实现的设计、方程与已实测的行为——不是开发计划，也不是变更日志。

- **现行能力声明（实测边界）**单独维护在
  `phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md`（2026-10-11）。
  本文给出的包线、极限、误差数字是设计依据；一个能力能不能信、信到哪，
  以能力边界文书为准。
- 活代码里约 40 处 docstring 按节号引用本文（"design.md §2"、"§6.3"、
  "§12 risk 3" 等）。**§1–§13 的节号骨架是契约，不许重排**；§12 由
  `tests/F/test_F10_risk_gate_index.py` 机械解析（风险条数、守门指针、
  「无门」登记数都有锁），格式要求见该节。
- 本文前身（英文为主、维护至 phase 6 收口的活文档）归档在
  `phases/p6/docs/design.md`；两版节号一致。

---

## 1. 定位与适用范围

**目标问题。** 机翼与翼身组合体的定常跨声速外流：来流 Mach 0.3–0.87，
附着流或弱激波流（局部法向激波 Mach ≲ 1.3）。等熵、无旋、无粘。

★ 2026-10-10 使用者裁决（H11 缩包线）：上界 0.87 是**声明**；α = 0 / 低升力的
**实测**可信上界为 **M 0.78–0.80**（D14：`tests/D/test_D14_alpha0_transonic_nonuniqueness.py`），
§2/§12 登记的 Steinhoff–Jameson 带（0.82–0.85，低升力）为非唯一区，
不计入能力包线。

**全势（FP）相对 Euler 买到什么。** 每节点一个标量未知量 φ，而不是五个守恒
变量；椭圆主导的算子适配强隐式求解器（Newton + AMG 预条件的 Krylov）；单个
解便宜 1–2 个量级。适合做优化/代理模型循环里的气动引擎，以及粘性–无粘耦合
的无粘内核（2D 里 TSFOIL 扮演的角色——本求解器是 pyTSFoil TSD→FP 升级路线
的三维非结构类比，且摆脱了小扰动假设与笛卡尔网格限制）。

**FP 做不到什么。** 激波后有旋流（激波被替换为满足质量守恒跳跃的等熵压缩，
不是 Rankine–Hugoniot）；激波位置误差随激波强度增长；已知 FP 解在
M∞ ≈ 0.82–0.85 低升力窄带内对常规翼型不唯一（Steinhoff–Jameson）。这些是
已接受的模型限制，登记在 §12，不是实现目标。

**网格类别。** 非结构四面体（以后可选混合 prism/tet），顶点中心
（vertex-centered）。尾迹必须是网格中的**贴体内表面**（conforming internal
surface，见 §4、§12——这是对网格生成器施加的唯一最强约束）。

---

## 2. 控制方程与归一化

速度势 Φ，**V** = ∇Φ。长度用参考弦长 c、速度用 U∞、密度用 ρ∞ 归一化。
记 φ 为无量纲势，q = |∇φ|：

**守恒形全势方程**

    ∇·(ρ ∇φ) = 0                                   (2.1)

**等熵密度律**

    ρ(q²) = [ 1 + (γ−1)/2 · M∞² (1 − q²) ]^{1/(γ−1)}          (2.2)

**当地 Mach 数**

    M²(q²) = q² M∞² / [ 1 + (γ−1)/2 · M∞² (1 − q²) ]           (2.3)

**声速（无量纲，平方）**

    a² = 1/M∞² + (γ−1)/2 · (1 − q²)                            (2.4)

常用导出量（实现为纯 numba 函数，对手算值做单元测试）：

- 临界速度 q\*²（M = 1 处）：q\*² = [2 + (γ−1)M∞²] / [(γ+1)M∞²]。
- dρ/d(q²) = −(M∞²/2) ρ^{2−γ} = −ρ/(2a²) · …；实际实现直接由 (2.2) 写
  ρ′(q²) = −(M∞²/2) · ρ^{(2−γ)}。
- 沿流线的质量通量 ρq 恰在 M = 1 处取最大值。这一非单调性就是中心格式在
  超声速区失败、必须对 ρ 做迎风的数学原因（§3）。

**压强系数（精确等熵式，非线性化）：**

    Cp = 2/(γ M∞²) · [ ρ^γ − 1 ]                               (2.5)

其中 ρ 来自 (2.2)。永远不要用线性化的 Cp = −2u——本求解器的全部意义就在于
在线性化失效的地方仍然成立。

**(2.1) 的型性。** 拟线性二阶；非守恒形 (a²δij − uiuj) ∂²φ/∂xi∂xj = 0——
M < 1 处椭圆，M > 1 处双曲。必须做随型变化的离散；我们通过密度迎风来做，
而不是算子切换（非结构网格上不用 Murman–Cole 式开关逻辑）。

---

## 3. 跨声速处理：人工密度（通量偏置）

沿 Hafez–South–Murman / Holst 人工压缩性一族：超声区把 (2.1) 中的 ρ 换成
迎风密度

    ρ̃ = ρ − ν · Δℓ · (∂ρ/∂ℓ)_upwind                            (3.1)

ℓ 为当地流线方向，Δℓ 为当地网格长度，开关函数

    ν = C · max( 0 , 1 − M_c² / M² )                            (3.2)

截止 Mach M_c ≈ 0.95–1.0，常数 C ∈ [1, 2]（C、M_c 是求解器参数——与
pyTSFoil 参数集同一套标定哲学；从第一天起就暴露在配置里）。

**非结构实现（逐单元，迎风项一阶）：**

对每个单元 e，速度 **V**_e = ∇φ|_e（线性 tet 上为常数）：

1. 找**上游邻居**单元 u(e)：在 e 的面邻居中，面外法向与 **V**_e 最反平行
   的那个。面邻接一次性预计算；u(e) 每个非线性迭代重算（便宜：每 tet 4 个
   点积）。
2. 取 ρ̃_e = ρ_e − ν_e (ρ_e − ρ_{u(e)})。这就是 (3.1)，其中
   Δℓ·∂ρ/∂ℓ ≈ ρ_e − ρ_upstream，形式鲁棒且与网格无关。
3. 单元通量/残差（§6）中用 ρ̃_e。亚声速单元（ν = 0）恒等退化为中心
   Galerkin 格式——二阶精度。

必须保持的性质（为它们写断言/测试）：

- ν = 0 ⇒ 格式对称，Picard 矩阵 SPD（§8）。
- 激波在 2–3 个单元内单调捕捉；无膨胀激波（ρ 的迎风只在穿过声速的压缩中
  提供所需耗散）。
- 均匀流保持：φ = x 的均匀场在任意网格上——包括跨越尾迹切口——必须给出
  机器零残差。这是每次内核改动后的第一道回归测试。

### 3.1 现行 P4 算子及其两个不同的缺陷（P6 设计复查）

P4 实现（`kernels/upwind.py`）是 López 有限元人工压缩性的*逐单元迎风*实现
（学位论文 Eq. 3.19–3.27，§13）：

    ρ̃_e = ρ_e − μ_s (ρ_e − ρ_{u(e)}),   μ_s = C · max(0, μ_e, μ_{u(e)})     (3.3)
    μ = 1 − M_c²/M²          （López 3.21；UP3D 把分母钳为
                                max(M², M_c²)，亚临界时 μ ≡ 0，逐位精确）

其中 `C = upwind_c`（López 的 μ_c），`M_c = m_crit`。上游单元 u(e) 由
**离散多跳定向行走**选出（`upstream_elements`：穿过面邻居，最多 4 跳，每跳
取最逆流向的邻居，目标行程 0.8·单元流向尺度）。ρ̃_e 是**冻结的逐单元标量**，
同时加权残差 (6.1) 与 Picard 矩阵 (6.2)。两个硬钳制守正值性：M_cap = 3 的
q² 限速（`limit_q2_field`，López 的 M²_max ~ O(3)，§3.4）与 0.05 的 ρ̃ 下限
（`rho_tilde_sweep`）。

这个算子位于 roadmap 原先混为一谈的**两个独立问题**的中心——N1 调查
（2026-07-08）用测量证明它们**根因不同、修法不同**：

- **表面 Cp 锯齿（精度——门 G6.1）是壁面 Cp 梯度*恢复*的伪迹，不是通量的
  问题。** 超声速算例壁面 Cp 的 ≈2 单元锯齿由*逐三角形常值的 P1 梯度恢复*
  （`triangle_tangential_gradients`，无平均）在棱柱劈分细长三角形的壁面剖分
  上产生：相邻细长三角形即使对*光滑* φ 也给出振荡梯度。决定性证据（均在
  NACA0012 M0.80 粗网格解上）：(i) 对**同一**行走解换用法向门控的边邻居
  平滑梯度重算壁面 Cp，G6.1 锯齿指标降 **~330×**（0.0758 → 0.00023，斜率
  反转 39 → 1），激波位置 0.604→0.607，cl_p 变化 0.3% 以内——φ 场本来就
  光滑；(ii) **更光滑的人工密度通量**（§3.2 的核）对该指标**完全没有**
  改善（每种标定下持平或更差）。所以 G6.1 在后处理修（§9 `smooth_passes`），
  不动通量。*这推翻了先前（demo_report §P4）归因于 ρ̃ 选择翻转的说法——
  该归因已被否证。*
- **通量对 φ 的不可微性（挡住精确 Newton Jacobian）是真实且独立的缺陷——
  但它属于 P7，不属于 G6.1。** 整数 u(e)、`max(0, μ_e, μ_{u(e)})` 开关和
  两个硬钳制都对 φ 不光滑；精确 Newton Jacobian（§6.3）需要的是**冻结选择
  下** ∂ρ̃/∂φ 有定义——现行**行走**已经满足（见下段）。推导并 FD 验证该
  ∂ρ̃/∂φ 才是 **P7** 交付物与 **P8 Newton 前置**，而不是修锯齿；§3.2 的核
  只是一个*可选*的更光滑替代，不是必需品。

PDF 核实的洞察（López Appendix B，§6.3）进一步明确了 P7 通量要提供什么：
**它不要求*选择*可微。** López 在每个 Newton 步内冻结 u(e)，只对 ρ_{u(e)}、
μ、ρ *本身*求导（Eq. B.3–B.8）；他自己的 `μ_c·max(0,μ,μ_up)` 只是 C⁰ 也
照样二次收敛（不光滑轨迹是零测集）。两个钳制在收敛物理解处不激活
（P5 medium：**0 触底 / 0 限幅**），所以它们的不光滑性从不进入收敛
Jacobian——光滑化它们（§3.2）只是瞬态鲁棒性选项，不是正确性要求。

### 3.2 可微迎风密度（P7——P8 Newton 前置）

> 范围说明（2026-07-08 重圈定）：这个核是**可选**的平滑器 / 快 ~10× 的
> Picard 通量，**不是** P7 Newton 前置的强制项。P7 前置是现行**行走**的
> *冻结选择*可微性（§3.1、§6.3）——固定 u(e) 时其 ∂ρ̃/∂φ 本来就有定义且
> 稀疏（~+1 单元/行）。这个核只是 P8 Newton 通量的*候选*，前提是其更密的
> Jacobian（§6.3）实测净划算。它不是 G6.1 锯齿的修法（§9），保持 opt-in
> （`UpwindOperator(weighted=True, mode="kernel")`）；现行默认 = P4 行走。

朴素形式——把单个整数邻居换成**光滑的面邻居加权**

    ρ_up(e) = Σ_f w_f(e) · ρ_{nb(f)} / Σ_f w_f(e),
    w_f(e) = max(0, −V_e · n̂_f)^p         (n̂_f = 单位面外法向)     (3.4)

对 V_e = ∇φ 是 C¹，但在薄片 tet 上**瞬态不稳定**：N1 探针（从行走收敛粗场
出发的冻结-Γ 快速 Picard）对所有 p 与行程增益都在 ~60 次迭代内把它吹爆
（M_max 20–40），起点还是*好*场。原因**不是**耗散量级（收敛场上它给出
行走的 1.3×）——*近*环只达到 ≈0.3·extent，对它平均采样的是过局部、过
中心的区域，破坏了稳定双曲区所需的迎风特性。（途中排除两条死路：面法向
入流加权 `max(0,−V·n̂_f)` 在薄片上产生*负行程*的反耗散单元——要用形心
位移对齐 `max(0,−V̂·d̂_{e,nb})`；行程倍增 `G·(ρ_e−ρ_up)` 也救不了，近环
够不到上游时 `G·0 = 0`。）

**现行 opt-in 算子**因此是有真实行程的**多环流线-高斯核**：在预计算的
深度 3 BFS 邻域上，把每个单元 c 用中心位于上游一个流向尺度处的高斯加权：

    up_c = −V̂·(c_c − c_e),  perp²_c = |c_c − c_e|² − up_c²
    w_c  = exp(−½[(up_c − reach·ext)²/σ_s² + perp²_c/σ_p²]),   up_c > 0
    ρ_up = Σ_c w_c ρ[c] / Σ_c w_c                                       (3.4′)

它保住了行走的稳定**行程**（采样的是真上游单元，不是近环），在固定邻域下
对 V̂ = ∇φ/|∇φ| 是 C^∞（Newton 就绪），且——没有每迭代要重算的行走——
每次密度迭代比行走快 ~10×。实测稳定：reach ≈ 0.8–1.2·extent、
σ_s,σ_p ≈ 0.25–0.4·extent 时收敛到行走*自己的*物理状态（M_max 1.37，
0 触底/限幅），而近环混合发散；reach ≈ 0.6 太短。保留：ν ≡ 0 亚临界精确
无操作（G4.2）、正值性钳制、(M²−1)/M² 耗散下限。激波点算子保持对光滑
混合上游 Mach 的**硬** max(ν_e, ν_up)——亚临界下恰为 0，所以不用光滑
`max_ε`（max_ε(0,0)=ε≠0 会破坏无操作性质）。

**它并不降低表面 Cp 锯齿**（那是恢复伪迹，§3.1/§9）：按复现行走的粗网格
激波/cl 标定（reach 1.0、C≈2.5：激波 0.619 vs 0.604，cl_KJ 0.356 vs
0.364）后，G6.1 指标*持平或更差*。它的价值在于**可选**的 Picard 加速 /
更光滑通量；P7/P8 Newton 默认通量是稀疏的冻结**行走**（§6.3）。若哪天
升格为 Newton 通量，需要重标定 `upwind_c`（≈2.0–2.5，行走是 1.5），因为
核的有效耗散与行走不同。现行默认 = P4 行走。

**可选——光滑密度钳制（仅当 Newton 在触碰钳制的瞬态上卡死）。** 把硬
M_cap 钳制换成 σ = ½(1 + tanh((M² − M²_max)/ε))，
ρ_clamped = (1−σ)ρ(M²) + σ ρ(M²_max)，它跨 cap 是 C¹。López 用的是**硬**
钳制也照样二次收敛（钳制只在孤立的翼尖-TE 几何奇点处激活——学位论文
Fig. 3.5、§3.4 "only local parts of the boundary are affected"），且 UP3D
的收敛态不触碰它，所以列为兜底，不是前置。

---

## 4. 环量、尾迹模型、三维 Kutta 条件

有升力 ⇒ φ 多值；计算域必须用一张从后缘发出、延伸到远场的**尾迹面** W
切开。

**网格要求。** W 是贴体内表面：W 上的节点复制成两份（上侧拷贝 i⁺、下侧
拷贝 i⁻），上方单元挂 i⁺，下方挂 i⁻。Gmsh 里把尾迹建成嵌入面再用
plugin/crack 步骤复制，或者在求解器网格预处理器里从带标签的内表面生成
复制（推荐——网格生成保持 vanilla，节点劈开在 `mesh/wake_cut.py` 里做）。

**W 上的跳跃条件：**

    [φ] ≡ φ⁺ − φ⁻ = Γ(s)      （沿流向尾迹线常值）   (4.1)
    [ρ ∂φ/∂n] = 0             （质量通量连续）        (4.2)

s 参数化展向。(4.2) 自然满足——把尾迹单元当内部单元组装（两侧通量都加进
约束系统）；(4.1) 用主–从线性约束：

    φ_{i⁺} = φ_{i⁻} + Γ(s_i)                                            (4.3)

实现：消去 i⁺ 自由度（把它们的行/列折进 i⁻，外加一个正比于 Γ 的右端项）。
消去映射一次性预计算；Γ 只进右端项，所以 Picard 迭代里更新 Γ 不需要重新
组装矩阵。

**复制必须包含后缘节点。** 连续意义下 [φ](TE) = Γ——这本身就是 Kutta
条件——所以切口（和节点加倍）要延伸到壁面上。若 TE 节点保持单值，跳跃会
在第一个尾迹单元上被压成零，等价于在 TE 停了一个强度 Γ 的点涡：壁面吸力
~ (Γ/2πr)² 积出 ~ Γ²/h 的虚假力，且随加密*发散*（P2 实测：粗 NACA0012
网格上六个 TE 相邻三角形贡献 −0.27，而期望 cl ≈ 0.6；同一算例 TE 加倍后
cl = 0.6012，Kutta–Joukowski 0.6）。

**Kutta 条件。** TE 上/下压强相等 ⇒ 势形式下每个展向站位 j：

    Γ_j^{new} = φ_{TE,j}^{upper} − φ_{TE,j}^{lower}                     (4.4)

在 TE 偏一个节点处取值（或外推到 TE），在非线性外循环中更新，欠松弛
ω_Γ ≈ 0.7–1.0。线性（Laplace）驱动器上 Γ → Γ_target 映射是仿射的，实测
斜率 b ≈ 0.93（TE 偏一节点的光滑流跳跃几乎就是 Γ 本身），所以朴素松弛
收敛很慢；P2 驱动器从第二次更新起用逐站位割线（Aitken）步，2–3 次更新
收敛。欠松弛形式保留为兜底，也是非线性（P3+）外循环的基线。站位之间的
Γ(s)：展向参数上分段线性；尾迹线继承其发出站位的 Γ——尾迹网格预处理器
必须建这个站位→尾迹线映射。收敛监控：‖ΔΓ‖∞ 与残差范数并列。

> **实现注记（A3 2026-07-18，kimi 代码审查）。** 现行
> `mesh/wake_cut.py::node_station` 给每个尾迹节点赋**最近 z 站位**的
> Γ——分段**常值**，不是上文指定的分段线性。两者在分层/挤出片网格上完全
> 一致（每个尾迹节点都落在某站位的 z 上，所有已提交的 quasi-2D 族都属此
> 类），在完全非结构的 M1/M2 片上差 O(Δz)——与该片本来就有的 P1 插值
> 误差同阶，不值得盲改。记录在此，免得规格与代码悄悄不一致。

**(4.4) 自 P14（2026-07-17）起是 LEGACY（探针）估计器。** A2 实测了它的
两个缺陷：单节点探针采样*制造*展向 Γ(z) 抖动（固定-Γ 判别量
D = 7.33/25.70），且等*势*跳跃只是等*压强*的 O(h) 有偏代理（TE Cp 间隙
比 level-set 路径差 34×/133×）。P14 替代品（`constraints/te_pressure.py`，
在 `solve_laplace_lifting` 与耦合 Newton 驱动器上以
`kutta_estimator="pressure"` opt-in；默认 "probe" 逐位不变）是**贴壁控制体
压强相等估计器**——B4 的对象移植到 conforming 切口网格上：

    F_j = 站位 TE 节点上 (|q_u|² − |q_l|²) 的均值 = 0   (4.4′)

q_u/q_l 是 TE 节点上侧（从拷贝）/下侧（主拷贝）贴壁单元扇上体积加权的
P1 速度（精确壁面归属；整扇取会错 +11–15%，B4）。精确分解
(q_u+q_l)·(q_u−q_l) 冻结均值给出割线驱动器消费的仿射隐含目标；Newton
驱动器用精确二次残差，每步重建行 (K_p, D = ∂F/∂Γ) 与广义消去
K̃ = −D⁻¹K_p（probe：D = −I 精确退化为历史代码）。两条实测性质必须尊重：
(i) 两种闭包的逐点一致只到探针自己的 O(h) 读数偏差（交叉读数 M6 M0.5 上
粗 3.7% → 中 1.05%），接近 1 的映射斜率 b ≈ 0.93 经 1/(1−b) ≈ 14× 放大进
收敛 Γ——收敛升力合法地移动（M0.5 粗/中 +2.1%/+4.5%）；(ii) 二次行的
Newton 吸引盆比仿射探针行小——用探针解做压力解的种子（P14 配方；M6
medium 上冷 Picard-5 种子会游走并快速失败）。`solve_subsonic_lifting` 的
内层割线与 `continuation.py` 的 Γ-割线/抛光保持探针形式（P14 接线范围，
使用者裁决）。

尾迹几何：平面，与来流（或弦平面）对齐——标准 FP 做法；无力尾迹松弛超出
范围（附着流下误差是载荷的二阶量）。

**已设计并已删除的替代——Track B level-set 尾迹**（设计 2026-07-07；
phase 1 大体交付；phase 3 裁决 D5 整线删除，库文件 9 个 / −4624 行；归档
原件 `phases/p1/docs/design_track_b.md`、台账
`phases/p1/docs/roadmap/track_b.md`）。思路是 level-set + 多值单元
（CutFEM 式）尾迹表示，Kutta 隐式施加（Γ 从解中涌现——无 Γ 自由度、无
割线），替代上文的 conforming 片 + 主–从消去，使网格无需预嵌尾迹面
（α 扫描不重划网格；翼身多尾迹）。它留下两条仍然成立的设计结论：原始
罚 Kutta 降级为可选诊断量，真正的闭包是非线性 TE 压强相等 Kutta（B4——
这条结论以 (4.4′) 的形式活在 conforming 路径里）；弯曲/自由尾迹**搁置**
（2026-07-10：直尾迹载荷误差 O(θ²) ≈ 0.1%，几何更新迫使每步重建切口集/
自由度，与 Newton 冲突）。本节的 conforming 路径是唯一的现行路径。

### 4.1 翼尖/尾迹边缘奇点（P13）

上一节的刚性**平面**尾迹终止于一条自由边——三维里后掠翼尖的 TE 角点，
尾迹片从那里向下游延伸，但*越过*翼尖没有延续（M1：翼尖节点保持单值，
Γ(tip) = 0 离散成立）。这条自由边是**奇异的**，且奇点是尾迹**模型**的
性质，与具体离散或壁面几何无关。

**绑定 Γ vs 尾涡。** 翼尖的正则条件是 Γ(tip) = 0——上/下速度相等，那里
没有*绑定*涡。这是**必要不充分**。物理上卸进尾迹的对象是强度
γ(s) = −dΓ/ds 的**尾涡片**（s 为展向参数）：机翼朝翼尖卸载，在翼尖达到
峰值的是*卸载率* |dΓ/ds|——不是 Γ 本身（有限翼近椭圆载荷下连续意义
dΓ/ds → −∞；B7 level-set 解的光滑 Γ(z) 上实测翼尖处比半展向大 ~10×，
而 Γ 本身从 0.10 降到 0）。一张**平**的、到此为止的尾涡片诱导出的横流会
翻转自由边——与有迎角平板前缘同一机制——所以边缘速度行为是 **1/√r**
（r 为到边的距离），即**平板边缘（平方根）奇点**。它*不是*集中线涡的
1/r（那需要把片的涡量收到一条线上——正是刚性平面模型省略的卷绕）。加密
研究用峰值速度 ~ h^(−p) 的指数 p 区分两者：1/√r ⇒ p ≈ 0.5，1/r ⇒ p = 1。

**证据（P13/G13.1，2026-07-13）。** 亚声速探测（M∞ 0.5，无激波/人工
密度/限幅），ONERA M6 翼尖盒上，conforming 三点（粗/中/细）峰值当地
Mach 对 1/h 的 log-log 指数为 **p ≈ 0.59**——与 1/√r 一致，否证了早先
"1/r 型"的措辞。（**用这幅连续图像推理修法之前先读 §4.2。** 指数是对的，
但*求解器实际看到的机制*是离散的——最外 TE 站位把保留的 Γ 在一个单元上
卸掉——能预言哪种渐变有效的是那幅离散图像，不是这幅连续图像。）峰值在
conforming 和 Track-B level-set 两条路上都发散（换表示磨不掉它），而周围
流场保持网格收敛（*局部*边缘奇点）。加密足够时边缘过速即使在亚声速也
触发限速器/密度下限，于是细网格不再是离散解（P9/G9.1 跨声速发现的
M∞ 0.5 类比）。证据：`phases/p1/cases/demo/p13_tip_edge_singularity/`。

### 4.2 修法：展向载荷渐变（P13/G13.2）

曲面*壁面*单元（P11）去不掉尾迹片边缘，Track B 的 level-set 改的是尾迹
*表示*而不是模型。但解药最终**既不是卷绕也不是涡核**——本节早期草稿提过
这两条。实测机制是**离散的**：

**求解器永远看不到连续边缘。** 它看到的是**最外 TE 站位**：它保有环量
Γ_last，并把它作为*集中涡在最后一个单元上*卸掉——尾迹片自由边节点是
单值的，跳跃从 Γ_last 跨一个单元掉到 0。这个卸下的涡诱导速度
~ Γ_last / h。于是取 Γ_last ~ h^q，

    边缘峰值 ~ h^(q−1)        ⇒        p ≈ 1 − q                       (4.5)

正则化判据是 **q ≥ 1**。它精确预言了无渐变基线：Γ ~ √u 且 u_last ~ h 给出
Γ_last ~ √h，即 q = 0.44（实测），故 p ≈ 0.56——对实测 p = 0.52。要紧的
*不是*连续边缘指数，也*不是*粘性核，而是最后一个站位的环量随加密消失得
多快。

**修法随之直接。** 把接受的环量朝翼尖渐变：

    Γ_eff(z) = F(z) · Γ_Kutta(z),   F = smoothstep((z_tip − z)/r_c)      (4.6)

施加在逐站位 Kutta 目标上（`constraints/wake.py::tip_taper_factors`；
`solve_newton_lifting(tip_taper=…)`，默认 `None` = 逐位不变）。F 只依赖
几何（与 φ 无关），Newton Jacobian 结构不变——Kutta 行只是被同一个 F_j
缩放。现行模型：`vanish_smooth`，**r_c = 0.05 · b_semi**。

**F 必须紧支。** tanh 渐变（F(tip)=½）也能到 q≈1 也确实正则化，但 tanh
永远到不了 1：它的尾巴把机翼卸载到 η≈0.77，代价 −7.4% 升力、−4.9% 内段
环量，并在没有奇点要修的中–外展向破坏 TE 压强闭包。翼尖模型必须*局部*。

**渐变是被放大的，不是直接施加的。** Γ 是 Γ = F·Γ_Kutta(Γ) 的不动点，
Kutta 映射斜率 b ≈ 0.93（§4，P2 实测），所以

    Γ/Γ*  =  F (1 − b) / (1 − F b)                                       (4.7)

——F = 0.8 给出 0.21×，不是 0.8×。这就是为什么实测 q（≈3.3）远超 F 的
朴素指数，也是 r_c 必须保持小的原因。

**结果。** 边缘峰值抹平（粗/中/细上 p：0.592 → 0.009），M6 细网格成为
**真正的离散解**（0 限幅/触底）——P9/G9.1 没能做到的事——M∞0.84 下
medium 的 M_max 从 2.13 降到 1.725。**代价：** ~1–1.6% cl 的蓄意翼尖卸载
（随 r_c 缩放），局部在 η≈0.95 以内。证据：
`phases/p1/cases/demo/p13_tip_edge_singularity/`。

**警示（P13/G13.3）：第三个奇点，在壁面上。** 去掉翼尖奇点并**没有**交出
3D 网格收敛——M6 升力序列仍不在渐近区（增量随均匀加密*变大*，这是仍有
奇点在被分辨的特征）。自相似阶梯上的三区盒研究定位了剩下的东西，而且它
**不是**中展向/前缘分辨率（早先对展向"晃动"的一种读法，已撤回——那晃动
是比较两个不同网格族的伪迹，见下）：

| 区域 | coarse_ss → medium → fine | p | 判定 |
|---|---|---|---|
| **翼尖帽边缘（壁面）** | 0.662 → 0.824 → **1.015** | **+0.321** | **发散** |
| 尾迹自由边（§4.2 的修法） | 0.536 → 0.565 → 0.570 | +0.045 | 有界 |
| 机翼 p99（对照） | 0.642 → 0.628 → 0.629 | −0.014 | 有界 |

§4.2 的修法成立，普通机翼场收敛；发散的是 `meshgen/wing3d.py` 的**平翼尖
帽**与上/下表面相交的尖凸边。真实 ONERA M6 的帽是**圆**的——平帽是有
记录的蓄意简化——所以这正是 §4.1–§4.2 处理的那类边缘奇点，只是在**机体**
上而不是尾迹上。**它不是 P11（曲面壁面*单元*）问题：** 等参单元无法正则化
一个真正尖的几何*棱边*；错的是几何本身，所以修法是 Track M（把帽做圆）。
证据：`phases/p1/cases/demo/p13_tip_edge_singularity/run_g133_ladder.py`
（5/5）。

**网格阶梯前置（2026-07-13 修复）。** 任何 M6 三点研究必须跑在
`generate_onera_m6.RICHARDSON_LADDER = (coarse_ss, medium, fine)` 上。原
`coarse` **不在加密射线上**：`_level_params` 的所有长度都从 `h_wall` 导出，
唯独远场被钳（`h_far = min(2.5, 120·h_wall)`），而钳制只在 `coarse` 档
咬合——于是 coarse→medium 把远场加密 1.39× 而壁面是 2×。`coarse_ss` 是
去掉钳制重切的同级；阶梯现在每个长度尺度都精确加密 2.000。这作废所有更
早的 M6 三点 Richardson，**包括 P9/G9.1 的**。

---

## 5. 边界条件

| 边界 | 条件 | 实现 |
|---|---|---|
| 固壁（机翼、机身） | ρ ∂φ/∂n = 0 | Galerkin 弱形式中的自然（do-nothing）条件——面积分为零。无罚、无鬼单元。 |
| 远场 | φ → φ∞ + φ_vortex | 外节点 Dirichlet：φ∞ = x cosα cosβ + y sinβ + z sinα cosβ（迎角 α、侧滑 β 的均匀流），加一个总环量 ΣΓ_j Δs_j 的可压缩马蹄涡修正（Prandtl–Glauert 缩放）。R ≳ 25–50 弦长时修正很小，但带上它域可缩到 ~15 弦长。 |
| 尾迹 | (4.1)–(4.2) | 约束消去，§4。 |
| 对称面（半模） | ∂φ/∂n = 0 | 自然，免费。 |

一个细节：全 Neumann 壁面 + Dirichlet 远场下系统适定；若哪天用纯 Neumann
变体（如管道流），钉住一个节点。

**Track V transpiration BC（已实现，粘性耦合）。** VII 耦合（Drela IBL3；
phase 1 Track V 与 phase 4 GS4.1 交付，代码在 `pyfp3d/viscous/`）把齐次壁面
条件换成壁面面元上的吹气源 ρ ∂φ/∂n = d(ρ_e u_e δ*)/ds，并在尾迹片上加类似
的 δ*_wake 质量源——Picard 层只动右端项，不动网格，δ* = 0 与无粘路径逐位
一致。**现行能力警示：** 有粘耦合在举力机翼上跑不起来（D09/D10：第 0/第 1
次外迭代即死）——以能力边界文书为准。

### 5.1 壁面几何误差的边界通量修正（G1.6 候选修法；门 G1.3–G1.5 + DP1）

> **P11 勘误（2026-07-19，实测——§5.1.3 结局 + PROJECT_STRUCTURE.md
> "Known gaps" P11 块）：** 下文的 variational-crime 根因已被**推翻**——
> 该 crime 占 medium 网格误差的份额 ≈0.2 pp；11.6% 是 P1 场在 h=0.08 的
> 固有 max-范数能力，而次一阶的扫描特征是固定体网格地板伪迹。
> §§5.1–5.1.2 原样保留，作为当时的设计记录。

背景（门 G1.6，原 G1.2，不可压球）：medium 网格 Cp 误差 ~11.6%，目标
<2%，根因（见 PROJECT_STRUCTURE.md "Known gaps"）归为几何/变分不一致，
而不是表面梯度恢复。分部积分后 P1 Galerkin 壁面项 ⟨v, ∇φ·ñ⟩ 作为自然 BC
被丢掉——它强制的是「穿过**平直面元**（法向 ñ）的通量为零」，而物理条件是
「穿过**真实曲面**（法向 n）的通量为零」。面元法向与真法向差 O(h)，产生
壁面速度/Cp 的一阶几何误差；原始节点势确实呈现次一阶壁面收敛，与该机制
一致。先前记录的结论是「真修需要曲面/等参壁面单元，一个单独圈定的工作」。
下面的文献普查找到了足迹小得多的中间路线；应先验证它们，再决定曲面单元
的功夫是否还要。选项按实现足迹排序。（长期禁令不变：不要再提进一步
h-加密或恢复格式微调——两者都已有证据排除。）

> **验证结局（2026-07-06，见 §5.1.2）：** G1.3/G1.4 oracle 实验排除了
> Option A——贴体网格上它要修的边界数据缺陷（近）为零，实测天花板
> ≈11.3%，对 <2% 目标。DP1 走了 ">5%" 分支：Option C + 单独圈定的曲面
> 单元。§§ Option A/B 保留为设计记录。

**Option A（原推荐、最先实现——已被 G1.3/G1.4 oracle 排除，见 §5.1.2）：
真法向弱通量修正（滞后通量修正）。**

出处（两条独立工作线，同一实用配方）：

- Krivodonova & Berger（JCP 2006）：直边网格上用*物理*几何的法向而不是
  *计算*几何的法向施加固壁条件；对 Euler 方程这显著改善解质量而不需曲面
  网格。Ciallella, Gaburro, Lorini & Ricchiuto（Appl. Math. Comput. 2023）
  把同一多项式修正推广到一般 2D/3D 边界条件。
- Shifted Boundary Method（Main & Scovazzi 2018；Atallah, Canuto &
  Scovazzi 分析系列）与更早的 Bramble–Dupont–Thomée（1972）边界值修正：
  把边界条件的**数据**修正到近似边界上，以抵消几何一致性误差。

数学形式。在壁面面元上的点 x，记 p(x) 为 x 到真曲面（球）的最近点投影，
n = n(p(x)) 真单位外法向，ñ 面元单位外法向。分解

    ñ = (ñ·n) n + t,    t := ñ − (ñ·n) n    （t 在真切平面内，|t| = O(h)）

于是 ∇φ·ñ = (ñ·n)(∇φ·n) + ∇φ·t。物理边界条件 ∇φ·n(p(x)) = 0 消去第一
项；壁面项不再为零而是 ⟨v, ∇φ·t⟩。以**滞后（Picard）**方式移到右端：

    ∫_Ω ∇φ^{k+1}·∇v dV = ⟨ v,  ∇φ^k · ( ñ − (ñ·n) n ) ⟩_{Γ_h,wall}
    （其他边界条件不变）

实现注记（对最终实现有约束力）：

- 刚度矩阵**完全不变**：保持 SPD，AMG 层级与预条件子原样复用。唯一改动
  是壁面面元上的一个右端组装循环。
- 新几何数据：面元积分点的最近点投影与真法向。球有解析式
  n = (x − c)/|x − c|；日后机翼阶段来自 CAD/解析几何或细参考面——接口
  设计成可替换的 `closest_point_normal(x)` 回调，预计算进 SoA 数组（按
  agent-rules 的 Numba 硬约束：热循环里无 Python 对象操作）。
- 积分：每面元 3 点 Gauss（或现有面积分用的规则）；∇φ^k 是相邻单元的
  分片常值梯度。
- 不动点循环：修正是 O(h)，迭代收缩。纯 Laplace 验证跑固定 3–5 次外迭代
  （同矩阵重复解——便宜）；进入全势范围后折进现有 ρ-Picard 循环，零边际
  成本。
- 均匀流保持检查：均匀流 φ = U∞·x 不满足球壁面条件，所以带壁网格上出现
  非零修正项是预期的；但必须确认 V0 均匀流门
  （`tests/A/test_A01_freestream_preservation.py`，用无壁/全远场配置）永远
  不触发这条修正路径。
- 前置：壁面面元绕向一致——`_wall_vertex_normals` 里现有的绕向断言是本
  修正的前提，不许删。

#### 5.1.1 圆柱预研（门 G1.3，原 G1.2-a0）

`cases/meshes/cylinder_2.5d/` 是 Option A/B 路线的**指定快速试验台**，先于
球使用（门 G1.3 先于 G1.4）。理由：它呈现**同一个**曲壁 variational
crime，且已量化——max |Cp err| 0.091（粗）→ 0.045（中），~O(h)
（`tests/A/test_A14_meshgen_cylinder.py`）；Option A 需要的每个几何要素都有
闭式；网格便宜（6.9k / 17.3k tets）；诊断是一维曲线 Cp(θ)。

圆柱闭式（半径 a，轴 z，来流 U 沿 x，r_xy² = x² + y²）：

    φ_exact = U x (1 + a²/r²),    p(x) = a (x, y, 0)/r_xy,    n = (x, y, 0)/r_xy

§5.1 的 `closest_point_normal(x)` 回调对圆柱和球各有一个解析实现（球：
n = (x − c)/|x − c|）；接口设计本身不变。

警示（三条都有约束力；引用本预研时一条都不许丢）：

1. **必要不充分。** 圆柱单曲率，球双曲率。**G1.6 门只在球上关闭。** 圆柱
   结果只作为进入球实验 G1.4 的前置证据。
2. **展向噪声地板。** quasi-2D 网格带有 3-tet 棱柱劈分固有的 O(h) 展向
   噪声（粗网格 max |w|/U∞ ≈ 2.9e-2），以同样 O(h) 污染面内梯度恢复。
   圆柱验收判据因此是：修正后误差显著低于未修正**且** Cp 误差收敛阶从
   次一阶恢复到 ≈一阶——**不设绝对阈值**。oracle 运行顺带测这个地板；
   其量级反馈进 G2.5(b) 的重定规格。
3. 圆柱吸力峰 Cp = −3 比球的 −1.25 更强——更苛刻的压力测试；两个几何
   之间的百分比数字不许互相比较。

可视化：单层网格上切片退化为平凡操作——对称面自带 2D 三角剖分，场直接
tripcolor。辅助函数按 P2 `post/section_cut.py` 的最终接口原型定义（签名
预留 z = const 参数）；P2 再加一般 3D 插值路径。

与被否的 Nitsche/罚原型划清界限（明确记录，防止被误归类为重复实验）：
早先尝试改的是同一条件在同一错误（平直面元）几何上的**施加机制**；
Option A 修正的是边界条件的**数据本身**——最近点投影引入真法向。这才是
SBM 中 "shift" 的实质贡献，与施加机制无关。

理论预期与天花板（记录在案以管理预期）：

- SBM 分析表明对 P1 单元——分片常值梯度、零 Hessian、没有二阶 Taylor
  展开可用——朴素 shifted Neumann 条件在 L² 掉一阶，但 H¹ 半范保持最优
  一阶。Cp 是受 H¹ 控制的梯度量，所以这个天花板不挡本项目的路。
- 本项目用**贴体**网格（顶点在真曲面上）：几何缝隙厚 O(h²)，远好于
  unfitted SBM 的 O(h) 缝；一阶法向旋转修正因此是主导误差项。
- 预期：Cp 恢复到接近一阶收敛，medium 网格误差有望降到 2% 以下；确切
  天花板由 oracle 实验（门 G1.4）实测。

#### 5.1.2 预研结局（2026-07-06）：贴体网格上 Option A 被排除

当天完成的 G1.3 圆柱预研与 G1.4 球 oracle 否证了上面的预期（完整证据：
roadmap G1.3/G1.4 条目；`tests/C/test_C02_cylinder_wall_correction.py`；
oracle 结果在 `cases/demo/p1_laplace/results/`）。★ 原稿此处引用的
`artifacts/G1.3/` 已不存在：`artifacts/` 目录于 2026-08-24 删除——它是
gitignored 的本地输出、从未进 HEAD，按证据纪律它从来不算证据；门证据
（PNG/CSV）现写到被跟踪的 `cases/gates/<gate_id>/`，旧 artifacts 内容可用
对应测试本地重生成。事后看机制是初等的：

- 对**调和**势与**贴体**壁面顶点，平直面元与真曲面之间的区域恰由这两张
  面封闭；真曲面带零通量，由散度定理精确解**穿过每个面元的净通量恰为
  零**。自然 BC 的一致性缺陷 ⟨v, ∇φ_exact·ñ⟩ 只是一阶矩残差——粗圆柱上
  实测 ~2e-5 max，按 ~O(h⁴) 收缩。SBM/BDT 修正在 O(h) 缝的 *unfitted*
  边界上挣钱；O(h²) 的贴体缝上没有东西可修。
- Option A 的 t 形式 ⟨v, ∇φ_exact·t⟩ 在圆柱上额外组装到机器零（均匀分布
  圆上相邻面元精确抵消）；组装与手算单面元算例核对过，所以这个零是真的。
- 实测天花板：圆柱 Cp 误差不变（粗/中/细 max 9.10e-2/4.49e-2/2.22e-2，
  斜率 1.02，修正 ≡ 未修正）；球 medium 网格 max |Cp err| 0.1156 →
  0.1133（完整恢复一致性缺陷）——**11.3% 天花板**，对 <2% 目标。
- 改写圆柱案例定位的附带发现：其每个层级 ~76% 的 Cp 误差可以通过把*精确*
  势喂给表面恢复复现——quasi-2D 单层条带上细长三角形的壁面剖分主导
  恢复误差，而壁面节点 φ 以健康的 ~1.2 阶收敛。圆柱因此**并不**呈现球的
  次一阶病理，不再是该病理的指定试验台（§5.1.1 的「同一 variational
  crime」理由在这一点上作废；网格与解析 Cp 验证仍留作 M0 端到端检查）。

后果（DP1 ">5%" 分支，roadmap 已记录）：球的剩余误差归因于域扰动
（∫_Ωh vs ∫_Ω 缺的面元/曲面薄片）加 P1 梯度近似——正是 Gap-SBM 的缝隙
积分建模的项（下面的 Option B），所以 B 作为曲面单元设计复查的可选预研
材料存活；但 O(h²) 缝隙厚度下其收益预期是二阶小。精度路线是曲面/等参
壁面单元，单独圈定——**后来的 P9 阶段**（G9.1 球 Cp <2%；它还接手
V6 <1% 残余地板，G9.2），G1.6 门按 Option C 重定义。*（2026-07-19 执行为
phase P11——实测**负结果**，上文球的「域扰动」归因被 §5.1.3 结局取代：
主导项是朴素的 P1 近似，域项 ≈0.2 pp。）*

**Option B（Option A 不够时的升级——按 §5.1.2 现仅为可选预研材料）：
Gap-SBM 缝隙修正。**

Collins, Li, Lozinski & Scovazzi, "Gap-SBM"（arXiv:2508.09613, 2025）正好
针对上面的 P1 Neumann 次优性。它由真边界与计算边界之间的距离图构造近似
缝隙几何，把解与检验函数延拓进缝隙，对修正后的变分形式施加近似积分；无
额外自由度、无切割单元、无鬼罚项，Neumann 问题有最优 L² 与 H¹ 收敛证明。
具体足迹：壁面面元面积分乘一个缝隙厚度系数，加壁面节点处几个距离向量
修正项——可实现为 SoA 数组 + 一个面元循环，Numba 友好。警示：论文分析
是 2D（作者声明 3D 是直接推广），且贴体情形下缝隙系数是 O(h²)，这些额外
项本来就小——这正是「Option A 优先」的理由。

**Option C（务实兜底）：重定义门而不是改格式。**

若曲面单元最终仍需自己的功夫，把 G1.6 验收判据重定义为「相对**几何一致
参考解**（同一多面体域上的高精度参考，如 BEM 或超细网格）的 Cp 误差/
收敛阶」，把几何模型误差从代码正确性验证中剥离。这符合 §10 验证阶梯不
混淆模型误差与代码误差的原则。本选项不提高物理精度，只调整验证标尺。

#### 5.1.3 P11 实现设计：曲壁面相邻单元（2026-07-19）

DP1 ">5%" 分支的曲面单元路线，作为 phase P11 开（roadmap track_p.md §P11；
范围 = 先做球/Laplace 腿，见
`phases/p1/docs/analysis/next_phase_priorities_2026-07-19.md`）。设计取能移动
*积分域*——G1.3/G1.4 证明真正错的东西——的最小改动，同时不动自由度、
稀疏性、SPD 性与任何已提交求解路径：

- **几何：只在贴壁层用二次（tet10 式）映射。** 收集壁面三角剖分的边；把
  每条这种边的中点经现有可替换 `closest_point_normal(x)` 回调（§5.1，
  `solve/wall_correction.py`——球/圆柱解析）投影到真曲面。含 ≥1 条壁面边
  的每个 tet 成为*曲*单元：其 10 节点几何映射在壁面边上用投影中点、其余
  用直中点。任何共享一条投影边的 tet 按构造含壁面边，所以共享面从两侧
  拿到相同的二次几何——曲层几何上 conforming，无缝无叠。顶点本来就贴体；
  中点偏移 O(h²)（球上 ≈ h²/8R），于是域误差从 O(h²) 降到 O(h³)，且——
  G1.6 的真正机制——壁面法向误差从 O(h) 降到 O(h²)。
- **场：映射 P1（超参）。** 标量 φ 每 tet 保持 4 个顶点自由度；曲单元上
  基函数是参考线性函数与二次映射的复合。单元刚度用积分，
  K_ij = Σ_q w_q (J⁻ᵀĝ_i)·(J⁻ᵀĝ_j) |det J(ξ_q)|，ĝ_i 为常值参考梯度，J 为
  二次映射的（ξ 线性）Jacobian；默认规则 = 4 点 2 次 tet 规则，与 5 点
  3 次规则交叉核对（G11.4）。每个积分点断言 det J > 0。
- **组装：标准矩阵上的稀疏增量。** A_curved = A_P1 + ΔA，
  ΔA = Σ_curved (K_quadrature − K_flat)，同样的 4×4 节点块 ⇒ CSR 模式是
  A_P1 的子集，和保持对称且（带 Dirichlet 行）SPD，CG+AMG 求解路径原样
  复用。以 opt-in 的 `stiffness_delta` 参数暴露在 `solve_laplace` 上（仿
  `body_source_rhs`）；默认路径——因此每个已提交结果与 V0 均匀流门——
  逐位不变。
- **预登记风险——超参一致性。** 曲单元上的映射 P1 空间不再精确包含物理
  线性多项式（只差到 O(h²) 映射扰动）；曲层内均匀流精确性从机器零退化
  为 O(h²)。这只限 opt-in 路径（V0 配置无壁）。若它占主导——判别量 =
  门 G11.2 的 φ 收敛阶扫描加平投影零测试 G11.3（平面
  `closest_point_normal` ⇒ ΔA ≡ 0 逐位）——记录的兜底是等参 P2 壁面层
  （新边中点自由度），它自己是决策点，不是悄悄升级。
- **此处范围外**（记录在案，不是忘了）：尖几何棱边（§4.1——曲单元无法
  正则化它们，Track M 负责）；可压/跨声速组装路径与 Kutta/尾迹 ↔ 曲单元
  耦合（同类的翼身/机身应用，球腿落地后的独立阶段）。

**结局（2026-07-19 实测，P11 收口——完整记录 roadmap track_p.md §P11；
证据 `cases/demo/p11_curved_walls/`）。** 这条路线是有记录的**负结果**，
且预登记风险正是机制：映射 P1 的梯度一致性误差是 **O(h)**（粗球上线性
复现偏差 max 0.138），与曲化去掉的面元法向误差同阶，所以 medium 球 Cp
只移动 11.56% → 11.33%——与 G1.4 边界数据 oracle 天花板相同。对照组
重新归因了 G1.6 本身：同样平直面元的结构化壳层壁面-φ 以 ~2 阶收敛
（E6），h_min 扫描的阶塌缩是固定体网格污染地板（E8：h_min=0.03 只加密
远场，壁面 φ 误差降 3.17×，阶恢复到 1.89）。「medium 上 max 2%」判据
要求 h = 0.08 处壁面速度 O(h²) 精确——超出任何网格上任何 P1 场方法的
能力；记录的三岔（Option C 重定规格并附实测可通过形式 / 等参 P2 壁面层 /
接受限制）由使用者裁决（Option C 后于 phase 1 采纳）。不要再提超参
（映射 P1）曲单元。

---

## 6. 空间离散

**线性四面体上的 Galerkin FEM（P1）。** P1 tet 上它代数上等价于顶点中心
median-dual 有限体积格式，但 FEM 视角实现最干净：不用构造对偶网格度量，
没有面循环——只有单元循环。

### 6.1 弱形式与逐单元装配

求 φ（满足 Dirichlet + 尾迹约束），使对所有检验函数 N_i，

    R_i = Σ_e  ρ̃_e ( ∇φ_e · ∇N_i ) V_e  = 0                            (6.1)

∇φ_e = Σ_k φ_k ∇N_k|_e。每单元预计算：体积 V_e 与 4×3 形函数梯度矩阵
B_e（常值）。之后每次非线性求值、每单元：一次 4×3 gemv 得 ∇φ_e、密度律、
迎风查询、一次 3×4 gemv 散布——天生 numba 友好。

### 6.2 Picard（冻结密度）矩阵

    A_ij = Σ_e ρ̃_e (∇N_i · ∇N_j) V_e                                    (6.2)

——加权的刚度矩阵。ν ≡ 0 时 SPD。有迎风时 ρ̃_e 耦合到上游单元的自由度；
**Picard** 里我们蓄意在 A 中忽略该耦合（把 ρ̃_e 当冻结标量）——A 保持
对称，迎风物理通过残差/右端滞后进入。这是经典、极鲁棒的格式；其收敛速度
随激波强度退化，这是 §8 Newton 选项的动机。

### 6.3 Newton 精确 Jacobian（需要 P7 可微通量；Newton 求解本身是 P8）

(6.1) 对 φ_k 求导得 López 的两项形式（学位论文 Eq. 3.24，§13）：

    ∂R_i/∂φ_k = Σ_e V_e [ ρ̃_e ∇N_i·∇N_k                     ← Term 1（Picard，SPD）
                + (∂ρ̃_e/∂φ_k) (∇φ_e·∇N_i) ]                  ← Term 2（密度）   (6.3)

Term 1 正是冻结密度 Picard 矩阵 (6.2)——现成的着色装配器原样产出。Term 2
带密度敏感度；沿等熵律与开关函数展开（López Appendix B，已对 PDF 核实），
按流态分支，超声速加速情形为

    ∂ρ̃_e/∂φ_k     = (∂ρ_e/∂φ_k)(1 − μ) − (ρ_e − ρ_up)(∂μ/∂φ_k)          (B.3)
    ∂ρ̃_e/∂φ_k^up  = +μ (∂ρ_up/∂φ_k^up)                                   (B.4)

其中 ∂ρ/∂φ = (∂ρ/∂u²)(∂u²/∂φ)，∂ρ/∂u² = −(ρ∞/2a∞²)[…]^{(2−γ)/(γ−1)} (B.8)，
∂μ/∂φ = (**+**M_c²/M⁴)(∂M²/∂u²)(∂u²/∂φ)（B.9–B.12；**符号 2026-07-10 由
P7 G7.3 FD 门更正**——μ = 1 − M_c²/M² 有 dμ/dM² = +M_c²/M⁴，而
dM²/du² > 0，整条链为正（μ 随速度增大，理应如此）；此处早先的「−」是
转录笔误。现行导数（`kernels/upwind.py::rho_tilde_sensitivities_sweep`，
P7）用 + 号，在每种流态下与冻结选择上现行通量的中心差分吻合到 ~4e-10，
`tests/A/test_A03_differentiated_flux.py`）。**(B.4) 非零**：当前单元的残差
依赖*上游*单元的自由度，所以 Newton 模板**比 Picard 宽**——稀疏图与单元
着色必须由单元 + 上游图重建。

**模板宽度在 UP3D 不是「一层单元」（2026-07-08 更正）。** López 的*单跳*
上游给出宽一层的模板；UP3D 走得更远，因为薄片棱柱劈分 tet 需要多跳行程
（§3：单跳只达到流向尺度的 ~0.37）。**默认 Newton 通量是冻结行走**——
对现行迎风密度在固定选择上求导，不是换通量：

- **冻结行走（默认）**（`mode="walk"`，u(e) 是 ≤4 跳外的单个单元）：
  Term 3 把 e 耦合到恰*一个*上游单元——nnz 增量温和（~+1 单元/行，最
  接近 López），但图距离 ≤4，所以着色/CSR 图要加远距离 e→u(e) 边（色数
  更高）。冻结选择下本来就可微（等熵 ρ + 分支 ∂ν/∂φ，B.3–B.6）；`max_ε`
  只是可选，用来阻尼远离解处的活跃集抖动（López 的 C⁰ `max` 在解附近——
  选择与活跃集不再变化处——已能达到二次收敛）。这就是 P7 交付物。
- **流线-高斯核（可选）**（`mode="kernel"`）：ρ_up 是*整个深度 3 BFS
  邻域*的加权和，所以 Term 3 把 e 耦合到该邻域的*每个*单元——实质性
  **更密**的 Jacobian（行 nnz ~O(邻域大小)，不是 +1）；权重依赖
  w_f(∇φ_e) 只增加自耦合。更光滑、无冻结选择，但系统更密——一条 Picard
  加速路径，只有实测划算才升格为 Newton 通量。

旧的「~30% 内存 / 一层」估计已退役；在承诺更宽的着色/CSR 重建之前，先在
N2/N3 取**实测的** nnz/GMRES 迭代/AMG setup 数字。

来自 PDF 的两个塑造实现的事实：

- López 保留开关函数导数 ∂μ/∂φ（B.3、B.6）——这是*完整* Jacobian
  （「strategy A」）；丢掉它（只留 Term-1 + 局部密度）会把二次收敛降为
  超线性。
- 推导把**上游选择 u(e) 当冻结**，只对 ρ_up 求导；选择每个 Newton 步
  重估，但不求导。这正是为什么只需让通量在固定选择下可微（§3.1），而不
  必对选择本身求导。

精确 Jacobian（Term 1 + Term 2 + 上游耦合）**非对称**且在超声区不定——
用 GMRES + AMG（在对称部分上聚集）或 GMRES + ILU，不用 CG（§8）。

精度阶：亚声区二阶，捕捉激波处局部一阶——标准且可接受。

---

## 7. Numba 内核架构

设计规则（血汗换来的 numba 约束——已编入项目 skill）：

1. **只要 Struct-of-arrays。** 所有网格/解状态用扁平、连续、显式类型的
   `np.float64` / `np.int32` 数组。`@njit` 内不要 Python 对象、不要
   dict、不要 dataclass 实例。一个薄薄的 Python 层 `Mesh` dataclass
   *持有*数组；内核以参数收数组。
2. **`prange` 用两遍或着色装配。** 散布-加竞态是并行化的中心隐患：
   - 残差：贪心**单元着色**（同色任两单元不共节点）；外层串行循环颜色，
     内层对该色单元 `prange`。着色启动时算一次（tet 典型 ~8–12 色）。
   - 矩阵：从符号稀疏模式一次性预计算 `elem_to_csr[e, 4, 4] → nnz 索引`
     映射；装配写入随之落到不相交或同色的位置。
3. **所有内核 `@njit(cache=True, fastmath=True)`**；`parallel=True` 只用在
   profile 证明要紧的地方（单元循环；不用在小的 BC 循环上）。
4. **热内核内零分配。** 残差、∇φ 工作区（逐色或逐线程）、密度数组全部
   预分配传入。
5. **物理用纯函数。** `density(q2, minf2, gamma)`、`mach2(q2, ...)`、
   `cp(...)` 写成标量 njit 函数——可单独单元测试，后处理复用。
6. 线性代数留在 **SciPy/PyAMG 地界**（已编译，无需重写）：
   `scipy.sparse.csr_matrix` + `pyamg.smoothed_aggregation_solver` 作
   `scipy.sparse.linalg.cg`（Picard/亚声速）或 `gmres`（Newton/跨声速）
   的预条件子。可选的 matrix-free numba matvec 以后再说。

**模块布局（现行树，仓库根 `pyfp3d/`）**

```
pyfp3d/
├── _compat.py           # numpy/numba 版本兼容垫片
├── mesh/
│   ├── reader.py        # meshio 前端：gmsh/.su2/.cgns → 裸数组
│   ├── wake_cut.py      # 沿带标签尾迹面的节点复制；TE 站位→尾迹线映射
│   ├── metrics.py       # numba：B_e、V_e、面邻接、质量检查
│   ├── coloring.py      # 并行装配的贪心单元着色
│   └── manifest.py      # 网格 manifest 指纹（provenance，F01）
├── meshgen/             # 参数化网格生成器：extrude / planar / wing3d /
│                        # wingbody / fuselage / structured（M6 C-H）
├── physics/
│   └── isentropic.py    # numba 标量：ρ(q²)、M²(q²)、a²、Cp、q*²、ρ′
├── kernels/
│   ├── gradient.py      # 逐单元 ∇φ_e
│   ├── upwind.py        # 上游单元 id、ν_e、ρ̃_e（行走 + 可选核；P7 敏感度）
│   ├── entropy.py       # 熵修正（GS1b.3）
│   ├── residual.py      # R 的着色装配 (6.1)
│   └── jacobian.py      # Picard 矩阵 (6.2)；Newton (6.3)
├── constraints/
│   ├── dirichlet.py     # 远场值，含涡修正
│   ├── wake.py          # 主–从消去、Γ 处理、Kutta 更新 (4.4)、翼尖渐变
│   └── te_pressure.py   # P14 压强相等 Kutta 估计器 (4.4′)
├── solve/
│   ├── linear.py        # CG/GMRES + AMG/ILU 包装、容差
│   ├── picard.py        # 外循环：ρ̃ 更新 → Γ 更新 → 线解 → 松弛
│   ├── newton.py        # 精确 Jacobian、全耦合 (φ,Γ)、Mach 延拓
│   ├── continuation.py  # 难算例的 Mach/α 爬坡
│   ├── wall_correction.py  # §5.1 边界通量修正（oracle 路径）
│   ├── curved_wall.py   # §5.1.3 曲壁面单元刚度增量（P11，负结果保留）
│   └── timing.py        # 计时仪器
├── post/
│   ├── surface.py       # 节点 Cp、剖面 cl、力与力矩、Kutta–Joukowski 对照
│   ├── shock.py         # 激波位置诊断
│   ├── section_cut.py   # z = const 切片
│   └── vtk_out.py       # meshio VTK 导出：φ、V、M、ρ、Cp
├── viscous/             # Track V：IBL3、transpiration、松/紧耦合（§5）
└── （回归算例在仓库根 `cases/`，不在包内）
```

---

## 8. 非线性求解策略

**基线：松弛 Picard（phase 1–5 的主力，仍是冷启动/种子路径）。**

```
assemble sparsity, coloring, wake elimination map      # 一次性
φ ← farfield initial guess (uniform flow)
repeat:
    ∇φ_e, q²_e, ρ_e            # 单元扫描
    u(e), ν_e, ρ̃_e             # 迎风扫描
    A(ρ̃), b(Dirichlet, Γ)      # 着色装配（ρ̃ 再滞后时只装右端）
    solve A δφ = −R             # CG+AMG（对称）/ GMRES+ILU（加非对称件时）
    φ ← φ + ω δφ                # 亚声速 ω ≈ 1；跨声速 0.7–0.9
    Γ_j ← Γ_j + ω_Γ (Kutta target − Γ_j)
until ‖R‖₂/‖R₀‖₂ < 1e−10 and ‖ΔΓ‖∞ < tol
```

预期行为：亚临界算例 O(10–30) 次 Picard 迭代收敛；跨声速 O(100–300)。
AMG 每次线解 O(N)——工作站 ~1–3 M 节点总耗时分钟级，这是设计目标。

**加速项，按实现价值排序：**

1. Picard 迭代间复用 AMG 预条件子（每 k≈5–10 步重建）。
2. 非精确求解：早期松线性容差（相对 1e−2），接近收敛收紧
   （Eisenstat–Walker 式调度；`forcing` 参数）。
3. Mach 延拓：先收敛 M∞ − 0.05，再重启。FP 非唯一带附近与带激波的
   ONERA M6 必需。
4. 完整 Newton，(6.3) + 伪瞬态阻尼——终端二次收敛；只在 Picard 刀枪不入
   之后。

### 8.1 Γ 作为未知量的全耦合 Newton（P8）

Picard/割线架构嵌套三层循环（Mach 延拓 → 逐站位 Γ 割线 → 密度 Picard）。
P5 medium 失败把脆弱性根因定到**割线–密度耦合**：割线从尚未收敛的密度场
上读 Kutta 失配，坏的斜率估计过冲 Γ，迫使更多密度工作，再污染下一次失配
——一推就散的正反馈（P5 的修法是早停 + 固定-Γ 抛光把它正则化，不是
治愈）。López 完全避开这类失败：解一个 **Γ 与 φ 并列未知量**的 Newton
系统，φ 和 Γ 一起收敛，没有读滞后残差的割线。

UP3D 的主–从归约让这便宜，因为 Γ-Jacobian 块几乎已经在代码里。由
φ_full = T·φ_red + g(Γ)（尾迹跳跃），且远场 Dirichlet 数据也通过涡修正
带 Γ，耦合步为

    ⎡ TᵀJ T        ∂R_red/∂Γ ⎤ ⎡ δφ_red ⎤     ⎡ −Tᵀ R      ⎤
    ⎢                        ⎥ ⎢        ⎥  =  ⎢            ⎥
    ⎣ ∂F/∂φ_red    −I        ⎦ ⎣ δΓ     ⎦     ⎣ −F         ⎦

其中

- `TᵀJ T` 是归约完整 Jacobian (6.3)，Term 1 复用 `assemble_matrix`；
- `F_j(φ_red) = kutta_targets_j(φ_red) − Γ_j` 是 Kutta 闭包残差
  （`wake.py::kutta_targets`，原样复用），∂F_j/∂φ_red 是稀疏平均矩阵
  （该站位上/下 TE 探针节点上 ±1/n_j）；
- **∂R_red/∂Γ_j = TᵀJ g_j**（g_j = 站位 j 从属指示向量）。Picard 层
  TᵀA g_j = h_j，`wake.py` **已经物化**为 `self._h`；精确 Newton 列用
  完整 J，可以 lag 到 A 作为有记录的近似。
- **远场列。** Γ 以*第二种*方式进入 φ_full：远场涡 Dirichlet 值对 Γ 线性
  （`dirichlet.py`，`farfield_spanwise_gamma` 渐变或展向均匀均值）。所以
  ∂R_red/∂Γ_j 多一项 −A_coupling·(∂vals_red/∂Γ_j)。**这项容易漏**——它
  不在 Picard 代码路径里（在那里它被悄悄折进右端）——必须加进 Newton
  列。

全局化（López 核实，学位论文 §4.4/§4.8）：从均匀流直接裸 Newton 跨声速
发散，所以包进 **Mach 延拓 + 载荷步进**。对早先内部调度提案的两条更正，
已对 PDF 核实（Tables 4.7/4.8/4.13）：

- 单个算例内 López 保持 M_crit 与 μ_c **固定**，只爬 M∞；算例越难，整个
  爬坡选*越低*的 M_crit 和*越高*的 μ_c（Case 1：M_crit 0.99、μ_c 1.0；
  Case 4 α=2°/M∞0.75：M_crit 0.90、μ_c 1.1）。**没有**逐步的
  M_crit=0.99→0.90 扫描。
- ONERA M6（Table 4.13）：12 步，M_crit **恒定 0.95**，M∞ 爬 0.50→0.84
  （步 1–8）期间 μ_c 保持 2.0，然后 M∞ 固定、μ_c **降** 2.0→1.6
  （步 9–12）在场建立后锐化激波。所以 μ_c 调度是*到目标后的耗散削减*，
  不是爬坡中加大。

内层 GMRES 加 Eisenstat–Walker forcing 调度（早期松、近收敛紧——现有
`forcing` 参数是单参数替身），AMG 层级跨几个 Newton 步复用。López 数据
给出的预期成本：亚声速总共 5–10 次 Newton 迭代；跨声速每载荷步 4–9 次；
ONERA M6 ≈ 12 步 × 5–9 ≈ 60–110 次总 Newton 迭代（Tables 4.5/4.6/4.9，
§4.7）——对照当年的 ~10⁴ 次 Picard 迭代。

N-Γ 分裂（保留 Γ 割线外循环，只对密度内解做 Newton）保留**仅为兜底**：
它保留割线，因此保留 P5 耦合不稳定风险。

---

## 9. 后处理与力

- 节点 q²、M、Cp 由体积加权单元平均（以后可换超收敛片恢复）。
- 表面力：壁面三角形上积 Cp·n → C_L、C_D_pressure、C_m。**警告：** FP
  近场压阻被等熵激波模型与离散污染；照报，但可信的分解是 Γ(s) 的
  Trefftz 平面诱导阻力（需要时再加激波扫掠熵代理的波阻估计）。
- 展向 Γ 给出剖面 cl：cl = 2Γ/(c(η))——与表面积分的直接交叉核对；收敛
  解上两者必须 ~1% 相符（尾迹/Kutta 代码的极好 bug 探测器）。
- 每 k 次迭代 VTK 输出，监控激波位置。

### 9.1 壁面 Cp 梯度恢复——G6.1 锯齿的修法（P6）

壁面速度是 φ 的面内切向梯度（自然 BC，§5）。默认恢复是**逐三角形常值**
（`triangle_tangential_gradients` / `wall_cp_curve` /
`wall_force_coefficients`，无平均）——选它是为了尖 TE 折痕不被两侧切平面
混合污染。它在棱柱劈分薄片壁面三角剖分上的代价就是 ≈2 单元的**表面 Cp
锯齿**：相邻细长三角形即使对光滑 φ 也给出振荡的常值梯度。它曾被误归到
人工密度通量；它是**恢复伪迹**（N1，2026-07-08——§3.1：对*同一*解做
节点平均，G6.1 指标降 ~330×；更光滑的通量不降低它）。

**修法（门 G6.1）：`smooth_wall_tangential_gradients`（`post/surface.py`），
以 `smooth_passes` 穿进 Cp 曲线与力积分。** 几遍 Jacobi 把每个壁面三角形
的梯度与其边邻居平均，**按外法向对齐门控**
（`n̂_out·n̂_out(nb) > cos_thresh`），所以永不跨尖 TE 平均（两个 TE 三角形
的外法向几乎反平行）。NACA0012 M0.80 粗网格行走解上实测，一遍：G6.1
指标 0.0758 → 0.00023（330×，反转 39 → 1），激波 0.604 → 0.607，cl_p
0.3% 以内，|Cp|_max 不变（TE 保住）——且它是梯度上的固定线性算子
（可微）。`smooth_passes = 0` 是逐位不变的默认；近场 cd_pressure 会动
（粗网格 ~15%），但它是明确不被信任的 FP 量（阻力用 Trefftz）。

**用它修 Cp 曲线，不修载荷（G6.3，2026-07-08 实测）。** ONERA M6 粗网格
上 V6 一致性 |CL_p − CL_KJ|/CL_KJ 随平滑**变差**（`smooth_passes` 0/1/2：
2.40% → 3.35% → 3.88%）：±锯齿在面积分里大部抵消，而平均反而抹平真实的
LE 吸力峰，把 CL_p 再压低 ~1%，离可信的 CL_KJ 更远。所以整个 V6 地板是
尖 TE/LE 的 P1 壁面梯度问题，不是锯齿。（后续链条：尖棱边本就在曲单元
射程之外，§4.1；曲单元路线 P11 于 2026-07-19 以实测负结果关闭，§5.1.3；
M6 尾巴由 P13 重归因到翼尖帽几何，§4.2。）建议：报告的 Cp 曲线用
`smooth_passes > 0`；`wall_force_coefficients` 用 `smooth_passes = 0`
（参数保留为 opt-in，但 raw 载荷更准）。

---

## 10. 验证阶梯

每一项都是 pytest；粗网格上保持 CI 可跑。门证据（PNG/CSV）落到被跟踪的
`cases/gates/<gate_id>/`。各级现行通过状态以能力边界文书为准；V2 的
「<2% @ medium」判据实测超出 P1 场在 h=0.08 的固有能力（§5.1.3），
G1.6 已按 Option C 重定规格。

| # | 算例 | 检查 | 门限 |
|---|---|---|---|
| V0 | 均匀流保持，任意 tet 网格 + 尾迹切口 | ‖R(φ=x)‖∞ | < 1e−12 |
| V1 | Laplace MMS（ν=0, ρ=1），制造三角函数 φ（实现：sin πx · sin πy · sin πz） | L2 误差 vs h | 斜率 ≥ 1.9 |
| V2 | 不可压球 | Cp vs 1 − (9/4)sin²θ | max err < 2%（medium 网格） |
| V3 | 亚临界挤出 NACA0012，M 0.5，α 2°（展向周期或大展弦比 + 对称） | Cp 对 2D FP/面元参考；cl | Δcl < 2% |
| V4 | 跨声速 NACA0012 M 0.80 α 1.25°（挤出） | 激波 x/c 对已发表 FP（Holst AGARD FP workshop 结果——注意 FP ≠ Euler 激波位置） | 激波 Δx/c < 0.03 |
| V5 | ONERA M6，M 0.84，α 3.06° | λ 激波拓扑，η = 0.44/0.65/0.90 剖面 Cp 对 FP 参考（FLO 类 / TRANAIR）；CL | 定性 λ + CL 在 FP 文献散布内 |
| V6 | Γ 一致性 | **2.5D（G2.4）：** 剖面 cl = 2Γ/(Uc) 对表面压强积分。**3D（G5.2）：** 总 CL_p 对*展向积分*的 Kutta–Joukowski CL_KJ = 2∫Γ(z)dz/(U·S)。 | < 1% |

V6 警示（2D vs 3D）：*剖面*恒等式 cl(z) = 2Γ(z)/(U c) 只在 2D /
quasi-2D 流成立。有限 3D 机翼上剖面当地流不是 2D（展向流 + 诱导下洗），
所以压强积分的剖面 cl **不**逐站位等于 2Γ(z)/(U c)——只有*展向积分*的
总升力精确：L = ρ∞U∞∫Γ(z)dz（远场/Trefftz 平面 Kutta–Joukowski）。所以
V6 只对挤出 quasi-2D 算例逐剖面检查（G2.4，`sectional_cl_from_gamma`）；
ONERA M6 上是积分的 `cl_kj_3d`（`post/surface.py`）对总 CL_p。M6 的展向
`2Γ/c` 图是*载荷分布*，不是逐剖面一致性门。

参考数据警示：先对*全势*已发表结果验证（Holst 的 Progress in Aerospace
Sciences 2000 综述收集了它们），Euler 其次——否则模型误差与代码误差
混淆。

---

## 11. 开发路线（指针）

- **阶段历史与归档索引：** `phases/README.md`——phase 1–6 各阶段一句话
  总结、归档内容、保证能跑的历史树恢复点。
- **phase 1 的 track roadmap（冻结历史）：** `phases/p1/docs/roadmap.md`
  （track 索引 + 工作规则），逐 track 台账在 `phases/p1/docs/roadmap/`：
  Track P（求解器 P0–P14）、Track M（网格）、Track B（level-set 尾迹，
  phase 3 已删）、Track V（粘性耦合）、Track A（校验）。阶段 ↔ §10 验证
  阶梯的映射在各 track 文件的门条目里。
- **现行能力声明（实测边界）：**
  `phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md`
  （2026-10-11）。
- 本文前身（英文版，维护至 phase 6 收口）归档：`phases/p6/docs/design.md`。

---

## 12. 已知风险与缓解

1. **尾迹贴体网格是一号实践风险。** 其余都是教科书内容；让 Gmsh 嵌入尾迹
   面、让预处理器鲁棒劈节点（TE 缝、翼尖边、尾迹–远场相交）会吃掉真实
   工作量。缓解：先对一个平凡简单的挤出翼型网格建 `wake_cut.py`；写拓扑
   断言（每个尾迹面恰有一个 ⁺ 单元和一个 ⁻ 单元）。⚠ **勘误（P2，
   2026-07-06；文档 2026-07-19 改）：TE 节点*是*复制的**——原稿「TE
   节点不复制」的规格在 P2 期间带发散证据重定规格（单值 TE 定量上是错
   的；Kutta 需要 b≈0.93 的割线）。见 §4 与 P2 阶段条目。
   **守它的门**:`tests/A/test_A28_wingbody_topology.py`（翼身拓扑/结合部）、`tests/B/test_B04_vtk_io.py`（wake 切割 roundtrip）、`tests/F/test_F01_provenance.py`（网格 manifest 指纹）。〔2026-10-12,W3.2〕
2. **FP 非唯一**（M∞ ≈ 0.82–0.85，低 α，常规翼型）：延拓方向会选出不同
   分支。记录在案；永远从亚临界收敛态向上爬 Mach。
   **守它的门**:`tests/D/test_D14_alpha0_transonic_nonuniqueness.py`（α=0 非唯一带的能力边界门）；运行时缓解：`newton.py` 的 `_FP_NONUNIQUE_MACH` 冷启动提示（W2.4）。〔2026-10-12,W3.2〕
3. **强激波**（M_local > 1.3）：等熵跳跃低估激波压升，激波位置比 Euler
   靠后。不是 bug；在 README 写明有效包线，并在 max(M) 超过 1.35 时
   运行时告警。
   **守它的门**:**无门（登记）** —— 有运行时缓解（`newton.py` 的 `_FP_SHOCK_VALIDITY_MACH` 的 print 警告,W2.4）,没有任何测试断言它触发;登记：可加一条警告触发的机制测试。〔2026-10-12,W3.2〕
4. **带激波算例上 Picard 停滞**：症状——残差平台，激波振荡一个单元。
   缓解按序：(3.2) 中加大 C、降低 ω、延拓，然后 Newton+PTC。
   **守它的门**:`tests/A/test_A30_seed_fallback_contract.py`（冷启动种子回退合约）、`bench/failure_modes.py::classify_failure`（五模式分类,经 D11/D12 使用）。〔2026-10-12,W3.2〕
5. **vibe coding 中的 Numba 编译期摩擦**：内核保持小、参数列表显式；
   `cache=True`；`PYFP3D_NOJIT=1` 环境开关把 `njit` 换成恒等装饰器，
   使 pdb 调试可行。
   **守它的门**:`tests/A/test_A01_freestream_preservation.py::test_import_pyfp3d`（守 NOJIT 默认值语义）；机制通道本身 = `PYFP3D_NOJIT` 各模块级开关 + 暖态纪律（CLAUDE.md 全套基线段）。〔2026-10-12,W3.2〕

---

## 13. 核心参考文献

全势跨声速与人工密度：

- Holst, T.L., "Transonic flow computations using nonlinear potential
  methods," *Progress in Aerospace Sciences* 36 (2000) —— 标准综述；含
  人工密度一族与 workshop 验证数据。
- Hafez, M., South, J., Murman, E., "Artificial compressibility methods for
  numerical solutions of transonic full potential equation," *AIAA J.* 17(8),
  1979.
- Jameson, A., "Iterative solution of transonic flows over airfoils and
  wings," *CPAM* 27 (1974)；机翼 FP 基线见 FLO 系列报告。
- Steinhoff, J., Jameson, A., "Multiple solutions of the transonic potential
  flow equation," *AIAA J.* 20(11), 1982 —— 非唯一性。
- Neel, R.E., *Advances in the Computation of Transonic Full Potential Flows
  on Unstructured Grids*, PhD thesis, Virginia Tech, 1997 —— 与本求解器
  概念最接近的先驱（顶点型非结构 FP）。
- Caughey, D.A., Jameson, A., "Basic advances in the finite-volume method for
  transonic potential-flow calculations" —— FV/FEM 等价视角。
- Drela, M., *Flight Vehicle Aerodynamics*, MIT Press 2014 —— 势流方法中
  的尾迹/Kutta 与 Trefftz 平面处理。

平直面元网格上曲壁的边界通量修正（§5.1）：

- Krivodonova, L., Berger, M., "High-order accurate implementation of solid
  wall boundary conditions in curved geometries," *J. Comput. Phys.* 211
  (2006) 492–512.
- Ciallella, M., Gaburro, E., Lorini, M., Ricchiuto, M., "Shifted boundary
  polynomial corrections for compressible flows: high order on curved domains
  using linear meshes," *Appl. Math. Comput.* (2023).
- Main, A., Scovazzi, G., "The shifted boundary method for embedded domain
  computations. Part I: Poisson and Stokes problems," *J. Comput. Phys.* 372
  (2018) 972–995.
- Atallah, N.M., Canuto, C., Scovazzi, G., "The high-order Shifted Boundary
  Method and its analysis," *CMAME* 394 (2022) 114885.
- Collins, J.H., Li, C., Lozinski, A., Scovazzi, G., "Gap-SBM: A New
  Conceptualization of the Shifted Boundary Method with Optimal Convergence
  for the Neumann and Dirichlet Problems," arXiv:2508.09613 (2025).
- Bramble, J.H., Dupont, T., Thomée, V., "Projection methods for Dirichlet's
  problem in approximating polygonal domains with boundary-value corrections,"
  *Math. Comp.* 26 (1972) 869–879.
- Burman, E., Hansbo, P., Larson, M.G., "A cut finite element method with
  boundary value correction," *Math. Comp.* 87 (2018) 633–657.

FE 全势、嵌入尾迹与全耦合 Newton（§3.1–3.2、§6.3、§8.1）：

- López Canalejo, I., *A Finite-Element Transonic Potential Flow Solver with an
  Embedded Wake Approach for Aircraft Conceptual Design*, PhD dissertation,
  Technische Universität München, 2021 —— **未纳入仓库**（旧版写的
  `docs/references/Dissertation_Inigo_Lopez.pdf` 路径不存在；仓库根的
  `reference/` 只有两份 Drela PDF 与 XFOIL 6.99 源码包，见下）。FE 人工
  密度通量与全 Jacobian Newton 求解器最接近的先驱。§3.3 人工压缩性
  （Eq. 3.19–3.27，开关函数 μ_s = μ_c·max(0,μ,μ_up)）；§3.4 限速 /
  M²_max ~ O(3) 钳制；§3.5 嵌入尾迹（2D/3D 条件，TE 与小切口处理）；
  Ch.4 Newton 收敛（Tables 4.5/4.6/4.9，严格二次）与载荷步进调度
  （Tables 4.7/4.8 NACA0012，4.13 ONERA M6）；ONERA M6 载荷 Table 4.15
  （CL 0.288 KRATOS = Tranair）；Appendix B Eq. (B.1)–(B.17) 完整敏感度
  推导，含非零上游耦合 (B.4/B.6) 与开关函数导数 (B.3/B.6)。NACA0012
  验证用**尖** TE（Eq. 4.2 修改最后一个厚度系数以闭合原本钝的翼型），
  Newton 全二次收敛——即 2D 尖 TE 本身不是 Newton 障碍；只有 3D 翼尖-TE
  几何奇点需要密度钳制。
- López Canalejo, I., Núñez, M., Baiges, J., Rossi, R., "An embedded approach
  for the solution of the full potential equation with finite elements,"
  *CMAME* 388 (2022) 114244 —— 学位论文 Ch.3–4 的期刊浓缩版。
- Cai, X.-C., Keyes, D.E., Young, D.P., "Parallel Newton–Krylov–Schwarz
  algorithms for the transonic full potential equation," *SIAM J. Sci. Comput.*
  19 (1998) 246–265 —— FPE 的非精确 Newton + GMRES。
- Eisenstat, S.C., Walker, H.F., "Choosing the forcing terms in an inexact
  Newton method," *SIAM J. Sci. Comput.* 17 (1996) 16–32.

仓库内文献（`reference/`）：

- `reference/drela-2013-three-dimensional-integral-boundary-layer-formulation-for-general-configurations.pdf`
  —— Drela IBL3 原始文献，Track V 粘性耦合（§5）的方程来源。
- `reference/drela-giles-2012-viscous-inviscid-analysis-of-transonic-and-low-reynolds-number-airfoils.pdf`
  —— 跨声速/低 Re 翼型粘性–无粘分析（XFOIL 体系）。
- `reference/XFOIL6.99.zip` —— XFOIL 6.99 源码包（标定参照）。
