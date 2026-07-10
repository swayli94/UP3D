# 尾流表示方案设计：从 Conforming 面到 Level-Set 多值有限元
> Design note for UP3D / pyFP3D
> Date: 2026-07-07 (updated 2026-07-07 23:22 — baseline 同步 e3d0386)
> Code baseline: HEAD e3d0386
> References:
> - Davari et al. 2019, *Comput. Mech.* 63:821–833（body-fitted 几何 + level-set 尾流，UP3D 方案 B 的直接前作）
> - Núñez et al. 2022, *CMAME* 388:114244（全嵌入式几何 + level-set 尾流，方案 B 的理论泛化参考）
---
## 1. 问题背景
### 1.1 当前架构的局限
UP3D 当前的尾流处理采用 **conforming 网格尾流面 + 物理节点复制 + master-slave 消元** 方案：
```
Gmsh 生成网格（尾流面作为 embedded surface）
   → mesh/wake_cut.py 物理节点复制
   → constraints/wake.py T^T A T 消元 + Γ RHS
   → solve/picard.py 外层 Γ secant 迭代
   → solve/linear.py AMG-CG
```
这个方案在 P1–P4 阶段工作正常（P4 G4.1 medium 已于 2026-07-07 关闭），但有三个结构性限制：
1. **网格生成必须预先嵌入尾流面**——Gmsh 需要 `model.mesh.embed` 尾流线/面，网格拓扑与尾流几何绑定
2. **改攻角需要重新生成网格**——尾流面方向跟着来流方向变
3. **钝尾缘没有明确的尖缘来锚定尾流面**——尾流从钝体后缘的哪个位置出发？上下表面分离点在哪？
此外，当前方案在多物体场景（翼身组合体、多段翼）下，网格生成的拓扑复杂度急剧上升——每条尾流面需要单独 embed，尾流面之间的交叉关系需要在网格中处理。
### 1.2 设计需求
- 网格生成不再依赖尾流面
- 改攻角不需要重新生成网格
- 支持钝尾缘绕流
- 天然支持多物体（多段翼、翼身组合体）
- 为未来曲线尾流/自由尾迹打开可能
- 保持与 P1–P4 验证结果的一致性
---
## 2. 三个候选方案
### 2.1 方案 A：罚函数 Kutta + 保留网格尾流面（最小改动）
**思路**：保留当前的 conforming 网格尾流面和 `wake_cut.py` 节点复制，仅把 Kutta 条件从显式 $\Gamma$ secant 改为罚函数。
```
Gmsh 生成网格（尾流面作为 embedded surface）        ← 不变
   → mesh/wake_cut.py 物理节点复制                    ← 不变
   → constraints/wake.py T^T A T 消元                 ← 不变
   → wake/kutta_penalty.py 罚函数 Kutta               ← 新增（替代 Γ secant）
   → solve/picard.py 无 Γ 外层循环                    ← 修改
```
**优点**：
- 改动量最小，主要在 `solve/picard.py` 和新增 `wake/kutta_penalty.py`
- 可以快速验证罚函数 Kutta 是否比显式 $\Gamma$ secant 更好
- 去掉 $\Gamma$ 外层循环后，P4 transonic 的三层嵌套（Mach × $\Gamma$ × Picard）变成两层，可能减少 30-50% 总迭代数
**缺点**：
- **不解决核心问题**：网格生成仍然需要嵌入尾流面
- 不支持钝尾缘
- 不支持多物体的尾流面拓扑简化
- 改攻角仍需重新生成网格
**定位**：快速验证罚函数 Kutta 的数值行为，为方案 B 铺路。可以作为独立的小步改进先落地。
---
### 2.2 方案 C：网格不嵌入尾流 + 求解器内节点复制（中等改动）
**思路**：网格生成时不 embed 尾流面，求解器读入网格后用 level-set 函数标记尾流面位置，然后在求解器内部做节点复制（和当前 `wake_cut.py` 类似，但基于 level-set 而非网格 tag）。Kutta 仍用显式 $\Gamma$ secant。
```
Gmsh 生成网格（无尾流面，只需 wall + farfield）      ← 简化
   → wake/levelset.py 定义 φ_wake(x) = 0             ← 新增
   → wake/cut_and_duplicate.py 基于 level-set 做节点复制  ← 新增（替代 wake_cut.py）
   → constraints/wake.py T^T A T 消元 + Γ RHS         ← 基本不变
   → solve/picard.py Γ secant 外层循环                ← 不变
```
**优点**：
- 网格生成不需要嵌入尾流面
- 改攻角只需更新 level-set 函数，不重新生成网格
- 保留了现有的 Kutta 机制（$\Gamma$ secant）和 `constraints/wake.py` 的消元逻辑
- 实现量中等
**缺点**：
- **仍然是物理节点复制**——尾流面所有节点都被复制，进全局编号，矩阵规模 $\approx 2N$
- 尾流面变化时需要重新做节点复制（虽然不重新生成网格，但拓扑结构变了）
- **不支持曲线尾流**——物理节点复制假设尾流面与网格边对齐，弯曲的尾流面无法处理
- Kutta 仍用 $\Gamma$ secant 外层循环，P4 的三层嵌套问题不解决
- 多尾流场景下节点复制的拓扑管理仍然复杂
**定位**：比方案 A 前进一步，解耦了网格生成和尾流面，但保留了节点复制的底层机制。是一个可行的中间方案，但天花板有限。
---
### 2.3 方案 B：Level-Set 尾流 + 多值有限元（目标方案）
**思路**：用 level-set 隐式定义尾流面，被切单元用多值有限元（局部富集）处理势函数跳跃，Kutta 用罚函数。彻底消除物理节点复制和 $\Gamma$ 外层循环。
```
Gmsh 生成网格（无尾流面，只需 wall + farfield）      ← 简化
   → wake/levelset.py 定义 φ_wake(x) = 0             ← 新增
   → wake/cut_elements.py 识别被切单元 + 从 DOF        ← 新增
   → kernels/cut_assembly.py 多值 FE 装配             ← 新增
   → wake/kutta.py 罚函数 Kutta                       ← 新增
   → solve/picard.py 无 Γ 外层循环                    ← 修改
```
**优点**：
- 网格生成不需要嵌入尾流面
- 改攻角只需更新 level-set 函数
- **局部富集**——只有被切单元增加从 DOF，矩阵规模 $\approx N + N_{\mathrm{ext}}$（$N_{\mathrm{ext}} \ll N$），远小于方案 C 的 $2N$
- 罚函数 Kutta 去掉 $\Gamma$ 外层循环，简化求解结构
- **支持曲线尾流**——level-set 可以表示任意曲面，为自由尾迹打开门
- **天然支持多物体**——多个 level-set 函数各自独立
- 钝尾缘可通过罚函数 $\mathbf{n}_{\mathrm{kutta}}$ 的定义来处理
**缺点**：
- 实现量最大——需要重写尾流模块，新增切单元装配逻辑
- 多值 FE 是新的数值基础设施，需要仔细验证
- 罚函数 $k_{\mathrm{kutta}}$ 参数需要调试
- AMG 在扩展矩阵上的性能需要验证
**定位**：终极方案。一次性投入大，但解决了所有结构性限制，为后续扩展打开空间。
---
### 2.4 三个方案的对比
| 维度                   | 方案 A（罚函数 Kutta） | 方案 C（Level-set + 节点复制） | 方案 B（Level-set + 多值 FE） |
| ---------------------- | ---------------------- | ------------------------------ | ----------------------------- |
| 网格无需嵌入尾流面     | ❌                      | ✅                              | ✅                             |
| 攻角变化不重新生成网格 | ❌                      | ✅                              | ✅                             |
| 钝尾缘支持             | ❌                      | ❌                              | ✅                             |
| 多物体支持             | ❌                      | 部分（拓扑复杂）               | ✅（天然）                     |
| 曲线尾流/自由尾迹      | ❌                      | ❌                              | ✅（可扩展）                   |
| 矩阵规模               | $N$（消元后）          | $2N$                           | $N + N_{\mathrm{ext}}$        |
| $\Gamma$ 外层循环      | 去掉                   | 保留                           | 去掉                          |
| P4 三层嵌套简化        | ✅（→两层）             | ❌                              | ✅（→两层）                    |
| 实现量                 | 小                     | 中                             | 大                            |
| 风险                   | 低                     | 中                             | 中高                          |
| 天花板                 | 低                     | 中                             | 高                            |
> **注（2026-07-07 22:45 更新）**：P4 的三层嵌套（Mach × $\Gamma$ secant × Picard）在当前代码中仍然存在。P4 修复（local $\theta \cdot \mathrm{diag}(A_{\mathrm{free}})$ 阻尼 + $h_j$ 批处理化）解决了 medium mesh 发散问题，但未改变三层嵌套结构。方案 B 去掉 $\Gamma$ 外层循环的收益依然成立。
### 2.5 为什么选择方案 B
方案 A 是"快速止血"——它改善了 P4 的求解结构，但不解决网格依赖这个核心问题。
方案 C 是"半步前进"——它解耦了网格生成和尾流面，但底层仍然是物理节点复制。这意味着：
- 矩阵规模翻倍（$2N$），在线性求解器层面比方案 B 更贵
- 不支持曲线尾流——物理节点复制要求尾流面与网格边对齐
- 多物体场景下节点复制的拓扑管理仍然复杂
- $\Gamma$ secant 外层循环保留，P4 的迭代效率问题不解决
方案 B 是"根治"——它用 level-set + 多值 FE 彻底替代了物理节点复制，带来三个本质性优势：
1. **矩阵规模最优**：只有被切单元增加从 DOF（$N_{\mathrm{ext}} \ll N$），不像方案 C 的全尾流面节点复制
2. **表示能力最强**：level-set 可以表示任意曲面（包括弯曲尾流），物理节点复制做不到
3. **求解结构最简**：罚函数 Kutta 去掉 $\Gamma$ 外层循环，P4 从三层嵌套变两层
方案 B 的实现量虽然最大，但其核心数值工作（多值 FE 装配）有一个关键简化——**切单元不需要几何切割积分**，只需用不同 DOF 编号装配两次（详见 §3.3）。这个简化使得方案 B 的实现量比初看可控得多。
**推荐路径**：先落地方案 A（验证罚函数 Kutta），再推进方案 B（彻底替代）。方案 C 不单独实施——如果方案 B 的某个阶段需要中间验证，方案 C 的部分组件可以临时使用。
---
## 3. 方案 B 的理论基础
以下理论内容引用两篇核心文献：
- **Davari 2019**：body-fitted 几何 + level-set 尾流 + Cut FE（UP3D 方案 B 的直接前作）
- **Núñez 2022**：全嵌入式几何 + level-set 尾流 + 多值 FE（理论泛化参考）
两篇论文的尾流面处理技术一脉相承，区别在于几何是否也用 level-set。UP3D 采用 body-fitted 几何，对应 Davari 2019 的架构。
### 3.1 尾流面的 Level-Set 表示
尾流面 $\partial \Omega_W$ 用一个符号距离函数 $\varphi_{\mathrm{wake}}(\mathbf{x})$ 隐式定义：
$$\varphi_{\mathrm{wake}}(\mathbf{x}) = 0 \iff \mathbf{x} \in \partial \Omega_W \quad \text{（尾流面）}$$
$$\varphi_{\mathrm{wake}}(\mathbf{x}) > 0 \iff \mathbf{x} \in \Omega_W^+ \quad \text{（上翼面侧）}$$
$$\varphi_{\mathrm{wake}}(\mathbf{x}) < 0 \iff \mathbf{x} \in \Omega_W^- \quad \text{（下翼面侧）}$$
Davari 2019（§2.2）首次提出在 body-fitted 几何上用 level-set 定义尾流面：
> *The wake is embedded in the CFD mesh by defining a level set function $\varphi_{\mathrm{wake}}$ that identifies the wake surface.*
> Núñez 2022（§2.4）沿用了相同的尾流 level-set 定义，并进一步将几何也用 level-set 表示。Núñez 2022 §2.3 明确区分了两者：
> *While a level set function is used in [8] to define the wake, the airfoil shapes considered are still body-fitted. In this paper, both the geometry of study and the wake are modelled using embedded methods.*
> 其中 [8] 即 Davari 2019。
> 对于直尾流（当前需求），$\varphi_{\mathrm{wake}}$ 是从 TE 点出发、沿来流方向的半平面的符号距离函数。对于曲线尾流（未来扩展），$\varphi_{\mathrm{wake}}$ 可以是任意曲面的符号距离函数——**表示能力不受限制**。
> **多尾流**：每个尾流面一个独立的 level-set 函数。被多个尾流面同时切到的单元可以有多个从 DOF。
### 3.2 多值有限元空间分解
被尾流面切到的单元中，势函数 $\Phi$ 是多值的——尾流面上侧用 $\Phi^+$，下侧用 $\Phi^-$。
**Davari 2019 的做法**（Eq.28–32）：通过 Heaviside 函数 $H(\varphi_{\mathrm{wake}})$ 构造扩展形函数。被切单元中，每个节点的势函数表示为：
$$\Phi_h(\mathbf{x}) = \sum_i \left[ \Phi_i \cdot N_i(\mathbf{x}) + \Psi_i \cdot N_i(\mathbf{x}) \cdot H(\varphi_{\mathrm{wake}}(\mathbf{x})) \right]$$
其中 $\Phi_i$ 是主 DOF，$\Psi_i$ 是辅助 DOF。$H(\varphi_{\mathrm{wake}})$ 在 $+$ 侧为 $1$，$-$ 侧为 $0$。因此：
- $+$ 侧：$\Phi_h = \sum_i (\Phi_i + \Psi_i) N_i$
- $-$ 侧：$\Phi_h = \sum_i \Phi_i N_i$
跳跃量 $[\Phi] = \sum_i \Psi_i N_i$。
**Núñez 2022 的做法**（§2.4）：用 FE 空间分解来表述同一思想：
$$V_h = V_h^{\mathrm{main}} \oplus V_h^{\mathrm{ext}}$$
- **$V_h^{\mathrm{main}}$**：背景网格所有标准节点的 FE 空间。与当前 UP3D 的 P1 空间完全相同。
- **$V_h^{\mathrm{ext}}$**：仅在被切单元中存在的从节点（$\Phi_{\mathrm{ext}}$）的 FE 空间。**支撑局限于被切单元**，出了被切单元就是零。
两种表述等价——Davari 的 Heaviside 扩展形函数和 Núñez 的空间分解是同一数学结构的两种实现方式。UP3D 采用 Núñez 的空间分解表述，因为它与现有的 P1 Galerkin 装配框架更自然地对接。
被切单元中，节点的 DOF 分配遵循侧标记（Núñez 2022 Eq.(11)-(12)）：
$$\Phi^+_i = \begin{cases} \Phi_i & \text{if } \varphi_{\mathrm{wake}} > 0 \text{ （+侧节点用主 DOF）} \\ \Phi_{\mathrm{ext},i} & \text{if } \varphi_{\mathrm{wake}} < 0 \text{ （-侧节点用从 DOF）} \end{cases}$$
$$\Phi^-_i = \begin{cases} \Phi_{\mathrm{ext},i} & \text{if } \varphi_{\mathrm{wake}} > 0 \text{ （+侧节点用从 DOF）} \\ \Phi_i & \text{if } \varphi_{\mathrm{wake}} < 0 \text{ （-侧节点用主 DOF）} \end{cases}$$
即："$+$"侧的积分用一套 DOF 编号，"$-$"侧用另一套。同一个形状函数 $N_i$ 在两侧有不同的 DOF 值，从而实现了势函数的跳跃。
这是标准的 XFEM/GFEM 局部富集思路——从 DOF 不进入全局编号的常规部分，只在被切单元的局部装配中出现。矩阵规模 $\approx N_{\mathrm{global}} + N_{\mathrm{ext}}$，远小于物理节点复制的 $2 \times N_{\mathrm{global}}$。
### 3.3 切单元装配——不需要几何切割积分
这是 Davari/Núñez 方案的关键工程简化。
对被切单元 $e$，其残差贡献分为"$+$"侧和"$-$"侧两部分（Núñez 2022 Eq.(28)-(29)；Davari 2019 Eq.(33)–(34) 中等价的 Heaviside 形式）：
$$R^+_{i,W+}(\Phi^+) = \int_{\Omega_W^+} \rho \nabla N_i \cdot \nabla N_j \, d\Omega \cdot \Phi^+_j - \int_{\partial \Omega_W} N_i g \, d\partial\Omega$$
$$R^-_{i,W-}(\Phi^-) = \int_{\Omega_W^-} \rho \nabla N_i \cdot \nabla N_j \, d\Omega \cdot \Phi^-_j - \int_{\partial \Omega_W} N_i g \, d\partial\Omega$$
Núñez 2022（§2.4）：
> *The main potential flow terms are evaluated in the wake elements on both the upper and lower side of the domain... The residuals and Jacobians derived from these equations are rewritten for the elements that are intersected by the wake using its corresponding values of potential according to Eqs. (11) and (12).*
> **关键简化**：在 P1 四面体上，形状函数梯度 $B_e$ 是常数，体积 $V_e$ 也是常数。切单元被尾流面分成两个子区域，但两个子区域的积分可以近似为**整个单元积分的一个分配**——不需要做精确的几何切割子区域积分。
> 具体来说，切单元的矩阵贡献是：
> $$A^e = \tilde{\rho}_e \cdot V_e \cdot (B_e^+ \cdot B_e^{+T} + B_e^- \cdot B_e^{-T})$$
> 其中 $B_e^+$ 和 $B_e^-$ 是**相同的** $B_e$（同一个四面体的形状函数梯度），但 DOF 编号不同：
- $B_e^+$："$+$"侧节点用 ext DOF，"$-$"侧节点用 main DOF
- $B_e^-$：所有节点用 main DOF（或反之，取决于约定）
所以 $A^e$ 的结构是：**标准 P1 刚度矩阵 + 一份替换了 DOF 编号的副本**。不需要切体积积分，不需要子四面体分解。
**代价**：切单元被"装配两次"（一次给每侧），但切单元数 $\ll$ 总单元数，开销可忽略。
切单元被尾流面分成两个子体积 $V_e^+$ 和 $V_e^-$。在 P1 四面体上 $B_e$ 是常数，最简单的近似是各取 $V_e/2$。更精确的做法是计算 level-set 面将四面体切成的两个子体积比——但这只是个几何计算，不影响装配逻辑。Davari 2019 和 Núñez 2022 都采用了类似的体积分数近似。
### 3.4 尾流面边界条件
尾流面上需要施加两个条件（Davari 2019 Eq.(6)–(7)；Núñez 2022 Eq.(8)–(9)）：
**1. 质量通量连续**：
$$\mathbf{n} \cdot (\rho \nabla \Phi^+ - \rho \nabla \Phi^-) = 0 \quad \text{on } \partial \Omega_W$$
**2. 压力相等**：
$$|\nabla \Phi^+|^2 - |\nabla \Phi^-|^2 = 0 \quad \text{on } \partial \Omega_W$$
**关键实现细节**：这两个条件不是直接强施加的，而是通过**最小二乘泛函**弱施加。两篇论文的做法有差异：
**Davari 2019 的做法**（Eq.23–27）：将压力相等条件线性化（小扰动假设 $u \ll v_\infty$），得到对流方程：
$$\mathbf{u}_\infty \cdot \nabla(\Phi^+ - \Phi^-) = 0 \quad \text{on } \partial \Omega_W \quad \text{(Eq.23)}$$
然后构造最小二乘泛函：
$$\Pi(\Phi) = \frac{1}{2} \int_{\Gamma_W} \left[\mathbf{u}_\infty \cdot \nabla(\Phi^+ - \Phi^-)\right]^2 d\Gamma \quad \text{(Eq.24)}$$
对应的残差和 Jacobian 通过对 $\Phi$ 求导得到（Eq.25–27）。这个做法的优点是 Jacobian 结构简单，但依赖小扰动线性化。
**Núñez 2022 的做法**（Eq.30–31）：将两个条件合成一个矢量方程：
$$\rho^+ \mathbf{v}^+ = \rho^- \mathbf{v}^- \quad \text{on } \partial \Omega_W \quad \text{(Eq.30)}$$
即 $\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^- = 0$（同时包含法向质量通量连续和切向压力相等）。然后构造最小二乘泛函：
$$\Pi(\Phi^+, \Phi^-) = \frac{1}{2} \int_{\partial \Omega_W} \|\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^-\|^2 \, d\partial\Omega \quad \text{(Eq.31)}$$
对应的残差和 Jacobian 从泛函求导得到。关键性质：**约束项的 Jacobian 等于主问题 Jacobian 的差** (Eq.38: $J_B = J_W^+ - J_W^-$)，这是因为形函数在切单元的子区域上是分片线性的。
Núñez 2022 的矢量形式不依赖小扰动线性化，更适用于跨声速。UP3D 采用 Núñez 的矢量形式。
约束项的残差贡献加到切单元的对应 DOF 上：
$$R_B^{e,+}(\Phi^+, \Phi^-) = R_{W}^{e,+}(\Phi^+) - R_{W}^{e,-}(\Phi^-) \quad \text{(Eq.36)}$$
$$R_B^{e,-}(\Phi^+, \Phi^-) = R_{W}^{e,-}(\Phi^-) - R_{W}^{e,+}(\Phi^+) \quad \text{(Eq.37)}$$
**压力相等在远场的自动满足**：对于直尾流，远场处两侧速度趋近来流，压力自动相等。Kutta 条件在 TE 处通过罚函数施加（§3.5）。
**注意**：在多值 FE 中，质量通量连续条件**不会**"自然满足"——它需要通过上述最小二乘泛函显式施加。这是因为切单元的两侧各自做 Galerkin 积分时，内部分界面上的法向通量项只在各自侧的积分中出现，不会自动抵消。最小二乘泛函的作用就是强制两侧通量一致。
### 3.5 罚函数 Kutta 条件
Kutta 条件的本质是：流动平滑离开尾缘，速度沿尾缘角平分线方向（Núñez 2022 §2.6）。
Núñez 2022：
> *the Kutta condition is defined as the fact that the velocity leaves smoothly the trailing edge. This condition is needed to find a unique solution in the potential flow equation, which gives the right amount of circulation in the solution.*
> 数学表达（Núñez 2022 Eq.(44)）：
> $$\mathbf{n}_{\mathrm{kutta}} \cdot (\nabla \Phi)\big|_{\mathrm{TE}} = 0$$
> 其中 $\mathbf{n}_{\mathrm{kutta}}$ 是尾缘角平分线的法向。
> Núñez 2022 通过罚函数施加（Eq.(45)–(46)）：
> $$\Lambda(\Phi, k_{\mathrm{kutta}}) = \frac{1}{2} \rho \|\nabla \Phi\|^2 - \frac{k_{\mathrm{kutta}}}{2} \|\psi(\Phi)\|^2 \quad \text{(Eq.45)}$$
> 其中 $\psi(\Phi) = \mathbf{n}_{\mathrm{kutta}} \cdot (\nabla \Phi)_{\mathrm{TE}}$。
> 离散化后的罚项矩阵贡献：
> $$K_{\mathrm{kutta}}[i,j] = k_{\mathrm{kutta}} \cdot \int_{\Omega_{\mathrm{kutta}}} (\mathbf{n}_{\mathrm{kutta}} \cdot \nabla N_i)(\mathbf{n}_{\mathrm{kutta}} \cdot \nabla N_j) \, d\Omega$$
> 在 P1 四面体上这是一个秩-1 矩阵：
> $$K_{\mathrm{kutta}}^e = k_{\mathrm{kutta}} \cdot V_e \cdot (\mathbf{n}_{\mathrm{kutta}} \cdot B_e)^T \cdot (\mathbf{n}_{\mathrm{kutta}} \cdot B_e)$$
> 加到 $\Omega_{\mathrm{kutta}}$ 单元对应的 DOF 上。
> **$\Omega_{\mathrm{kutta}}$ 的识别——body-fitted vs 嵌入式几何的差异**：
- **Núñez 2022（全嵌入式几何）**：需要 Algorithm 1 来识别 $\Omega_{\mathrm{kutta}}$——找同时被几何 level-set $\varphi$ 和尾流 level-set $\varphi_{\mathrm{wake}}$ 切到的单元及其邻居。这是因为后缘是 level-set 交叉点，不对应明确的网格节点。Núñez 2022 还指出，在 $\Omega_{\mathrm{kutta}}$ 区域内，几何 level-set 的切单元分割需要特殊处理——因为后缘是尖缘，level-set 在此处的描述不够准确，无穿透条件 可能给出错误的速度方向。因此这些单元忽略几何分割，只保留尾流分割和 Kutta 罚函数。
- **Davari 2019 / UP3D（body-fitted 几何）**：后缘节点明确存在，尾流从后缘出发。**不需要 Algorithm 1**——$\Omega_{\mathrm{kutta}}$ 就是 TE 节点所在的单元及其邻居。Kutta 条件直接在 TE 附近的单元上施加罚函数。
这是 body-fitted 几何的一个显著优势：Kutta 条件的实现更简单直接。
**钝尾缘的处理**：$\mathbf{n}_{\mathrm{kutta}}$ 不要求尖缘——钝尾缘的上下表面在后缘处有切向，取其角平分线法向即可。Núñez 2022：
> *In an embedded setting, one cannot rely on the no-penetration condition coming from the adjacent elements, as the shape of the trailing edge will not be properly defined in general. ... the equations above are used in conjunction with a penalty term, which ensures that the Kutta condition is satisfied even if the trailing edge is not perfectly described by $\varphi$.*
> 注意：这段论述针对的是全嵌入式几何。在 body-fitted 几何（UP3D）中，后缘形状由网格精确描述，钝尾缘的 $\mathbf{n}_{\mathrm{kutta}}$ 可以直接从网格几何计算。
### 3.6 $\Gamma$ 的涌现
在方案 B 中，**$\Gamma$ 不再是显式参数**——没有外层 secant 迭代。环量从解中自然涌现：
$$\Gamma = [\Phi] = \Phi^+_{\mathrm{TE}} - \Phi^-_{\mathrm{TE}}$$
这是因为多值 FE 允许 $\Phi$ 在尾流面两侧取不同值，罚函数 Kutta 约束了 TE 处的速度方向，整个系统自洽地确定 $\Gamma$。
对比当前方案：$\Gamma$ 是显式参数，需要外层 secant 迭代来找正确的 $\Gamma$ 值。在 P4 transonic 中，这个外层循环（每 Mach level 最多 12 次 gamma eval × 800 Picard its）是迭代量爆炸的主要原因之一。
---
## 4. 方案 B 的模块设计
### 4.1 新增模块
```
pyfp3d/wake/
├── __init__.py
├── levelset.py          # φ_wake(x) 定义与求值
├── cut_elements.py      # 被切单元识别 + 从 DOF 编号 + 侧标记
├── kutta.py             # Ω_kutta 识别 + 罚函数项
└── multivalued.py       # 多值 FE 空间管理（DOF 映射、CSR 模式）
pyfp3d/kernels/
├── cut_assembly.py      # 切单元的双侧残差 + Jacobian 装配（NEW）
```
### 4.2 修改模块
| 模块                       | 改动                                               |
| -------------------------- | -------------------------------------------------- |
| `mesh/reader.py`           | 不再需要 "wake" boundary tag                       |
| `mesh/wake_cut.py`         | **废弃**，被 `wake/` 替代                          |
| `constraints/wake.py`      | **废弃**，master-slave 消元不再需要                |
| `constraints/dirichlet.py` | 远场涡修正改用 level-set 尾流方向                  |
| `kernels/residual.py`      | 增加 `assemble_residual_cut()` 路径                |
| `kernels/jacobian.py`      | `PicardOperator` 扩展支持从 DOF                    |
| `solve/picard.py`          | 去掉 $\Gamma$ secant 外层循环                      |
| `solve/continuation.py`    | 去掉 $\Gamma$ 外层循环，简化为纯 Mach continuation |
| `solve/linear.py`          | AMG 预处理器需要处理扩展矩阵（主 DOF + 从 DOF）    |
| `post/vtk_out.py`          | 可视化需要处理双侧 $\varphi$                       |
| `post/surface.py`          | 壁面 Cp 提取不受影响（壁面单元不被切）             |
### 4.3 不变模块
| 模块                    | 原因                                      |
| ----------------------- | ----------------------------------------- |
| `physics/isentropic.py` | 物理关系与尾流表示无关                    |
| `kernels/upwind.py`     | 人工密度在单元层面操作，与 DOF 结构无关   |
| `mesh/metrics.py`       | $B_e$, $V_e$ 计算不变（切单元用相同几何） |
| `mesh/coloring.py`      | 着色基于单元邻接图，与 DOF 无关           |
### 4.4 核心数据结构
#### WakeLevelSet（`wake/levelset.py`）
```python
class WakeLevelSet:
    """尾流面的 level-set 表示：φ_wake(x) = 0 定义尾流面。
    直尾流：从 TE 点出发、沿来流方向的半平面的符号距离函数。
    多尾流：多个 WakeLevelSet 实例，各自独立。
    曲线尾流（未来）：更新 φ_wake 的几何以对齐流场。
    """
    def __init__(
        self,
        te_point: np.ndarray,       # 尾缘点 (3,)
        direction: np.ndarray,      # 尾流面方向 (3,)（来流方向）
        normal: np.ndarray,         # 尾流面法向 (3,)（⊥ direction, 翼面法向）
        extent: float = 50.0,       # 尾流延伸长度（到远场）
    ):
        ...
    def evaluate(self, points: np.ndarray) -> np.ndarray:
        """在节点上求值 φ_wake，返回 符号距离。
        φ_wake > 0: "+" 侧（上翼面侧）
        φ_wake < 0: "−" 侧（下翼面侧）
        φ_wake = 0: 在尾流面上
        """
        ...
    def update_direction(self, new_direction: np.ndarray):
        """更新尾流方向（改攻角时调用，网格不变）。"""
        ...
```
#### CutElementMap（`wake/cut_elements.py`）
```python
class CutElementMap:
    """被尾流 level-set 切到的单元的识别与 DOF 管理。
    核心数据：
    - cut_elems: 被切单元的 element id 列表
    - node_side: 每个节点在切单元中的侧标记 (+1 / -1 / 0=面上)
    - ext_dof_offset: 从 DOF 的全局编号偏移
    - n_ext_dofs: 从 DOF 总数
    """
    def __init__(
        self,
        elements: np.ndarray,       # (n_tets, 4) 全局单元连接
        phi_wake: np.ndarray,       # (n_nodes,) level-set 值
    ):
        # 1. 找被切单元：四面体 4 个节点中 φ_wake 既有正又有负
        # 2. 为每个被切单元的节点标记侧
        # 3. 为从 DOF 编号：被切单元中的每个节点获得一个 ext DOF
        #    编号方式：ext_dof[i] = n_main_dofs + i
        ...
        self.cut_elems = ...        # (n_cut,) 被切单元 id
        self.n_ext_dofs = ...       # 从 DOF 数量
        self.total_dofs = n_main + self.n_ext_dofs
    def get_dof_assignment(self, e: int) -> np.ndarray:
        """返回单元 e 的 DOF 分配表。
        非切单元: [n0, n1, n2, n3]（标准主 DOF）
        切单元: [(n0或ext0), (n1或ext1), (n2或ext2), (n3或ext3)]
            "+"侧节点 → ext DOF
            "−"侧节点 → 主 DOF
        """
        ...
```
**从 DOF 粒度的设计决策**：
建议用**逐节点从 DOF**（每个被切节点一个从 DOF），而非 Núñez 2022 原文的逐单元从 DOF。如果一个节点出现在多个切单元中，它只有一个从 DOF——这对应于"$\Phi_{\mathrm{ext}}$ 是全局连续的"这个物理假设，与主 DOF 的连续性一致。
|             | 逐单元从 DOF（Núñez 2022 原文）  | 逐节点从 DOF（建议）                                         |
| ----------- | -------------------------------- | ------------------------------------------------------------ |
| 从 DOF 数量 | $\sim 4 \times n_{\mathrm{cut}}$ | $\sim n_{\mathrm{cut,nodes}} \approx 1.5 \times n_{\mathrm{cut}}$ |
| 矩阵稀疏性  | 更稀疏（纯局部）                 | 稍密（跨切单元耦合）                                         |
| 实现复杂度  | 更简单（纯局部）                 | 需要全局节点去重                                             |
| 物理正确性  | ✓                                | ✓（从 DOF 共享 = 连续性）                                    |
#### MultivaluedFE（`wake/multivalued.py`）
```python
class MultivaluedFE:
    """多值有限元空间 V_h = V_h^main ⊕ V_h^ext 的管理器。
    职责：
    - 构建扩展 DOF 空间的 CSR 稀疏模式
    - 提供 DOF 映射：给定单元和侧，返回 4 个 DOF 编号
    - 管理边界条件：Dirichlet 只作用于主 DOF
    """
    def __init__(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        cut_map: CutElementMap,
    ):
        # 构建扩展 CSR 模式：
        # - 非切单元：标准 P1 Galerkin 模式（主 DOF × 主 DOF）
        # - 切单元 "+" 侧：ext DOF 在 + 侧节点, main DOF 在 - 侧节点
        # - 切单元 "−" 侧：main DOF 在所有节点
        ...
        self.n_total = n_main + cut_map.n_ext_dofs
        self.pattern = ...          # (n_total, n_total) CSR
    def assemble_matrix(self, rho_tilde, B, V, cut_map) -> sp.csr_matrix:
        """装配 Picard 矩阵。
        非切单元：A_ij = ρ̃_e (∇N_i·∇N_j) V_e（与现在相同）
        切单元 "+" 侧贡献：用 + 侧 DOF 编号装配
        切单元 "−" 侧贡献：用 − 侧 DOF 编号装配
        → 切单元的矩阵贡献 = A^+ + A^-
        """
        ...
```
**切单元装配不需要几何切割**（§3.3）。一个切单元 $e$ 的矩阵贡献是：
$$A^e = \tilde{\rho}_e \cdot V_e \cdot (B_e^+ \cdot B_e^{+T} + B_e^- \cdot B_e^{-T})$$
$B_e^+$ 和 $B_e^-$ 是**相同的** $B_e$，但 DOF 编号不同。不需要切体积积分。
#### KuttaPenalty（`wake/kutta.py`）
```python
class KuttaPenalty:
    """罚函数 Kutta 条件（Núñez 2022 Eq.(44)-(45)）。
    在 Ω_kutta 单元上施加：n_kutta · ∇Φ|_TE = 0
    罚项：k_kutta * ∫_Ωkutta (n_kutta · ∇Φ)^2 dΩ
    其中 n_kutta 是尾缘角平分线的法向。
    """
    def __init__(
        self,
        te_point: np.ndarray,
        te_bisector_normal: np.ndarray,    # 尾缘角平分线法向
        elements: np.ndarray,
        nodes: np.ndarray,
        cut_map: CutElementMap,
        geom_levelset: np.ndarray = None,  # 几何 level-set（嵌入式几何）
        k_kutta: float = 100.0,
    ):
        # 识别 Ω_kutta：同时被尾流 level-set 和几何 level-set
        # 切到的单元及其邻居（Núñez 2022 Algorithm 1）
        # body-fitted 几何时：TE 节点所在单元及其邻居
        ...
        self.kutta_elems = ...      # Ω_kutta 单元列表
        self.n_kutta = te_bisector_normal
    def assemble_penalty(self, B, V, dof_map) -> (data, rows, cols):
        """装配罚项矩阵和残差贡献。
        罚项矩阵（López Eq.(46)）：
        K_kutta[i,j] = k_kutta * V_e * (n·B_e[i]) * (n·B_e[j])
        这是一个秩-1 修正，加到全局矩阵的 Ω_kutta 相关 DOF 上。
        """
        ...
```
---
## 5. 装配公式
### 5.1 残差（切单元）
对被切单元 $e$，设节点的侧标记为 $s_k \in \{+1, -1\}$：
**主 DOF 的残差贡献（"$-$" 侧）**：
$$R^-_i = \tilde{\rho}_e \cdot V_e \cdot (B_e \cdot B_e^T)_{ij} \cdot \Phi_j^{\mathrm{main}} \quad \forall j$$
**从 DOF 的残差贡献（"$+$" 侧）**：
$$R^+_i = \tilde{\rho}_e \cdot V_e \cdot (B_e \cdot B_e^T)_{ij} \cdot \Phi_j^+ \quad \forall j$$
其中：
$$\Phi_j^+ = \begin{cases} \Phi_{\mathrm{ext},j} & \text{if } s_j = +1 \\ \Phi_{\mathrm{main},j} & \text{otherwise} \end{cases}$$
等价地，切单元的残差 = 两个标准 P1 残差之和，区别只在 DOF 编号。
### 5.2 Picard 矩阵（切单元）
**"$-$" 侧贡献**（所有节点用 main DOF）：
$$A^-[i_{\mathrm{main}}, j_{\mathrm{main}}] \mathrel{+}= \tilde{\rho}_e \cdot V_e \cdot (B_e \cdot B_e^T)_{ij}$$
**"$+$" 侧贡献**（$+$ 侧节点用 ext DOF, $-$ 侧节点用 main DOF）：
$$A^+[i_{\mathrm{ext}}, j_{\mathrm{ext}}] \mathrel{+}= \tilde{\rho}_e \cdot V_e \cdot (B_e \cdot B_e^T)_{ij} \quad (\text{both } +)$$
$$A^+[i_{\mathrm{ext}}, j_{\mathrm{main}}] \mathrel{+}= \tilde{\rho}_e \cdot V_e \cdot (B_e \cdot B_e^T)_{ij} \quad (i+, j-)$$
$$A^+[i_{\mathrm{main}}, j_{\mathrm{ext}}] \mathrel{+}= \tilde{\rho}_e \cdot V_e \cdot (B_e \cdot B_e^T)_{ij} \quad (i-, j+)$$
### 5.3 尾流面边界条件
**质量通量连续 + 压力相等**：通过 Núñez 2022 Eq.(31) 的最小二乘泛函施加：
$$\Pi(\Phi^+, \Phi^-) = \frac{1}{2} \int_{\partial \Omega_W} \|\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^-\|^2 \, d\partial\Omega$$
对应的残差和 Jacobian 加到切单元的对应 DOF 上。关键性质：Jacobian 等于主问题 Jacobian 的差（$J_B = J_W^+ - J_W^-$）。
**在 TE 处**：通过罚函数施加 Kutta 条件（§3.5）。远场部分对直尾流，两侧速度趋近来流，压力自动相等。
### 5.4 Kutta 罚函数
在 $\Omega_{\mathrm{kutta}}$ 单元上（Núñez 2022 Eq.(45)-(46)）：
$$\Pi_{\mathrm{kutta}} = \frac{k_{\mathrm{kutta}}}{2} \sum_{e \in \Omega_{\mathrm{kutta}}} V_e \cdot (\mathbf{n}_{\mathrm{kutta}} \cdot \nabla \Phi_e)^2$$
$$\nabla \Phi_e = \sum_k \Phi_k \cdot B_e[k] \quad \text{(DOF 编号取决于切单元侧)}$$
对矩阵的贡献：
$$K_{\mathrm{kutta}}[i,j] = k_{\mathrm{kutta}} \cdot V_e \cdot (\mathbf{n}_{\mathrm{kutta}} \cdot B_e[i]) \cdot (\mathbf{n}_{\mathrm{kutta}} \cdot B_e[j])$$
秩-1 矩阵，加到 $\Omega_{\mathrm{kutta}}$ 单元对应的 DOF 上。
### 5.5 $\Gamma$ 的提取
不再有显式 $\Gamma$ 参数。解收敛后，$\Gamma$ 从解中提取：
```python
gamma = phi[upper_te_node] - phi[lower_te_node]
# 或更精确地：在 TE 附近的切单元上积分 [φ]
```
---
## 6. 求解流程
### 6.1 新的 Picard 循环
```python
# 预处理（一次性）
mesh = read_mesh("naca0012.msh")           # 无尾流面
wls = WakeLevelSet(te_point, direction, normal)
phi_wake = wls.evaluate(mesh.nodes)
cut_map = CutElementMap(mesh.elements, phi_wake)
mfe = MultivaluedFE(mesh.nodes, mesh.elements, cut_map)
kutta = KuttaPenalty(te_point, te_bisector, mesh, cut_map, k_kutta=100.0)
op = PicardOperator(mesh.nodes, mesh.elements)  # B_e, V_e, coloring 不变
# 初始化
phi = np.zeros(mfe.n_total)                # 主 DOF + 从 DOF
phi[main_dofs[farfield]] = farfield_values  # Dirichlet 只作用于主 DOF
# Picard 循环（无 Γ 外层！）
for outer in range(n_picard_max):
    grad, q2 = mfe.velocities(phi, op)      # 扩展的单元速度计算
    rho = density_field(q2 / u_inf**2, m_inf)
    rho_tilde = upwind(rho, q2, ...)
    # 装配（非切单元 + 切单元 + 罚项）
    A = mfe.assemble_matrix(rho_tilde, op.B, op.V, cut_map)
    A += kutta.assemble_penalty(op.B, op.V)
    # Dirichlet 消元（只作用于主 DOF）
    A_free, b_free = apply_dirichlet(A, b, main_dirichlet_nodes, values)
    # 线性求解
    phi_free = solve_cg_amg(A_free, b_free)
    # 收敛判断
    if ||residual|| < tol:
        break
```
### 6.2 与当前流程的对比
| 步骤          | 当前                                                         | 新方案                                             |
| ------------- | ------------------------------------------------------------ | -------------------------------------------------- |
| 网格          | Gmsh embed wake                                              | Gmsh 标准（无 wake）                               |
| 尾流预处理    | wake_cut.py 物理复制                                         | CutElementMap 逻辑标记                             |
| $\Gamma$ 参数 | 显式，外层 secant 迭代                                       | 隐式，从解中涌现                                   |
| Kutta         | $\Gamma_{\mathrm{target}} = \varphi_{\mathrm{TE}}^+ - \varphi_{\mathrm{TE}}^-$ | 罚函数 $\mathbf{n} \cdot \nabla \Phi = 0$          |
| 矩阵大小      | $2 \times N$（master+slave 消元后 $\approx N$）              | $N + N_{\mathrm{ext}}$（$N_{\mathrm{ext}} \ll N$） |
| 外层循环      | Mach × $\Gamma$ secant × Picard                              | Mach × Picard（少一层）                            |
---
## 7. 线性求解器适配
扩展矩阵的结构：
$$A = \begin{bmatrix} A_{mm} & A_{me} \\ A_{em} & A_{ee} \end{bmatrix} \quad \begin{aligned} & m: \text{主 DOF } (N) \\ & e: \text{从 DOF } (N_{\mathrm{ext}} \ll N) \end{aligned}$$
- $A_{mm}$：主-主块，与当前 P1 Galerkin 矩阵结构相同
- $A_{me}$ / $A_{em}$：主-从耦合块，只在切单元处有非零（稀疏）
- $A_{ee}$：从-从块，只在切单元处有非零（很小）
**AMG 策略**：
- 对 $A_{mm}$ 用标准 PyAMG smoothed aggregation（与现在相同）
- 对 $A_{ee}$ 块用简单 Jacobi 或 ILU（很小，不需要 AMG）
- 整体用块预处理：
$$M = \begin{bmatrix} \mathrm{AMG}(A_{mm}) & 0 \\ -A_{em} \cdot \mathrm{AMG}(A_{mm}) & \mathrm{ILU}(A_{ee}) \end{bmatrix}$$
或者更简单——直接对整个 $A$ 做 PyAMG（$N + N_{\mathrm{ext}}$ 不比 $N$ 大多少），让 AMG 自动处理。
---
## 8. 实现阶段
### Phase B1：Level-set 尾流 + 切单元识别（验证基础设施）
**目标**：给定网格和 level-set 函数，正确识别切单元并建立 DOF 映射。
**交付物**：
- `wake/levelset.py`：直尾流 level-set（从 TE 点沿来流方向的半平面）
- `wake/cut_elements.py`：切单元识别 + 逐节点从 DOF 编号
- 测试：在现有 NACA 2.5D 网格上（不嵌入尾流面），验证切单元识别正确
**验证标准**：
- 切单元数量合理（$\approx$ 尾流面穿过的单元数）
- 从 DOF 数量 $\approx$ 切单元节点数去重
- 侧标记与几何一致（$+$ 侧在上翼面，$-$ 侧在下翼面）
### Phase B2：多值 FE 装配（核心数值验证）
**目标**：在 Laplace（$\rho \equiv 1$）问题上验证多值 FE 装配的正确性。
**交付物**：
- `wake/multivalued.py`：扩展 DOF 空间 + CSR 模式
- `kernels/cut_assembly.py`：切单元双侧残差 + 矩阵装配
- 修改 `kernels/jacobian.py`：PicardOperator 支持扩展 DOF
- 修改 `solve/linear.py`：扩展矩阵的 Dirichlet 消元
**验证标准**：
- V0 自由流保持：$\|R(\varphi = U_\infty \cdot \mathbf{x})\|_\infty < 10^{-12}$
- V1 MMS 收敛：L2 误差 slope $\geq 1.9$（与当前一致）
- 关键对比：在新方案下用 Laplace 求解 NACA0012 非升力（$\alpha = 0$），$c_l \approx 0$
### Phase B3：罚函数 Kutta + 升力求解
**目标**：实现罚函数 Kutta，恢复升力求解能力。
**交付物**：
- `wake/kutta.py`：$\Omega_{\mathrm{kutta}}$ 识别 + 罚项装配
- 修改 `solve/picard.py`：无 $\Gamma$ 外层的 Picard 循环
- 修改 `constraints/dirichlet.py`：远场涡修正改用 level-set 尾流方向
**验证标准**：
- V3 亚声速 NACA0012 $M = 0.5$, $\alpha = 2°$：$c_l$ 在 [PG, KT] 区间内（与 P3 一致）
- $\Gamma$ 从解中提取的值与当前 secant 收敛值一致（$< 1\%$ 差异）
- 无 $\Gamma$ 外层循环，Picard 迭代数 $\leq$ 当前
### Phase B4：跨声速 + Mach continuation
**目标**：恢复 P4 跨声速能力。
**交付物**：
- 修改 `solve/continuation.py`：去掉 $\Gamma$ 外层，纯 Mach continuation + Picard
- 适配 `kernels/upwind.py`：确认人工密度在切单元上正确工作
- 继承 `damping_theta`（local $\theta \cdot \mathrm{diag}(A_{\mathrm{free}})$）阻尼机制
**验证标准**：
- V4 NACA0012 $M = 0.80$, $\alpha = 1.25°$：激波位置 $x/c \approx 0.60$（与 P4 一致）
- 无 $\Gamma$ 外层循环后，总 Picard 迭代数减少（预期 30-50%）
- $c_l$ 与 P4 结果一致（$< 2\%$ 差异）
> **注（2026-07-07 22:45 更新）**：P4 当前性能基线——coarse 174s/10464 Picard its，medium 16m39s。`damping_theta=0.2` 是默认值，方案 B4 应直接复用。去掉 $\Gamma$ 外层后，`solve/continuation.py` 的 `max_gamma_evals=12` × `n_picard_eval=800` 嵌套消除，预期 medium 从 16m39s 降到 $\sim 8$-10min。
### Phase B5：多尾流验证
**目标**：验证多物体绕流。
**交付物**：
- 支持多个 WakeLevelSet 实例
- 多段翼网格生成（不需要嵌入任何尾流面）
- 多 $\Omega_{\mathrm{kutta}}$ 识别
**验证标准**：
- 两段翼（如 NLR 7301）：各自 $c_l$ 合理
- 翼身组合体：机身无升力，机翼升力不受影响
### Phase B6（搁置）：曲线尾流 / 自由尾迹

> **决策（2026-07-10）：搁置，不实施。** 详见 `20260707_2118_ibl_viscous_coupling_design.md` §4.5.6。原因：(1) 精度收益微小（$O(\theta_{\mathrm{wake}}^2) \sim O(0.1\%)$）；(2) 每次方向变化需重建 CutElementMap + enrichment DOF + 约束系统，工程代价不成比例；(3) 被切单元集离散跳变与 Newton 紧耦合结构性冲突；(4) López 2021 先例；(5) 优先级排序。保留方案 B `update_direction()` 的接口能力以备未来需要。

**目标**：尾流面几何迭代对齐流场。
**交付物**：
- `wake/levelset.py` 增加 `update_from_velocity()` 方法
- 外层尾流几何迭代循环
**验证标准**：
- 大攻角下尾流弯曲，$c_l$ 与直尾流有可见差异
- 尾流面收敛（几何变化 $<$ 容差）
---
## 9. 风险与缓解
| 风险                                                         | 严重度 | 缓解                                                         |
| ------------------------------------------------------------ | ------ | ------------------------------------------------------------ |
| 罚函数 $k_{\mathrm{kutta}}$ 的值敏感（太大病态，太小不满足 Kutta） | 中     | 从 $k = 100$ 开始，参考 López 的参数；用 V3 gate 扫描        |
| AMG 在扩展矩阵上性能下降                                     | 中     | 先尝试整体 AMG；如果不行用块预处理                           |
| 切单元的双侧装配导致矩阵条件数变差                           | 低     | 切单元数 $\ll$ 总单元数，影响有限                            |
| 从 DOF 的连续性假设可能不适用于卷起尾流                      | 低     | B6 阶段才需要考虑，B1-B5 不受影响                            |
| 现有 P1-P4 测试套件需要适配                                  | 中     | 保留旧的 wake_cut.py 作为 reference path，新方案通过后切换   |
| 切单元的侧体积分数近似（$V_e/2$）精度不够                    | 低     | P1 四面体上 $B_e$ 是常数，近似误差可控；如需精确可计算子体积 |
---
## 10. 与现有代码的共存策略
**不采用大爆炸式重写。** 两条路径并存：
1. `solve/picard.py` 保留 `solve_subsonic_lifting`（当前路径，用 wake_cut + $\Gamma$ secant + `damping_theta`）
2. 新增 `solve/picard_ls.py`（level-set 路径，用 multivalued FE + 罚函数 Kutta）
3. 每个 Phase 验证通过后，逐步切换默认路径
4. P5（ONERA M6）可以直接用新路径
> **注（2026-07-07 22:45 更新）**：当前 `solve_subsonic_lifting` 已新增 `damping_theta` 参数（local $\theta \cdot \mathrm{diag}(A_{\mathrm{free}})$ 阻尼），`constraints/wake.py::WakeConstraint.update_matrix` 已实现 $h_j$ 批处理化。方案 B 的 `picard_ls.py` 应继承这两个改进——`damping_theta` 直接复用，$h_j$ 批处理化在 level-set 路径下无对应物（无 $\Gamma$ RHS），不影响。
> 测试套件通过参数化同时跑两条路径，确保一致性。
---
## 11. 参考文献
### 核心文献（已下载到 `/home/z/my-project/references/`）
- **Davari, M., Rossi, R., Dadvand, P., López, I., Wüchner, R.** "A cut finite element method for the solution of the full-potential equation with an embedded wake." *Computational Mechanics* 63(5):821–833, 2019. DOI: 10.1007/s00466-018-1624-3.
  - §2.2: body-fitted 几何 + level-set 尾流面定义
  - Eq.5/10: 壁面 Neumann BC ($g=0$)
  - Eq.6–7: 尾流面质量通量连续 + 压力相等
  - Eq.23–27: 压力相等线性化 + 最小二乘泛函（小扰动假设）
  - Eq.28–32: Heaviside 扩展形函数 + DOF 复制（$\Phi/\Psi$）
  - Eq.33–34: 切单元残差（Heaviside 形式）
  - §3.2 + [33]: 与 boundary layer solver 耦合的方向
  - Appendix A: 边界层理论证明压力沿法向不变
- **Núñez, M., López, I., Baiges, J., Rossi, R.** "An embedded approach for the solution of the full potential equation with finite elements." *Computer Methods in Applied Mechanics and Engineering* 388:114244, 2022.
  - Abstract: 全嵌入式——几何和尾流都用 level-set
  - §2.3: 区分 body-fitted (Davari [8]) 与全嵌入式
  - §2.4: level-set 尾流 + 多值 FE 空间分解
  - Eq.11–12: 切单元 DOF 分配
  - Eq.24–25: 切单元无穿透条件（几何 level-set）
  - Eq.28–29: 切单元残差
  - Eq.30: 矢量尾流约束 $\rho^+ \mathbf{v}^+ = \rho^- \mathbf{v}^-$
  - Eq.31: 尾流面 BC 最小二乘泛函
  - Eq.34–39: 约束残差和 Jacobian
  - Eq.40–43: 切单元矩阵组装
  - Eq.44–46: Kutta 罚函数
  - Algorithm 1: $\Omega_{\mathrm{kutta}}$ 识别（仅全嵌入式需要）
### 补充文献
- López Canalejo, I., Baez Jones, E., Núñez, M., Zorrilla, R., Rossi, R., Bletzinger, K.-U., Wüchner, R. "A transonic potential solver with an embedded wake approach using multivalued finite elements." *Coupled Problems 2021*.
  - 多值 FE 的工程实现细节
  - Kratos MultiPhysics 框架中的实现
- Rodriguez, D., Sturdza, P., Suzuki, Y., Martins-Rivas, H., Peronto, A. "A rapid, robust, and accurate coupled boundary-layer method for CART3D." AIAA 2012-302. — Davari 2019 [33]，全速势 + IBL 耦合的工程实例
