# IBL 粘性修正耦合方案设计（修订版）

> Design note for UP3D / pyFP3D
> Date: 2026-07-07（updated 2026-07-09 02:28 — Drela IBL3 路线确定 + BLWF/TRANAIR 对比）
> Code baseline: HEAD 170eaa9
> Prerequisite: P6（smooth wall gradient）; P8（Newton Jacobian，V3 紧耦合需要）
> **本文档替代原版中的 Green's lag entrainment 方案，改用 Drela 2013 IBL3 6 方程形式**
> References:
> - **Drela 2013** (AIAA 2013-2437) — 3D IBL3 完整推导 — `references/Drela_2013_IBL3_general_configurations.pdf`
> - **BLWF58 方法文档** (TsAGI, 108pp) — quasi-simultaneous coupling 参考 — `references/BLWF58_method_algorithms.doc`
> - López Canalejo 2021 博士论文 (TUM) — §4.1 粘性偏差, §5.2 Outlook
> - Davari et al. 2019 (Comput. Mech. 63:821–833) — body-fitted FPE + 嵌入尾流

## 0. 文档目的
讨论在全速势方程框架下，如何通过积分边界层法（IBL）引入粘性修正，实现无粘-有粘迭代（VII）。
核心问题：
1. IBL 的输出如何影响全速势方程的边界条件？
2. IBL 的输出如何影响尾流面？
3. 是否需要自适应更新网格？
4. Level-set 能否用于引入 IBL 的位移厚度贡献？还是修改边界条件更合理？

## 1. 文献澄清：几何与尾流的 level-set 使用
### 1.1 两篇论文的引用关系
| 论文                                                         | 简称        | 几何                 | 尾流               | UP3D 对应                  |
| ------------------------------------------------------------ | ----------- | -------------------- | ------------------ | -------------------------- |
| Davari et al. 2019, *Comput. Mech.* 63:821–833               | Davari 2019 | body-fitted          | level-set (Cut FE) | **UP3D 方案 B 的直接前作** |
| Núñez et al. 2022, *CMAME* 388:114244                        | Núñez 2022  | level-set (全嵌入式) | level-set          | 不采用                     |
| 注意：Núñez 2022 的第一作者是 Marc Núñez（CIMNE/UPC），López 是第三作者。之前文档中将其称为"López 论文"不够准确，以下用"Núñez 2022"指代。 |             |                      |                    |                            |
### 1.2 Núñez 2022 的完整做法（全嵌入式）
Núñez 2022 论文中，**几何物面和尾流面都用了 level-set**：
> *A fully embedded approach to solve the full-potential equation is presented, where **both the geometry and the wake are defined implicitly with a level set function**.* (Abstract)
> 论文用了**两个 level-set 函数**：
> | Level-set                                                    | 定义             | 作用                                                         |
> | ------------------------------------------------------------ | ---------------- | ------------------------------------------------------------ |
> | $\varphi(\mathbf{x})$                                        | 几何符号距离函数 | $\varphi=0$ 定义物面，$\varphi \geq 0$ 流体域，$\varphi < 0$ 物体内部 (Eq.10) |
> | $\varphi_{\mathrm{wake}}(\mathbf{x})$                        | 尾流符号距离函数 | $\varphi_{\mathrm{wake}}=0$ 定义尾流面                       |
> | 物面边界条件通过切单元积分（只积分 $\varphi \geq 0$ 侧，即 $\Omega_{\varphi,\mathrm{in}}$）隐式施加，无穿透条件 $\mathbf{n} \cdot \nabla \Phi = 0$ 通过"切单元内边界项置零"实现 (Eq.24–25)。$\Omega_{\mathrm{kutta}}$ 识别（Algorithm 1）就是找**同时被 $\varphi$ 和 $\varphi_{\mathrm{wake}}$ 切到的单元**——这只有在几何也是 level-set 时才需要。 |                  |                                                              |
### 1.3 Davari 2019 的做法（body-fitted 几何 + 嵌入式尾流）
Davari 2019 是 Núñez 2022 的前作（Núñez 2022 参考文献 [8]）。Núñez 2022 §2.3 明确区分了两者：
> *While a level set function is used in [8] to define the wake, **the airfoil shapes considered are still body-fitted**. In this paper, **both the geometry of study and the wake are modelled using embedded methods**.*
> Davari 2019 的技术方案：
- **几何**：body-fitted 网格，物面由网格面贴合，后缘节点明确存在
- **尾流**：用 level-set 函数 $\varphi_{\mathrm{wake}}$ 隐式定义，尾流面切割背景网格单元
- **DOF 复制**：被尾流切割的单元，节点被复制为原 DOF ($\Phi$) 和辅助 DOF ($\Psi$)，通过 Heaviside 函数构造扩展形函数 (Eq.28–32)
- **尾流 BC**：质量守恒 (Eq.6) + 压力相等 (Eq.7)，通过最小二乘泛函施加
- **Kutta 条件**：后缘节点明确存在，尾流从后缘出发，不需要特殊识别算法
### 1.4 Davari 2019 对粘性修正的讨论
Davari 2019 已经明确提到了与边界层求解器耦合的可能性：
> *An improved accuracy of the method can be achieved with the coupling of the full potential solver with a boundary layer solver as shown in [33].* (§3.2)
> 其中 [33] 是 Rodriguez et al. 2012 (AIAA 2012-302)，一个快速 coupled boundary-layer 方法用于 CART3D。
> 此外，Davari 2019 Appendix A 用边界层理论证明了"压力沿物面法向不变"——这正是势流方法能准确预测升力的理论基础，也是 VII 耦合的理论依据。
### 1.5 UP3D 的定位
**UP3D 采用 body-fitted 几何**——物面由网格面贴合，不使用 level-set 表示几何。这对应 Davari 2019 的架构，而非 Núñez 2022 的全嵌入式架构。
UP3D 的 level-set 仅用于**尾流面**（方案 B），不用于物面。这是一个明确的设计决策。

## 2. IBL 方程选择：Drela IBL3（6 方程形式）

### 2.1 为什么选 Drela IBL3

原版文档推荐 Green's lag entrainment（3 方程）。基于 Drela 2013 (AIAA 2013-2437) 的分析和 BLWF58/TRANAIR 的对比（见 `20260709_0145_3d_vii_implementation_analysis.md` §15），改为 **Drela IBL3 6 方程**：

| 方面 | Green's 3-eq | **Drela IBL3 6-eq** | BLWF Keller box |
|------|-------------|---------------------|-----------------|
| 分离处理 | ❌ 不能 | ✅ 4 方程自然允许 | ✅ inverse mode 切换 |
| crossflow | ❌ 无 | ✅ Ψ 方程 | ✅ 直接 3D PDE |
| 计算量 | 基准 | ~2× Green | ~100× IBL3 |
| 闭合关系验证 | XFOIL | XFOIL/MSES 30 年 | 不需要 |
| 表面网格 | 流线积分 | **Galerkin FE 弱形式** | 结构化贴体 |

### 2.2 IBL3 的 6 个未知量和方程

参考 Drela 2013 §II，在壁面表面网格的每个节点上求解：

| 未知量 | 含义 | 方程 | Drela 2013 Eq. |
|--------|------|------|----------------|
| $\delta$ | 位移厚度 | 质量 | Eq.(21) |
| $A$ | 动能厚度形参 | 动能 | Eq.(28) |
| $B$ | 横向动量厚度 | 动量 | Eq.(24) |
| $\Psi$ | profile twist（crossover） | 横向曲率 | Eq.(29) |
| $C_{\tau1}$ | Reynolds stress lag（流向） | stress transport | §II.F |
| $C_{\tau2}$ | Reynolds stress lag（横向） | stress transport | §II.F |

**层流区**：$C_{\tau1}, C_{\tau2}$ 退化为不稳定波包络增长率（$e^N$ 转捩），与 XFOIL 一致。

**与经典 BL 的一致性**（Drela 2013 §II.I）：小曲率近似下退化为 von Kármán 动量积分方程。

### 2.3 Drela 的表面 FE 方法——不沿流线积分

**关键优势**：Drela 2013 §III 不沿流线积分——在表面网格上用 Galerkin 弱形式求解 IBL3 方程。这完全避开了非结构网格上的流线追踪问题。

每个壁面三角形上：
1. 定义局部 Cartesian 基 $(\hat{x}_i, \hat{y}_i, \hat{z}_i)$（Drela 2013 §III.1, Eq.(64)）
2. 速度和积分厚度用基函数插值（Drela 用 bilinear 四边形，UP3D 用 P1 三角形——退化情况）
3. 方程残差在节点上组装——与 UP3D 的体 FE 方法一致

### 2.4 物面与尾迹的统一处理

IBL3 在物面和尾流面上求解**同一套 6 方程**——不是两套独立方程。Drela 2013 §IV.C 原文：

> "the surface and wake paneling also serves as the boundary layer grid. Each grid node has the six primary viscous unknowns (δ A B Ψ Cτ1 Cτ2), governed by the six integral boundary layer equations"

DCV/DDCV 框架（§II.H, Figure 2-3）的设计让同一套积分方程自然覆盖物面 BL 和尾迹——DCV 是"所有有粘性缺陷的区域"。TE 处不需要特殊切换——6 方程在 TE 节点处自动从壁面边界条件过渡到尾迹边界条件。

物面与尾迹的区别只在闭合关系（速度剖面形状）：

| 区域 | 速度剖面边界条件 | 闭合关系 |
|------|----------------|---------|
| **物面** | no-slip: $\mathbf{q}_w = 0$ | wall profile（Coles wake + sublayer，Drela 2013 §II.C） |
| **尾迹** | 两侧速度相等: $\mathbf{q}^+ = \mathbf{q}^-$ | wake profile（中心对称，无壁面） |

**IBL3 的输出**：
- $\delta^*_{\mathrm{wall}}(s)$ — 物面位移厚度
- $\theta(s)$, $H(s)$, $C_f(s)$ — 动量厚度、形状因子、摩擦系数
- $\delta^*_{\mathrm{wake}}(s)$ — 尾迹位移厚度（从 TE 向下游衰减，同一套方程自然延续）
- $A, B, \Psi, C_{\tau1}, C_{\tau2}$ — 其他积分量（动能、横向动量、profile twist、stress lag）
- 分离线位置（4 方程自然捕获，Drela 2013 §IV.C "analogous to captured shocks"）

**对 UP3D 的意义**：尾流面 IBL 修正（V4）不是独立于 V1 的额外工作——它是 IBL3 方程的自然延续。只需在尾流面网格上也定义 6 个未知量，用尾迹闭合关系替代壁面闭合关系。V1 和 V4 可以在同一个求解器中实现，区别只是闭合关系的切换。

## 3. $\delta^*$ 影响全速势方程的两种传统方法
### 3.1 方法 A：等效物体（几何修正）
把有效物面从真实物面外推 $\delta^*(s)$：
$$y_{\mathrm{eff}}(s) = y_{\mathrm{wall}}(s) + \delta^*(s)$$
势流方程在加厚的物面上施加 $\mathbf{n} \cdot \nabla \varphi = 0$。
**特点**：
- 物理直观——等效物体就是"边界层把物面撑开了"
- **需要网格变形或重新生成**——有效几何变了，网格必须跟上
- 传统工业代码（FLO 系列早期版本）用过这个方法
### 3.2 方法 B：Transpiration 边界条件
保持网格不动，修改壁面 Neumann 条件。从边界层连续性方程在薄层近似下推出：
| 条件                                                         | 公式                                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| 无粘                                                         | $\mathbf{n} \cdot (\rho \nabla \Phi) = 0$（不可穿透，$g=0$ in Davari Eq.5） |
| 有粘修正后                                                   | $\mathbf{n} \cdot (\rho \nabla \Phi) = \frac{d}{ds}[\rho_e U_e \delta^*(s)]$（质量穿透） |
| 物理含义：边界层的排挤效应等效为在壁面上吹/吸质量，流量等于位移厚度的沿程变化率。 |                                                              |
| 在 FEM 框架中，Davari 2019 Eq.10 的弱形式右端项为 $\int_{\Gamma_N} N_i \, g \, d\Gamma$，其中 $g=0$ 对应壁面无穿透。Transpiration 修正就是让 $g \neq 0$： |                                                              |
| $$R_{\mathrm{wall},i} \mathrel{+}= \int_{\partial \Omega_{\mathrm{wall}}} N_i \cdot \dot{m}(s) \, d\partial\Omega$$ |                                                              |
| 其中 $\dot{m}(s) = \frac{d}{ds}[\rho_e U_e \delta^*(s)]$。   |                                                              |
| **可压缩性说明**：$\rho_e$ 和 $U_e$ 由势流解在壁面处的边界层外缘给出。对可压缩流，$\rho_e$ 通过 Bernoulli 关系 (Davari Eq.2 / Núñez Eq.2) 从 $\nabla \Phi$ 求得。$\delta^*$ 的定义在可压缩流中需考虑密度变化： |                                                              |
| $$\delta^*_{\mathrm{compressible}} = \int_0^\infty \left(1 - \frac{\rho u}{\rho_e U_e}\right) dn$$ |                                                              |
| IBL 求解器需要使用可压缩动量积分方程（如 Green's lag entrainment method）。 |                                                              |
| **特点**：                                                   |                                                              |
- **不需要动网格**——只修改边界条件
- 工程实现简单——在壁面 Neumann 项中增加一个源
- 传统工业代码（FLO 系列、EPPLER、XFOIL 全速势模式）几乎都选这个方法
- Davari 2019 §3.2 已明确指出全速势 + boundary layer solver 耦合是改进精度的方向
### 3.3 两种方法的等价性
在薄边界层近似下（$\delta^*/c \ll 1$），两种方法等价到一阶。差异在于：
|                                                              | 等效物体           | Transpiration |
| ------------------------------------------------------------ | ------------------ | ------------- |
| 网格                                                         | 需要变形/重生成    | 不动          |
| 壁面曲率                                                     | 改变了（几何偏移） | 保留原始曲率  |
| 大曲率区域（前缘）                                           | 有二阶差异         | 更准确        |
| 实现复杂度                                                   | 高（网格变形）     | 低（改 BC）   |
| 对于 UP3D 的典型工况（亚声速/跨声速薄翼），$\delta^*/c \sim 0.01\text{--}0.03$，两者差异可忽略。但分离附近 $\delta^*$ 急剧增长时差异变大——不过那时 IBL 本身已经不可靠了。 |                    |               |

## 4. $\delta^*$ 如何影响尾流面
### 4.1 Davari 2019 / Núñez 2022 的尾流面条件（纯无粘）
两篇论文的尾流面边界条件相同，均为两个条件：
1. **质量通量连续** (Davari Eq.6, Núñez Eq.8)：
$$\mathbf{n} \cdot (\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^-) = 0 \quad \text{on } \partial \Omega_W$$
2. **压力相等** (Davari Eq.7, Núñez Eq.9)：
$$|\nabla \Phi^+|^2 - |\nabla \Phi^-|^2 = 0 \quad \text{on } \partial \Omega_W$$
**关键实现细节**：这两个条件不是直接强施加的，而是通过**最小二乘泛函**弱施加。
Davari 2019 的做法（Eq.24–27）：将压力相等条件线性化（小扰动假设 $u \ll 1$），得到对流方程 $\mathbf{u}_\infty \cdot \nabla(\Phi^+ - \Phi^-) = 0$ (Eq.23)，然后构造最小二乘泛函：
$$\Pi(\Phi) = \frac{1}{2} \int_{\Gamma_W} \left[\mathbf{u}_\infty \cdot \nabla(\Phi^+ - \Phi^-)\right]^2 d\Gamma$$
Núñez 2022 的做法（Eq.30–31）：将两个条件合成一个矢量方程 $\rho^+ \mathbf{v}^+ = \rho^- \mathbf{v}^-$ (Eq.30)，然后构造最小二乘泛函：
$$\Pi(\Phi^+, \Phi^-) = \frac{1}{2} \int_{\partial \Omega_W} \|\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^-\|^2 \, d\partial\Omega$$
对应的残差 (Eq.34–35) 和 Jacobian (Eq.38–39) 都是从这个泛函求导得到的。关键性质：**约束项的 Jacobian 等于主问题 Jacobian 的差** (Eq.38: $J_B = J_W^+ - J_W^-$)，这是因为形函数在切单元的子区域上是分片线性的。
### 4.2 IBL 尾迹修正后的尾流面条件
边界层在后缘并没有结束，它延续为尾迹。尾迹也有位移厚度 $\delta^*_{\mathrm{wake}}(s)$，从后缘向下游衰减。
加入 IBL 尾迹修正后，尾流面的质量通量连续条件不再为零——尾迹的位移厚度贡献了一个沿尾流面的质量源：
| 条件                                                         | 公式                                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| 纯无粘                                                       | $\mathbf{n} \cdot (\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^-) = 0$ |
| IBL 修正                                                     | $\mathbf{n} \cdot (\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^-) = \frac{d}{ds}[\rho_e U_e \delta^*_{\mathrm{wake}}(s)]$ |
| 物理含义：尾迹的位移厚度使得尾流面两侧不再只是势函数跳跃，还有一个沿尾流面的质量源项，代表尾迹的排挤效应。 |                                                              |
### 4.3 在最小二乘泛函框架中的实现
**关键点**：不能简单地把 IBL 源项直接加到残差上。Davari/Núñez 的尾流 BC 是通过最小二乘泛函 $\Pi$ 弱施加的，IBL 修正需要修改泛函本身。
有两种实现路径：
**路径 ①：修改泛函约束（推荐）**
将 Núñez 2022 Eq.30 的约束修改为：
| 约束                                                         | 公式                                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| 原始                                                         | $\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^- = 0$ (Núñez Eq.30) |
| 修正                                                         | $\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^- = \mathbf{m}_{\mathrm{wake}}(s)$（引入源项） |
| 其中 $\mathbf{m}_{\mathrm{wake}}(s) = \frac{d}{ds}[\rho_e U_e \delta^*_{\mathrm{wake}}(s)]$ 是沿尾流面坐标的已知函数（由 IBL 给出）。 |                                                              |
| 最小二乘泛函变为：                                           |                                                              |
| $$\Pi_{\mathrm{modified}}(\Phi^+, \Phi^-) = \frac{1}{2} \int_{\partial \Omega_W} \|\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^- - \mathbf{m}_{\mathrm{wake}}(s)\|^2 \, d\partial\Omega$$ |                                                              |
| 对应的残差增加一项（与 $\mathbf{m}_{\mathrm{wake}}$ 相关的线性项），Jacobian 不变（因为 $\mathbf{m}_{\mathrm{wake}}$ 是已知函数，不依赖于 $\Phi$）。 |                                                              |
| **路径 ②：分解为法向 + 切向**                                |                                                              |
| 将矢量约束 $\rho^+ \mathbf{v}^+ - \rho^- \mathbf{v}^- = \mathbf{m}_{\mathrm{wake}}$ 分解： |                                                              |
- 法向分量：质量通量连续 + IBL 源 → 修改 Eq.8
- 切向分量：压力相等 → 不变（IBL 不直接修正压力跳跃）
路径 ① 更自然地嵌入到 Núñez 2022 的泛函框架中，实现改动量最小。
### 4.4 在 body-fitted 几何（UP3D/方案 B）中的特殊性
UP3D 采用 body-fitted 几何 + level-set 尾流（Davari 2019 架构），与 Núñez 2022 的全嵌入式有一个重要区别：
- **Davari 2019 / UP3D**：后缘节点明确存在，尾流从后缘出发。Kutta 条件通过后缘节点的 DOF 处理自然实现，不需要 $\Omega_{\mathrm{kutta}}$ 识别算法。
- **Núñez 2022**：后缘是 level-set 交叉点，需要 Algorithm 1 识别 $\Omega_{\mathrm{kutta}}$，并通过罚函数 (Eq.20) 施加 Kutta 条件。
因此，UP3D 的尾流面 IBL 修正**不需要处理 $\Omega_{\mathrm{kutta}}$ 区域的特殊性**——在 Kutta 区域（后缘附近），尾流面的 DOF 复制和 BC 施加方式与 Davari 2019 一致，IBL 源项只需在尾流面的常规切单元中加入。

### 4.5 TE 拐折与有攻角尾流偏折的处理

有攻角来流条件下，TE 处发生**三重叠加**的几何/物理效应：(1) 翼型 TE 楔角导致的剪切层几何折转；(2) 上下翼面 BL 在 TE 合并为尾迹；(3) 尾流受环量诱导下洗而偏折。以下分别分析 Drela IBL3 和 López 嵌入式尾流框架如何处理这些效应，并给出 UP3D 的分阶段方案。

#### 4.5.1 Drela IBL3 对 TE 拐折的处理：局部基自适应

Drela 2013 的核心设计是**每个节点拥有独立的局部 Cartesian 基** $(\hat{x}, \hat{y}, \hat{z})$（§III.B, Eq.(3), Figure 4）：

- $\hat{x}, \hat{z}$ → 沿剪切层表面切向
- $\hat{y}$ → 沿剪切层法向（TSL 近似下 $\hat{n}_w \simeq \hat{y}$）

关键性质（§III.B, Eq.(5) 后原文）：

> "the integral boundary layer equations derived later are **invariant to such in-plane coordinate rotation**, any such offset in $\psi$ has no effect."

这意味着即使 TE 处剪切层表面有几何折转（上翼面 → 尾流面、下翼面 → 尾流面的方向突变），每个节点用自己的局部基，**方程形式不变**。"拐折"被局部基的逐点旋转吸收——不需要额外的 TE 修正方程。

#### 4.5.2 TSL 近似在 TE 处的局部失效

薄剪切层近似（§II.E, Eq.(11)）丢弃了曲率项：

$$O\!\left(\frac{v}{q_i},\; \kappa(y_e - y_w)\right) \ll 1$$

其中 $\kappa$ 是剪切层表面的曲率。在尖锐 TE 处 $\kappa \to \infty$，TSL 近似**局部失效**。Drela 的处理方式是不做特殊处理——方程照常运行，依赖 6 方程在 TE 节点处的连续演化来"过渡"过去。实际验证（Drela 2013 §IV.C 双锥翼案例）表明这种处理在工程精度下可接受。

#### 4.5.3 上下翼面 BL 在 TE 的合并

DCV 框架下，TE 处的匹配是隐式的——6 个未知量 $(\delta, A, B, \Psi, C_{\tau1}, C_{\tau2})$ 在 TE 节点处连续过渡：

| 物理量 | TE 匹配条件 |
|--------|------------|
| $\delta$ | $\delta_{\mathrm{wake}}(0) = \delta^*_{\mathrm{upper}}(TE) + \delta^*_{\mathrm{lower}}(TE)$（厚度叠加） |
| $A, B$ | 从上下翼面的壁面流线方向过渡到尾迹中心线方向 |
| $\Psi$ | profile twist 连续过渡 |
| $C_{\tau1}, C_{\tau2}$ | stress lag 从两侧值过渡 |

关键点：Drela **不在 TE 做显式匹配**——6 方程在 TE 节点处自然演化。因为 DCV 框架下物面和尾流是同一个计算域，TE 只是闭合关系的切换点，不是方程的边界。

对 UP3D 的意义：master-slave TE 节点复制结构天然支持这一合并——上下面 BL 的 6 个未知量在 TE slave 节点处汇合，切换到尾迹闭合关系后继续沿尾流面积分。

#### 4.5.4 尾流的几何偏折

有攻角时尾迹中心线不沿来流方向——它受到机翼环量诱导的下洗（downwash）偏折。偏折角的主要来源：

1. **下洗角** $\epsilon_i \approx C_L / (\pi AR)$（有限翼展），2D 中更小
2. **TE 后流动偏转** ≈ 0（Kutta 条件要求流动平滑离开 TE）

对薄翼小攻角：$\theta_{\mathrm{wake}} \approx \alpha + \epsilon_i$，偏折量级 ~1–3°。

**Drela 2013 的处理**：面板法中尾流面几何是预先给定的（$D_{ij}$ 只依赖面板+尾流几何，§IV.C Eq.(79)），论文中尾流面是直的。Drela 承认这是一个近似——$D_{ij}$ 只算一次，不做迭代松弛。IBL3 方程在给定尾流面几何上求解，局部基跟随尾流面方向。如果尾流面是直的（沿来流），但实际流动有下洗偏折，IBL3 解的是"几何上略有偏差的剪切层"。

**López 2021 的处理**（§2.4, Eq.(2.53)）：显式选择直线/平面尾流沿来流方向：

> "The wake is modeled as a straight surface in the freestream direction. This assumption neglects the roll-up and downwash effects but avoids iteratively computing the wake's geometry."

松弛手段：允许质量穿透尾流面——

$$\hat{n} \cdot (\rho_u \mathbf{u}_u - \rho_l \mathbf{u}_l) = 0 \quad \text{across the wake}$$

即尾流不要求无穿透，只要求上下质量通量连续 + 压力相等。这在势流层面部分补偿了尾流几何不准的问题。López 同时承认（§4.1）："it has been shown that the wake's geometry affects the solution"。

#### 4.5.5 误差量级评估

尾流偏折对势流解的影响是 $O(\theta_{\mathrm{wake}}^2)$ 量级——López §4.1 的 $C_L$ 误差 ~0.17% 印证了这个量级。

对 IBL3 在尾流上的影响：
- $\delta^*_{\mathrm{wake}}$ 沿尾流衰减率依赖流向梯度 $d/dx$——如果 $x$ 方向偏了 $\theta_{\mathrm{wake}}$，梯度计算有 $O(\theta_{\mathrm{wake}}^2)$ 误差
- 对工程精度（$C_L$ 误差 < 1%），直尾流 + 质量穿透松弛是可接受的

**"拐折"三重叠加的完整图景**：

```
           实际流动方向（受环量诱导偏折）
            ↗  θ_wake = α + ε_i（下洗偏折）
           ↗
TE ●────────  尾流面（López: 直线，沿来流 α）
 ╲            ↕ 几何折转角 = θ_TE_upper vs α
  ╲          ↕ 几何折转角 = θ_TE_lower vs α
   ╲
    上下翼面在 TE 的楔角 2ε
```

三个效应的处理方式汇总：

| 效应 | Drela IBL3 的处理 | López 框架的影响 |
|------|-------------------|-----------------|
| TE 局部基折转 | 局部基逐点旋转，方程不变 | 无影响——尾流面网格在 TE 处自然有不同的切向 |
| TSL 近似在 TE | $\kappa(y_e-y_w) \ll 1$ 局部失效，Drela 不做特殊处理 | 与 Drela 一致 |
| 上下 BL 合并 | 6 方程自然过渡，$\delta_{\mathrm{wake}}(0) = \delta_{\mathrm{up}} + \delta_{\mathrm{low}}$ | master-slave TE 节点复制为合并提供天然框架 |
| **尾流几何偏折** | IBL3 在给定几何上求解——几何不对则解有偏 | **直尾流 ≠ 实际偏折尾流——这是主要误差源** |
| 攻角变化时尾流方向 | 面板法需重算 $D_{ij}$ | 嵌入式尾流不需要重构网格——$\alpha$ 变化只改来流方向 |

#### 4.5.6 UP3D 的处理方案：直尾流 + 质量穿透松弛，不做几何松弛

**方案：直尾流 + 质量穿透松弛（López 路线），不实施尾流几何松弛。**

- 尾流面沿来流方向，TE 处局部基折转由 Drela 框架自动处理
- 上下 BL 在 TE 的合并由 master-slave 节点结构支持（方案 B 中由 CutElementMap + enrichment DOF 支持）
- 闭合关系从 wall → wake 切换
- 势流侧：尾流面允许质量穿透（López Eq.(2.53)），不要求无穿透
- IBL3 在直尾流面几何上求解，$\delta^*_{\mathrm{wake}}$ 沿来流方向积分
- 精度：$O(\theta_{\mathrm{wake}}^2) \sim O(1\%)$，对工程验证可接受

对 UP3D 的典型验证工况（NACA 0012 $\alpha \leq 5°$、ONERA M6 $\alpha \leq 4°$），$\theta_{\mathrm{wake}} < 3°$，直尾流的精度已经足够。

**为什么不做尾流几何松弛——即使方案 B 的 level-set 框架在技术上支持它**

方案 B（`20260707_1505_levelset_wake_design.md`）的 `WakeLevelSet.update_direction()` 天然支持尾流几何变化——不动网格，只改 level-set 函数，对应 §8 Phase B6（曲线尾流/自由尾迹）。VSPAERO VLM 等 Lagrangian 面板法也支持尾流自适应。但 UP3D 决定不实施尾流几何松弛，原因如下：

1. **精度收益微小**：尾流几何偏折误差是 $O(\theta_{\mathrm{wake}}^2) \sim O(0.1\%)$。López §4.1 验证表明直尾流 + 质量穿透 vs XFOIL（隐式精确尾流）的 $C_L$ 误差 ~0.17%，远小于其他误差源（可压缩性、BL 耦合缺失、人工密度）。对 UP3D 的目标工况（小到中等攻角），这个误差量级不构成瓶颈。

2. **工程代价不成比例**：即使 level-set 改方向不动网格，每次方向变化仍需重建 `CutElementMap`、重新分配 enrichment DOF、重新组装约束后系统、重新识别 small-cut 配置。这是一套非平凡的预处理流程，需要额外的外层迭代循环控制、收敛判据和鲁棒性处理。

3. **与 Newton 紧耦合的结构性冲突**：尾流方向 $\theta_{\mathrm{wake}}$ 变化时，被切单元集离散跳变（某个单元从"不被切"变为"被切"），这是 **离散不连续性**，不能直接放进 Newton Jacobian。只能用外层固定点迭代 × 内层 Newton 的嵌套结构，增加实现复杂度。

4. **López 的先例**：López 2021 在 level-set 嵌入式尾流框架（正是方案 B 的直接前作）中明确选择直尾流 + 质量穿透松弛，并在 §5.2 Outlook 中未将尾流几何松弛列为未来方向。这表明在该精度等级下，直尾流是合理的工程选择。

5. **优先级排序**：UP3D 的核心目标是全速势方程 + IBL 粘性耦合的紧耦合 Newton 求解。尾流几何松弛是"锦上添花"而非"必需"——在 V1 松耦合和 V3 紧耦合验证完成之前，投入尾流几何松弛会分散精力。

## 5. 是否需要自适应更新网格？
**取决于方法选择**：
| 方法                                                         | 要动网格吗？ | 说明                                  |
| ------------------------------------------------------------ | ------------ | ------------------------------------- |
| 等效物体（几何修正）                                         | **要**       | 每次迭代变形网格                      |
| Transpiration BC                                             | **不要**     | 只改边界条件                          |
| Level-set 等效物体                                           | **不要**     | 更新 level-set 函数（需要嵌入式几何） |
| UP3D 是 body-fitted 网格，如果用 Transpiration BC，**完全不需要动网格**。VII 迭代只在边界条件层面做。 |              |                                       |

## 6. Level-set 能否用于 IBL 粘性修正？
### 6.1 上一轮讨论的修正
上一轮讨论中提出了"用 level-set 表示等效物面（$\varphi_{\mathrm{eff}} = \varphi_{\mathrm{geom}} - \delta^*$）来隐式施加 $\delta^*$"的方案。这个方案在技术上是成立的，但有一个前提条件当时没有说清楚：
**Level-set 等效物体方案需要几何本身也是 level-set 表示的。**
UP3D 采用 body-fitted 几何，物面不靠 level-set 表示。因此"用 level-set 做物面粘性修正"这个方案**不适用于 UP3D 当前架构**。它只适用于像 Núñez 2022 那样的全嵌入式几何框架（$\varphi$ 表示真实物面，理论上可以改为 $\varphi_{\mathrm{eff}} = \varphi - \delta^*$）。
### 6.2 UP3D 中 level-set 的角色
在 UP3D 的架构中，level-set 的职责是明确的：
| 对象                                                         | 表示方式                                | 是否用 level-set |
| ------------------------------------------------------------ | --------------------------------------- | ---------------- |
| **物面**                                                     | body-fitted 网格面                      | ❌                |
| **尾流面**                                                   | level-set $\varphi_{\mathrm{wake}} = 0$ | ✅                |
| 这个分工是设计决策，不是临时方案——body-fitted 几何在物面精度、壁面 Cp 提取、边界层后处理方面有天然优势，不需要为了 IBL 耦合而改变。 |                                         |                  |
### 6.3 Level-set 用于尾流面的 IBL 修正
虽然物面不适用 level-set 修正，但**尾流面的 IBL 修正可以在 level-set 框架内实现**。
在方案 B 的多值 FE 中，尾流面已经有切单元装配机制和最小二乘泛函。IBL 尾迹修正就是在最小二乘泛函的约束中增加一个位移厚度贡献的源项（见 §4.3 路径 ①）：
| 约束                                                         | 公式                                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| 原始                                                         | $\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^- = 0$ (Núñez Eq.30) |
| IBL 修正                                                     | $\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^- = \mathbf{m}_{\mathrm{wake}}(s)$（引入源项） |
| 这个实现不需要用 level-set 来表示 $\delta^*_{\mathrm{wake}}$ 的几何——$\delta^*_{\mathrm{wake}}$ 是一个沿尾流面坐标 $s$ 的标量函数，直接作为已知源项注入最小二乘泛函即可。level-set $\varphi_{\mathrm{wake}}$ 的作用是定义尾流面的位置和切单元识别，与 IBL 修正正交。 |                                                              |

## 7. UP3D 的 IBL 耦合方案选择
### 7.1 方案对比
|                                    | 物面修正                                                     | 尾流修正                                                     | 需要动网格？ | 需要嵌入式几何？        |
| ---------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------ | ----------------------- |
| **方案 VI-1：纯 Transpiration**    | Transpiration BC ($g \neq 0$ in Eq.5)                        | 最小二乘泛函增加 $\mathbf{m}_{\mathrm{wake}}$ 源 (Eq.30 修改) | ❌            | ❌                       |
| **方案 VI-2：Level-set 等效物体**  | $\varphi_{\mathrm{eff}} = \varphi_{\mathrm{geom}} - \delta^*$ | 最小二乘泛函增加 $\mathbf{m}_{\mathrm{wake}}$ 源             | ❌            | ✅（需 Núñez 2022 架构） |
| **方案 VI-3：等效物体 + 网格变形** | 几何偏移 + remesh                                            | 几何偏移                                                     | ✅            | ❌                       |
### 7.2 推荐：方案 VI-1（纯 Transpiration）
**理由**：
1. **与 UP3D 架构一致**：UP3D 是 body-fitted 网格（Davari 2019 架构），物面 BC 已经在网格面上施加。Transpiration 只是把齐次 Neumann ($g=0$) 改成非齐次 ($g=\dot{m}$)，不引入新的几何表示。
2. **不需要动网格**：VII 迭代只更新边界条件值，网格不变。
3. **不需要嵌入式几何**：不需要把物面转成 level-set 表示。
4. **实现量最小**：壁面 Neumann 装配增加源项 + 尾流面最小二乘泛函修改约束。
5. **精度足够**：对薄翼 $\delta^*/c \sim 0.01\text{--}0.03$，Transpiration 与等效物体一阶等价。大曲率区域（前缘）Transpiration 更准确。
6. **工业代码标准做法**：FLO 系列、XFOIL 等都用 Transpiration。
7. **文献支持**：Davari 2019 §3.2 明确提到全速势 + boundary layer solver 耦合是改进精度的方向，其 body-fitted 架构与 Transpiration BC 天然兼容。
### 7.3 方案 VI-2 的定位
Level-set 等效物体方案是一个**优雅的理论方案**，但它需要 UP3D 先转向全嵌入式几何（像 Núñez 2022 那样，$\varphi$ 表示真实物面）。这在当前路线图中没有规划。如果未来 UP3D 考虑全嵌入式几何（比如为了处理复杂几何或多物体），那么 VI-2 可以作为那时的 IBL 耦合方案——届时只需把 $\varphi$ 替换为 $\varphi_{\mathrm{eff}} = \varphi - \delta^*$。
**当前不采用 VI-2。**
### 7.4 方案 VI-3 的定位
等效物体 + 网格变形是传统做法，但在 body-fitted FEM 框架中代价高——每次 VII 迭代都要变形网格、重算 $B_e$ 和 $V_e$。相比 Transpiration（只改 BC 值），代价不成比例。**不采用。**

## 8. VII 耦合策略：松→紧渐进路线

### 8.1 路线总览

基于 BLWF58/TRANAIR/Drela 的对比分析（见 `20260709_0145_3d_vii_implementation_analysis.md` §15），采用三阶段渐进路线：

| 阶段 | 耦合方式 | 依赖 | 收敛阶 | 典型迭代数 |
|------|---------|------|--------|-----------|
| **V1 松耦合** | 固定点迭代（势流→IBL→transpiration→重解势流） | P6 | 线性 | 5-10（附着流） |
| **V2 quasi-simultaneous**（可选） | Hilbert 积分预估势流响应 | V1 | 线性（更快） | 3-5 |
| **V3 紧耦合 Newton** | 增广 $(φ, Γ, δ, A, B, Ψ, C_{τ1}, C_{τ2})$ 同时求解 | P8 | **二次** | 5-10 |

**TRANAIR 的教训**（Boeing 2003）：从 indirect coupling 升级到 direct coupling 是为了收敛可靠性。松耦合在接近分离时不够——V3 是最终目标。

**Drela 2013 的论断**（§I.A）："simultaneous solution is highly desirable, if not essential"。

### 8.2 V1 松耦合迭代循环

```text
VII 松耦合循环：
1. 势流求解（当前 φ, master-slave 尾流 + Kutta）
   → 得到壁面 q_e(s), ρ_e(s) 和尾流面 q_e,wake(s), ρ_e,wake(s)

2. IBL3 求解（Drela 2013 §III 表面 Galerkin FE）
   → 物面：δ*(s), θ(s), H(s), Cf(s), A, B, Ψ, C_τ1, C_τ2
   → 尾迹：δ*_wake(s) 从 TE 向下游

3. 更新边界条件（只改 RHS，矩阵不变）：
   壁面 transpiration：ṁ_wall(s) = d/ds [ρ_e U_e δ*_wall(s)]
   尾流面源项：       ṁ_wake(s) = d/ds [ρ_e U_e δ*_wake(s)]

4. 回到 1，直到 δ* 收敛
```

**势流侧的改动**（松耦合只改 RHS）：

| 模块 | 改动 | Drela 2013 依据 |
|------|------|----------------|
| `kernels/residual.py` | 壁面 Neumann 项增加 transpiration 源 | Eq.(76) 前两项 |
| `constraints/wake.py` | RHS 增加尾流面 δ*_wake 源项 | $T^T b - \sum_j \Gamma_j h_j + T^T b_{\text{wake}}(\delta^*_{\text{wake}})$ |
| `physics/ibl3.py`（新建） | Drela 6 方程表面 FE 求解器 | §II Eq.(21)-(29), §III |

### 8.3 V3 紧耦合 Newton

增广未知量 $\mathbf{x} = (\varphi, \Gamma, \delta, A, B, \Psi, C_{\tau1}, C_{\tau2})$，增广残差：

$$\mathbf{R} = \begin{pmatrix} R_\varphi \\ R_\Gamma \\ R_\delta \\ R_A \\ R_B \\ R_\Psi \\ R_{C\tau1} \\ R_{C\tau2} \end{pmatrix}, \qquad \mathbf{J} = \begin{pmatrix} J_{\varphi\varphi} & J_{\varphi\Gamma} & J_{\varphi,\delta^*} & \cdots \\ J_{\Gamma\varphi} & J_{\Gamma\Gamma} & 0 & \cdots \\ J_{\delta^*,\varphi} & 0 & J_{\delta^*\delta^*} & \cdots \\ \vdots & & & \ddots \end{pmatrix}$$

其中：
- $J_{\varphi,\delta^*} = \partial \dot{m}_{\text{wall}} / \partial \delta^*$ — transpiration 对势流的耦合（Drela Eq.(76) 前两项的 Jacobian）
- $J_{\delta^*,\varphi} = \partial R_\delta / \partial \varphi$ — 势流速度对 BL 的耦合（通过 $q_e$）
- 对角块 $J_{\delta^*\delta^*}$ 等 — IBL3 方程自身的 Jacobian

**GMRES + 块预处理**：势流块用 AMG，BL 块用 ILU——两个块的预处理特性不同。

**与 P8 Newton 的关系**：V3 是 P8 增广 $(φ, Γ)$ 系统的进一步扩展——加入 6 个 BL 变量。P8 的 GMRES+AMG 框架可以复用，只需扩展增广 Jacobian 的维度。

### 8.4 收敛判据

**V1 松耦合**：
$$\frac{\|\delta^{*,n+1} - \delta^{*,n}\|_2}{\|\delta^{*,n+1}\|_2} < \varepsilon_{\mathrm{vii}}, \qquad \varepsilon_{\mathrm{vii}} \sim 10^{-4}$$

**V3 紧耦合**：Newton 残差 $\|\mathbf{R}\|_2 < 10^{-10}$，二次收敛验证：$\|\mathbf{R}_{k+1}\| / \|\mathbf{R}_k\|^2 \to \text{常数}$。

## 9. VII 的适用范围与限制
### 9.1 适用
- 亚声速附着流（无分离）
- 跨声速弱激波（激波前附着，激波后可能轻微分离但 IBL 仍可处理）
- 薄翼、小到中等攻角
### 9.2 不适用
- 大攻角分离流——IBL 在分离点附近失效
- 强激波诱导分离——$\delta^*$ 急剧增长，IBL 不可靠
- 钝体绕流——边界层假设不成立
### 9.3 与 UP3D 验证矩阵的关系
| 验证项        | 是否需要 IBL | 说明                 |
| ------------- | ------------ | -------------------- |
| V0 自由流     | ❌            | 纯无粘               |
| V1 MMS        | ❌            | 纯无粘               |
| V3 亚声速升力 | 可选         | IBL 修正量很小       |
| V4 跨声速     | 推荐         | IBL 对激波位置有修正 |
| V5 ONERA M6   | 推荐         | 工业级验证           |

## 10. 实现阶段

### Phase V1：IBL3 求解器 + 松耦合（独立验证）

**目标**：实现 Drela 6 方程 IBL3 表面 FE 求解器，松耦合到势流。

**交付物**：
- `viscous/ibl3.py`：Drela 2013 §II-III 的 6 方程表面 Galerkin FE 求解器
  - 局部 Cartesian 基 + P1 三角形元（Drela 用 bilinear 四边形，UP3D 适配三角形）
  - 6 个未知量 $(\delta, A, B, \Psi, C_{\tau1}, C_{\tau2})$ per node
  - 闭合关系：streamwise + crossflow profile shapes（Drela 2013 §II.C）
  - 转捩：$C_{\tau}$ 方程的 $e^N$ 形式
- `viscous/transpiration.py`：$\delta^* \to \dot{m}$ 转换（Drela Eq.(76)）
- `viscous/coupling.py`：VII 松耦合迭代控制器
- 修改 `kernels/residual.py`：壁面 Neumann 增加 transpiration 源
- 修改 `constraints/wake.py`：RHS 增加尾流面 $\delta^*_{\text{wake}}$ 源项

**验证标准**：
- 2D NACA 0012 亚声速：$\delta^*$ 与 XFOIL 对比
- VII 松耦合 5-10 次迭代收敛（附着流）
- 势流矩阵不变（只改 RHS），G4.2 bitwise no-op 当 $\delta^*=0$

**前置条件**：P6（smooth wall gradient 作为 IBL 输入）

### Phase V2：quasi-simultaneous coupling（可选中间步骤）

**目标**：加入 Hilbert 积分预估，提高松耦合收敛速度。

**交付物**：
- `viscous/hilbert.py`：Hilbert 积分算子（BLWF58 §3.2 路线）
- 在非结构网格上实现——需要在表面网格的局部邻域上积分

**验证标准**：
- 松耦合迭代数减少 30-50%
- 接近分离时收敛性改善

**前置条件**：V1

> **V2 是可选的**。如果 V1 已经足够快（5-10 次迭代），可以跳过 V2 直接做 V3。BLWF58 用 quasi-simultaneous 是因为它不做全局 Newton——UP3D 的 P8 Newton 路线使得 V3 直接可做。

### Phase V3：紧耦合 Newton

**目标**：增广 $(φ, Γ, \delta, A, B, \Psi, C_{\tau1}, C_{\tau2})$ 同时 Newton 求解。

**交付物**：
- `viscous/jacobian.py`：IBL3 方程的 Jacobian（6×6 block per node）
- 增广 Newton 驱动器：扩展 P8 的 $(φ, Γ)$ 增广到 $(φ, Γ, \text{BL})$
- GMRES + 块预处理（势流 AMG + BL ILU）

**验证标准**：
- Newton 二次收敛：$\|\mathbf{R}_{k+1}\| / \|\mathbf{R}_k\|^2 \to$ 常数
- 2D NACA 0012 跨声速：VII 修正后激波位置偏移与实验对比
- 3D ONERA M6：$C_L$ 从无粘 0.245 向实验值 0.288 靠近

**前置条件**：P8（Newton Jacobian 框架）+ V1

### Phase V4：尾流面 IBL 修正（V1 的自然延续，非独立阶段）

**目标**：将 IBL3 方程从物面延续到尾流面，输出 $\delta^*_{\text{wake}}$ 作为势流尾流面 RHS 源项。

**关键认识**：这不是独立于 V1 的工作。Drela IBL3 的 6 方程在物面和尾流面上是**同一套方程**，区别只在闭合关系（§2.4）。V1 实现时就应预留尾流面网格的 6 个未知量，V4 只是切换闭合关系。

**交付物**：
- 尾迹闭合关系：wake profile（中心对称，无壁面），替代物面的 wall profile
- `constraints/wake.py` 修改：RHS 增加 $T^T b_{\text{wake}}(\delta^*_{\text{wake}})$
- TE 过渡：6 方程在 TE 节点处自动从壁面闭合切换到尾迹闭合——无需特殊处理
  - UP3D 的 master-slave TE 节点复制为 TE 处上下面 BL 汇合提供自然框架
  - 与 BLWF58 §3.3.2 的 TE 特殊处理不同——Drela 的 DCV 框架不需要 TE 切换逻辑

**验证标准**：
- NACA 0012 跨声速：尾流修正后 $C_D$ 改善
- 尾流面 $\Delta\varphi$ 沿流向保持常数（V6 一致性）
- $\delta^*_{\text{wake}}$ 从 TE 处的 $\delta^*_{\text{wall,upper}} + \delta^*_{\text{wall,lower}}$ 向下游衰减

**前置条件**：V1（同一求解器，切换闭合关系）

## 11. 参考文献
### IBL3 理论核心（已下载）
- **Drela, M. 2013. "Three-Dimensional Integral Boundary Layer Formulation for General Configurations." AIAA Paper 2013-2437.** — `references/Drela_2013_IBL3_general_configurations.pdf`
  - §II Eq.(21)-(29): 6 方程 IBL3 完整推导（δ, A, B, Ψ, C_τ1, C_τ2）
  - §II.H: DCV/DDCV 框架（尾迹方程）
  - §II.C: 闭合关系（streamwise + crossflow profile shapes）
  - §III Eq.(76)-(77): 壁面 transpiration BC 的 Galerkin 弱形式
  - §III.1-3: 表面 FE 方法（局部 Cartesian 基 + bilinear 元）
  - §IV: 测试案例（torpedo body、低展弦比机翼带分离）
  - **本文档的 IBL 方程、表面 FE 方法、transpiration BC 公式均来源于此**

### BLWF/TRANAIR 对比参考（已下载）
- **Karas, O.V., Kovalev, V.E. *BLWF58 Computational Method and Algorithms.* TsAGI, 108pp.** — `references/BLWF58_method_algorithms.doc`
  - §3.1-3.2: 四种耦合方案对比（direct/semi-inverse/quasi-simultaneous/inverse）
  - §3.3.2: TE 处 BL→wake 过渡处理
  - quasi-simultaneous coupling（Hilbert 积分）——V2 可选中间步骤的参考
- **BLWF58 User's Guide.** TsAGI, 166pp. — `references/BLWF58_user_guide.doc`
- **Zhang, K., Hepperle, M. 2010. *Evaluation of the BLWF Code.* DLR IB 124-2010/3.** — `references/DLR_2010_BLWF_evaluation.pdf`
- **Johnson, F.T. et al. 2003. "Thirty Years of Development and Application of CFD at Boeing." AIAA 2003-3439.** — `references/Boeing_2003_AIAA_30yr_CFD.pdf`
  - TRANAIR 从 indirect→direct coupling 的演进——验证 V3 紧耦合是最终目标

### 势流 + 嵌入尾流框架（已下载）
- **Davari, M. et al.** "A cut finite element method for the solution of the full-potential equation with an embedded wake." *Comput. Mech.* 63:821–833, 2019.
  - §2: body-fitted 几何 + level-set 尾流 + Cut FE 方案
  - Eq.5/10: 壁面 Neumann BC，Transpiration 修正的起点
  - §3.2 + [33]: 明确提出与 boundary layer solver 耦合的方向
- **Núñez, M. et al.** "An embedded approach for the solution of the full potential equation with finite elements." *CMAME* 388:114244, 2022.
  - Eq.30–31: 矢量尾流约束 + 最小二乘泛函（尾流面 IBL 修正的框架）

### VII 方法论文献
- Nishida, B., Drela, M. 1995. "Fully simultaneous coupling for three-dimensional viscous/inviscid flows." AIAA 1995-1806.
  - 3D simultaneous Newton coupling；López 引用 [23]
- Lock, R.C. "The determination of the viscous correction for the lift of a thin wing." RAE TR 6709 (1967).
- Rodríguez, D. et al. 2012. "A rapid, robust, and accurate coupled boundary-layer method for CART3D." AIAA 2012-302.

## 附录 A：代码适配确认（2026-07-09，HEAD 170eaa9）

| 检查项 | 状态 | 说明 |
|--------|------|------|
| VII 迭代接口 | ✅ 兼容 | VII 在势流求解外层，P6/P7/P8 的改动不影响 VII 接口 |
| Transpiration BC 实现 | ✅ 兼容 | 壁面 Neumann 源项只进 RHS；V1 松耦合不改矩阵 |
| 跨声速稳定性 | ✅ 正面 | P6 smooth wall gradient 使 IBL 输入更稳定 |
| `constraints/wake.py` $h_j$ 批处理 | ✅ 兼容 | 尾流面 $\delta^*_{\text{wake}}$ 源项加到 RHS，与 $h_j$ 无关 |
| P6 smooth wall gradient | ✅ 前置 | IBL3 的壁面速度输入需要光滑梯度——P6 已提供 |
| P8 Newton 框架 | ✅ V3 前置 | 增广 $(φ, Γ, \text{BL})$ 扩展 P8 的 $(φ, Γ)$ 框架 |
| master-slave TE 节点复制 | ✅ TE 耦合 | V4 尾流面 IBL 的 TE 过渡可利用 TE 节点复制 |

## 附录 B：与 `20260709_0145_3d_vii_implementation_analysis.md` 的关系

本文档是 IBL 耦合方案的**设计文档**（做什么、怎么做）。
`20260709_0145_3d_vii_implementation_analysis.md` 是**分析报告**（为什么、与谁对比、风险在哪），其中 §15 包含 BLWF58/TRANAIR/Drela 的详细对比。

两份文档共同构成 UP3D VII 的完整技术路线：
- 本文档 §2（IBL3 方程选择）← 分析报告 §15.4-15.5（对比和修正）
- 本文档 §8（耦合策略）← 分析报告 §15.2-15.3（BLWF/TRANAIR 经验）
- 本文档 §10（实现阶段）← 分析报告 §14.2（依赖关系）
