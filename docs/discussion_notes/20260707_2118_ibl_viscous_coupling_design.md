# IBL 粘性修正耦合方案设计
> Design note for UP3D / pyFP3D
> Date: 2026-07-07 (updated 2026-07-07 23:22 — baseline 同步 e3d0386)
> Code baseline: HEAD e3d0386
> Prerequisite: `20260707_1505_levelset_wake_design.md`（方案 B：Level-Set 尾流 + 多值 FE）
> References: Núñez et al. 2022 (CMAME 114244); Davari et al. 2019 (Comput. Mech. 63:821–833)

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

## 2. IBL 的输出
IBL（积分边界层法）在物面上求解动量积分方程，核心输出沿物面坐标 $s$：
- $\delta^*(s)$ —— 位移厚度（最关键）
- $\theta(s)$ —— 动量厚度
- $H = \delta^*/\theta$ —— 形状因子
- $C_f(s)$ —— 壁面摩擦系数
- 分离点位置（$H \to \sim 3.5$ 时）
同时，IBL 有**尾迹方程**——从后缘出发继续积分，给出 $\delta^*_{\mathrm{wake}}(s)$ 沿尾迹的演化。
因此 IBL 提供两套位移厚度：
- $\delta^*_{\mathrm{wall}}(s)$：物面上的位移厚度
- $\delta^*_{\mathrm{wake}}(s)$：尾迹中的位移厚度，从 TE 向下游衰减

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

## 8. 方案 VI-1 的 VII 迭代循环
```text
VII 循环：
1. 势流求解（当前 φ_wake level-set + 多值 FE + 罚函数 Kutta）
   → 得到壁面 U_e(s), ρ_e(s) 和尾流面 U_e,wake(s), ρ_e,wake(s)
2. IBL 求解
   → 物面积分：δ*_wall(s), θ(s), H(s), Cf(s)
   → 尾迹积分：δ*_wake(s) 从 TE 向下游
3. 更新边界条件：
   壁面 Transpiration：m_dot_wall(s) = d/ds [ρ_e U_e δ*_wall(s)]
   尾流面源项： m_dot_wake(s) = d/ds [ρ_e U_e δ*_wake(s)]
4. 回到 1，直到 δ* 收敛
```
对应的数学表达：
**势流求解输出**：
- 壁面：$U_e(s), \rho_e(s)$
- 尾流面：$U_{e,\mathrm{wake}}(s), \rho_{e,\mathrm{wake}}(s)$
**IBL 求解输出**：
- 物面积分：$\delta^*_{\mathrm{wall}}(s), \theta(s), H(s), C_f(s)$
- 尾迹积分：$\delta^*_{\mathrm{wake}}(s)$ 从 TE 向下游
**更新边界条件**：
壁面 Transpiration：
$$\dot{m}_{\mathrm{wall}}(s) = \frac{d}{ds} \left[ \rho_e U_e \delta^*_{\mathrm{wall}}(s) \right]$$
尾流面源项：
$$\dot{m}_{\mathrm{wake}}(s) = \frac{d}{ds} \left[ \rho_e U_e \delta^*_{\mathrm{wake}}(s) \right]$$
**收敛判据**：$\delta^*$ 的相对变化
### 8.1 势流侧的改动
| 模块                                     | 改动                                                 | 文献依据                                                     |
| ---------------------------------------- | ---------------------------------------------------- | ------------------------------------------------------------ |
| `kernels/residual.py`                    | 壁面 Neumann 项增加 transpiration 源 ($g = \dot{m}$) | Davari Eq.5/10, $g$ 从 0 变为 $\dot{m}$                      |
| `kernels/cut_assembly.py`（方案 B 新增） | 最小二乘泛函约束增加 $\mathbf{m}_{\mathrm{wake}}$ 源 | Núñez Eq.30 修改：$\rho^+ \nabla \Phi^+ - \rho^- \nabla \Phi^- = \mathbf{m}_{\mathrm{wake}}$ |
| `physics/ibl.py`（新增）                 | IBL 求解器：可压缩动量积分方程                       | Davari §3.2, 引用 Rodriguez 2012 [33]                        |
### 8.2 IBL 侧的模块
```text
pyfp3d/viscous/
├── __init__.py
├── ibl.py              # 积分边界层求解器（物面 + 尾迹）
├── transpiration.py    # 将 δ* 转换为壁面/尾流面 Neumann 源项
└── coupling.py         # VII 迭代控制器
```
### 8.3 收敛判据
VII 迭代的收敛判据是 $\delta^*$ 的相对变化：
$$\frac{\|\delta^{*,n+1} - \delta^{*,n}\|_2}{\|\delta^{*,n+1}\|_2} < \varepsilon_{\mathrm{vii}}$$
典型 $\varepsilon_{\mathrm{vii}} \sim 10^{-4}$，迭代 3–5 次收敛（弱耦合）。

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

## 10. 实现阶段（与方案 B 的阶段对齐）
### Phase V1：IBL 求解器（独立验证）
**目标**：实现二维 IBL 求解器，在给定 $U_e(s)$ 下求解 $\delta^*(s)$。
**交付物**：
- `viscous/ibl.py`：可压缩动量积分方程（Green's lag entrainment 或 Thwaites + 可压缩修正）
- 测试：与 XFOIL 的 $\delta^*$ 结果对比
**前置条件**：无（可与方案 B 并行开发）
### Phase V2：壁面 Transpiration BC
**目标**：在 body-fitted 网格上实现壁面 transpiration BC。
**交付物**：
- `viscous/transpiration.py`：$\delta^* \to \dot{m}$ 转换
- 修改 `kernels/residual.py`：壁面 Neumann 增加源项
- 测试：NACA0012 无升力 transpiration ON，验证 $c_l \approx 0$（$\delta^*$ 小）
**前置条件**：Phase V1
### Phase V3：VII 迭代
**目标**：实现完整的 VII 循环。
**交付物**：
- `viscous/coupling.py`：VII 迭代控制器
- 测试：NACA0012 亚声速，VII 收敛（3–5 次迭代）
**前置条件**：Phase V2
### Phase V4：尾流面 IBL 修正
**目标**：在尾流面切单元装配中增加 $\delta^*_{\mathrm{wake}}$ 源项。
**交付物**：
- 修改 `kernels/cut_assembly.py`：切单元残差增加尾流面源项
- `viscous/ibl.py` 增加尾迹积分
- 测试：NACA0012 跨声速，VII + 尾流修正后激波位置偏移
**前置条件**：Phase V3 + 方案 B 的 Phase B3

## 11. 参考文献
### 核心文献（已下载）
- **Davari, M., Rossi, R., Dadvand, P., López, I., Wüchner, R.** "A cut finite element method for the solution of the full-potential equation with an embedded wake." *Computational Mechanics* 63(5):821–833, 2019. DOI: 10.1007/s00466-018-1624-3.
  - §2: body-fitted 几何 + level-set 尾流 + Cut FE 方案
  - Eq.5/10: 壁面 Neumann BC ($g=0$)，Transpiration 修正的起点
  - Eq.6–7: 尾流面质量通量连续 + 压力相等
  - Eq.23–27: 压力相等线性化 + 最小二乘泛函
  - Eq.28–32: Heaviside 扩展形函数 + DOF 复制
  - §3.2 + [33]: 明确提出与 boundary layer solver 耦合的方向
  - Appendix A: 边界层理论证明压力沿法向不变
- **Núñez, M., López, I., Baiges, J., Rossi, R.** "An embedded approach for the solution of the full potential equation with finite elements." *Computer Methods in Applied Mechanics and Engineering* 388:114244, 2022.
  - §2.3: 区分 body-fitted (Davari [8]) 与全嵌入式
  - Eq.24–25: 切单元无穿透条件
  - Eq.30–31: 矢量尾流约束 + 最小二乘泛函
  - Eq.34–39: 约束残差和 Jacobian
  - Eq.40–43: 切单元矩阵组装
  - Algorithm 1: $\Omega_{\mathrm{kutta}}$ 识别（仅全嵌入式需要）
### VII 方法论文献
- Lock, R.C. "The determination of the viscous correction for the lift of a thin wing." RAE TR 6709 (1967). — 经典 VII + transpiration 理论
- Yoshihara, H., Spee, B.M. "Viscous-inviscid interactions." AGARD AR-138 (1979). — VII 综述
- Holst, T.L. "Viscous transonic airfoil workshop compendium of results." AIAA J. 1987. — 跨声速 VII 验证
- Rodriguez, D., Sturdza, P., Suzuki, Y., Martins-Rivas, H., Peronto, A. "A rapid, robust, and accurate coupled boundary-layer method for CART3D." AIAA 2012-302. — Davari 2019 [33]，全速势 + IBL 耦合的工程实例

## 附录：代码适配确认（2026-07-07 22:45，HEAD e0353d9）
P4 修复（local $\theta \cdot \mathrm{diag}(A_{\mathrm{free}})$ 阻尼 + $h_j$ 批处理化）对本方案的影响：
| 检查项                                         | 影响   | 说明                                                         |
| ---------------------------------------------- | ------ | ------------------------------------------------------------ |
| VII 迭代接口                                   | 无     | VII 在势流求解外层，`damping_theta` 不影响接口               |
| Transpiration BC 实现                          | 无     | 壁面 Neumann 源项与阻尼机制独立                              |
| 跨声速稳定性                                   | ✅ 正面 | local damping 使 transonic 更稳定，VII 在 transonic 场景下更容易收敛 |
| `constraints/wake.py` $h_j$ 批处理化           | 无     | V1（IBL 求解器）独立于尾流表示；V4（尾流面 IBL 修正）依赖方案 B 的 B3，与 $h_j$ 无关 |
| **结论：本文档无需修改，与当前代码完全兼容。** |        |                                                              |
