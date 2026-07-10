# 三维无粘-有粘迭代（VII）实现路径分析

> Design note for UP3D / pyFP3D
> Date: 2026-07-09 01:45
> Code baseline: HEAD 170eaa9
> 前置文档：`20260707_2118_ibl_viscous_coupling_design.md`（IBL 方案 VI-1 选型）
> References:
> - López Canalejo 2021 博士论文 (TUM) — Ch.2 §2.4 气动力分析, §5.2 Outlook
> - Nishida & Drela 1995 (AIAA 1995-1806) — 3D fully simultaneous viscous/inviscid coupling
> - Drela 2014, *Flight Vehicle Aerodynamics* — IBL 理论与 2D VII
> - Davari et al. 2019 (Comput. Mech. 63:821–833) — body-fitted FPE + 嵌入尾流, Appendix A BL 理论
> - UP3D `docs/design.md` §5 (BC), §9 (Post-processing), §11 P10 (VII hook)
> - UP3D `docs/roadmap.md` P10 (backlog: VII transpiration BC)

---

## 目录

1. [问题定位：为什么需要 VII](#1-问题定位)
2. [理论框架：IBL 方程与无粘-有粘耦合](#2-理论框架)
3. [López 论文的定位与局限](#3-lópez-论文的定位与局限)
4. [2D VII 经典方法（XFOIL/MSES 路线）](#4-2d-vii-经典方法)
5. [3D VII 的核心困难](#5-3d-vii-的核心困难)
6. [Nishida-Drela 1995 的 3D simultaneous coupling](#6-nishida-drela-1995)
7. [UP3D 的 3D VII 实现路径](#7-up3d-的-3d-vii-实现路径)
8. [壁面 Transpiration BC 的具体实现](#8-壁面-transpiration-bc)
9. [尾流面 IBL 修正](#9-尾流面-ibl-修正)
10. [IBL 求解器的设计](#10-ibl-求解器的设计)
11. [迭代策略与收敛性](#11-迭代策略与收敛性)
12. [与 P8 Newton 的交互](#12-与-p8-newton-的交互)
13. [验证策略](#13-验证策略)
14. [风险与未解决问题](#14-风险与未解决问题)

---

## 1. 问题定位

### 1.1 纯全速势方程的局限

全速势方程假设等熵、无旋、无粘。对于附着边界层的薄翼绕流，这在天[SYSTEM_NOTE: Content compressed. Read the full version's theory if needed.]r $\theta$ —— 动量厚度
- $H = \delta^*/\theta$ —— 形状因子
- $C_f(s)$ —— 壁面摩擦系数
- 分离点位置（$H \to \sim 3.5$）

在 3D 中，IBL 需要在物面的二维参数空间 $(s_1, s_2)$ 上求解，输出是面场 $\delta^*(s_1, s_2)$ 等。

### 2.2 可压缩 BL 积分方程

对于跨声速流动，BL 方程必须考虑可压缩性。二维可压缩动量积分方程（ equilibrium / Green's lag entrainment 方法）：

$$\frac{d\theta}{ds} + (2+H-M_e^2)\frac{\theta}{U_e}\frac{dU_e}{ds} = \frac{C_f}{2}$$

其中 $M_e$ 是边界层外缘 Mach 数（从势流解读取），$U_e$ 是外缘速度。

3D 形式需要增加横向动量方程（crossflow）：

$$\frac{\partial \theta_{11}}{\partial s_1} + \frac{\partial \theta_{12}}{\partial s_2} + \text{(geometric terms)} = \frac{C_{f1}}{2}$$

$$\frac{\partial \theta_{21}}{\partial s_1} + \frac{\partial \theta_{22}}{\partial s_2} + \text{(geometric terms)} = \frac{C_{f2}}{2}$$

这些方程在物面网格上离散，需要知道物面的曲面度量（第一、第二基本形式）。

### 2.3 边界层外缘量的提取

IBL 方程的输入是边界层外缘的速度分布 $U_e(s)$ 和压力梯度 $dU_e/ds$。在势流框架中：

$$U_e = |\mathbf{u}_\infty + \nabla\varphi|_{\text{wall}}$$

$$\frac{dU_e}{ds} = \frac{\partial U_e}{\partial s} = \hat{s} \cdot \nabla U_e|_{\text{wall}}$$

UP3D 的壁面速度恢复在 `post/surface.py`——P6 已经实现了 `smooth_wall_tangential_gradients`，可以直接给 IBL 提供光滑的壁面速度梯度。这是 VII 与 P6 的一个自然接口。

### 2.4 Transpiration BC 的推导

从边界层连续性方程在薄层近似下积分，壁面处的等效质量通量为：

$$\dot{m}_{\text{wall}}(s) = \frac{d}{ds}[\rho_e U_e \delta^*(s)] \quad \text{(2D)}$$

在 3D 中，沿物面坐标 $(s_1, s_2)$：

$$\dot{m}_{\text{wall}}(s_1, s_2) = \nabla_\Gamma \cdot (\rho_e U_e \delta^* \hat{e}_\text{streamline})$$

其中 $\nabla_\Gamma$ 是曲面散度，$\hat{e}_\text{streamline}$ 是壁面流线方向。这是 UP3D 壁面 Neumann 条件中 $g=0$ 变为 $g \neq 0$ 的来源。

---

## 3. López 论文的定位与局限

### 3.1 López 的粘性修正意识

López 在 §4.1（NACA 0012 不可压缩）明确讨论了粘性效应的缺失：

> *"To account for viscous effects, instead of setting the mass flow at the airfoil's surface to zero, it can be set to a known transpiration value (e.g., to model displacement or wall transpiration effects)."*

他在 lift polar 对比中观察到：大 $\alpha$ 时吸力峰被低估、升力被高估——这正是忽略边界层位移效应的经典症状。

### 3.2 López 的展望

§5.2 Outlook 最后一句话：

> *"Finally, the proposed FPS can be coupled to a boundary layer solver. Although this approach still would not be as general as a RANS solver, it would allow accounting for viscous effects to compute more accurate solutions than the single FPS at a lower cost [23]."*

[23] = Nishida & Drela 1995。

### 3.3 López 没有做的

López **没有实现** VII——他的求解器是纯无粘的。他只是：
- 意识到了粘性修正的需要（§4.1 的升力极曲线偏差）
- 引用了 Nishida-Drela 1995 作为 3D 耦合的先例
- 在 Outlook 中指出耦合 BL solver 是改进方向

所以 López 论文不是 VII 的实现参考——它是 **UP3D 的无粘基线**的参考。VII 的实现需要参考 Nishida-Drela 1995 和 Drela 2014。

### 3.4 López 的壁面 BC 结构对 VII 的支持

López 的壁面 BC 是 slip BC（$q=0$），在 Galerkin 弱形式中是 natural BC（do-nothing）。这恰好是 transpiration BC 最容易嵌入的结构——不需要改矩阵结构，只改 RHS：

$$R_i = \sum_e \int_{\Omega_e} \tilde{\rho} \nabla\varphi \cdot \nabla N_i \, d\Omega - \int_{\Gamma_{\text{wall}}} N_i \cdot \dot{m}(s) \, d\Gamma$$

UP3D 的 `kernels/residual.py` 已经有壁面积分的框架（自然 BC 的 surface term），只是 $g=0$ 时这项为零。VII 只需让 $g \neq 0$。

---

## 4. 2D VII 经典方法

### 4.1 XFOIL 的做法（Drela 1989）

XFOIL 是 2D VII 的黄金标准。耦合策略：

1. **同时求解**（fully simultaneous）：势流 φ 和 BL 变量 $(\delta^*, \theta, H)$ 在同一个 Newton 系统中求解
2. **势流离散**：panel method（高阶面元），不是 FE
3. **BL 方程**：2D 可压缩动量积分 + lag entrainment 方程，沿弧长 $s$ 离散
4. **耦合**：BL 的 $\delta^*$ 通过 transpiration BC 进入势流；势流的 $U_e(s)$ 进入 BL
5. **后缘到尾迹**：BL 在 TE 平滑过渡到 wake 方程，$\delta^*$ 连续

### 4.2 MSES 的做法（Drela-Giles 1990）

MSES 把 XFOIL 的思想推广到多段翼和欧拉方程：

1. **势流→Euler**：用 Euler 方程代替全速势，但 BL 耦合方式不变
2. **同时求解**：Euler + BL 在同一个 Newton 系统中
3. **Newton-Jacobian**：BL 方程的 Jacobian $\partial \delta^*/\partial U_e$ 和势流的 Jacobian $\partial R/\partial \varphi$ 交叉耦合

### 4.3 对 UP3D 的启示

2D 经典方法的核心设计决策：
- **同时求解 > 迭代耦合**：同时求解消除了 lag 误差，收敛更快更鲁棒
- **Newton 统一框架**：势流和 BL 的残差在同一个增广 Jacobian 中
- **Transpiration > 几何修正**：不动机架网格

但这些都是在 2D 中，1D BL 方程沿弧长积分。3D 的情况完全不同。

---

## 5. 3D VII 的核心困难

### 5.1 BL 方程从 1D 变 2D

2D VII 中，BL 方程是沿弧长 $s$ 的 1D ODE——离散后每个壁面节点几个未知量。3D 中：

- BL 方程在物面的 2D 参数空间上求解——是 PDE 不是 ODE
- 需要处理 **crossflow**（横向流动）——后掠翼上前缘附近显著
- 需要物面的 **曲面几何信息**（第一基本形式 $E, F, G$；主曲率方向）
- 附面线（streamline tracing on surface）——BL 方程沿壁面流线积分，需要先追踪流线

### 5.2 壁面流线追踪

3D IBL 的最基础步骤是确定壁面流线方向：

$$\hat{e}_\text{streamline} = \frac{\mathbf{u}_{e,\text{tangent}}}{|\mathbf{u}_{e,\text{tangent}}|}$$

然后沿流线积分 BL 方程。但在分离附近、后掠翼前缘附近，壁面流线方向可能急剧变化。流线追踪在非结构三角形网格上是一个数值敏感操作——需要稳健的 streamline tracing 算法。

### 5.3 分离判定与 VII 的有效域

IBL 方法在附着流中可靠。分离附近（$H \gtrsim 2.6$），动量积分方程的假设开始失效。3D 分离更复杂——不是一条线而是一个区域，且可能有无后缘分离的后掠翼分离泡。

VII 的有效域：
- ✅ 全附着流（design point）
- ✅ 轻度激波-BL 相互作用（shock-induced transition 但未分离）
- ⚠️ 激波后分离（IBL 可预测分离点但不计算分离区）
- ❌ 大范围分离（需要 RANS/DES）

### 5.4 尾流面修正的 3D 复杂性

2D 中，尾迹是从 TE 出发的一条线，$\delta^*_{\text{wake}}(s)$ 是 1D 函数。3D 中：

- 尾流面是 2D 曲面，$\delta^*_{\text{wake}}(s_1, s_2)$ 是面场
- 尾流面的流线方向可能与来流方向不同（翼尖涡卷起）
- 翼尖附近 $\delta^*$ 演化复杂（弦向流→展向流的过渡）

---

## 6. Nishida-Drela 1995 的 3D simultaneous coupling

López 引用的 [23] 是 3D VII 的关键文献。其核心设计：

### 6.1 方法概要

Nishida & Drela 1995 提出了 3D 全同时耦合（fully simultaneous coupling）：

1. **势流**：3D 全速势方程，FE 离散
2. **BL 方程**：3D 积分 BL 方程（含 crossflow），在物面网格上离散
3. **同时 Newton 求解**：增广系统

$$\begin{pmatrix} J_{\varphi\varphi} & J_{\varphi\text{BL}} \\ J_{\text{BL}\varphi} & J_{\text{BL,BL}} \end{pmatrix} \begin{pmatrix} \Delta\varphi \\ \Delta\text{BL} \end{pmatrix} = -\begin{pmatrix} R_\varphi \\ R_\text{BL} \end{pmatrix}$$

4. **BL 变量**：每个壁面节点 $(\delta^*, \theta, H, C_{f1}, C_{f2})$——5 个未知量/节点
5. **耦合项**：
   - $J_{\varphi\text{BL}} = \partial R_\varphi / \partial \text{BL}$：transpiration BC 中 $\dot{m}$ 对 $\delta^*$ 的导数
   - $J_{\text{BL}\varphi} = \partial R_\text{BL} / \partial \varphi$：BL 方程中 $U_e$ 对 $\varphi$ 的导数

### 6.2 关键技术

- **流线追踪**：每个 Newton 迭代重新追踪壁面流线，但流线上的 BL 变量在 Newton 步间连续
- **曲面坐标系**：用物面网格的参数化（$(u,v)$ 坐标），BL 方程在参数空间中离散
- **后缘到尾迹的过渡**：TE 处 BL 方程切换为 wake 方程，$\delta^*$ 连续但 $C_f \to 0$
- **crossflow 模型**：用 Cooke's 假设或 Green's 3D entrainment 方法

### 6.3 验证

Nishida & Drela 在 Boeing 747 wing-body 上验证，与 N-S 解和实验数据对比。结论：
- 升力线斜率与实验一致到 2%
- 激波位置与 N-S 解一致到 3% chord
- 分离预测在 $C_L > 0.5$ 后开始偏离

---

## 7. UP3D 的 3D VII 实现路径

### 7.1 推荐架构：松耦合 → 紧耦合渐进路线

UP3D 的 VII 实现应分两个阶段：

| 阶段 | 策略 | 描述 | 复杂度 |
|------|------|------|--------|
| **V1-V2: 松耦合** | 迭代耦合 | 势流→IBL→修 BC→重解势流→… | 低 |
| **V3: 紧耦合** | simultaneous Newton | 增广 Jacobian（需 P8 Newton 基础） | 高 |

理由：
- 松耦合可以先验证 IBL 求解器本身的正确性和 transpiration BC 的实现
- 紧耦合需要 P8 Newton 的 Jacobian 基础——VII 紧耦合本质上是增广 Newton
- 松耦合在附着流中收敛足够快（通常 5-10 次外迭代）
- 紧耦合的收益主要在接近分离时（lag 误差大）

### 7.2 松耦合的具体流程

```
1. 求解纯势流 → φ⁰, U_e⁰(s)
2. 追踪壁面流线 → streamline mesh on wall
3. 沿流线积分 IBL → δ*¹(s), θ¹(s), Cf¹(s)
4. 计算 transpiration → ṁ_wall¹(s) = d/ds[ρ_e U_e δ*¹]
5. 修正壁面 BC → g = ṁ_wall¹
6. 重新求解势流 → φ¹, U_e¹(s)
7. 检查收敛 → ||δ*^{k+1} - δ*^k|| / ||δ*^k|| < tol
8. 未收敛 → 回到 2
```

### 7.3 紧耦合的具体流程（V3，需 P8）

```
组装增广残差:
  R_φ(φ, δ*) = 势流残差 + transpiration 项
  R_BL(δ*, φ) = IBL 残差（含 U_e(φ)）

组装增广 Jacobian:
  | J_φφ    J_φ,δ*   |
  | J_δ*,φ  J_δ*,δ*  |

Newton 步:
  | J_φφ    J_φ,δ*   | | Δφ   |   | -R_φ  |
  | J_δ*,φ  J_δ*,δ*  | | Δδ*  | = | -R_BL |

更新: φ ← φ + Δφ, δ* ← δ* + Δδ*
```

### 7.4 与 UP3D 现有代码的接口

| UP3D 模块 | VII 修改 | 阶段 |
|-----------|---------|------|
| `post/surface.py` | 输出 $U_e(s_1, s_2)$, $dU_e/ds$ 给 IBL | V1 |
| `kernels/residual.py` | 壁面 Neumann 项 $g \neq 0$：$\int_{\Gamma_w} N_i \dot{m} \, d\Gamma$ | V2 |
| `kernels/jacobian.py` | Picard 矩阵不变（transpiration 只进 RHS） | V2 |
| `solve/picard.py` | 外层加 VII 迭代循环 | V2 |
| `constraints/wake.py` | 尾流面 RHS 加 $\dot{m}_{\text{wake}}$ | V2 |
| `kernels/jacobian_newton.py`（P8 新建） | 增广 Jacobian 含 $J_{\varphi,\delta^*}$ | V3 |
| 新建 `pyfp3d/viscous/` | IBL 求解器 | V1 |

---

## 8. 壁面 Transpiration BC 的具体实现

### 8.1 当前壁面 BC

UP3D 当前壁面 BC 是 slip BC：

$$\rho \frac{\partial\varphi}{\partial n} = 0 \quad \text{on } \Gamma_{\text{wall}}$$

在 Galerkin 弱形式中，积分分部后壁面项自然为零（natural BC / do-nothing）：

$$R_i = \sum_e \int_{\Omega_e} \tilde{\rho} (\nabla\varphi \cdot \nabla N_i) \, d\Omega - \underbrace{\int_{\Gamma_{\text{wall}}} N_i \cdot g \, d\Gamma}_{=0 \text{ when } g=0}$$

### 8.2 加入 transpiration 后

$$g(s_1, s_2) = \dot{m}_{\text{wall}}(s_1, s_2) = \nabla_\Gamma \cdot (\rho_e U_e \delta^* \hat{e}_{\text{sl}})$$

其中 $\nabla_\Gamma$ 是曲面散度，$\hat{e}_{\text{sl}}$ 是壁面流线方向。

残差变为：

$$R_i = \sum_e \int_{\Omega_e} \tilde{\rho} (\nabla\varphi \cdot \nabla N_i) \, d\Omega - \int_{\Gamma_{\text{wall}}} N_i \cdot \dot{m}_{\text{wall}}(s_1, s_2) \, d\Gamma$$

### 8.3 离散实现

壁面积分在壁面三角形上做：

```python
for tri in wall_triangles:
    # 顶点 i, j, k
    # 形函数 N_i, N_j, N_k（面积坐标）
    # ṁ_wall 在三角形内线性插值（节点值）
    
    # 积分 ∫ N_i * ṁ_wall dΓ（面积 = A_tri）
    # 对线性 N_i 和线性 ṁ_wall:
    #   ∫ N_i * ṁ_wall dΓ = A_tri/3 * (ṁ_i + ṁ_j + ṁ_k)/3
    # 简化: A_tri/3 * mean(ṁ) 分配到三个节点
    
    R[i] -= A_tri/3 * (ṁ_i + ṁ_j + ṁ_k)/3
    # 同理 R[j], R[k]
```

这是标准的线性载荷离散——`kernels/residual.py` 中壁面项的位置已经有占位（$g=0$ 时为零），只需填入 $\dot{m}$。

### 8.4 可压缩 transpiration 的密度

可压缩流中，transpiration 公式中的 $\rho_e$ 需要从壁面速度计算：

$$\rho_e = \rho_\infty \left[1 + \frac{\gamma-1}{2} M_\infty^2 \left(1 - \frac{U_e^2}{U_\infty^2}\right)\right]^{1/(\gamma-1)}$$

这与 UP3D 的等熵密度公式（`physics/isentropic.py::density_isentropic`）完全一致——只是用在壁面外缘速度上。

### 8.5 对 Picard 矩阵的影响

关键：transpiration 只进 RHS，**不进矩阵**。

Picard 矩阵 $A_{ij} = \sum_e \tilde{\rho}_e (\nabla N_i \cdot \nabla N_j) V_e$ 不含 $\dot{m}$——$\dot{m}$ 是已知函数（从上一次 IBL 解读出），不是 $\varphi$ 的函数。

所以在松耦合中：
- 矩阵不变（$g$ 是 lagged 的）
- 只改 RHS
- `solve/linear.py` 不动

在紧耦合中：
- $\dot{m}$ 依赖 $\delta^*$，$\delta^*$ 是未知量
- $\delta^*$ 依赖 $U_e(\varphi)$，$\varphi$ 也是未知量
- 需要增广 Jacobian 的交叉项 $J_{\varphi,\delta^*}$

---

## 9. 尾流面 IBL 修正

### 9.1 尾迹方程

BL 在后缘不结束——它变为尾迹。尾迹方程与 BL 方程类似，但 $C_f = 0$（无壁面）：

$$\frac{d\theta_{\text{wake}}}{ds} + (2+H_{\text{wake}}-M_e^2)\frac{\theta_{\text{wake}}}{U_e}\frac{dU_e}{ds} = 0$$

初始条件来自 TE 处的 BL 解：$\theta_{\text{wake}}(s_{\text{TE}}) = \theta_{\text{upper}} + \theta_{\text{lower}}$。

### 9.2 尾流面 transpiration

尾流面的质量通量条件从纯无粘的 $\rho^+ u^+ - \rho^- u^- = 0$ 变为：

$$\rho^+ u^+ - \rho^- u^- = \dot{m}_{\text{wake}}(s_1, s_2)$$

其中 $\dot{m}_{\text{wake}} = \nabla_W \cdot (\rho_e U_e \delta^*_{\text{wake}} \hat{e}_{\text{streamline}})$。

### 9.3 UP3D 的尾流面实现

UP3D 用 master-slave 消元处理尾流面。当前 RHS 是：

$$b_{\text{red}}(\Gamma) = T^T b - \sum_j \Gamma_j h_j$$

加入尾流 transpiration 后，RHS 增加一项：

$$b_{\text{red}}(\Gamma, \delta^*_{\text{wake}}) = T^T b - \sum_j \Gamma_j h_j + T^T b_{\text{wake}}(\delta^*_{\text{wake}})$$

其中 $b_{\text{wake},i} = \int_{\Gamma_W} N_i \dot{m}_{\text{wake}} \, d\Gamma$。

这与当前代码的接口完全兼容——`constraints/wake.py` 的 `reduce` 方法只需接受一个额外的 RHS 贡献。

### 9.4 尾流面流线追踪的困难

3D 尾流面上的流线方向在翼尖附近可能与来流方向偏差很大（翼尖涡卷起）。$\delta^*_{\text{wake}}$ 沿尾流面流线的演化需要在尾流面网格上做流线追踪——这是一个曲面上的 2D 流线追踪问题。

简化策略：假设尾流面流线方向 = 来流方向（López 的尾流面线性化假设 Eq.3.44）。这在翼尖附近不准确，但对升力/阻力的影响是二阶的（$\delta^*_{\text{wake}}$ 在翼尖附近趋于零）。

---

## 10. IBL 求解器的设计

### 10.1 模块结构

```
pyfp3d/viscous/
├── __init__.py
├── streamline.py      # 壁面流线追踪
├── bl_equations.py    # 3D 可压缩 BL 积分方程
├── wake_equations.py  # 尾迹方程
├── transpiration.py   # δ* → ṁ_wall, ṁ_wake
└── solver.py          # IBL 求解器（沿流线积分 + Newton 局部解）
```

### 10.2 IBL 方程的选择

| 方法 | 模型 | 优点 | 缺点 |
|------|------|------|------|
| **Thwaites + Head** | 2D 不可压 | 最简单 | 不可压；无 crossflow |
| **Green's lag entrainment** | 2D 可压缩 | 跨声速适用；成熟 | 无 crossflow |
| **3D integral (Nishida-Drela)** | 3D 可压缩 | 有 crossflow | 复杂；需曲面几何 |
| **MSES-style (Drela-Giles)** | 2D 可压缩 | 同时求解 | 2D only |

推荐：**V1 用 Green's lag entrainment（2D 流线积分），V2 扩展到 3D crossflow**。

理由：
- Green's lag entrainment 是跨声速 IBL 的标准方法，XFOIL/MSES 都用
- V1 沿壁面流线逐条积分，每条流线是 1D 问题——可以用 2D 成熟方法
- crossflow 效应在后掠翼前缘附近重要，但在 UP3D 当前的 NACA 0012 extruded 准 2D case 中不显著
- ONERA M6 有 30° 后掠，V2 需要 crossflow

### 10.3 壁面流线追踪

在壁面三角形网格上追踪流线：

```python
def trace_streamlines(nodes, triangles, u_tangent, seed_points):
    """
    从 seed points 出发，沿 u_tangent 方向追踪壁面流线。
    
    Algorithm:
    1. 在每个三角形内，u_tangent 是常数（P1 恒定梯度）
    2. 从种子点出发，沿 u_tangent 方向走到三角形边界
    3. 在相邻三角形继续追踪
    4. 终止条件：到达 TE 或远场
    """
```

种子点策略：在翼面前缘均匀布点，或用壁面节点本身作为种子。

### 10.4 BL 方程的离散

沿流线用有限差分：

$$\frac{d\theta}{ds} = f(\theta, H, U_e, dU_e/ds, M_e, ...)$$

用向后差分（upwind，因为 BL 方程沿流向）或 Newton 迭代（同时求解 $\theta, H$）。

每条流线是 1D 问题，可以独立求解——天然并行。

---

## 11. 迭代策略与收敛性

### 11.1 松耦合的收敛行为

松耦合 VII 的收敛取决于：
- 边界层对势流的反馈强度（$\delta^*$ 的大小）
- 势流对边界层的敏感度（$dU_e/d\delta^*$ 的大小）

在附着流中：
- $\delta^*/c \sim 0.01$——小扰动
- 每次修正后势流变化小
- 通常 5-10 次外迭代收敛

在接近分离时：
- $\delta^*$ 急剧增长
- 反馈变强
- 可能需要 under-relaxation：$\delta^{*}_{k+1} = \omega \delta^{*}_{\text{IBL}} + (1-\omega) \delta^{*}_k$，$\omega \sim 0.3-0.5$

### 11.2 与 Mach continuation 的交互

UP3D 的 Mach continuation 在势流层做。VII 需要在每个 Mach 步做完整的 VII 迭代：

```
for M in Mach_steps:
    solve potential at M (with current δ*)
    VII iterate:
        solve IBL at M
        update transpiration
        solve potential at M (with new δ*)
        check convergence
    proceed to next M
```

第一个 Mach 步（$M_\infty - \Delta M$）可能需要更多 VII 迭代（冷启动）。后续步可以用上一步的 $\delta^*$ 作为初值。

### 11.3 与 P8 Newton 的交互

P8 的全耦合 Newton 可以自然扩展为 VII 紧耦合：

增广残差：

$$\mathbf{R}_{\text{aug}} = \begin{pmatrix} R_\varphi(\varphi, \Gamma, \delta^*) \\ R_\Gamma(\varphi, \Gamma) \\ R_{\text{BL}}(\delta^*, \varphi) \end{pmatrix}$$

增广 Jacobian：

$$\mathbf{J}_{\text{aug}} = \begin{pmatrix} J_{\varphi\varphi} & J_{\varphi\Gamma} & J_{\varphi,\delta^*} \\ J_{\Gamma\varphi} & J_{\Gamma\Gamma} & 0 \\ J_{\delta^*,\varphi} & 0 & J_{\delta^*,\delta^*} \end{pmatrix}$$

其中：
- $J_{\varphi,\delta^*} = \partial \dot{m} / \partial \delta^*$：transpiration 对 $\delta^*$ 的导数——解析可推导
- $J_{\delta^*,\varphi} = \partial R_{\text{BL}} / \partial \varphi = \partial R_{\text{BL}} / \partial U_e \cdot \partial U_e / \partial \varphi$：BL 残差对势流的导数

这使得 VII 紧耦合的额外代价主要是：
- 增广矩阵的 size 增加（壁面节点 × 2-3 BL 变量）
- 交叉 Jacobian 的计算

但 GMRES+AMG 可以处理增广系统——只是矩阵更大更不对称。

### 11.4 推荐路线

```
P8 Newton (无粘) → V1 松耦合 VII → V2 3D crossflow → V3 紧耦合
```

V1 可以与 P8 并行开发（松耦合不需要 Newton Jacobian）。V3 需要 P8 完成后才能做。

---

## 12. 与 P8 Newton 的交互

### 12.1 P8 为 VII 提供什么

P8 的全耦合 Newton Jacobian 是 VII 紧耦合的基础：

- $J_{\varphi\varphi}$：势流 Jacobian（已有）
- $J_{\varphi\Gamma}$：Kutta 残差对 $\varphi$ 的导数（已有）
- $J_{\Gamma\varphi}$：$\varphi$ 对 Kutta 残差的导数（已有）

VII 只需添加：
- $J_{\varphi,\delta^*}$：transpiration 对 $\delta^*$ 的导数
- $J_{\delta^*,\varphi}$：BL 残差对 $\varphi$ 的导数
- $J_{\delta^*,\delta^*}$：BL 残差对 $\delta^*$ 的导数

### 12.2 $J_{\varphi,\delta^*}$ 的推导

壁面残差项：

$$R_{\text{wall},i} = -\int_{\Gamma_w} N_i \dot{m} \, d\Gamma, \quad \dot{m} = \frac{d}{ds}[\rho_e U_e \delta^*]$$

对 $\delta^*$ 求导（固定 $\varphi$，即 $U_e$, $\rho_e$ 固定）：

$$\frac{\partial R_{\text{wall},i}}{\partial \delta^*_j} = -\int_{\Gamma_w} N_i \frac{\partial \dot{m}}{\partial \delta^*_j} \, d\Gamma = -\int_{\Gamma_w} N_i \frac{d}{ds}[\rho_e U_e N_j] \, d\Gamma$$

这是一个壁面三角形上的积分，解析可计算。

### 12.3 $J_{\delta^*,\varphi}$ 的推导

BL 残差（动量积分方程）：

$$R_{\text{BL}} = \frac{d\theta}{ds} + (2+H-M_e^2)\frac{\theta}{U_e}\frac{dU_e}{ds} - \frac{C_f}{2}$$

$U_e$ 和 $dU_e/ds$ 依赖 $\varphi$（壁面速度梯度）。所以：

$$\frac{\partial R_{\text{BL}}}{\partial \varphi_j} = \frac{\partial R_{\text{BL}}}{\partial U_e} \frac{\partial U_e}{\partial \varphi_j} + \frac{\partial R_{\text{BL}}}{\partial (dU_e/ds)} \frac{\partial (dU_e/ds)}{\partial \varphi_j}$$

$\partial U_e / \partial \varphi_j$ 就是壁面速度对 $\varphi$ 的导数——已经在 P8 的 Jacobian 中（壁面梯度恢复的导数）。

---

## 13. 验证策略

### 13.1 验证梯级

| 级别 | 测试 | 参考数据 | 预期 |
|------|------|---------|------|
| **V1.1** | 2D NACA 0012, M=0.5, α=5°, Re=6e6 | XFOIL VII 解 | $C_L$ 偏差 < 3% |
| **V1.2** | 2D NACA 0012, M=0.7, α=2° | XFOIL VII 解 | 激波位置偏差 < 5% chord |
| **V2.1** | 3D ONERA M6, M=0.84, α=3.06° | 实验 + Nishida-Drela 1995 | $C_L$ 偏差 < 5%（vs 无粘改善 > 3%） |
| **V2.2** | 3D 矩形翼 NACA 0012, AR=4 | 无粘解 vs VII 解 | VII 应降低 $C_L$ 过预测 |
| **V3.1** | 2D NACA 0012 接近失速, α=12° | XFOIL | 松耦合发散但紧耦合收敛 |

### 13.2 关键验证指标

1. **粘性修正的方向性**：VII 应该降低 $C_L$（无粘过预测），尤其在吸力峰处
2. **激波前移**：BL 位移使有效翼型变厚，激波位置应前移
3. **$C_L$ vs α 极曲线**：VII 应延迟失速角（vs 无粘过早失速）
4. **壁面 Cp**：VII 的吸力峰应低于无粘解（BL 位移降低了有效攻角）

### 13.3 与 López 数据的对比

López 报告了纯无粘解的 $C_L$ 偏差：
- NACA 0012 $\alpha=10°$：吸力峰低估，$C_L$ 高估——这是无粘 FP 的典型症状
- ONERA M6：$C_L = 0.288$（KRATOS）vs 0.286（SU2 Euler）——都是无粘，偏差小

VII 的效果应该让 $C_L$ 从无粘值降低，向实验值靠近。López 的 ONERA M6 CL=0.288（无粘），实验约 0.26-0.27（含粘性），VII 应给出 0.27-0.28。

---

## 14. 风险与未解决问题

### 14.1 技术风险

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| **壁面流线追踪在非结构网格上不稳定** | 高 | 高 | 用 P6 的 smooth_wall_tangential_gradients 先光滑速度场；稳健的 RK4 追踪；seed 点加密 |
| **IBL 在 3D 后掠翼前缘 crossflow 效应** | 中 | 中 | V1 先用 2D 流线积分（忽略 crossflow），V2 加 3D correction |
| **松耦合在接近分离时不收敛** | 高 | 中 | under-relaxation；紧耦合（V3）作为 fallback |
| **尾流面 $\delta^*_{\text{wake}}$ 在翼尖附近不准确** | 中 | 低 | 翼尖 $\delta^* \to 0$，误差影响小；用来流方向近似流线方向 |
| **可压缩 BL 方程在激波后不适定** | 中 | 高 | 限制 VII 有效域到 $M_\text{shock} < 1.3$；激波后分离用经验修正 |
| **紧耦合增广 Jacobian 的 GMRES 收敛慢** | 中 | 中 | 块预处理（BL 块用 ILU，势流块用 AMG） |

### 14.2 与 UP3D 其他 phase 的依赖关系

```
P5 (已关闭) → P6 (已关闭) → P7 (differentiable flux) → P8 (Newton)
                                                          ↓
                                                     V3 (紧耦合)
                                                          ↑
P6 (smooth wall gradient) → V1 (松耦合) → V2 (3D crossflow) ↑
```

V1 只依赖 P6（光滑壁面梯度作为 IBL 输入）和 IBL 求解器自身。
V3 依赖 P8（Newton Jacobian 作为增广系统的基础）。

### 14.3 开放问题

1. **3D BL 方程的 crossflow 模型选择**：Cooke's 假设（简单但精度有限）vs Green's 3D entrainment（复杂但更准）vs 直接求解 3D parabolized BL 方程（最准但最贵）
2. **壁面流线追踪的网格依赖性**：非结构三角形上的流线追踪可能产生网格敏感的结果——需要量化
3. **TE 处的 BL→wake 过渡**：3D 中 TE 是一条线，上下翼面的 BL 在此汇合——过渡条件的 3D 推广需要小心
4. **transpiration BC 在 Newton Jacobian 中的远场影响**：壁面 transpiration 改变了远场涡修正中的 $\Gamma$——$J_{\varphi,\delta^*}$ 是否需要包含远场项？需要分析

---

## 15. 成熟工具对比分析（BLWF / TRANAIR / Drela IBL3）

> 基于 BLWF58 方法文档 (108pp)、DLR 评估报告、Boeing 30 年 CFD 回顾、Drela 2013 的交叉对比

### 15.1 三个工具的架构对比

| 维度 | **BLWF58** (TsAGI) | **TRANAIR** (Boeing) | **Drela IBL3** (MIT) | **UP3D** (目标) |
|------|-----------|-----------|-----------|---------|
| **无粘求解器** | 全速势, 结构化 Chimera 网格, FVM | 全速势, 自适应矩形网格, FEM | 面涡法 (低阶 panel) / 可换 FE 势流 | 全速势, 非结构四面体, Galerkin P1 |
| **BL 方法** | 3D 有限差分 (Keller box), 可压缩 | 2D IBL + 后掠/锥度修正 (Drela 原始) | 3D 积分 (6 方程: δ, A, B, Ψ, C_τ1, C_τ2) | 待定 (推荐 Drela IBL3) |
| **尾流 BL** | 2D 积分 或 3D 有限差分 | 2D IBL | DDCV 积分方程 | 待定 |
| **耦合方案** | quasi-simultaneous (Hilbert 积分) | 直接耦合 (direct solver) | simultaneous Newton | 松→紧渐进 |
| **壁面 BC** | transpiration (位移厚度) | transpiration (位移厚度) | transpiration Eq.(76) | transpiration (VI-1 方案) |
| **分离处理** | inverse mode BL + quasi-simultaneous | 2D IBL 可处理分离 | 4 方程可表示分离 + crossover | 限于附着流 (V1) |
| **求解器** | AF (approximate factorization) + GMRES | sparse direct solver | Newton + direct (小问题) / GMRES (大问题) | Picard→Newton+GMRES |
| **网格类型** | 结构化贴体 (翼/身/短舱/尾各一块) | 自适应矩形 (embedded geometry) | 双线性四边形 (表面) | 非结构三角形 (表面) + 四面体 (体) |

### 15.2 BLWF58 的关键实现细节

BLWF58 是目前工程上最成熟的 3D FPE+BL 耦合工具之一。从方法文档 (108pp) 提取的关键技术：

**3D BL 求解器——Keller box 有限差分**：

- 在翼面贴体坐标 $(x, z)$ 上建立 3D 可压缩 BL 方程
- **Keller box scheme**：二阶精度隐式差分，自然处理 3D crossflow
- **predictor-corrector**：沿流向步进，predictor 用线性外推，corrector 用 Newton 迭代
- **inverse mode**：在分离区切换为规定 $C_p$ 求 $\delta^*$（而非规定 $u_e$ 求 $\delta^*$），避免 Goldstein 奇异性

**四种耦合方案的分析**（BLWF58 方法文档 §3.1-3.2）：

| 方案 | BL 解法 | 势流 BC | 适用域 | UP3D 对应 |
|------|--------|---------|--------|----------|
| **Direct** | $u_e$ 给定 → 求 $\delta^*$ | transpiration from $\delta^*$ | 附着流，弱相互作用 | V1 松耦合 |
| **Semi-inverse** | inverse ($C_p$ 给定) | $u_e$ 给定 | 分离区 | — |
| **Quasi-simultaneous** | direct + Hilbert 积分预估 $u_e$ 变化 | mixed BC (transpiration + interaction law) | 中等相互作用 | V2 |
| **Inverse** | inverse | $C_p$ 给定 | 强分离 | — |

BLWF58 选择 **quasi-simultaneous** 作为默认方案——它在 direct 和 inverse 之间取平衡：BL 计算时用 Hilbert 积分预估势流对 $u_e$ 变化的响应，势流计算时用 mixed BC。这使得耦合不需要 Newton 全局求解就能处理中等强度的粘性-无粘相互作用。

**关键工程细节**：
- TE 处的特殊处理（§3.3.2）：在翼面 TE 上下面的 BL 汇合处，需同时满足上下面的 transpiration BC + 尾流面连续性——BLWF58 在 TE 节点上直接求解耦合方程
- 弹性变形耦合（Appendix A）：简单梁理论，VII 迭代中同时更新气动和结构

### 15.3 TRANAIR 的关键实现细节

从 Boeing 2003 论文提取：

- **早期用 A411 BL 代码 + indirect coupling**——收敛可靠性不够，激波-BL 相互作用模型可疑
- **后来改用 Drela 2D IBL + 直接耦合**：TRANAIR 团队在 Drela 的 2D IBL 基础上加入后掠和锥度修正，直接耦合到全速势 FEM
- **关键优势**：直接求解器（sparse direct）使耦合矩阵的组装和求解很可靠
- **精度**："agreement with experiment was considered by project users to be quite remarkable"

**TRANAIR 的路线选择**：从 indirect coupling（= 松耦合）→ direct coupling（= 紧耦合 Newton），原因是松耦合的收敛可靠性不够。这与 Drela 2013 的论断一致——"simultaneous solution is highly desirable, if not essential"。

### 15.4 Drela IBL3 vs BLWF 的关键差异

| 方面 | Drela IBL3 (6 方程) | BLWF (有限差分 3D BL) | 对 UP3D 的启示 |
|------|---------------------|----------------------|---------------|
| **BL 表征** | 积分厚度 (δ, A, B, Ψ) — 少量未知量 | 全速度剖面 (Keller box) — 大量未知量 | IBL3 更适合 UP3D 的快速定位 |
| **分离处理** | 4 方程自然允许分离 + crossover | inverse mode 切换 | IBL3 更简洁 |
| **crossflow** | Ψ 方程（lateral curvature） | 有限差分直接解 3D 方程 | IBL3 的 crossflow 是近似的（薄层假设）；有限差分更精确但更贵 |
| **闭合关系** | 需要经验 profile 形状（streamwise + crossflow） | 不需要（直接解 PDE） | IBL3 的闭合关系是关键不确定性 |
| **网格要求** | 表面网格（bilinear 四边形） | 3D body-fitted 网格 | IBL3 只需表面网格——与 UP3D 的壁面三角形兼容 |
| **计算量** | ~100× 少于有限差分（Drela 2013 原话） | 基准 | IBL3 与 UP3D 的快速定位一致 |

### 15.5 对 UP3D VII 文档的影响

基于以上对比，确认和修正 `20260707_2118_ibl_viscous_coupling_design.md` 和 `20260709_0145_3d_vii_implementation_analysis.md` 中的决策：

**确认的决策**：

1. **Transpiration BC（VI-1 方案）**：BLWF58 和 TRANAIR 都用 transpiration（位移厚度），不用几何修正。Drela Eq.(76) 给出精确公式。✅ 确认。

2. **松→紧渐进路线**：BLWF58 用 quasi-simultaneous（介于松和紧之间）作为默认；TRANAIR 从松耦合升级到紧耦合。UP3D 的 V1→V3 渐进路线与此一致。✅ 确认。

3. **松耦合在分离附近不够**：BLWF58 和 TRANAIR 的经验都确认这一点。UP3D V1 限于附着流是合理的。✅ 确认。

**需要修正的决策**：

4. **IBL 方程选择**：原文档推荐 Green's lag entrainment（3 方程）→ **修正为 Drela IBL3（6 方程）**。理由：
   - Drela 2013 提供完整的 3D 理论推导和闭合关系
   - 4 方程可处理分离和 crossover（Green 不行）
   - 计算量远小于有限差分（BLWF 路线）
   - Drela 2013 的 FE 方法（局部 Cartesian 基 + bilinear 元）可适配 UP3D 的三角形壁面网格

5. **表面网格处理**：原文档说"非结构三角形上的流线追踪"→ **修正为 Drela 的表面 FE 方法**。Drela 2013 不沿流线积分——而是在表面网格上用 Galerkin 弱形式求解 IBL3 方程。这完全避开了流线追踪问题。具体：
   - 每个壁面三角形上定义局部 Cartesian 基 $(\hat{x}_i, \hat{y}_i, \hat{z}_i)$
   - 速度和积分厚度用 bilinear 基函数插值
   - 方程残差在节点上组装——与 UP3D 的体 FE 方法一致

6. **quasi-simultaneous vs simultaneous**：原文档推荐松耦合（V1）→ 紧耦合（V3）。BLWF58 的经验表明 quasi-simultaneous（Hilbert 积分）是更好的中间步骤。**修正**：V2 应考虑 quasi-simultaneous 而非纯松耦合。Hilbert 积分预估势流响应，只需局部操作，不需要全局 Newton。

7. **TE 处理**：BLWF58 §3.3.2 专门处理 TE 上下面 BL 汇合——在 TE 节点上直接求解耦合方程。**新增**：UP3D 的 V4（尾流面 IBL 修正）需要考虑 TE 处的 BL→wake 过渡。UP3D 的 master-slave 消元在 TE 处复制了节点，这为 TE 耦合提供了自然的框架。

### 15.6 UP3D 相对 BLWF/TRANAIR 的优势和劣势

**优势**：
- **非结构网格**：处理复杂几何比 BLWF 的结构化 Chimera 更灵活
- **嵌入尾流**：master-slave 消元比 BLWF 的 conforming wake 更自动化
- **Python+Numba**：原型迭代快，适合 VII 方案探索

**劣势**：
- **BL 求解器不存在**：BLWF 有 30 年成熟的 3D Keller box BL 求解器；UP3D 从零开始
- **结构化网格上的 BL 更稳定**：BLWF 的结构化网格天然适合 Keller box；UP3D 的非结构三角形需要 Drela 的表面 FE 方法
- **quasi-simultaneous 的 Hilbert 积分在非结构网格上更复杂**：BLWF 在结构化网格上用沿流向的 Hilbert 积分；UP3D 需要在非结构网格上实现等价操作

---

## 参考文献

- López Canalejo, I. 2021. *A Finite-Element Transonic Potential Flow Solver with an Embedded Wake Approach for Aircraft Conceptual Design.* PhD diss., TUM.
  - §2.4: 三种力计算方法（近场/远场/势跳跃）
  - §4.1: 无粘解的粘性偏差讨论
  - §5.2: Outlook — BL 耦合展望
- Nishida, B., Drela, M. 1995. "Fully simultaneous coupling for three-dimensional viscous/inviscid flows." AIAA Paper 1995-1806.
  - 3D VII 的核心文献；FE 势流 + 3D IBL + simultaneous Newton
- Drela, M. 2014. *Flight Vehicle Aerodynamics.* MIT Press.
  - IBL 理论、2D VII（XFOIL/MSES 路线）、动量积分方程推导
- Davari, A. et al. 2019. "A CutFEM approach for the solution of the full potential equation with embedded wakes." *Comput. Mech.* 63:821–833.
  - body-fitted FPE + 嵌入尾流；Appendix A BL 理论；§3.2 BL 耦合展望
- Rodríguez, D. et al. 2012. "A fast coupled boundary-layer method for Cartesian grids." AIAA Paper 2012-302.
  - Davari 2019 [33]；CART3D + BL 耦合
- Green, J.E. et al. 1972. "Prediction of turbulent boundary layers and wakes in compressible flow by a lag-entrainment method." RAE TR 72231.
  - Green's lag entrainment 方法的原始文献
- Stock, H.W., Haase, W. 1999. "Feasibility study of e^t transition prediction in Navier-Stokes methods for airfoils." *Aerospace Science and Technology* 3(6): 405-412.
  - 3D BL 方程的工程实现参考
- **Drela, M. 2013. "Three-Dimensional Integral Boundary Layer Formulation for General Configurations." AIAA Paper 2013-2437.** — `references/Drela_2013_IBL3_general_configurations.pdf`
  - 3D IBL3 完整推导：6 方程 (δ, A, B, Ψ, C_τ1, C_τ2)；DCV/DDCV 框架；§III Eq.(76) transpiration BC；FE 方法；分离流测试
- **Karas, O.V., Kovalev, V.E. *BLWF58 Computational Method and Algorithms.* TsAGI, 108pp.** — `references/BLWF58_method_algorithms.doc`
  - 3D Keller box BL 求解器；quasi-simultaneous coupling (Hilbert 积分)；四种耦合方案分析；TE 处理；结构化 Chimera 网格
- **Karas, O.V., Kovalev, V.E. *BLWF58 User's Guide.* TsAGI, 166pp.** — `references/BLWF58_user_guide.doc`
  - 输入格式；网格生成；结果可视化
- **Zhang, K., Hepperle, M. 2010. *Evaluation of the BLWF Code.* DLR IB 124-2010/3.** — `references/DLR_2010_BLWF_evaluation.pdf`
  - DLR 对 BLWF 的独立评估；DLR-F4 翼身组合体验证
- **Johnson, F.T., Tinoco, E.N., Yu, N.J. 2003. "Thirty Years of Development and Application of CFD at Boeing Commercial Airplanes, Seattle." AIAA Paper 2003-3439.** — `references/Boeing_2003_AIAA_30yr_CFD.pdf`
  - TRANAIR 发展史；从 indirect coupling 到 direct coupling 的演进
- **Bolsunovsky, A.L. et al. 2010. "An Experience in Aerodynamic Design of Transport Aircraft." ICAS 2010.** — `references/TsAGI_2010_ICAS_aircraft_design_experience.pdf`
  - TsAGI 运输机设计经验；BLWF 在设计流程中的应用
