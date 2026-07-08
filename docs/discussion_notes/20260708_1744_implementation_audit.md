# UP3D 实现路径全面分析报告

> 基于 López 2021 博士论文 (TUM) + Kratos 框架架构的交叉验证
> Date: 2026-07-08 17:44
> Code baseline: HEAD 15c908e
> References:
> - López Canalejo 2021 博士论文 — `references/Lopez_2021_dissertation_transonic_FPE_embedded_wake.pdf`
> - Kratos Multiphysics — https://github.com/KratosMultiphysics/Kratos
> - UP3D `docs/design.md`, `docs/roadmap.md`

---

> **⚠️ 后续修正（2026-07-08，roadmap/design P7·P8 复核）。** 本审计两处判断已修正，
> 读者以 `docs/roadmap.md` P7/P8 + `design.md` §6.3 为准：
> 1. **§4.1/§4.2 "UP3D 缺 max(μ, μ_up)" 不准确。** P4 walk 已实现激波点算子
>    ν=max(ν_e, ν_up)（roadmap P4 G4.1 scheme-hardening (3)），P7 kernel 设计亦在
>    光滑混合的上游 Mach 上取 max(ν_e, ν_up)（design.md §218）。差异是"当前单元 vs
>    当前+上游"的施加范围，不是完全缺失。
> 2. **§11.2/§12（stencil "宽一层"、"+30% 内存"）对 UP3D 不成立。** 那是 López 单跳
>    结论；UP3D 上游为多跳（sliver tet）。**P7 定案（re-scoped 2026-07-08）**：默认
>    Newton 通量 = 冻结选择下的 walk（∂ρ̃/∂φ 稀疏，~+1 上游单元/行，但图距 ≤4）；
>    kernel（耦合整个 depth-3 邻域，明显更稠密）降级为可选，仅当 N2 实测划算才作 P8
>    通量。非既定"一层 +30%"。

## 目录

1. [总体架构对比](#1-总体架构对比)
2. [理论：控制方程与物理模型](#2-理论控制方程与物理模型)
3. [数值格式：空间离散](#3-数值格式空间离散)
4. [数值格式：人工可压缩性](#4-数值格式人工可压缩性)
5. [数值格式：尾流与 Kutta 条件](#5-数值格式尾流与-kutta-条件)
6. [数值格式：非线性求解策略](#6-数值格式非线性求解策略)
7. [数值实现：线性代数与并行](#7-数值实现线性代数与并行)
8. [数值实现：边界条件](#8-数值实现边界条件)
9. [数值实现：网格与预处理](#9-数值实现网格与预处理)
10. [数值实现：后处理与验证](#10-数值实现后处理与验证)
11. [P6/P7 实现路径技术细节](#11-p6p7-实现路径技术细节)
12. [风险与未解决问题](#12-风险与未解决问题)
13. [总结：实现路径合理性评估](#13-总结实现路径合理性评估)

---

## 1. 总体架构对比

### 1.1 UP3D vs Kratos/López 架构定位

| 维度 | UP3D | López/Kratos |
|------|------|-------------|
| **语言** | Python + Numba | C++ (Kratos kernel) + Python (driver) |
| **FE 框架** | 自研（PicardOperator + colored prange） | Kratos Multiphysics（OOP 层次架构） |
| **线性代数** | SciPy sparse + PyAMG | AMGCL (C++, Demidov 2020) |
| **网格** | 非结构四面体 P1，body-fitted | 非结构四面体 P1，embedded wake |
| **尾流** | conforming 面 + master-slave 消元 | embedded 面 + auxiliary DOF + least-squares BC |
| **非线性求解** | Picard（固定矩阵）+ secant Γ + Mach continuation | Newton-Raphson（完整 Jacobian）+ load stepping |
| **并行** | Numba prange（共享内存，colored） | Kratos MPI + OpenMP |
| **适用阶段** | 概念设计/优化（与 López 相同） | 同 |

### 1.2 合理性评估

UP3D 的架构选择对 Python 生态中的快速原型开发是**合理的**：

- **Numba + SciPy/PyAMG** 避免了 C++ 编译开销，适合 vibe-coding 迭代
- **body-fitted + conforming wake** 比 López 的 embedded approach 更简单（不需要 cut-cell 处理和小切单元全积分），代价是网格生成更受约束
- **master-slave 消元** 比 López 的 least-squares wake BC 更简洁（Γ 只进 RHS，矩阵保持 SPD），且对全耦合 Newton 更友好（∂R/∂Γ = -h_j 已知）

**关键差异的影响**：

1. **尾流处理**：UP3D 的 master-slave 消元在 Picard 路径上比 López 的 least-squares 更精确（jump 恰好等于 Γ，不是最小二乘近似），且不需要 small-cut 全积分（§3.5.5）。代价是网格必须在尾流面 conforming，不能像 López 那样隐式嵌入。
2. **非线性求解**：UP3D 的三层嵌套（Mach continuation → secant Γ → density Picard）是**临时方案**，P5 暴露的 H3 耦合失稳证明此架构在跨声速不可靠。P7 全耦合 Newton 是正确方向。
3. **性能**：Python + Numba 的 overhead 约 3-5×（vs C++），但 PyAMG 的 AMG setup 是编译路径，性能差距主要在装配和 matvec。对 medium mesh（~350k tets）可接受。

---

## 2. 理论：控制方程与物理模型

### 2.1 全速势方程

| 项 | UP3D design.md §2 | López §2.3 | 一致性 |
|---|---|---|---|
| 守恒形式 | $\nabla \cdot (\rho \nabla \varphi) = 0$ (2.1) | $\nabla \cdot (\rho \mathbf{u}) = 0$ (3.1) | ✅ 一致 |
| 等熵密度 | $\rho = [1 + \frac{\gamma-1}{2} M_\infty^2 (1 - q^2)]^{1/(\gamma-1)}$ (2.2) | 同 (Eq.2.41) | ✅ 一致 |
| 局部 Mach | $M^2 = q^2 M_\infty^2 / [1 + \frac{\gamma-1}{2} M_\infty^2 (1 - q^2)]$ (2.3) | 同 (Eq.2.43) | ✅ 一致 |
| Cp | $C_p = \frac{2}{\gamma M_\infty^2}[\rho^\gamma - 1]$ (2.5) | 同 (Eq.2.46) | ✅ 一致 |

**代码验证** (`physics/isentropic.py`)：
- `density_isentropic(q2, m_inf, gamma)` — 实现了 (2.2)，经单元测试验证 ✅
- `mach_number_squared(q2, m_inf, gamma)` — 实现了 (2.3) ✅
- `pressure_coefficient(q2, m_inf, gamma)` — 实现了 (2.5) ✅
- `critical_speed_squared(m_inf, gamma)` — 实现 $q^{*2}$，López 也在 §2.3 推导 ✅

### 2.2 假设范围

| 假设 | UP3D | López | 备注 |
|---|---|---|---|
| 等熵 | ✅ | ✅ | 激波用等熵压缩近似，非 Rankine-Hugoniot |
| 无旋 | ✅ | ✅ | 势流假设 |
| 无粘 | ✅ | ✅ | 高 Re 外流 |
| 定常 | ✅ | ✅ | 无时间导数 |
| M∞ 范围 | 0.3–0.87 | 0.3–0.85 | UP3D 略宽 |
| 激波强度 | 局部法向激波 Mach ≲ 1.3 | 同 | 超出此范围 FP 模型误差增大 |

### 2.3 极限速度 / 密度截断

| 项 | UP3D | López §3.4 | 一致性 |
|---|---|---|---|
| 极限速度 | $u_{\mathrm{lim}} = u_\infty \sqrt{1 + \frac{2}{(\gamma-1)M_\infty^2}}$ | 同 Eq.(3.29) | ✅ |
| M² 截断 | `limit_q2_field`: $M^2_{\max} = m_{\mathrm{cap}}^2 = 9.0$ (m_cap=3) | $M^2_{\max} \sim O(3.0)$ §3.4 | ✅ 一致 |
| 截断方式 | hard clamp | hard clamp | ✅ |
| 激活区域 | 翼尖 TE 几何奇异性 | 同（Fig.3.5, 2888 m/s） | ✅ |
| 收敛态 | 0 floored / 0 limited (P5 medium) | 0 (Table 4.9 收敛后) | ✅ |

**代码验证** (`physics/isentropic.py::limit_q2_field`)：
- hard clamp 实现：$M^2 > M^2_{\max}$ 时 $q^2$ 被截断到 $q^2(M^2_{\max})$
- López 用同样方式，且在收敛态截断非激活——对 Newton 可微性无阻碍

---

## 3. 数值格式：空间离散

### 3.1 Galerkin P1 有限元

| 项 | UP3D design.md §6 | López §3.2 | 一致性 |
|---|---|---|---|
| 单元 | 线性四面体 P1 | 线性三角/四面体 P1 | ✅ |
| 弱形式 | $R_i = \sum_e \tilde{\rho}_e (\nabla\varphi_e \cdot \nabla N_i) V_e = 0$ (6.1) | $R_i = \sum_e \int \tilde{\rho} \mathbf{u} \cdot \nabla N_i \, d\Omega_e - \sum_c \int N_i q \, d\Gamma_c$ (3.9) | ⚠️ |
| 预计算 | $B_e$ (4×3 shape gradient), $V_e$ — 一次性 | 同 | ✅ |
| 速度 | $\nabla\varphi_e = \sum_k \varphi_k \nabla N_k\|_e$ | 同 | ✅ |

**差异分析**：

UP3D 的残差 (6.1) 用 $\tilde{\rho}_e$ 作为**冻结的逐元素标量**乘以 $\nabla\varphi_e \cdot \nabla N_i \cdot V_e$。López 的 (3.9) 形式上相同，但多了一个面通量项 $\sum_c \int N_i q \, d\Gamma_c$——这是 Neumann BC 的面积分。

**UP3D 的处理**：design.md §5 说壁面用 natural BC（do-nothing），即面通量项为零。远场用 Dirichlet 消元。因此 UP3D 不显式组装面通量项，而是通过 Dirichlet 消元 + natural BC 隐式处理。这在数学上等价，实现更简洁。

**代码验证** (`kernels/residual.py::assemble_residual`)：
- 逐元素循环：`R[tet[i]] += V * dotp(grad_phi, grads[i])` — 纯体积积分，无面通量 ✅
- 与 (6.1) 一致 ✅
- 参考路径（serial P1）保留作为回归基准 ✅

### 3.2 Picard 矩阵

| 项 | UP3D (6.2) | López (3.16) | 一致性 |
|---|---|---|---|
| 矩阵 | $A_{ij} = \sum_e \tilde{\rho}_e (\nabla N_i \cdot \nabla N_j) V_e$ | 同（不可压缩形式 $\rho_\infty = 1$） | ✅ |
| 对称性 | SPD when $\nu \equiv 0$ | 同 | ✅ |
| 冻结密度 | $\tilde{\rho}_e$ 作为标量，忽略 $\partial\tilde{\rho}/\partial\varphi$ 耦合 | Newton 时包含；Picard 时忽略 | ✅ |

**代码验证** (`kernels/jacobian.py::PicardOperator`)：
- `assemble_matrix(rho_tilde)` — 用预计算 sparsity + colored prange 装配 ✅
- `assemble_residual(phi, rho_tilde)` — 残差装配 ✅
- `velocities(phi)` — 逐元素 $\nabla\varphi_e$ + $q^2_e$ ✅
- 零分配热路径（预分配 buffer）✅
- elem_to_csr scatter map 一次性预计算 ✅
- colored prange 保证无竞争写入 ✅

### 3.3 精度阶

| 区域 | UP3D | López | 备注 |
|---|---|---|---|
| 亚声速 | 2 阶 | 2 阶 | P1 Galerkin 标准精度 |
| 激波附近 | 1 阶 | 1 阶 | 人工密度迎风一阶 |
| 壁面 Cp | ~1 阶（P1 几何误差） | ~2 阶（embedded，不需要壁面拟合） | UP3D 的 G1.6 问题 |

---

## 4. 数值格式：人工可压缩性

### 4.1 人工密度公式

| 项 | UP3D (3.1')-(3.2) | López Eq.(3.19)-(3.20) | 一致性 |
|---|---|---|---|
| 迎风密度 | $\tilde{\rho}_e = \rho_e - \nu_e(\rho_e - \rho_{u(e)})$ | $\tilde{\rho} = \rho - \mu_s(\rho - \rho_{\mathrm{up}})$ | ✅ 符号不同，形式相同 |
| 切换函数 | $\nu_e = C \cdot \max(0, 1 - M_c^2/M_e^2)$ | $\mu_s = C \cdot \max(0, \mu_e, \mu_{\mathrm{up}})$ | ⚠️ **不同** |
| 上游选择 | multi-hop directional walk (≤4 hops) | 逐元素上游（integer-walk，每步更新） | ⚠️ UP3D 多跳 |

**关键差异 1：切换函数**

López 的 $\mu_s = C \cdot \max(0, \mu, \mu_{\mathrm{up}})$ 取**当前元素和上游元素的最大值**。UP3D 的 $\nu_e = C \cdot \max(0, 1 - M_c^2/M_e^2)$ 只取**当前元素**的值。

- López 形式：当当前元素或上游元素任一超临界时激活迎风——更保守
- UP3D 形式：只有当前元素超临界时激活——可能漏掉上游超临界但当前亚临界的情形

**影响**：在激波附近，López 的形式更早激活人工耗散，可能更稳定。UP3D 的 design.md §3.2 (P6 target) 已计划改为 smooth weighting，可同时修正此差异。

**关键差异 2：multi-hop walk**

UP3D 的 multi-hop walk（≤4 hops，目标 0.8 元素流向长度）是 P4 硬化路径的证据驱动扩展。López 的上游选择是单步 integer-walk。

- UP3D 的 rationale：M0 准 2D 网格的 sliver tet 单跳只达 0.25-0.39 元素流向长度，需要多跳才能达到有效耗散
- López 的 KRATOS 网格可能没有 sliver 问题，单跳足够
- P6 consistent flux (design.md §3.2 (3.4)) 用 smooth face-neighbor weighting 替代 integer-walk，消除了多跳需求

**代码验证** (`kernels/upwind.py`)：
- `UpwindOperator` 类：预计算 face_neighbors, centroids ✅
- `upstream_elements()` — Numba kernel，multi-hop walk ✅
- `rho_tilde_sweep()` — $\tilde{\rho}_e = \rho_e - \nu_e(\rho_e - \rho_{u(e)})$ ✅
- $\nu_e$ 的 $\max(0, ...)$ guard 保证亚临界时 $\nu \equiv 0$ bitwise（G4.2）✅
- `rho_floor=0.05` 防止负密度 ✅
- 监控：`nu_max`, `n_supersonic`, `n_floored` ✅

### 4.2 切换函数的三种流态

López §3.3.1 区分三种流态，UP3D 当前未显式区分：

| 流态 | López | UP3D | 备注 |
|---|---|---|---|
| 亚声速 $M < M_c$ | $\mu_s = 0 \Rightarrow \tilde{\rho} = \rho$ | $\nu = 0 \Rightarrow \tilde{\rho} = \rho$ | ✅ 一致 |
| 超音速加速 $M > M_c, M > M_{\mathrm{up}}$ | $\mu_s = \mu$ | $\nu = C(1 - M_c^2/M^2)$ | ⚠️ 缺 max(μ, μ_up) |
| 超音速减速 $M > M_c, M < M_{\mathrm{up}}$ | $\mu_s = \mu_{\mathrm{up}}$ | 同上（不区分） | ⚠️ López 取上游 μ |

**影响**：减速流态（激波后侧）López 用上游 μ，UP3D 用当前 μ。在激波后减速区，上游 μ 更大（上游更超临界），所以 López 引入更多耗散。UP3D 在激波后侧的耗散可能不足，导致振荡。

**P6 修正**：design.md §3.2 的 smooth max_ε 同时处理了此问题。

### 4.3 Jacobian 敏感度（P7 预备）

| 项 | López Appendix B | UP3D design.md §6.3 | 一致性 |
|---|---|---|---|
| Term 1 (Picard) | $\rho \nabla N_i \cdot \nabla N_j$ | 同 | ✅ |
| Term 2 (密度) | $(\partial\tilde{\rho}/\partial\varphi_j) u_a \partial N_i/\partial x_a$ | 同 | ✅ |
| Term 3 (上游) | $\partial\tilde{\rho}/\partial\varphi_j^{\mathrm{up}} \neq 0$ (B.4) | 同（design.md 已修正） | ✅ |
| 上游选择可微性 | 冻结 $u(e)$，只对 $\rho_{\mathrm{up}}$ 求导 | 同（design.md §3.1 已修正） | ✅ |
| 切换函数导数 | $\partial\mu/\partial\varphi$ 纳入 (B.3, B.6) | 同（design.md §6.3 已修正） | ✅ |

---

## 5. 数值格式：尾流与 Kutta 条件

### 5.1 尾流几何与表示

| 项 | UP3D | López §3.5 | 差异 |
|---|---|---|---|
| 尾流面 | conforming 内部面（Gmsh 建模） | embedded 面（隐式定义，不 conforming） | **架构差异** |
| 节点复制 | solver-side (`wake_cut.py`) | solver-side (auxiliary DOF) | 类似 |
| TE 节点 | 复制（jump reaches wall） | 复制（aux DOF 作为下方元素普通 DOF） | ✅ 目标一致 |
| 自由边（翼尖） | 不复制，pin $\Gamma_{\mathrm{tip}} = 0$ | N/A（embedded 自动处理） | UP3D 显式 |
| 网格生成约束 | 尾流必须 conforming | 无约束（隐式嵌入） | UP3D 更受约束 |

### 5.2 尾流边界条件

| 项 | UP3D | López | 差异 |
|---|---|---|---|
| 质量连续 | master-slave 消元自动满足 | $g_1 = \hat{n} \cdot (\mathbf{u}_u - \mathbf{u}_l) = 0$ (3.43)，least-squares 弱形式 | **方法不同** |
| 压力相等 | Kutta 条件间接保证 | $g_2 = \hat{u}_\infty \cdot (\mathbf{u}_u - \mathbf{u}_l) = 0$ (3.44)，least-squares | **方法不同** |
| 线性化 | 不需要（master-slave 是精确约束） | $\rho_u, \rho_l \to \rho_\infty$; $\bar{u} \to \hat{u}_\infty$ | UP3D 更精确 |
| Jacobian block | $T^T J T$（投影，无额外 block） | $J^{uu}, J^{ul}, J^{lu}, J^{ll}$ (3.48-3.51) | UP3D 更简单 |
| Γ 进入方式 | RHS only（$b_{\mathrm{red}} = T^T b - \sum_j \Gamma_j h_j$） | least-squares 残差中显式 | UP3D 更简洁 |

**合理性评估**：

UP3D 的 master-slave 消元在数学上比 López 的 least-squares 更精确：
- jump 恰好等于 Γ（不是最小二乘近似）
- 矩阵保持 $T^T A T$（SPD when $\nu=0$），CG+AMG 可用
- Γ 只进 RHS，Picard 路径不需要重新装配
- 对全耦合 Newton 更友好：$\partial R_{\mathrm{red}}/\partial \Gamma_j = -h_j$ 已在 `wake.py` 中预算

代价是网格必须在尾流面 conforming，不能像 López 那样隐式嵌入（对 aeroelastic optimization 不友好，因为攻角改变时尾流方向变）。

**代码验证** (`constraints/wake.py`)：
- `WakeConstraint` 类：`T` (投影矩阵), `A_reduced = T^T A T`, `h_j = T^T A g_j` ✅
- `reduced_rhs(gamma)` = $T^T b - \sum_j \Gamma_j h_j$ ✅
- `expand(phi_red, gamma)` = $\varphi_{\mathrm{full}} = T \varphi_{\mathrm{red}} + g(\Gamma)$ ✅
- `kutta_targets(phi, wc)` = 探针处 $\varphi_{\mathrm{upper}} - \varphi_{\mathrm{lower}}$ 的 per-station 均值 ✅

### 5.3 Kutta 条件

| 项 | UP3D | López | 差异 |
|---|---|---|---|
| Kutta 条件 | $\Gamma_j = \varphi_{\mathrm{TE,upper}} - \varphi_{\mathrm{TE,lower}}$ (4.4) | 通过压力相等 $g_2$ 隐式施加 | **方法不同** |
| 探针 | 壁面 TE 一阶邻居（`kutta_upper/lower`） | N/A（least-squares 自动处理） | UP3D 显式 |
| 更新方式 | secant (Aitken) 外层 | Newton 内层（Γ 作为未知数） | UP3D 临时方案 |
| 收敛标准 | $\|\Gamma_{\mathrm{target}} - \Gamma\|_\infty < \mathrm{tol}$ | Newton 残差 $\|F\| \to 0$ | 不同的收敛量 |

**H3 耦合失稳分析**：

UP3D 的三层嵌套（Mach continuation → secant Γ → density Picard）在跨声速暴露的 H3 失稳：
- secant 读到的 $F = \mathrm{kutta\_targets}(\varphi(\Gamma, \rho(\Gamma))) - \Gamma$，$\rho$ 未完全收敛时 $F$ 被密度滞后污染
- 坏 slope → Γ 过冲 → $\varphi$ 大变 → $\rho$ 需更多步 → $F$ 污染更重 → 正反馈

López 的 Newton 架构不存在此问题——残差是 $R_i(\varphi) = 0$，每步同时更新 $\varphi$ 和 Jacobian，不存在 secant 在未收敛 $\rho$ 上读污染 $F$ 的结构。

**P7 全耦合 Newton** 是消除 H3 的正确方向（design.md §8.1 已规划）。

### 5.4 小切单元 (small-cut) 问题

| 项 | UP3D | López §3.5.5 | 备注 |
|---|---|---|---|
| 小切单元 | N/A（conforming 网格不存在切单元） | 需要全积分技术（§3.5.5） | UP3D 优势 |
| 条件数 | 由网格质量决定 | 小切单元导致病态条件数 | UP3D 无此问题 |

---

## 6. 数值格式：非线性求解策略

### 6.1 当前 UP3D（Picard + secant + continuation）

```
UP3D:  Mach continuation → secant(Γ) ↘ Picard(ρ) ↗ 交叉耦合
       (三层嵌套，H3 耦合失稳风险)
```

**代码验证** (`solve/continuation.py`, `solve/picard.py`)：

`continuation.py`：
- Mach levels: $M_0, M_0 + \Delta M, ..., M_{\mathrm{target}}$ ✅
- 每个 Mach level 调用 `solve_subsonic_lifting` ✅
- secant Γ root-find：per-station，`n_gamma_evals` 次 ✅
- `TRANSONIC_DEFAULTS`: evals=12, n_picard_eval=800 等 ✅

`picard.py::solve_subsonic_lifting`：
- 外层 density Picard：`n_picard_max` 次迭代 ✅
- 内层 Kutta secant：`max_kutta_updates` 次（或 `kutta_per_outer=1`）✅
- AMG 复用：`amg_rebuild_every` 次 ✅
- inexact forcing：`forcing` 参数 ✅
- 伪时间：`pseudo_dt`（全局质量集中）或 `damping_theta`（局部对角）✅
- 密度松弛：`omega_rho` ✅

**问题**：
1. 三层嵌套的迭代代价高（P4 transonic coarse: 10,464 Picard its）
2. secant-density 耦合失稳（H3）
3. `kutta_per_outer=1` 是权宜之计（避免 secant 过冲，但不治愈根因）

### 6.2 López（Newton-Raphson + load stepping）

```
López:  load step → Newton-Raphson(φ, ρ̃(φ)) → 收敛 → next step
       (单层 Newton，完整 Jacobian，无 secant)
```

| 项 | López | UP3D (P7 目标) | 一致性 |
|---|---|---|---|
| 非线性方法 | Newton-Raphson | 全耦合 Newton | ✅ |
| Jacobian | 完整（Term 1+2+3, Appendix B） | 策略 A（同） | ✅ |
| Γ 处理 | 隐式（least-squares 残差进入 Newton） | Γ 作为未知数（增广 Jacobian） | 方法不同，目标一致 |
| 全局化 | load stepping（Mach ramp） | Mach continuation + load stepping + damping_theta | ✅ |
| 线性求解 | GMRES + AMGCL | GMRES + PyAMG | ✅ |
| 收敛验证 | Table 4.9: 严格二次 | G7.1: 严格二次 | ✅ |

### 6.3 Load stepping 参数调度

**López 实际做法**（PDF 核对修正后）：

| Case | M∞ ramp | Mcrit | μc | 步数 |
|---|---|---|---|---|
| NACA (1) | 0.70→0.72 | 0.99 (固定) | 1.0 (固定) | 3 |
| NACA (4) | 0.70→0.75 | 0.90 (固定) | 1.1 (固定) | 6 |
| ONERA M6 | 0.50→0.84, 然后 μc 2.0→1.6 | 0.95 (固定) | 2.0→1.6 (到达后降) | 12 |

**UP3D 当前**：`TRANSONIC_DEFAULTS` 固定 `upwind_c=1.5, m_crit=0.95`，Mach continuation 用固定 0.05 步长。

**P7 需修正**：
- Mcrit 在 load stepping 过程中固定
- μc 调度模式：ramp 阶段高值，到达目标 M∞ 后降低
- 步长随难度调整

---

## 7. 数值实现：线性代数与并行

### 7.1 线性求解器

| 项 | UP3D | López/Kratos | 差异 |
|---|---|---|---|
| Picard 路径 | CG + PyAMG (SPD) | CG + AMGCL | ✅ 等价 |
| Newton 路径 (P7) | GMRES + PyAMG/ILU (非对称) | GMRES + AMGCL | ✅ 等价 |
| AMG setup | smoothed_aggregation_solver | AMGCL (Demidov 2020) | ✅ 算法同类 |
| AMG 复用 | 每 `amg_rebuild_every` 步 | Kratos strategy 管理 | ✅ |
| 预处理种子 | `np.random.seed(0)` 固定（bit-reproducible） | N/A | UP3D 特色 |

**代码验证** (`solve/linear.py`)：
- `apply_dirichlet()` — principal submatrix reduction（SPD 保持）✅
- `build_amg_preconditioner()` — 固定 RNG seed 保证 bit-reproducible ✅
- `solve_cg_amg()` — CG + AMG preconditioner ✅
- Newton 路径需要新增 `solve_gmres_amg()` — **P7 交付物** ⚠️

### 7.2 装配并行化

| 项 | UP3D | López/Kratos | 差异 |
|---|---|---|---|
| 策略 | 元素着色 + colored prange | Kratos OpenMP |  |
| 着色 | greedy，~8-12 colors (tets) | N/A | UP3D 自研 |
| 竞争安全 | 同色元素不共享节点 → prange 无竞争 | Kratos 保证 | ✅ |
| bit-deterministic | 是（着色序列固定） |  | UP3D 特色 |
| 矩阵装配 | elem_to_csr scatter map（一次性预计算） | Kratos builder/allocator | ✅ |

**代码验证** (`mesh/coloring.py`, `kernels/jacobian.py`)：
- greedy element coloring ✅
- `elem_to_csr` scatter map 预计算 ✅
- colored prange kernel（零分配）✅
- `PicardOperator` 封装了预计算 B_e, V_e, sparsity, coloring ✅

### 7.3 Numba 架构合规性

design.md §7 的 6 条硬规则：

| 规则 | 代码合规性 | 验证 |
|---|---|---|
| 1. SoA only（flat float64/int32 arrays） | ✅ | `PicardOperator` 持有 arrays，kernels 接收 arrays |
| 2. colored prange for scatter | ✅ | `assemble_matrix_colored`, `assemble_residual_colored` |
| 3. `@njit(cache=True, fastmath=True)` | ✅ | 所有 kernels |
| 4. 无热路径分配 | ✅ | 预分配 `_R`, `_rho_tilde`, `_upstream` |
| 5. 纯物理函数 | ✅ | `physics/isentropic.py` 全部 njit scalar |
| 6. SciPy/PyAMG 线性代数 | ✅ | `solve/linear.py` |

---

## 8. 数值实现：边界条件

### 8.1 壁面

| 项 | UP3D | López | 一致性 |
|---|---|---|---|
| BC 类型 | 自然 BC（do-nothing，零通量） | 同 | ✅ |
| 实现 | Galerkin 弱形式中省略面积分 | 同 | ✅ |
| 几何误差 | P1 flat facet vs 真实曲面 O(h) | embedded 不存在此问题 | UP3D 的 G1.6 |

**UP3D 的 G1.6 问题**：body-fitted P1 网格的 flat facet 法向与真实曲面法向偏差 O(h)，导致壁面 Cp 精度 ~1 阶。López 的 embedded approach 不存在此问题（几何不依赖壁面拟合）。design.md §5.1 已详细分析，结论是需曲面/等参单元。

### 8.2 远场

| 项 | UP3D | López | 一致性 |
|---|---|---|---|
| BC 类型 | Dirichlet | Dirichlet | ✅ |
| 均匀流 | $\varphi_\infty = U_\infty(x\cos\alpha + y\sin\alpha)$ | 同 | ✅ |
| 涡修正 | 2D point vortex + Prandtl-Glauert scaling | López 用horseshoe vortex | ⚠️ |
| 3D 扩展 | `farfield_spanwise_gamma=True`: Γ(z) 渐缩 | horseshoe vortex | 方法不同 |
| 域大小 | ~15 chords（含涡修正） | ~50 chords | UP3D 更小域 |

**差异分析**：

UP3D 的 2D point vortex + PG scaling 在准 2D（M0）情况完全正确。3D 扩展用 `spanwise_gamma=True` 的 Γ(z) 线性插值 + 翼尖渐缩到 0，这是 horseshoe vortex 的一阶近似。López 未详细描述远场处理（论文中用 $\varphi \to \varphi_\infty$ at far field，未提涡修正）。

P5 暴露的远场 2D 点涡分支射线伪影（翼尖外 M~5 cell 簇）已通过 `spanwise_gamma=True` 修复。

**代码验证** (`constraints/dirichlet.py`)：
- `freestream_phi()` — 均匀流 ✅
- `vortex_phi_2d()` — 2D 点涡 + PG scaling + branch cut on wake ✅
- `farfield_dirichlet()` — 组合均匀流 + 涡修正，支持 spanwise_gamma ✅
- branch cut 与 master-slave 一致（θ_w=0 上方 / 2π 下方）✅

### 8.3 对称面

| 项 | UP3D | López | 一致性 |
|---|---|---|---|
| BC 类型 | 自然 BC（零通量） | 同 | ✅ |

---

## 9. 数值实现：网格与预处理

### 9.1 网格类型

| 项 | UP3D | López | 差异 |
|---|---|---|---|
| 体积网格 | 非结构四面体 | 非结构三角/四面体 | ✅ |
| 尾流面 | conforming 内部面 | embedded（cut cells） | **架构差异** |
| 网格生成 | Gmsh + solver-side duplication | Gmsh + Kratos cut |  |
| 翼型 | NACA 0012 (sharp TE), ONERA M6 | NACA 0012 (sharp TE, Eq.4.2) | ✅ |

### 9.2 尾流切割预处理

**代码验证** (`mesh/wake_cut.py`)：
- `cut_wake()` — 节点复制 + 元素重指向 + 边界面重指向 ✅
- TE 节点复制（jump reaches wall）✅ — 有量化证据（P2: -0.27 cl error → 0.6012）
- 自由边节点不复制（pin Γ_tip=0）✅
- spanwise station 分组（(x,y) 位置）✅
- Kutta 探针选取（壁面 TE 一阶邻居）✅
- `assert_wake_topology()` — 5 条拓扑断言 ✅
- flood fill 分离 +/- 侧（不需要全局平面假设）✅

**对比 López**：
- López 不需要 `wake_cut.py`（embedded，不 conforming）
- López 的 §3.5.3 "nodes lying on the wake" 用容差距离处理——UP3D 不存在此问题（conforming）
- López 的 §3.5.4 "TE node treatment"——auxiliary DOF 作为体下方元素普通 DOF。UP3D 的 TE 节点复制后 slave 作为上方元素 DOF，机制不同但目标一致
- López 的 §3.5.5 "small-cut full integration"——UP3D 不存在此问题（conforming）

### 9.3 网格生成

**代码验证** (`meshgen/`)：
- `planar.py` — 2D 翼型网格生成
- `extrude.py` — 准 2D 单层拉伸
- `wing3d.py` — 3D 翼面网格生成

M0 spec: 单层拉伸准 2D（55.5k / 350.7k tets for coarse/medium）
M1 spec: 3D 翼面（ONERA M6）

---

## 10. 数值实现：后处理与验证

### 10.1 气动力计算

| 方法 | UP3D design.md §9 | López §2.4 | 一致性 |
|---|---|---|---|
| 近场 | $\int C_p \cdot \mathbf{n} \, dS$ over wall | 同 §2.4.1 | ✅ |
| Kutta-Joukowski | $c_l = 2\Gamma/(U_\infty c)$ | 同 §2.4.3 | ✅ |
| 远场 (Trefftz) | $\Gamma(s)$ → induced drag | 同 §2.4.2 | ✅ |
| 交叉验证 | 近场 vs KJ vs 远场（< 1%） | 同 | ✅ |

**代码验证** (`post/`)：
- `surface.py` — 壁面 Cp, 力/矩积分
- `section_cut.py` — 展向剖面
- `shock.py` — 激波检测
- `vtk_out.py` — VTK 导出

### 10.2 验证矩阵

| Case | UP3D P 状态 | López Ch.4 | 一致性 |
|---|---|---|---|
| NACA 0012 不可压 | P2 ✅ | §4.1 | ✅ |
| NACA 0012 亚声速 | P3 ✅ | §4.3 | ✅ |
| NACA 0012 跨声速 | P4 ✅ | §4.4 | ✅ |
| Korn 翼型 | ❌ | §4.5 | UP3D 未做 |
| 阻力发散 | ❌ | §4.6 | UP3D 未做 |
| 3D 矩形翼 | ❌ | §4.7 | UP3D 未做 |
| ONERA M6 | P5 ✅ (closed) | §4.8 | ✅ |
| NASA CRM | ❌ | §4.9 | UP3D 未做 |

---

## 11. P6/P7 实现路径技术细节

### 11.1 P6：可微人工密度通量

**两个独立缺陷**（design.md §3.1）：

| 缺陷 | 描述 | 修复 | 与 Newton 关系 |
|---|---|---|---|
| A (sawtooth) | 相邻超音速单元选择不同上游 → ρ̃ 携带 ~2h checkerboard | consistent flux (3.4): smooth face-neighbor weighting | 独立（Picard 中也存在） |
| B (不可微性) | max(0, μ, μ_up) + hard clamp 不可微 | smooth max_ε + 可选 smooth clamp | P7 前置 |

**López PDF 核对结论**：
- López 冻结 u(e)，只对 ρ_up, μ, ρ 求导（Appendix B B.3-B.8）
- P6 只需在**固定选择下可微**，无需可微化整数选择
- switching function 导数 ∂μ/∂φ 纳入（B.3, B.6）
- hard clamp 在收敛态非激活，光滑化为可选

**代码改动点**：

| 文件 | 改动 | 阶段 |
|---|---|---|
| `kernels/upwind.py` | smooth max_ε switching function + consistent flux weighting | P6 |
| `kernels/jacobian_newton.py` | **新建**：Term 2+3 装配（含上游耦合 B.4） | P7 |
| `kernels/elem_to_csr.py` | 扩展：上游 DOF scatter（stencil 宽一层） | P7 |
| `mesh/coloring.py` | 重新着色（宽 stencil → 更多颜色） | P7 |
| `physics/isentropic.py` | 可选：smooth clamp | P7 (可选) |

### 11.2 P7：全耦合 Newton

**增广 Jacobian**（design.md §8.1）：

$$\begin{pmatrix} T^T J T & \partial R_{\mathrm{red}}/\partial \Gamma \\ \partial F/\partial \varphi_{\mathrm{red}} & -I \end{pmatrix} \begin{pmatrix} \delta\varphi_{\mathrm{red}} \\ \delta\Gamma \end{pmatrix} = \begin{pmatrix} -T^T R \\ -F \end{pmatrix}$$

**已有的代码组件**：

| 组件 | 代码位置 | 状态 |
|---|---|---|
| $T^T J T$ (reduced Jacobian) | `wake.py::WakeConstraint.update_matrix()` | ✅ 已有 Picard 版（$T^T A T$），Newton 版需替换 $A \to J$ |
| $\partial R_{\mathrm{red}}/\partial \Gamma_j = -h_j$ | `wake.py::self._h` | ✅ 已预算 |
| $F_j = \mathrm{kutta\_targets}_j - \Gamma_j$ | `wake.py::kutta_targets()` | ✅ 可复用 |
| $\partial F_j/\partial \varphi_{\mathrm{red}}$ | 探针处 ±1/n_j 稀疏平均矩阵 | ⚠️ **需新建** |
| $\partial R_{\mathrm{red}}/\partial \Gamma_j$ (far-field column) | $-A_{\mathrm{coupling}} \cdot \partial \mathrm{vals}_{\mathrm{red}}/\partial \Gamma_j$ | ⚠️ **易遗漏**，design.md §8.1 明确标注 |

**关键风险**：远场列 $\partial R_{\mathrm{red}}/\partial \Gamma_j$ 的第二项（far-field vortex correction 对 Γ 的依赖）在 Picard 代码中被折叠进 RHS，Newton 路径必须显式提取。**如果遗漏此项，Newton 残差不会到机器精度**。

### 11.3 线性求解器

| 路径 | 当前 | P7 需要 | 改动 |
|---|---|---|---|
| Picard (SPD) | CG + PyAMG | 保留 | 无 |
| Newton (非对称) | N/A | GMRES + PyAMG/ILU | **新建** `solve/linear.py::solve_gmres_amg()` |
| AMG setup | `smoothed_aggregation_solver(A)` | `smoothed_aggregation_solver((A + A^T)/2)` (对称部分) | 修改 |
| Eisenstat-Walker | `forcing` 参数（简化版） | 完整 E-W schedule | 修改 |

---

## 12. 风险与未解决问题

### 12.1 高风险

| # | 风险 | 影响 | 缓解 | 状态 |
|---|---|---|---|---|
| 1 | **远场列遗漏** | Newton 不收敛到机器精度 | design.md §8.1 已标注；P7 实现时显式提取 | ⚠️ 待实现 |
| 2 | **stencil 宽一层的着色/内存** | 装配代价 +30%，着色需重算 | 预估内存增加 ~30%；重新着色一次性开销 | ⚠️ 待实现 |
| 3 | **GMRES + PyAMG 对非对称矩阵的效率** | GMRES 迭代数可能多 | PyAMG 对对称部分做 aggregation；ILU 作为备选 | ⚠️ 待验证 |
| 4 | **surface-Cp sawtooth (Defect A)** | V6 < 1% 精度目标 | P6 consistent flux | ⚠️ P6 待实现 |

### 12.2 中风险

| # | 风险 | 影响 | 缓解 | 状态 |
|---|---|---|---|---|
| 5 | **切换函数差异**（UP3D 缺 max(μ, μ_up)） | 激波后侧耗散不足 | P6 smooth max_ε 自动修正 | ⚠️ P6 |
| 6 | **G1.6 壁面 Cp 精度**（P1 flat facet O(h)） | V6 < 1% 阻塞 | 曲面/等参单元（独立 effort） | ⚠️ 已知 |
| 7 | **Python GMRES 性能** | medium mesh 可能 > 5min | 短期接受；中期迁 PETSc | ⚠️ P7 验证 |
| 8 | **H3 耦合失稳在 N-Γ 分裂 fallback 中仍存在** | fallback 不可靠 | 全耦合 Newton 是主路线；N-Γ 分裂仅 fallback | ✅ design.md 已定 |

### 12.3 低风险

| # | 风险 | 影响 | 缓解 | 状态 |
|---|---|---|---|---|
| 9 | **TE 奇异性 + Newton 不可微性** | 3D 翼尖 TE 可能发散 | López 用 hard clamp + 尖 TE 仍二次收敛；N0 光滑化可选 | ✅ PDF 核对 |
| 10 | **sawtooth 与 Newton 可微性混淆** | P6 scope 误判 | design.md §3.1 已分离两个缺陷 | ✅ 已修正 |

---

## 13. 总结：实现路径合理性评估

### 13.1 核心判断

| 维度 | 评估 | 理由 |
|---|---|---|
| **物理模型** | ✅ 正确 | 全速势方程 + 等熵密度 + 极限速度截断，与 López 完全一致 |
| **空间离散** | ✅ 正确 | Galerkin P1 + 预计算 B_e/V_e + colored prange，实现高效且正确 |
| **人工可压缩性** | ⚠️ 基本正确 | 迎风密度公式正确；切换函数缺 max(μ, μ_up)；P6 将修正 |
| **尾流处理** | ✅ 优秀 | master-slave 消元比 López 的 least-squares 更简洁、更精确，对 Newton 更友好 |
| **非线性求解（当前）** | ⚠️ 临时方案 | Picard + secant 三层嵌套在跨声速有 H3 耦合失稳；P7 全耦合 Newton 是正确方向 |
| **非线性求解（P7 目标）** | ✅ 正确 | 全耦合 Newton + Γ 作为未知数 + 增广 Jacobian，与 López 目标一致 |
| **线性代数** | ✅ 合理 | CG+AMG (Picard) / GMRES+AMG (Newton)，PyAMG 与 AMGCL 算法同类 |
| **并行** | ✅ 合理 | colored prange + elem_to_csr，零分配热路径，bit-deterministic |
| **边界条件** | ✅ 正确 | 远场 Dirichlet + 涡修正 + PG scaling；壁面 natural BC；对称面 natural |
| **网格预处理** | ✅ 优秀 | wake_cut.py 拓扑断言完善，TE 节点复制有量化证据，自由边处理正确 |
| **后处理** | ✅ 正确 | 近场/KJ/远场三路交叉验证 |

### 13.2 与 López 的关键架构差异

| 差异 | UP3D 选择 | López 选择 | 评估 |
|---|---|---|---|
| 尾流 | conforming + master-slave | embedded + least-squares | UP3D 更简洁精确；López 更灵活（aeroelastic） |
| 非线性 | Picard → Newton (P7) | Newton from start | UP3D 的 Picard 是合理原型路径 |
| 网格 | body-fitted | embedded | UP3D 有 G1.6 壁面精度问题；López 有 small-cut 问题 |
| 语言 | Python + Numba | C++ (Kratos) | UP3D 适合快速迭代；性能 3-5× overhead |

### 13.3 P6/P7 实现路径技术清单

**P6（可微通量）**：
1. ✅ design.md §3.1-3.2 已定义两个独立缺陷和修复方案
2. ✅ López Appendix B 已验证策略 A（完整 Jacobian 含上游耦合）
3. ✅ 上游选择保持 integer-walk（López 冻结 u(e)）
4. ⚠️ 需实现 smooth max_ε switching function
5. ⚠️ 需实现 consistent flux weighting (3.4)
6. ⚠️ 需验证 G4.2 bitwise no-op

**P7（全耦合 Newton）**：
1. ✅ design.md §8.1 已定义增广 Jacobian
2. ✅ `wake.py` 的 $h_j$ 和 `kutta_targets()` 可复用
3. ⚠️ 需新建 $\partial F/\partial \varphi_{\mathrm{red}}$ 稀疏矩阵
4. ⚠️ 需显式提取远场列 $\partial R_{\mathrm{red}}/\partial \Gamma_j$（far-field vortex correction 对 Γ 的依赖）
5. ⚠️ 需新建 GMRES + AMG 线性求解路径
6. ⚠️ 需重新着色（宽 stencil）
7. ⚠️ 需实现 Eisenstat-Walker inexact Newton schedule
8. ⚠️ 需实现 load stepping 参数调度（Mcrit 固定，μc 先高后降）

### 13.4 深度开发前的 action items

1. **P6 实现**：按 design.md §3.2 (3.4) 实现 consistent flux + smooth max_ε
2. **P7 实现**：按 design.md §8.1 实现全耦合 Newton 增广 Jacobian
3. **关键验证**：远场列的完整性（action item #4 in §11.2）
4. **性能验证**：GMRES + PyAMG 对非对称 Jacobian 的效率（medium mesh）
5. **回归验证**：G4.2 bitwise no-op（P6 后）、G7.1 严格二次收敛（P7 后）

---

## 参考文献

- López Canalejo, I.P. "A Finite-Element Transonic Potential Flow Solver with an Embedded Wake Approach for Aircraft Conceptual Design." Dissertation, TUM, 2021.
  - Ch.2: 全速势方程推导
  - Ch.3 §3.2: FE 离散 + Jacobian Eq.(3.8)-(3.16)
  - Ch.3 §3.3: 人工可压缩性 Eq.(3.19)-(3.27)
  - Ch.3 §3.4: 极限速度 / 密度截断
  - Ch.3 §3.5: 嵌入尾流 2D/3D + TE + small-cut
  - Ch.4: 验证（NACA 0012, Korn, ONERA M6, NASA CRM）
  - Appendix B: 完整敏感度推导 Eq.(B.1)-(B.17)
- Kratos Multiphysics — https://github.com/KratosMultiphysics/Kratos
- UP3D `docs/design.md` §1-§13
- UP3D `docs/roadmap.md` P1-P8
- Demidov, D. "AMGCL — A C++ library for efficient solution of large sparse linear systems." *Software Impacts* 6 (2020) 100037.
- Cai, X.-C., Keyes, D.E., Young, D.P. "Parallel Newton-Krylov-Schwarz Algorithms for the Transonic Full Potential Equation." *SIAM J. Sci. Comput.* 19 (1998) 246–265.
- Eisenstat, S.C., Walker, H.F. "Choosing the forcing terms in an inexact Newton method." *SIAM J. Sci. Comput.* 17 (1996) 16–32.
- Hafez, M., South, J., Murman, E. "Artificial Compressibility Methods for Numerical Solutions of Transonic Full Potential Equation." *AIAA J.* 17(8) (1979) 838–844.
- Holst, T.L. "Transonic flow computations using nonlinear potential methods." *Progress in Aerospace Sciences* 36 (2000) 1–61.
