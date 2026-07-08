# ❌ 已弃用 — Newton 求解器实现策略：从 Picard 到 Newton-Krylov

> **本文档已于 2026-07-08 16:46 弃用**，替代文档为 `20260708_1646_newton_solver_strategy.md`
>
> **弃用原因**：基于 López 2021 博士论文 Appendix B 的交叉验证，发现本文档的三个核心判断有误：
> 1. López 实际用策略 A（完整 Jacobian 含上游耦合），不是策略 B——Appendix B Eq.(B.3)-(B.6) 明确推导了 $\partial\tilde{\rho}/\partial\varphi^{\mathrm{up}}$
> 2. P6 可微通量确实是 P7 的硬前置依赖（策略 A 需要）
> 3. N-$\Gamma$ 分裂保留 secant 外层，H3 耦合失稳风险未消除——应改为全耦合 Newton
>
> 以下内容保留供历史参考，不再作为实现指导。

---

# Newton 求解器实现策略：从 Picard 到 Newton-Krylov
> Design note for UP3D / pyFP3D
> Date: 2026-07-07 23:13（updated 2026-07-08 01:05 — P6 前置依赖修正 + Jacobian 策略分析）
> Code baseline: HEAD e3d0386
> Prerequisite: P4 已关闭。**P6（可微人工密度通量）不是 P7 的硬前置依赖**——取决于 Jacobian 构建策略，见 §1.5
> References:
> - design.md §6.3 (Newton Jacobian), §8 (Nonlinear solution strategy)
> - roadmap.md P7 (deliverables + gates，原 P6)
> - López et al. 2022 (CMAME 114244) — Newton-Raphson + AMGCL
> - Cai, Keyes, Young 1996-1998 (SIAM J. Sci. Comput.) — Newton-Krylov-Schwarz for transonic FPE

## 1. 为什么需要这个文档
roadmap.md **P7**（原 P6，2026-07-07 重编号）只写了 deliverables 和 gate：
```
Deliverables: exact Newton Jacobian (6.3) with widened sparsity, pseudo-transient continuation, GMRES + ILU path, AMG setup reuse, Eisenstat–Walker inexact-solve schedule, profiling report. Consumes the P6 differentiable flux.
Gates:
  G7.1 Newton terminal quadratic convergence on G4.1 case
  G7.2 ONERA M6 medium mesh < 5 min single node, end to end
  G7.3 full regression suite runtime < 10 min (CI budget)
```
这不够指导实现。本文档补充：算法选择、模块改动、分阶段路线、线性求解器策略、与现有代码的对接。

> **注（2026-07-08 修正）**：roadmap 中"Consumes the P6 differentiable flux"过于保守。López 2022 的实践证明，不含上游耦合的近似 Jacobian 也能让 Newton 在 13 次迭代内收敛。P6 和 P7 可以并行或反序。详见 §1.5。

## 1.5 P6 可微通量与 P7 的关系

### 1.5.1 问题根源

当前 P4 的人工密度（`kernels/upwind.py`）有两处不可微操作：

1. **Integer-walk 上游搜索**（`upstream_elements`）：通过面邻居 `argmin` 选出"最负流向位移"的邻居——离散选择，$\nabla\varphi$ 方向微小变化可导致 $u(e)$ 整数跳变
2. **激波点算子**：$\nu = \max(\nu_e, \nu_u)$——`max` 在 $\nu_e = \nu_u$ 处不可微

实测效果：supersonic pocket 的表面 Cp 有 $\approx 2h$ 的 odd-even 锯齿振荡。

### 1.5.2 三种 Jacobian 策略

Newton Jacobian 的完整形式（design.md (6.3)）有三项：

$$J_{ij} = \underbrace{\sum_e V_e \, \tilde{\rho}_e \, \nabla N_i \cdot \nabla N_k}_{\text{项 1: Picard（对称）}} + \underbrace{\sum_e V_e \, \frac{\partial \tilde{\rho}_e}{\partial q_e^2} \cdot 2(\nabla\varphi_e \cdot \nabla N_k)(\nabla\varphi_e \cdot \nabla N_i)}_{\text{项 2: 局部密度导数（非对称）}} + \underbrace{\sum_e V_e \, \frac{\partial \tilde{\rho}_e}{\partial \rho_{u(e)}} \cdot \frac{\partial \rho_{u(e)}}{\partial \varphi} \cdot \ldots}_{\text{项 3: 上游耦合（非对称，stencil 宽一层）}}$$

| 策略 | 包含项 | 需要可微通量？ | 收敛阶 | 文献先例 |
|------|--------|--------------|--------|---------|
| **A：精确 Jacobian** | 1 + 2 + 3 | **是**——项 3 的 $\partial\tilde{\rho}_e/\partial\rho_{u(e)}$ 要求 $u(e)$ 可微 | 严格二次 | design.md (6.3) 完整版；Jameson-Caughey（结构网格，上游方向固定） |
| **B：近似 Jacobian** | 1 + 2 | **否**——项 2 的 $\max$ 用光滑近似（softplus 或 $\max_\varepsilon$）即可 | 超线性 | **López 2022 实际做法** |
| **C：JFNK** | 有限差分 matvec | 否 | 超线性 | Cai-Keyes-Young NKS |

### 1.5.3 López 2022 的实际做法

López 论文的 Jacobian（CMAME Eq.26 区域）**明确包含了密度导数** $\partial\rho/\partial|v|^2$：

$$J_{ij}^e = \int_{\Omega_e} \rho \nabla N_i \cdot \nabla N_j + 2\frac{\partial\rho}{\partial|v|^2} (\nabla N_j \cdot \nabla\Phi)(\nabla N_i \cdot \nabla\Phi) \, d\Omega$$

但**不包含**上游耦合项 $\partial\tilde{\rho}_e/\partial\rho_{u(e)}$——论文中没有讨论上游元素的导数。这对应策略 B。

他们只在**尾流面 BC 的最小二乘泛函**（Eq.31）中忽略了密度导数（理由是 "small velocity perturbations"），主方程的 Jacobian 包含局部密度导数。

López 的数据：72k DOF，13 Newton iterations，7s（4 核 3.6GHz）。证明策略 B 足够有效。

### 1.5.4 结论：P6 不是 P7 的硬前置依赖

| 问题 | 回答 |
|------|------|
| P6 对 P7 Newton 必要吗？ | **不必要**——P7 走策略 B（López 路线），不需要可微通量 |
| P6 的价值是什么？ | 消除 Cp 锯齿振荡（精度问题），与 P7（性能问题）独立 |
| P6 和 P7 的关系？ | 可并行或反序。P6 只在追求策略 A（完整精确 Jacobian）时才是 P7 前置 |
| roadmap 的串行依赖合理吗？ | 过于保守，应修正为"P6 gates P5's section-Cp acceptance; P7 can proceed independently with strategy B" |

### 1.5.5 对 P7 实现的影响

P7 走策略 B（近似 Jacobian），具体做法：

- **项 1**（Picard 部分）：复用当前 `assemble_matrix_data_colored`
- **项 2**（局部密度导数）：新增 kernel，per-element 计算 $\frac{\partial\tilde{\rho}_e}{\partial q_e^2} \cdot 2(\nabla\varphi \cdot \nabla N_k)(\nabla\varphi \cdot \nabla N_i) \cdot V_e$
  - $\max(\nu_e, \nu_u)$ 用光滑近似 $\max_\varepsilon(a, b) = \frac{a+b}{2} + \sqrt{\frac{(a-b)^2}{4} + \varepsilon}$（$\varepsilon \sim 10^{-8}$）
  - 上游元素选择 $u(e)$ 仍用 integer-walk（不可微），但 $\partial\tilde{\rho}/\partial q^2$ 只依赖 $\nu_e$ 和 $\rho_e$，不涉及 $u(e)$ 的导数
- **项 3**（上游耦合）：**不做**——López 也没做
- **stencil 不变**：因为不含项 3，sparsity pattern 与 Picard 矩阵相同，**不需要重新着色或宽 stencil**

> 这意味着 §4.2 和 §4.3 中关于"stencil 宽一层"和"重新着色"的描述**不适用于策略 B**。N1 阶段的实现比原计划更简单。

## 2. Picard vs Newton 本质对比
### 2.1 算法结构
**Picard（当前）**：
$$\varphi_k \to \nabla \varphi,\; q^2,\; \rho(\varphi_k) \to \text{冻结 } \tilde{\rho} \to A(\tilde{\rho})\, \delta\varphi = -R(\varphi_k) \to \varphi_{k+1} = \varphi_k + \omega \cdot \delta\varphi$$
- $A(\tilde{\rho})$ 对称 SPD（$\tilde{\rho}$ 作为标量冻结）
- 线性收敛：$\|e_{k+1}\| \leq C \cdot \|e_k\|$，$C < 1$ 但接近 1（强激波时）
- CG + AMG 求解
**Newton**：
$$\varphi_k \to \nabla \varphi,\; q^2,\; \rho(\varphi_k) \to J(\varphi_k)\, \delta\varphi = -R(\varphi_k) \to \varphi_{k+1} = \varphi_k + \delta\varphi$$
- $J$ 非对称（密度耦合项，策略 B）
- 超线性收敛（非严格二次，因不含上游耦合）
- GMRES + 预处理求解
### 2.2 Jacobian 结构
策略 B 的 Jacobian（推荐路线）：
$$J_{ij} = \sum_e V_e \left[ \tilde{\rho}_e \, \nabla N_i \cdot \nabla N_k + \frac{\partial \tilde{\rho}_e}{\partial q^2_e} \cdot 2(\nabla \varphi_e \cdot \nabla N_k)(\nabla \varphi_e \cdot \nabla N_i) \right]$$
第一项 = Picard 矩阵（对称），sparsity pattern 不变。
第二项 = 密度-速度耦合（**非对称**），在超声速区激活（$\nu > 0$ 时 $\partial \tilde{\rho} / \partial q^2 \neq 0$），sparsity pattern 与第一项相同（不宽 stencil）。

> **与 design.md (6.3) 的差异**：design.md (6.3) 描述的是策略 A（含上游耦合项 3）。本文档推荐 P7 走策略 B（不含项 3），与 López 实践一致。design.md (6.3) 需要后续更新，标注策略 B 作为 P7 的推荐路线。

### 2.3 当前 Picard 性能基线
| 场景                                                         | 网格                       | 外层迭代                                | 时间   |
| ------------------------------------------------------------ | -------------------------- | --------------------------------------- | ------ |
| P3 subsonic                                                  | NACA medium (20,888 nodes) | 15 Picard                               | 35s    |
| P4 transonic coarse                                          | NACA coarse (5,610 nodes)  | 10,464 Picard (含 Mach × $\Gamma$ 嵌套) | 174s   |
| P4 transonic medium                                          | NACA medium (20,888 nodes) | ~12,931 Picard                          | 16m39s |
| P4 的 10,000+ 迭代来自三层嵌套：Mach continuation (3 levels) × $\Gamma$ secant ($\leq 12$ evals/level) × density Picard ($\leq 800$ its/eval)。Newton 去掉 $\Gamma$ secant 层后预期 20-40 次 Newton 步。 |                            |                                         |        |
### 2.4 参考文献的求解器选择
| 来源                                                         | 非线性求解器         | 线性求解器         | 预处理          | 并行        | Jacobian 策略 |
| ------------------------------------------------------------ | -------------------- | ------------------ | --------------- | ----------- | ------------- |
| López 2022 (CMAME)                                           | Newton-Raphson       | AMGCL (AMG)        | AMG             | 4 核 OpenMP | **B（近似）** |
| Cai-Keyes-Young 1996 (NKS)                                   | Inexact Newton       | GMRES              | Schwarz 域分解  | MPI 并行    | C（JFNK）     |
| Holst 1999 (NASA 综述)                                       | Newton（TRANAIR 等） | 多重网格 / GMRES   | AMG / ILU       | —           | A 或 B        |
| UP3D design.md §8                                            | Newton（P7）         | GMRES + ILU 或 AMG | AMG setup reuse | 未规划      | **B（推荐）** |
| López 的数据：72k DOF，13 次 Newton 迭代，7s（4 核 3.6GHz）。 |                      |                    |                 |             |               |

## 3. 实现策略
### 3.1 核心决策
| 决策          | 选择                                      | 理由                                                      |
| ------------- | ----------------------------------------- | --------------------------------------------------------- |
| 非线性方法    | Inexact Newton                            | 允许早期步线性求解不精确（Eisenstat-Walker），节省计算    |
| Jacobian 策略 | **B（近似 Jacobian）**                    | 不需要 P6 可微通量；与 López 实践一致；stencil 不变       |
| 线性求解器    | GMRES                                     | Jacobian 非对称，CG 不可用                                |
| 预处理        | AMG (PyAMG) + 可选 ILU                    | AMG 对椭圆算子最优；ILU 作为非对称部分的补充              |
| $\Gamma$ 处理 | 保留 conforming 尾流 + master-slave 消元  | 不依赖 level-set 方案 B；Newton 对 $\Gamma$ 的处理见 §3.5 |
| 伪时间        | 保留 damping_theta 作为 Newton 全局化策略 | Newton 从远场初值可能发散；damping_theta 提供稳定下降     |
| 并行          | 短期 Numba prange（装配）；中期 petsc4py  | 短期改动最小                                              |
### 3.2 不做什么
- **不做无矩阵 Newton**（`scipy.optimize.newton_krylov`）：不利用已有的显式 Jacobian 装配，效率低
- **不做策略 A（精确 Jacobian 含上游耦合）**：需要 P6 可微通量 + stencil 宽一层 + 重新着色，远期目标
- **不一次性替换 Picard**：Newton 和 Picard 并存，Newton 用于 transonic，Picard 保留作为 subsonic 的快速路径和 fallback
- **不做 Schwarz 域分解**（短期）：实现复杂，Numba 不友好；P7 先用单机 AMG，后续再考虑 petsc4py
### 3.3 $\Gamma$ 在 Newton 中的处理
当前 conforming 尾流的 $\Gamma$ 是通过 master-slave 消元 + secant 外层循环处理的。Newton 有两种选择：
| 方案                                                         | 描述                                                         | 优点                                  | 缺点                                                     |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------- | -------------------------------------------------------- |
| **N-$\Gamma$ 分裂**（推荐）                                  | Newton 只对 $\varphi$ 求 Jacobian；$\Gamma$ 仍用 secant 外层 | 改动最小；$\Gamma$ 的 Jacobian 不需要 | 保留 $\Gamma$ 外层循环（但 Newton 步数少，外层迭代也少） |
| 全耦合 Newton                                                | $\Gamma$ 也作为未知数，Jacobian 包含 $\partial R / \partial \Gamma$ | 完全消除外层循环                      | 需要推导 $\Gamma$ 的 Jacobian；尾流约束的线性化复杂      |
| **推荐 N-$\Gamma$ 分裂**：先用 Newton 替代密度 Picard 内层，保留 $\Gamma$ secant 外层。三层嵌套变成两层（Mach × Newton），已经大幅减少迭代数。全耦合 Newton 留作后续优化。 |                                                              |                                       |                                                          |
### 3.4 Inexact Newton + Eisenstat-Walker
Eisenstat-Walker 策略：第 $k$ 步的线性容差 $\eta_k$ 自适应：
$$\eta_k = \max\left(\eta_{\min}, \min\left(\eta_{\max}, \gamma \cdot \frac{\|R(\varphi_k)\|^\alpha}{\|R(\varphi_{k-1})\|^{\alpha-1}}\right)\right)$$
典型参数：$\eta_{\min} = 10^{-4}$, $\eta_{\max} = 0.9$, $\gamma = 0.9$, $\alpha = 1$（线性）或 $2$（二次区）。
效果：
- 早期 Newton 步：线性容差 $10^{-2} \sim 10^{-3}$，GMRES 5-10 次迭代
- 终端 Newton 步：线性容差 $10^{-6} \sim 10^{-8}$，GMRES 20-50 次迭代
- 总体减少 50-70% 的线性求解开销
当前代码的 `forcing` 参数已经实现了 Eisenstat-Walker 的简化版（`solve_subsonic_lifting` 的 `forcing=0.01`），但 P4 transonic path 设为 0.0（精确求解每步）。Newton 路径需要打开它。
### 3.5 伪时间全局化
纯 Newton 从远场初值在 transonic 下会发散。全局化策略：
| 策略                                                         | 描述                                                         | 当前状态                                             |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ---------------------------------------------------- |
| Mach continuation                                            | 从 $M_{0.70}$ 收敛解重启 $M_{0.80}$                          | ✅ 已实现                                             |
| damping_theta                                                | $J + \theta \cdot \mathrm{diag}(J)$ 提供阻尼                 | ✅ 已实现（Picard path）                              |
| 伪时间步进                                                   | $J + \mathrm{diag}(m/\Delta \tau)$                           | ⚠️ 已有 `pseudo_dt`（Picard path），Newton 路径可复用 |
| 线搜索                                                       | $\varphi_{k+1} = \varphi_k + \lambda \cdot \delta \varphi$，$\lambda \in (0,1]$ 使 $\|R\|$ 下降 | ❌ 未实现，可选                                       |
| **推荐**：Mach continuation + damping_theta 作为全局化。Newton 从 Mach continuation 的上一级收敛解重启，damping_theta 在前几步提供阻尼，收敛后自动退化为纯 Newton（$\theta \cdot \mathrm{diag}(J)$ 相对 $J$ 趋于零）。 |                                                              |                                                      |

## 4. 模块改动
### 4.1 新增模块
```
pyfp3d/
├── kernels/
│   ├── jacobian_newton.py      # 新增：Newton Jacobian 装配（策略 B：项 1 + 项 2）
│   └── jacobian.py             # 不变（Picard 路径保留）
├── solve/
│   ├── newton.py               # 新增：Newton 驱动器
│   ├── picard.py               # 不变（Picard 路径保留）
│   └── linear.py               # 扩展：新增 GMRES + AMG 求解
```
### 4.2 改动清单（策略 B——stencil 不变）

> **注（2026-07-08 修正）**：策略 B 不含上游耦合项，sparsity pattern 与 Picard 矩阵相同。原方案中"stencil 宽一层"和"重新着色"只适用于策略 A，已删除。

| 文件                         | 改动                                                         | 改动量 |
| ---------------------------- | ------------------------------------------------------------ | ------ |
| `kernels/jacobian_newton.py` | **新建**：Newton Jacobian 装配 kernel（项 1 + 项 2，sparsity 同 Picard） | 中     |
| `solve/newton.py`            | **新建**：Newton 驱动器（Inexact Newton + Eisenstat-Walker + Mach continuation + $\Gamma$ secant） | 大     |
| `solve/linear.py`            | **扩展**：新增 `solve_gmres_amg()` 函数                      | 小     |
| `solve/continuation.py`      | **修改**：`solve_transonic_lifting` 新增 `solver="newton"` 选项 | 中     |
| `mesh/coloring.py`           | **不变**（策略 B stencil 不变，无需重新着色）                | —      |
| `kernels/jacobian.py`        | **不变**（复用现有 sparsity pattern 和 elem_to_csr）         | —      |
| `solve/picard.py`            | **不变**                                                     | —      |
### 4.3 Newton Jacobian 装配
策略 B 的 Jacobian 包含两项，sparsity pattern 与 Picard 矩阵相同：
1. **项 1**（Picard 部分）：复用当前 `assemble_matrix_data_colored`
2. **项 2**（局部密度导数）：新增 kernel，per-element 计算 $\frac{\partial\tilde{\rho}_e}{\partial q_e^2} \cdot 2(\nabla\varphi \cdot \nabla N_k)(\nabla\varphi \cdot \nabla N_i) \cdot V_e$

$\max(\nu_e, \nu_u)$ 的光滑近似：
$$\max_\varepsilon(a, b) = \frac{a+b}{2} + \sqrt{\frac{(a-b)^2}{4} + \varepsilon}, \quad \varepsilon \sim 10^{-8}$$

在 $\nu = 0$（亚临界）时，$\partial\tilde{\rho}/\partial q^2 = 0$ exactly（项 2 消失），Jacobian = Picard 矩阵，G4.2 bitwise no-op 保持。

### 4.4 GMRES + AMG 求解
```python
# solve/linear.py 新增
def solve_gmres_amg(A, b, rtol=1e-6, maxiter=500, restart=30):
    """GMRES with AMG preconditioner for nonsymmetric Newton Jacobian."""
    _, M = build_amg_preconditioner(A)  # 复用现有 PyAMG setup
    x, info = spla.gmres(A, b, M=M, rtol=rtol, maxiter=maxiter, restart=restart)
    if info != 0:
        raise RuntimeError(f"GMRES did not converge (info={info})")
    return x
```
**预处理策略**：
- $J$ 非对称，但 PyAMG 的 `smoothed_aggregation_solver` 对非对称矩阵也工作（作为预处理子，取对称部分 $J + J^T$ 做 AMG setup）
- 如果 AMG 预处理效果不够，加一层 ILU(0) 作为后处理
- AMG setup 复用：Newton 步之间 Jacobian 变化不大，AMG hierarchy 可复用 5-10 步再重建

## 5. 分阶段实现路线
### Phase N1：Newton Jacobian 装配 + 验证
**目标**：实现并验证 Newton Jacobian（策略 B）的正确性。
**交付物**：
- `kernels/jacobian_newton.py`：Newton Jacobian 装配 kernel（项 1 + 项 2）
- 单元测试：有限差分验证 $J \cdot \mathbf{v} \approx (R(\varphi + \varepsilon \mathbf{v}) - R(\varphi)) / \varepsilon$
**验证标准**：
- 有限差分 Jacobian 与解析 Jacobian 一致（相对误差 $< 10^{-6}$）
- 亚临界流（$\nu = 0$）：Newton Jacobian = Picard 矩阵 + 零（项 2 在 $\nu = 0$ 时消失）
- 与 G4.2 bitwise no-op 一致性
**不做**：Newton 驱动器、线性求解。
### Phase N2：GMRES + AMG 线性求解
**目标**：实现非对称线性求解路径。
**交付物**：
- `solve/linear.py` 新增 `solve_gmres_amg()`
- 验证：在 N1 的 Jacobian 上求解 $J \cdot \delta \varphi = -R$，$\delta \varphi$ 方向正确
**验证标准**：
- GMRES 收敛（20-50 次迭代，给定好的预处理）
- AMG setup 复用 5 步后重建，不影响收敛
- 与 CG+AMG 在对称情形（$\nu = 0$）下结果一致
### Phase N3：Newton 驱动器 + 亚声速验证
**目标**：实现 Inexact Newton 驱动器，在亚声速上验证。
**交付物**：
- `solve/newton.py`：Newton 驱动器
  - Eisenstat-Walker 自适应线性容差
  - damping_theta 全局化
  - Mach continuation 支持
  - N-$\Gamma$ 分裂（$\Gamma$ secant 外层保留）
- `solve/continuation.py` 新增 `solver="newton"` 选项
**验证标准**：
- V3 NACA0012 $M_\infty = 0.5$, $\alpha = 2°$：$c_l$ 与 P3 一致（$< 0.5\%$ 差异）
- Newton 收敛在 5-10 步（亚声速，非线性弱）
- $\|R\|$ 单调下降
- G4.2 bitwise no-op 仍然成立
### Phase N4：跨声速验证
**目标**：在 transonic 上验证 Newton，与 P4 对比。
**交付物**：
- `solve/continuation.py` 的 `solve_transonic_lifting` 默认走 Newton 路径
- Newton + Mach continuation + damping_theta 全局化
**验证标准**：
- V4 NACA0012 $M_\infty = 0.80$, $\alpha = 1.25°$ coarse：激波 $x/c \approx 0.60$（与 P4 一致）
- **G7.1**：Newton 终端超线性收敛——$\|e_{k+1}\| / \|e_k\|^{1.5} \to$ 减小趋势
- 总 Newton 步 $\leq 30$（vs P4 的 10,464 Picard its）
- coarse 求解时间 $< 30$s（vs P4 的 174s）

> **注**：策略 B 的收敛阶是超线性而非严格二次，因为忽略了上游耦合项。G7.1 的验收标准应相应调整为"超线性收敛"而非"二次收敛"。如果后续加上 P6 + 策略 A，可追求严格二次。
### Phase N5：ONERA M6 + 性能 gate
**目标**：在 ONERA M6 上验证，达到 G7.2 性能目标。
**交付物**：
- ONERA M6 medium mesh 求解
- 性能 profiling 报告（装配 / 线性求解 / 后处理时间分解）
**验证标准**：
- **G7.2**：ONERA M6 medium mesh $< 5$ min single node（当前 P4 估计 $> 30$ min）
- **G7.3**：full regression suite $< 10$ min（当前 $\sim 5$ min，Newton 不应拖慢 subsonic 路径）
- P5 V5 gate（$\lambda$激波拓扑、剖面 Cp）通过

## 6. 预期性能
### 6.1 迭代次数对比
| 场景                | Picard 迭代             | Newton 预期 | 每步代价    | 总代价比                            |
| ------------------- | ----------------------- | ----------- | ----------- | ----------------------------------- |
| P3 subsonic         | 15                      | 5-8         | 2-3× Picard | $\sim 1$-2×（亚声速 Newton 收益小） |
| P4 transonic coarse | 10,464                  | 20-30       | 2-3× Picard | **10-30× 更快**                     |
| P4 transonic medium | $\sim 12,931$           | 25-40       | 2-3× Picard | **5-15× 更快**                      |
| ONERA M6 medium     | 未知（估计 $> 30,000$） | 30-50       | 2-3×        | **10-30× 更快**                     |

> **注（2026-07-08 修正）**：每步代价比从 3-5× 降为 2-3×，因为策略 B 不含上游耦合项，stencil 不变，装配代价与 Picard 相当。线性求解用 GMRES（非对称）比 CG 贵约 2×，但 AMG 预处理后迭代数相近。
### 6.2 时间预估
| 场景            | 当前时间 | Newton 预估 | 依据                                                     |
| --------------- | -------- | ----------- | -------------------------------------------------------- |
| P4 coarse       | 174s     | 5-15s       | 30 Newton steps × (装配 $\sim 0.5$s + GMRES $\sim 0.3$s) |
| P4 medium       | 16m39s   | 1-3min      | 40 steps × (装配 $\sim 2$s + GMRES $\sim 1.5$s)          |
| ONERA M6 medium | 未知     | 2-5min      | López 72k DOF 7s × 规模系数 × Python overhead            |
| 目标 G7.2       | —        | $< 5$ min   | ✅ 可达                                                   |
### 6.3 López 数据校准
López CMAME：72k DOF，13 Newton iterations，7s（4 核 3.6GHz，C++ + AMGCL）。
UP3D 的差异：
- Python + Numba（装配快，但 SciPy GMRES 比 C++ 慢 2-5×）
- 单线程线性求解（短期）
- 规模系数：ONERA M6 medium 63k nodes $\approx$ López 的 72k DOF
预估：López 7s × (Python overhead 3-5×) × (单线程 2×) $\approx$ **42-70s**。加上 Mach continuation 和 $\Gamma$ secant 外层，$\sim 2$-5min。在 G7.2 的 5min 目标内。

## 7. 开源库选择
### 7.1 短期（P7 初始版）
| 库                                                           | 用途           | 接口                                                |
| ------------------------------------------------------------ | -------------- | --------------------------------------------------- |
| `scipy.sparse.linalg.gmres`                                  | 非对称线性求解 | `gmres(A, b, M=M, rtol=η, maxiter=500, restart=30)` |
| `pyamg.smoothed_aggregation_solver`                          | AMG 预处理     | 当前已在用，改为预处理非对称 $J$                    |
| `scipy.sparse.linalg.spilu`                                  | ILU 后备预处理 | `spilu(J.tocsc(), drop_tol=1e-4)`                   |
| **优点**：零安装成本（已在依赖中），接口与当前 CG 调用几乎一样。 |                |                                                     |
| **缺点**：单线程；SciPy GMRES 比 PETSc 慢 2-5×。             |                |                                                     |
### 7.2 中期（性能优化）
| 库                                                           | 用途                           | 优势                              |
| ------------------------------------------------------------ | ------------------------------ | --------------------------------- |
| `petsc4py`                                                   | 并行 GMRES + 并行 AMG + 域分解 | MPI 并行，多核扩展性好            |
| `amgcl` (pybind11)                                           | López 用的库                   | C++ 高性能 AMG，Python 绑定需自建 |
| **迁移策略**：`solve/linear.py` 的接口保持不变（`solve_gmres_amg(A, b, ...)`），内部实现从 SciPy 切换到 PETSc。 |                                |                                   |
### 7.3 不推荐
| 库                             | 原因                                                    |
| ------------------------------ | ------------------------------------------------------- |
| `scipy.optimize.newton_krylov` | 无矩阵 Newton-Krylov（策略 C），不利用显式 Jacobian，效率低 |
| `scipy.sparse.linalg.lgmres`   | LGMRES 适合无预处理场景，有 AMG 预处理时标准 GMRES 更好 |
| `scipy.sparse.linalg.bicgstab` | BiCGStab 不保证下降，GMRES 更稳健                       |

## 8. 并行策略
### 8.1 短期（单机多核）
| 层次                                                         | 方案                                                         | 加速比           |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ---------------- |
| Jacobian 装配                                                | Numba prange colored assembly（已有，stencil 不变）         | 4-8×（核数限制） |
| GMRES 线性求解                                               | SciPy 单线程                                                 | 1×（瓶颈）       |
| AMG setup                                                    | PyAMG 单线程                                                 | 1×               |
| AMG solve（预处理）                                          | PyAMG 单线程                                                 | 1×               |
| **瓶颈在线性求解**。装配已经并行，但 GMRES 单线程会成为限制。 |                                                              |                  |
### 8.2 中期（PETSc 并行）
| 层次                                                         | 方案                    | 加速比                 |
| ------------------------------------------------------------ | ----------------------- | ---------------------- |
| Jacobian 装配                                                | PETSc MatSetValues 并行 | $N_{\mathrm{procs}}$ × |
| GMRES                                                        | PETSc KSPGMRES          | $N_{\mathrm{procs}}$ × |
| AMG                                                          | PETSc GAMG              | $N_{\mathrm{procs}}$ × |
| 域分解                                                       | PETSc PCASM (Schwarz)   | $N_{\mathrm{procs}}$ × |
| PETSc 的 Newton-Krylov-Schwarz 正是 Cai-Keyes-Young 1996 的方案。 |                         |                        |

## 9. 风险与缓解
| 风险                               | 概率 | 影响 | 缓解                                                         |
| ---------------------------------- | ---- | ---- | ------------------------------------------------------------ |
| Newton 从远场初值在 transonic 发散 | 高   | 高   | Mach continuation + damping_theta 全局化；保留 Picard 作为 fallback |
| GMRES 收敛慢（Jacobian 不适定）    | 中   | 高   | AMG 预处理；如果不够加 ILU；monitor GMRES 迭代数             |
| AMG 对非对称矩阵效果差             | 中   | 中   | PyAMG 对非对称矩阵取对称部分做 aggregation；如果不行用 ILU   |
| Python GMRES 太慢                  | 中   | 中   | 短期接受；中期迁 PETSc                                       |
| Newton Jacobian 装配 bug           | 中   | 高   | 有限差分验证（N1 阶段）；bitwise no-op 测试（$\nu = 0$ 时）  |
| 策略 B 收敛阶不足（不含上游耦合）  | 低   | 低   | López 实践证明足够；如需严格二次可后加 P6 + 策略 A           |

## 10. 与现有代码的共存
- `solve/picard.py` 保留，作为 subsonic 快速路径和 transonic fallback
- `solve/continuation.py` 新增 `solver="newton"` 参数，默认仍为 `"picard"`
- 测试套件参数化：`@pytest.mark.parametrize("solver", ["picard", "newton"])`
- N3 验证通过后，`solve_transonic_lifting` 默认切换为 `"newton"`
- P5（ONERA M6）可以直接用 Newton 路径

## 11. 与其他 design notes 的关系
| 文档                                           | 关系                                                         |
| ---------------------------------------------- | ------------------------------------------------------------ |
| `20260707_1505_levelset_wake_design.md`        | 独立——Newton 和 level-set 尾流互不依赖。方案 B 的 B4 阶段去掉 $\Gamma$ secant 外层，与本文档的 N-$\Gamma$ 分裂是替代关系 |
| `20260707_2118_ibl_viscous_coupling_design.md` | 依赖——VII 迭代在势流求解外层，Newton 求解器使 transonic 更稳定，VII 更容易收敛 |

## 12. 参考文献
- López Canalejo, I., Núñez, M., Baiges, J., Rossi, R. "An embedded approach for the solution of the full potential equation with finite elements." *CMAME* 388 (2022) 114244.
  - Newton-Raphson + AMGCL，72k DOF 13 iterations 7s（4 核）
  - Jacobian 策略 B（含局部密度导数，不含上游耦合）
  - Eq.26: Jacobian 推导（含 $\partial\rho/\partial|v|^2$）
  - Eq.31: 尾流面 BC 忽略密度导数（"small velocity perturbations"）
- Cai, X.-C., Keyes, D.E., Young, D.P. "Parallel Newton-Krylov-Schwarz Algorithms for the Transonic Full Potential Equation." *SIAM J. Sci. Comput.* 19 (1998) 246–265.
  - Inexact Newton + GMRES + Schwarz 域分解，FPE 领域最成熟的并行 Newton 方案
  - Jacobian 策略 C（JFNK，有限差分 matvec）
- Holst, T.L. "Transonic flow computations using nonlinear potential methods." *Progress in Aerospace Sciences* 36 (2000) 1–61.
  - FPE 求解器演化谱系：SLOR → AF → 多重网格 → Newton
- Eisenstat, S.C., Walker, H.F. "Choosing the forcing terms in an inexact Newton method." *SIAM J. Sci. Comput.* 17 (1996) 16–32.
  - Eisenstat-Walker 自适应线性容差策略
- Demidov, D. "AMGCL — A C++ library for efficient solution of large sparse linear systems." *Software Impacts* 6 (2020) 100037.
  - López 使用的 AMG 库
- Hafez, M., South, J., Murman, E. "Artificial Compressibility Methods for Numerical Solutions of Transonic Full Potential Equation." *AIAA J.* 17(8) (1979) 838–844.
  - 人工密度方法原始文献，López 引用为 [4]
