# Newton 求解器实现策略：全耦合 Newton-Krylov（修订版）

> Design note for UP3D / pyFP3D
> Date: 2026-07-08 16:46
> Code baseline: HEAD a3f75ef
> **本文档替代 `20260707_2313_newton_solver_strategy.md`**（已标记为 ❌ 弃用，见 §1.1）
> Prerequisite: P4 已关闭；P5 medium 失败暴露 TE 奇异性 + secant-密度耦合失稳（H3）
> References:
> - **López Canalejo 2021 博士论文** (TUM) — `references/Lopez_2021_dissertation_transonic_FPE_embedded_wake.pdf`
>   - Ch.3 §3.2: FE 离散 + residual/Jacobian Eq.(3.8)-(3.16)
>   - Ch.3 §3.3: 人工可压缩性 Eq.(3.19)-(3.27)，switching function §3.3.1
>   - Ch.3 §3.4: 极限速度 / 密度截断 M²max~O(3.0)
>   - Ch.3 §3.5: 嵌入尾流 2D/3D，TE 节点处理 §3.5.4，small-cut 全积分 §3.5.5
>   - Ch.4 Table 4.6/4.9: Newton 收敛验证（严格二次）
>   - Ch.4 Table 4.7-4.8, 4.13: load stepping 策略（含 Mcrit/μc 调度）
>   - Appendix B Eq.(B.1)-(B.17): 完整敏感度推导（**含上游元素耦合**）
> - design.md §6.3 (Newton Jacobian), §8 (Nonlinear solution strategy)
> - roadmap.md P5-P7
> - Cai, Keyes, Young 1998 (SIAM J. Sci. Comput.) — Newton-Krylov-Schwarz for transonic FPE

---

## 1. 为什么需要修订

### 1.1 旧文档的三个核心错误

旧文档 `20260707_2313_newton_solver_strategy.md` 基于一个错误前提——"López 用策略 B（近似 Jacobian 不含上游耦合）"。这个前提在三个层面导致了错误决策：

| # | 旧文档声称 | 博士论文实际内容 | 后果 |
|---|-----------|----------------|------|
| 1 | López 不含上游耦合项（§1.5.3） | Appendix B Eq.(B.3)-(B.6) 明确推导了 $\partial\tilde{\rho}/\partial\varphi_j^{\mathrm{up}}$ | 策略 B 的文献依据不成立 |
| 2 | P6 不是 P7 的硬前置依赖（§1.5.4） | López 用策略 A（完整 Jacobian），P6 可微通量确实是前置 | 依赖关系链错误 |
| 3 | N-$\Gamma$ 分裂推荐（§3.3） | López 用全耦合 Newton（无 secant 外层） | H3 耦合失稳风险未消除 |

### 1.2 H3 根因的新认识（2026-07-08）

P5 medium mesh 失败后的调试确认了 H3 根因：**secant-密度耦合失稳**。

- secant 的 $F = \text{kutta\_targets}(\varphi(\Gamma, \rho(\Gamma))) - \Gamma$
- $\rho$ 未完全收敛时 $F$ 被密度滞后污染
- secant slope 估计失准 → $\Gamma$ 过冲 → $\varphi$ 大变 → $\rho$ 需更多步 → $F$ 污染更重 → 正反馈发散
- 全量 continuation（非缓存）M_max 29.3 vs 缓存 5.2，证明"加预算+降松弛"方向错误

**López 的架构不存在此问题**：他用 Newton-Raphson，残差是 $R_i(\varphi) = 0$（Eq.3.9），Jacobian 含完整 $\tilde{\rho}$ 敏感度（Appendix B），每步同时更新 $\varphi$ 和矩阵。不存在 secant 在未收敛 $\rho$ 上读污染 $F$ 的结构。

### 1.3 P5 暴露的 TE 奇异性问题

P5 medium mesh 失败的另一个根因是**尖锐 TE 的 P1 离散奇异性**：TE 单元缩小加剧梯度，密度截断 $M^2_{\max}\sim O(3.0)$ 产生不可微区域。这对 Newton 有直接影响：Jacobian 在密度截断区不可微，可能违反 Newton 的光滑性假设。

López 在 NACA 0012 验证中使用了**钝 TE** 修改（Eq.4.2），ONERA M6 同样避免尖锐 TE。UP3D 的 M1 spec 要求 sharp zero-thickness TE，需要特别处理。

---

## 2. López 博士论文的关键技术要点

### 2.1 Jacobian 结构（§3.2, Eq.3.11-3.14）

López 的完整 Jacobian 有两项：

$$J_{ij} = \underbrace{\sum_e \int_{\Omega_e} \tilde{\rho} \frac{\partial N_i}{\partial x_a} \frac{\partial N_j}{\partial x_a} \, d\Omega_e}_{\text{Term 1: Picard（对称）}} + \underbrace{\sum_e \int_{\Omega_e} \frac{\partial \tilde{\rho}}{\partial \varphi_j} u_a \frac{\partial N_i}{\partial x_a} \, d\Omega_e}_{\text{Term 2: 密度耦合（非对称）}}$$

其中 Term 2 展开后含上游元素 DOF（Appendix B）：

$$\frac{\partial \tilde{\rho}}{\partial \varphi_j} = \frac{\partial \rho}{\partial \varphi_j} - \mu_s \left(\frac{\partial \rho}{\partial \varphi_j} - \frac{\partial \rho_{\mathrm{up}}}{\partial \varphi_j}\right) + \text{(switching function derivatives)}$$

以超音速加速流为例（Eq.B.3-B.4）：

$$\frac{\partial \tilde{\rho}}{\partial \varphi_j} = \frac{\partial \rho}{\partial \varphi_j}(1-\mu_s) - (\rho - \rho_{\mathrm{up}})\frac{\partial \mu}{\partial \varphi_j} + \mu_s \frac{\partial \rho_{\mathrm{up}}}{\partial \varphi_j} \quad \text{(B.3)}$$

$$\frac{\partial \tilde{\rho}}{\partial \varphi_j^{\mathrm{up}}} = -\mu_s \frac{\partial \rho_{\mathrm{up}}}{\partial \varphi_j^{\mathrm{up}}} \quad \text{(B.4)}$$

**(B.4) 是非零的上游耦合项**——当前元素对上游元素 DOF 的 Jacobian 贡献。这意味着 stencil 比纯 Picard 宽一层。

### 2.2 收敛验证（§4.3-4.4, Table 4.6/4.9）

**亚声速**（Table 4.6, $M_\infty = 0.72$, $\alpha = 0°$）：

| k | $\|R\|_{L_2}$ |
|---|---------------|
| 1 | $5 \times 10^{-1}$ |
| 2 | $2 \times 10^{-1}$ |
| 3 | $5 \times 10^{-2}$ |
| 4 | $2 \times 10^{-3}$ |
| 5 | $8 \times 10^{-7}$ |
| 6 | $9 \times 10^{-12}$ |

**跨声速**（Table 4.9, NACA 0012 $M_\infty = 0.71$, load step 2）：

| k | Case 1-2 | Case 3 | Case 4 |
|---|----------|--------|--------|
| 1 | $1 \times 10^{-2}$ | $1 \times 10^{-2}$ | $5 \times 10^{-1}$ |
| 2 | $4 \times 10^{-5}$ | $1 \times 10^{-6}$ | $9 \times 10^{-2}$ |
| 3 | $4 \times 10^{-9}$ | $9 \times 10^{-12}$ | $4 \times 10^{-3}$ |
| 4 | $9 \times 10^{-12}$ | — | $5 \times 10^{-4}$ |
| 5 | — | — | $2 \times 10^{-2}$ |
| 6 | — | — | $7 \times 10^{-5}$ |
| 7 | — | — | $8 \times 10^{-10}$ |

Case 1-3 展示了**严格二次收敛**（残差每次平方下降：$10^{-2} \to 10^{-5} \to 10^{-9} \to 10^{-12}$）。这只有在完整 Jacobian（含上游耦合）时才可能——近似 Jacobian 只能达到超线性。

Case 4（$\alpha = 2°$, $M_\infty = 0.75$, $M_{\mathrm{crit}} = 0.90$, $\mu_c = 1.1$）在第 5 步有反弹（$5 \times 10^{-4} \to 2 \times 10^{-2}$），然后恢复二次收敛到 $8 \times 10^{-10}$。这是强激波 case，反弹可能来自激波位置微调。

### 2.3 Load Stepping 策略（§4.4, Table 4.7-4.8, 4.13）

**NACA 0012 跨声速**：

| Case | $\alpha$ | $M_\infty$ | $M_{\mathrm{crit}}$ | $\mu_c$ | Load steps |
|------|----------|-----------|---------------------|---------|------------|
| (1) | 1° | 0.72 | 0.99 | 1.0 | 3 步: 0.70→0.71→0.72 |
| (2) | 1° | 0.73 | 0.99 | 1.0 | 4 步: 0.70→0.71→0.72→0.73 |
| (3) | 1° | 0.75 | 0.95 | 1.1 | 5 步: 0.70→...→0.75 |
| (4) | 2° | 0.75 | 0.90 | 1.1 | 6 步: 0.70→...→0.75 |

**ONERA M6**（Table 4.13）：12 步，$\mu_c$ 从 2.0 逐步降到 1.6。

**关键观察**：
- $M_{\mathrm{crit}}$ 随 $M_\infty$ 和 $\alpha$ 增加而降低（0.99→0.90）——人工耗散在更低的 M 就开始激活
- $\mu_c$ 随非线性程度增加——更强激波需要更多耗散
- 步数随难度增加（3→4→5→6→12）
- 每步 5-9 Newton-Raphson 迭代（§4.7：transonic solver needs eight load steps, each requiring from five to nine Newton Raphson iterations）

### 2.4 尾流处理（§3.5, Eq.3.35-3.51）

López 的尾流用 least-squares 弱形式：

**2D**（Eq.3.35-3.42）：单个矢量方程 $\mathbf{u}_u - \mathbf{u}_l = 0$，least-squares 泛函 $\Pi = \frac{1}{2}\int (\mathbf{u}_u - \mathbf{u}_l)^2 \, d\Omega_W$。

**3D**（Eq.3.43-3.51）：两个线性化条件：
- $g_1 = \hat{n} \cdot (\mathbf{u}_u - \mathbf{u}_l) = 0$（质量连续，线性化 $\rho_u, \rho_l \to \rho_\infty$）
- $g_2 = \hat{u}_\infty \cdot (\mathbf{u}_u - \mathbf{u}_l) = 0$（压力相等，线性化 $\bar{u} \to \hat{u}_\infty$）

Jacobian 有四个 block：$J^{uu}, J^{ul}, J^{lu}, J^{ll}$（Eq.3.48-3.51）。这些 block **不含密度导数**（因为线性化用了 $\rho_\infty$），但含上游元素 DOF 耦合。

**UP3D 的尾流**（master-slave 消元）与 López 的对比见 §4.3。

### 2.5 TE 节点处理（§3.5.4）

López 的 TE 节点处理：
- TE 节点同时属于体面和尾流面
- 对体面施加 slip BC（$q = 0$），**不施加尾流 BC**
- 辅助 DOF $\varphi_i^{\mathrm{aux}}$ 仅用于断开 TE 上下面的单元，作为体下方单元的普通势 DOF

**UP3D 的 TE 处理**（`mesh/wake_cut.py`）：尾流面节点复制 + TE 节点复制（jump reaches wall）+ Kutta 探针选取（一阶壁面邻居）。机制不同但目标一致。

### 2.6 极限速度 / 密度截断（§3.4, Eq.3.28-3.29）

López 的密度截断：$M^2_{\max} \sim O(3.0)$，超过此值密度保持常数。Figure 3.5 显示翼尖 TE 处速度达 2888 m/s，密度截断在此激活。

**关键问题**：截断区 $\partial \rho / \partial u^2 = 0$，但跨过截断阈值时有跳变。这违反 Newton 的光滑性假设。

---

## 3. UP3D 与 López 的架构对比

### 3.1 求解器架构

```
López:   load step → Newton-Raphson(φ, ρ̃(φ)) → 收敛 → next step
         (每步组装完整 Jacobian 含上游耦合，4-12 步 ramp)

UP3D:    Mach continuation → secant(Γ) ↘ Picard(ρ) ↗ 交叉耦合 → 失稳
         (固定矩阵 + 外层 secant + 内层 density Picard，三层嵌套)
```

### 3.2 尾流处理

| 特征 | López | UP3D |
|------|-------|------|
| 方法 | least-squares 弱形式 | master-slave 消元 |
| 未知数 | $\varphi^u, \varphi^l$（上下分开） | $\varphi_{\mathrm{red}} + \Gamma$（消元后） |
| Jacobian block | $J^{uu}, J^{ul}, J^{lu}, J^{ll}$ | $T^T J T$（投影） |
| $\Gamma$ 处理 | 隐含在 least-squares 泛函中 | 显式分离，RHS-only |
| 压力相等约束 | 显式（$g_2 = 0$） | 隐式（Kutta 条件） |
| 密度导数 | 不含（线性化用 $\rho_\infty$） | 不含（master-slave 不涉及密度） |

**UP3D 的方式更简单**——Jacobian 只是 $T^T J T$，不需要额外的尾流 block。而且 $\partial R / \partial \Gamma = -h_j$ 已经在代码中（`wake.py` 的 `_h` 数组），全耦合 Newton 的 $\Gamma$ Jacobian 几乎是现成的。

### 3.3 代码现状

| 模块 | 文件 | 当前状态 | Newton 需要 |
|------|------|---------|------------|
| Picard 求解器 | `solve/picard.py` | 三层嵌套 | 保留作 fallback |
| Continuation | `solve/continuation.py` | Mach continuation + secant | 替换为 Newton + load stepping |
| Wake 消元 | `constraints/wake.py` | master-slave $T^T A T$ + $h_j$ | **直接复用** |
| Kutta 探针 | `constraints/wake.py` | `kutta_targets()` | **直接复用** |
| 人工密度 | `kernels/upwind.py` | integer-walk 上游搜索 | P6 需可微化 |
| 密度截断 | `kernels/density.py` | $M^2_{\max}$ clamp | 需光滑化 |
| 线性求解 | `solve/linear.py` | CG + AMG | 新增 GMRES + AMG |

---

## 4. 推荐策略：全耦合 Newton（策略 A）

### 4.1 为什么选策略 A 而非 B

| 因素 | 策略 A（完整 Jacobian） | 策略 B（近似 Jacobian） |
|------|------------------------|------------------------|
| 上游耦合 | 含（Eq.B.4） | 不含 |
| 收敛阶 | **严格二次** | 超线性 |
| López 先例 | ✅ 实际做法 | ❌ 旧文档误读 |
| P6 依赖 | 硬前置 | 不前置 |
| Stencil | 宽一层 | 不变 |
| H3 风险 | **消除**（无 secant 外层） | **仍存在**（secant 保留） |
| TE 奇异性 | 需密度截断光滑化 | 同左 |

**核心论据**：H3 根因是 secant-密度耦合失稳。策略 B 保留 secant 外层，只是用 Newton 替代 Picard 内层——如果 Newton 不完全收敛（策略 B 的近似 Jacobian 导致），$F$ 仍有 lag noise，secant 仍可能失稳。策略 A 用全耦合 Newton 直接消除 secant，根除 H3。

### 4.2 全耦合 Newton 的 Jacobian

UP3D 的 reduced 系统（master-slave 消元后）：

$$\underbrace{(T^T J T)}_{A_{\mathrm{red}}} \, \delta\varphi_{\mathrm{red}} - \sum_j h_j \, \delta\Gamma_j = -T^T R(\varphi_{\mathrm{red}}, \Gamma)$$

加上 Kutta 条件 $F_j(\varphi_{\mathrm{red}}, \Gamma) = \text{kutta\_targets}_j - \Gamma_j = 0$：

$$\frac{\partial F_j}{\partial \varphi_{\mathrm{red}}} \, \delta\varphi_{\mathrm{red}} + \frac{\partial F_j}{\partial \Gamma_j} \, \delta\Gamma_j = -F_j$$

全耦合 Newton 系统：

$$\begin{pmatrix} T^T J T & -h_1 & -h_2 & \cdots & -h_{n_s} \\ \nabla_\varphi F_1 & -1 & 0 & \cdots & 0 \\ \nabla_\varphi F_2 & 0 & -1 & \cdots & 0 \\ \vdots & & & \ddots & \\ \nabla_\varphi F_{n_s} & 0 & 0 & \cdots & -1 \end{pmatrix} \begin{pmatrix} \delta\varphi_{\mathrm{red}} \\ \delta\Gamma_1 \\ \delta\Gamma_2 \\ \vdots \\ \delta\Gamma_{n_s} \end{pmatrix} = -\begin{pmatrix} T^T R \\ F_1 \\ F_2 \\ \vdots \\ F_{n_s} \end{pmatrix}$$

其中：
- $T^T J T$：体积 Jacobian 投影到 reduced 空间（含上游耦合，策略 A）
- $h_j = T^T A g_j$：**已在代码中**（`wake.py` 的 `self._h`）
- $\nabla_\varphi F_j = \partial(\text{kutta\_targets}_j) / \partial\varphi_{\mathrm{red}}$：需推导

### 4.3 Kutta 目标的 Jacobian 推导

`kutta_targets` 的定义（`wake.py` line 114-130）：

$$\Gamma_j^{\mathrm{target}} = \frac{1}{n_j} \sum_{k \in \mathrm{station}_j} \left(\varphi_{\mathrm{upper}_k} - \varphi_{\mathrm{lower}_k}\right)$$

其中 upper/lower 探针是壁面自由 DOF（reduced 空间中的原始节点）。因此：

$$\frac{\partial \Gamma_j^{\mathrm{target}}}{\partial \varphi_{\mathrm{red},m}} = \frac{1}{n_j} \sum_{k \in \mathrm{station}_j} \left(\delta_{\mathrm{upper}_k, m} - \delta_{\mathrm{lower}_k, m}\right)$$

这是一个非常稀疏的矩阵——每行最多 $2 n_j$ 个非零（$n_j$ = station 的 TE 节点数，通常 2-4）。**推导已完成，实现简单**。

### 4.4 N-$\Gamma$ 分裂作为 fallback

如果全耦合 Newton 在某些 case 不稳定（例如 $\Gamma$ 初值太差），回退到 N-$\Gamma$ 分裂：

- 外层：固定点迭代 $\Gamma \leftarrow \Gamma + \omega_\Gamma (\Gamma^{\mathrm{target}} - \Gamma)$
- 内层：Newton 只对 $\varphi_{\mathrm{red}}$ 求解，$\Gamma$ 固定

这等同于当前 Picard 架构但用 Newton 替代 density Picard。**仅作为 fallback，不作为主路线**，因为保留 secant 外层仍有 H3 风险。

---

## 5. TE 奇异性与 Newton 的交互

### 5.1 问题

P5 medium mesh 暴露的问题：18 个翼面 TE 处 $M > 2$ 的单元，$M_{\max} = 5.2$。根因是尖锐 TE 的 P1 离散奇异性——TE 单元缩小加剧梯度。

密度截断 $M^2_{\max} \sim O(3.0)$ 在此区域激活，$\partial \rho / \partial u^2$ 在截断阈值处有跳变：
- $M^2 < M^2_{\max}$：$\partial \rho / \partial u^2 \neq 0$（正常等熵导数）
- $M^2 \geq M^2_{\max}$：$\partial \rho / \partial u^2 = 0$（密度冻结）

这个跳变让 Jacobian 不光滑，Newton 在含截断单元的区域可能：
1. 收敛变慢（超线性而非二次）
2. 在截断/非截断边界振荡
3. 极端情况下不收敛

### 5.2 López 的处理

López 用 $M^2_{\max} \sim O(3.0)$ 同样的截断（§3.4），但他的验证 case 用了钝 TE（NACA 0012 Eq.4.2 修改了最后一个系数）。Figure 3.5 显示翼尖 TE 处速度达 2888 m/s，但这是 3D 翼尖的几何奇异性，不是 2D TE 的 P1 离散奇异性。

**López 没有在尖锐 TE + 细网格上验证 Newton。**

### 5.3 UP3D 的缓解方案

| 方案 | 描述 | 优点 | 缺点 | 工作量 |
|------|------|------|------|--------|
| **A: 密度截断光滑化** | 用 smooth clamp 替代 hard clamp：$\rho_{\mathrm{clamped}} = \rho_\infty (1 - \sigma) + \rho(M^2_{\max}) \cdot \sigma$，其中 $\sigma = \mathrm{sigmoid}((M^2 - M^2_{\max})/\epsilon)$ | Jacobian 可微；最小改动 | 引入 $\epsilon$ 参数 | 小 |
| **B: 钝 TE** | 修改几何使 TE 有小厚度 | 消除奇异性；与 López 一致 | 改变几何 | 中 |
| **C: 局部加密 + 曲面单元** | TE 处用高阶单元或局部加密 | 最物理 | 工作量大 | 大 |
| **D: 排除截断单元的 Term 2** | 在 $M^2 > M^2_{\max}$ 的单元，Term 2 设为 0（退化为 Picard） | 简单 | Jacobian 不一致 | 小 |

**推荐**：方案 A（密度截断光滑化）作为 N1 阶段的一部分。如果光滑化后 Newton 仍在 TE 处收敛差，再考虑方案 B。

### 5.4 密度截断光滑化的具体形式

当前 hard clamp（`kernels/density.py`）：

$$\rho_{\mathrm{clamped}} = \begin{cases} \rho(M^2) & M^2 < M^2_{\max} \\ \rho(M^2_{\max}) & M^2 \geq M^2_{\max} \end{cases}$$

光滑化版本：

$$\sigma = \frac{1}{2}\left(1 + \tanh\frac{M^2 - M^2_{\max}}{\epsilon}\right), \quad \epsilon \sim 0.1$$

$$\rho_{\mathrm{clamped}} = (1 - \sigma) \cdot \rho(M^2) + \sigma \cdot \rho(M^2_{\max})$$

导数：

$$\frac{\partial \rho_{\mathrm{clamped}}}{\partial u^2} = (1-\sigma) \frac{\partial \rho}{\partial u^2} + \frac{\partial \sigma}{\partial u^2} \left(\rho(M^2_{\max}) - \rho(M^2)\right)$$

其中 $\partial \sigma / \partial u^2$ 通过链式法则得到，$\partial M^2 / \partial u^2$ 已在 Appendix B Eq.(B.12) 给出。

---

## 6. Jacobian 策略对比（修正版）

### 6.1 三种策略

| 策略 | 包含项 | 需要可微通量？ | 收敛阶 | Stencil | López 先例 | H3 风险 |
|------|--------|--------------|--------|---------|-----------|---------|
| **A: 精确 Jacobian** | 1 + 2 + 3 | **是** | 严格二次 | 宽一层 | ✅ 实际做法 | **消除** |
| B: 近似 Jacobian | 1 + 2 | 否 | 超线性 | 不变 | ❌ 旧文档误读 | 仍存在 |
| C: JFNK | 有限差分 matvec | 否 | 超线性 | 不变 | Cai-Keyes-Young | 消除（无 secant） |

Term 1 = Picard（对称），Term 2 = 局部密度导数（非对称），Term 3 = 上游耦合（非对称，stencil 宽一层）。

### 6.2 依赖关系（修正版）

```
策略 A（推荐）:
  P5（P1 TE 奇异性缓解）→ P6（可微通量）→ P7（全耦合 Newton）
  严格二次收敛；直接消除 H3

策略 B（fallback）:
  P7（Newton + secant 外层）
  超线性收敛；H3 风险仍存在

策略 C（备选）:
  P6 非前置；P7（JFNK + 全耦合 Γ）
  超线性收敛；无显式 Jacobian，GMRES 迭代多
```

### 6.3 P6 可微通量的具体需求

策略 A 的 Term 3 要求 $u(e)$（上游元素选择）可微。当前 `kernels/upwind.py` 的 integer-walk 不可微。

P6 需要：
1. **可微的上游选择**：用 smooth max/argmax 替代离散选择，或用 weighting 方案（权重对所有邻居连续可微）
2. **可微的 switching function**：$\max(\nu_e, \nu_u)$ 用光滑近似 $\max_\varepsilon(a, b) = \frac{a+b}{2} + \sqrt{\frac{(a-b)^2}{4} + \varepsilon}$，$\varepsilon \sim 10^{-8}$
3. **可微的密度截断**：见 §5.4 光滑化方案

---

## 7. 实现路线

### 7.1 修订后的阶段

| 阶段 | 目标 | 前置依赖 | 交付物 |
|------|------|---------|--------|
| **N0** | 密度截断光滑化 + 验证 | 无 | 光滑 clamp + 单元测试 |
| **N1** | P6 可微通量 | N0 | 可微上游选择 + switching function |
| **N2** | Newton Jacobian 装配 + 验证 | N1 | 策略 A Jacobian + 有限差分验证 |
| **N3** | GMRES + AMG 线性求解 | N2 | 非对称线性求解路径 |
| **N4** | 全耦合 Newton 驱动器 + 亚声速验证 | N3 | Newton 驱动器 + $\Gamma$ Jacobian |
| **N5** | 跨声速验证 + load stepping | N4 | 参数调度 + G7.1 收敛验证 |
| **N6** | ONERA M6 + 性能 gate | N5 | G7.2/G7.3 性能目标 |

### 7.2 Phase N0：密度截断光滑化

**目标**：消除密度截断处的不可微性。

**交付物**：
- `kernels/density.py` 修改：`limit_q2_field` → `smooth_clamp_q2`，用 $\tanh$ 光滑化
- 单元测试：验证 $M^2 < M^2_{\max}$ 时与 hard clamp 一致（误差 $< 10^{-8}$）
- 验证 $\partial \rho / \partial u^2$ 在截断区连续

**验证标准**：
- P3/P4 回归测试无退化（截断区外 bitwise 一致）
- P5 medium TE 处 $M_{\max}$ 降低（光滑化不硬化截断，但 Jacobian 可微）

### 7.3 Phase N1：P6 可微通量

**目标**：使人工密度 $\tilde{\rho}_e$ 对 $\varphi$ 完全可微。

**交付物**：
- `kernels/upwind.py` 修改：
  - 上游选择：用 smooth weighting 替代 integer-walk，或保留 integer-walk 但在 Jacobian 中用 finite-difference 近似 $\partial u(e) / \partial \varphi$（López 的 Appendix B 推导假设上游固定，实际实现中可能也是 integer-walk + frozen upstream）
  - Switching function：$\max \to \max_\varepsilon$
- `kernels/jacobian_terms.py` 新增：Term 2 + Term 3 的 element-level kernel
- 单元测试：有限差分验证 $\partial \tilde{\rho} / \partial \varphi$ 在亚声速/加速/减速三种流态下正确

**验证标准**：
- 有限差分 Jacobian 与解析 Jacobian 一致（相对误差 $< 10^{-6}$）
- 亚临界流（$\nu = 0$）：Term 2 = Term 3 = 0，Jacobian = Picard 矩阵
- G4.2 bitwise no-op 一致性

> **López 实现的开放问题**：Appendix B 的推导假设上游元素固定（$\partial \rho_{\mathrm{up}} / \partial \varphi_j^{\mathrm{up}}$ 非零，但 $u(e)$ 的选择本身不可微）。López 可能也是 integer-walk + frozen upstream selection（每步重新选上游，但 Jacobian 用当前选择计算导数）。这与 P6 "可微通量" 的定义需要对齐——可能 P6 只需要 switching function 和密度截断可微，上游选择保持 integer-walk 但每步更新。N1 阶段需实验验证。

### 7.4 Phase N2：Newton Jacobian 装配

**目标**：实现策略 A 的完整 Jacobian 装配（Term 1 + 2 + 3）。

**交付物**：
- `kernels/jacobian_newton.py` 新建：
  - Term 1：复用现有 `assemble_matrix_data_colored`
  - Term 2：新增 kernel，per-element $\frac{\partial \tilde{\rho}_e}{\partial q_e^2} \cdot 2(\nabla\varphi \cdot \nabla N_k)(\nabla\varphi \cdot \nabla N_i) \cdot V_e$
  - Term 3：新增 kernel，per-element $\frac{\partial \tilde{\rho}_e}{\partial \rho_{\mathrm{up}}} \cdot \frac{\partial \rho_{\mathrm{up}}}{\partial \varphi_j^{\mathrm{up}}}$ 对上游 DOF 的贡献
- Stencil 扩展：`mesh/coloring.py` 需要重新着色（Term 3 宽一层）
- `kernels/elem_to_csr.py` 扩展：scatter map 支持上游邻居 DOF

**验证标准**：
- 有限差分 Jacobian 与解析 Jacobian 一致（相对误差 $< 10^{-6}$），包括超音速区
- 亚临界流（$\nu = 0$）：Jacobian = Picard 矩阵 + 零
- G4.2 bitwise no-op 一致性

### 7.5 Phase N3：GMRES + AMG 线性求解

**目标**：实现非对称线性求解路径。

**交付物**：
- `solve/linear.py` 新增 `solve_gmres_amg()`
- AMG 预处理：PyAMG `smoothed_aggregation_solver`，对 $J + J^T$ 做 setup（非对称预处理的标准做法）
- AMG setup 复用：Newton 步之间 Jacobian 变化不大，hierarchy 可复用 5-10 步再重建

**验证标准**：
- GMRES 收敛（20-50 次迭代，给定好的预处理）
- AMG setup 复用 5 步后重建，不影响收敛
- 与 CG+AMG 在对称情形（$\nu = 0$）下结果一致

### 7.6 Phase N4：全耦合 Newton 驱动器

**目标**：实现 Inexact Newton 驱动器，$\Gamma$ 作为未知数。

**交付物**：
- `solve/newton.py` 新建：
  - 全耦合 Jacobian 组装（§4.2）
  - Kutta 目标 Jacobian $\nabla_\varphi F$（§4.3）
  - Eisenstat-Walker 自适应线性容差
  - damping_theta 全局化（前几步提供阻尼）
  - Mach continuation 支持
- `solve/continuation.py` 新增 `solver="newton"` 选项

**验证标准**：
- NACA0012 $M_\infty = 0.5$, $\alpha = 2°$：$c_l$ 与 P3 一致（$< 0.5\%$ 差异）
- Newton 收敛在 5-10 步（亚声速，非线性弱）
- $\|R\|$ 单调下降
- G4.2 bitwise no-op 仍然成立
- **无 secant 外层**：$\Gamma$ 随 $\varphi$ 同时收敛

### 7.7 Phase N5：跨声速验证 + load stepping

**目标**：在 transonic 上验证 Newton，实现参数调度。

**交付物**：
- `solve/continuation.py` 的 `solve_transonic_lifting` 默认走 Newton 路径
- 参数调度模块：每步根据 $M_\infty$ 和 $\alpha$ 调整 $M_{\mathrm{crit}}$ 和 $\mu_c$
  - 参考 López Table 4.7-4.8 的调度表
  - 初始 $M_{\mathrm{crit}} = 0.99$，随 $M_\infty$ 增加逐步降低到 $\sim 0.90$
  - $\mu_c$ 从 1.0 逐步增加到 $\sim 1.6$（ONERA M6 到 $\sim 1.6$）
- Newton + Mach continuation + damping_theta 全局化

**验证标准**：
- NACA0012 $M_\infty = 0.80$, $\alpha = 1.25°$ coarse：激波 $x/c \approx 0.60$（与 P4 一致）
- **G7.1**：Newton 终端**严格二次收敛**——$\|R_{k+1}\| / \|R_k\|^2 \to$ 常数（参考 López Table 4.9 Case 1-3 的 $10^{-2} \to 10^{-5} \to 10^{-9} \to 10^{-12}$ 模式）
- 总 Newton 步 $\leq 30$ per load step（参考 López 每步 5-9 迭代）
- coarse 求解时间 $< 30$s（vs P4 的 174s）

### 7.8 Phase N6：ONERA M6 + 性能 gate

**目标**：在 ONERA M6 上验证，达到 G7.2 性能目标。

**交付物**：
- ONERA M6 medium mesh 求解
- 参数调度表（参考 López Table 4.13 的 12 步 $\mu_c$ 调度）
- 性能 profiling 报告

**验证标准**：
- **G7.2**：ONERA M6 medium mesh $< 5$ min single node
- **G7.3**：full regression suite $< 10$ min
- P5 V5 gate（$\lambda$ 激波拓扑、剖面 Cp）通过

---

## 8. Inexact Newton + Eisenstat-Walker

### 8.1 策略

Eisenstat-Walker 自适应线性容差：第 $k$ 步的线性容差 $\eta_k$：

$$\eta_k = \max\left(\eta_{\min}, \min\left(\eta_{\max}, \gamma \cdot \frac{\|R_k\|^\alpha}{\|R_{k-1}\|^{\alpha-1}}\right)\right)$$

典型参数：$\eta_{\min} = 10^{-4}$, $\eta_{\max} = 0.9$, $\gamma = 0.9$, $\alpha = 1$（线性区）或 $2$（二次区）。

### 8.2 效果

- 早期 Newton 步：线性容差 $10^{-2} \sim 10^{-3}$，GMRES 5-10 次迭代
- 终端 Newton 步：线性容差 $10^{-6} \sim 10^{-8}$，GMRES 20-50 次迭代
- 总体减少 50-70% 的线性求解开销

当前代码的 `forcing` 参数已实现简化版（`solve_subsonic_lifting` 的 `forcing=0.01`），Newton 路径需要用完整 Eisenstat-Walker 替换。

---

## 9. 全局化策略

### 9.1 Mach continuation + load stepping

纯 Newton 从远场初值在 transonic 下会发散。全局化策略：

| 策略 | 描述 | 当前状态 |
|------|------|---------|
| **Mach continuation** | 从 $M_{0.70}$ 收敛解重启 $M_{0.80}$ | ✅ 已实现 |
| **Load stepping** | 参考 López：每步调整 $M_\infty, M_{\mathrm{crit}}, \mu_c$ | ❌ 新增 |
| damping_theta | $J + \theta \cdot \mathrm{diag}(J)$ 提供阻尼 | ✅ 已实现（Picard path） |
| 伪时间步进 | $J + \mathrm{diag}(m/\Delta \tau)$ | ⚠️ 已有 `pseudo_dt`，Newton 路径可复用 |
| 线搜索 | $\varphi_{k+1} = \varphi_k + \lambda \cdot \delta\varphi$，$\lambda \in (0,1]$ 使 $\|R\|$ 下降 | ❌ 可选 |

**推荐**：Mach continuation + load stepping（参数调度）+ damping_theta。Newton 从上一级 Mach 收敛解重启，damping_theta 在前几步提供阻尼，收敛后自动退化为纯 Newton。

### 9.2 Load stepping 参数调度

参考 López Table 4.7-4.8, 4.13，建议 UP3D 的调度表：

| $M_\infty$ 范围 | $M_{\mathrm{crit}}$ | $\mu_c$ | 步长 $\Delta M$ |
|-----------------|---------------------|---------|-----------------|
| $< 0.70$ | 0.99 | 1.0 | 0.05 |
| $0.70 - 0.75$ | 0.95 | 1.0-1.1 | 0.01-0.02 |
| $0.75 - 0.80$ | 0.90 | 1.1-1.3 | 0.01 |
| $0.80 - 0.85$ | 0.85 | 1.3-1.6 | 0.005-0.01 |
| $> 0.85$ | 0.80 | 1.6-2.0 | 0.005 |

> 具体参数需在 N5 阶段实验标定。López 的值是 2D NACA 0012 和 3D ONERA M6 的，UP3D 的网格和几何可能需要不同的值。

---

## 10. 预期性能

### 10.1 迭代次数对比

| 场景 | Picard 迭代 | Newton 预期 | 每步代价 | 总代价比 |
|------|------------|-------------|---------|---------|
| P3 subsonic | 15 | 5-8 | 3-4× Picard | ~1-2×（亚声速 Newton 收益小） |
| P4 transonic coarse | 10,464 | 20-40 per load step × 3-6 steps = 60-240 total | 3-4× Picard | **5-15× 更快** |
| P4 transonic medium | ~12,931 | 25-50 per load step × 3-6 steps = 75-300 total | 3-4× Picard | **3-10× 更快** |
| ONERA M6 medium | 未知（估计 >30,000） | 5-9 per load step × 12 steps = 60-108 total | 3-4× Picard | **10-30× 更快** |

> **修正旧文档的错误**：旧文档引用"López 72k DOF, 13 iterations"——这是单个 load step 的迭代数，不是总数。López 的 ONERA M6 用 12 load steps × 5-9 Newton/step = 60-108 总迭代。

### 10.2 时间预估

| 场景 | 当前时间 | Newton 预估 | 依据 |
|------|---------|-------------|------|
| P4 coarse | 174s | 10-30s | ~80 Newton steps × (装配 ~0.5s + GMRES ~0.3s) |
| P4 medium | 16m39s | 2-5min | ~150 Newton steps × (装配 ~2s + GMRES ~1.5s) |
| ONERA M6 medium | 未知 | 3-6min | ~84 Newton steps × (装配 ~2s + GMRES ~1.5s) + 12 load steps overhead |
| 目标 G7.2 | — | < 5 min | ⚠️ 边缘——需要 AMG setup 复用 + Eisenstat-Walker 才能达标 |

> **比旧文档保守**：旧文档预估 42-70s 过于乐观。修正后 ONERA M6 预估 3-6min，在 G7.2 的 5min 目标边缘。需要 Eisenstat-Walker + AMG setup 复用才能确保达标。

### 10.3 López 数据校准

López CMAME 期刊论文的数据（72k DOF, 13 iterations, 7s）是**单个 load step** 的。博士论文 §4.7 的完整数据：

- NACA 0012 transonic：3-6 load steps × 5-9 Newton/step = 15-54 total
- ONERA M6：12 load steps × 5-9 Newton/step = 60-108 total

UP3D 的差异：
- Python + Numba（装配快，但 SciPy GMRES 比 C++ 慢 2-5×）
- 单线程线性求解（短期）
- 规模系数：ONERA M6 medium ~63k nodes ≈ López 的 72k DOF

---

## 11. 模块改动

### 11.1 新增/修改模块

```
pyfp3d/
├── kernels/
│   ├── density.py              # 修改：smooth clamp（N0）
│   ├── upwind.py               # 修改：可微 switching function（N1）
│   ├── jacobian_newton.py      # 新建：策略 A Jacobian 装配（N2）
│   └── jacobian.py             # 不变（Picard 路径保留）
├── mesh/
│   └── coloring.py             # 修改：宽 stencil 重新着色（N2）
├── solve/
│   ├── newton.py               # 新建：全耦合 Newton 驱动器（N4）
│   ├── picard.py               # 不变
│   ├── linear.py               # 扩展：GMRES + AMG（N3）
│   └── continuation.py         # 修改：solver="newton" + 参数调度（N5）
```

### 11.2 改动清单

| 文件 | 改动 | 阶段 | 改动量 |
|------|------|------|--------|
| `kernels/density.py` | smooth clamp 替代 hard clamp | N0 | 小 |
| `kernels/upwind.py` | 可微 switching function | N1 | 中 |
| `kernels/jacobian_newton.py` | **新建**：Term 1+2+3 装配 | N2 | 大 |
| `kernels/elem_to_csr.py` | 扩展：上游 DOF scatter | N2 | 中 |
| `mesh/coloring.py` | 重新着色（宽 stencil） | N2 | 中 |
| `solve/linear.py` | 新增 `solve_gmres_amg()` | N3 | 小 |
| `solve/newton.py` | **新建**：Newton 驱动器 + $\Gamma$ Jacobian | N4 | 大 |
| `solve/continuation.py` | `solver="newton"` + 参数调度 | N5 | 中 |

---

## 12. 风险与缓解

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| Newton 从远场初值在 transonic 发散 | 高 | 高 | Mach continuation + load stepping + damping_theta |
| TE 奇异性导致 Jacobian 不可微 | 高 | 高 | N0 密度截断光滑化；如果不够用钝 TE |
| GMRES 收敛慢（Jacobian 不适定） | 中 | 高 | AMG 预处理；如果不够加 ILU；monitor GMRES 迭代数 |
| AMG 对非对称矩阵效果差 | 中 | 中 | PyAMG 取对称部分做 aggregation；如果不行用 ILU |
| 上游选择不可微（integer-walk） | 中 | 中 | N1 阶段实验：可能 frozen upstream + 每步更新足够（López 可能也是这么做） |
| Python GMRES 太慢 | 中 | 中 | 短期接受；中期迁 PETSc |
| Newton Jacobian 装配 bug | 中 | 高 | 有限差分验证（N2 阶段）；bitwise no-op 测试 |
| Stencil 宽一层导致内存/着色问题 | 低 | 中 | 预估内存增加 ~30%；重新着色是一次性开销 |
| G7.2 性能目标不达标 | 中 | 中 | Eisenstat-Walker + AMG setup 复用；如果不够迁 PETSc |

---

## 13. 与其他 design notes 的关系

| 文档 | 关系 |
|------|------|
| `20260707_1505_levelset_wake_design.md` | 独立——Newton 和 level-set 尾流互不依赖。全耦合 Newton 消除 $\Gamma$ secant 后，方案 B 的 B4 阶段（去掉 secant）自动实现 |
| `20260707_2118_ibl_viscous_coupling_design.md` | 依赖——VII 迭代在势流求解外层，Newton 求解器使 transonic 更稳定，VII 更容易收敛 |
| `20260707_2313_newton_solver_strategy.md` | **❌ 弃用**——本文档替代。旧文档的策略 B + N-$\Gamma$ 分裂路线基于错误前提 |

---

## 14. 与现有代码的共存

- `solve/picard.py` 保留，作为 subsonic 快速路径和 transonic fallback
- `solve/continuation.py` 新增 `solver="newton"` 参数，默认仍为 `"picard"`
- 测试套件参数化：`@pytest.mark.parametrize("solver", ["picard", "newton"])`
- N4 验证通过后，`solve_transonic_lifting` 默认切换为 `"newton"`
- P5（ONERA M6）可以直接用 Newton 路径

---

## 15. 参考文献

- **López Canalejo, I.** (2021). *A Finite-Element Transonic Potential Flow Solver with an Embedded Wake Approach for Aircraft Conceptual Design.* PhD dissertation, Technische Universität München.
  - `references/Lopez_2021_dissertation_transonic_FPE_embedded_wake.pdf`
  - Ch.3 §3.2: FE 离散 + residual/Jacobian Eq.(3.8)-(3.16)
  - Ch.3 §3.3: 人工可压缩性 Eq.(3.19)-(3.27)，switching function 三种流态
  - Ch.3 §3.4: 极限速度 / 密度截断 M²max~O(3.0)
  - Ch.3 §3.5: 嵌入尾流 2D Eq.(3.35-3.42) / 3D Eq.(3.43-3.51)，TE 节点 §3.5.4，small-cut §3.5.5
  - Ch.4 Table 4.6/4.9: Newton 严格二次收敛验证
  - Ch.4 Table 4.7-4.8, 4.13: load stepping + Mcrit/μc 调度
  - Appendix B Eq.(B.1)-(B.17): 完整敏感度推导，**含上游元素耦合 (B.4, B.6)**
- López Canalejo, I., Núñez, M., Baiges, J., Rossi, R. "An embedded approach for the solution of the full potential equation with finite elements." *CMAME* 388 (2022) 114244.
  - 期刊版，内容是博士论文 Ch.3-4 的浓缩
- Cai, X.-C., Keyes, D.E., Young, D.P. "Parallel Newton-Krylov-Schwarz Algorithms for the Transonic Full Potential Equation." *SIAM J. Sci. Comput.* 19 (1998) 246–265.
  - Inexact Newton + GMRES + Schwarz 域分解
- Holst, T.L. "Transonic flow computations using nonlinear potential methods." *Progress in Aerospace Sciences* 36 (2000) 1–61.
  - FPE 求解器演化谱系：SLOR → AF → 多重网格 → Newton
- Eisenstat, S.C., Walker, H.F. "Choosing the forcing terms in an inexact Newton method." *SIAM J. Sci. Comput.* 17 (1996) 16–32.
- Hafez, M., South, J., Murman, E. "Artificial Compressibility Methods for Numerical Solutions of Transonic Full Potential Equation." *AIAA J.* 17(8) (1979) 838–844.
