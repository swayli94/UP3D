# GB31.4 — LS step 语义评估（④ 配套，证据性关闭）

> 日期：2026-07-21 · 预注册：`../PRE_REGISTRATION.md` §3 GB31.4
> 结论：**④ 关闭——step 接受/线搜索不是两翼天花板的约束（证据如下）**；
> "LS line search 移植"框架正式退役（LS 已有回退线搜索）。
> 复检钩：GB31.3 若改变 LS step 动力学（见 §4），本文件追加一节。

## 1. 事实：LS Newton 已有 safety backtrack，"移植"是误述

`solve_multivalued_newton` 的步长控制（newton_ls.py:899-922）：
best-of-tried 回退搜索——λ 从 1 起每次减半至 `lam_min=0.05`（实为 5
次试步），度量为自由 DOF 残差的 **inf-norm**、strict 减小判据、
全部试步非有限即 RuntimeError；与 conforming 的 merit 线搜索
（newton.py:938-972，|R|²+|F|² L2 merit、非严格减小、5 次减半无地板）
同属"安全网型回退"。LS 侧非精确 Krylov 步即由该网吸收
（newton_ls.py:378-379, 847-848）。因此 ④ 的"移植"只能是一次
**度量/语义升级**（L2 merit、Armijo、λ 地板、开关），不是新增能力。

## 2. 证据：两翼濒死级上 step 接受均非约束

**LS（B30 G2/G3 traces，committed `b30_transonic_ceiling/results/`）：**

| 级 | c | 末态残差 | 钳制 | 死法 |
|---|---|---|---|---|
| 0.7875 | 1.5 | 1.08e-13 | 5 lim + 1 flr（末 20 步计数 range=0，**卡死**） | freeze 窗外、80 步耗尽 |
| 0.80 | 2.0 | 1.93e-13 | 2+7 = 9 > 8 窗 | 同上 |
| 0.7875 | 1.6（降阶） | 3.32e-05 | 4+5 = 9 > 8 窗 | 同上 |

残差已到机器零附近而 Newton 在钳制态上**静止**——这是 active set
卡死，任何步长控制都无法把收敛态"走回"可接受态；死亡由验收门禁
（freeze 窗 8）判决，不由步长机制判决。

**CONF（同上 traces）：** 0.83 在 c=1.5 残差停滞 8.1e-6、c=2.0 停滞
1.4e-4（锯齿），0 limited、1–2 floored toggling——conforming 的
merit 线搜索**在位且无效**：停滞是光滑态上的真 Newton 失速
（翼尖奇点经物理本身污染雅可比/残差结构），不是步长过冲。

**结论**：两翼天花板的约束都是翼尖物理（LS = 钳制 active set 卡死 +
窗语义；CONF = 光滑停滞），step 接受机制在两侧均非约束 → ④ 关闭。

## 3. 存档：若未来要做 merit 语义升级的接口面与风险（探查记录）

- 触及面：newton_ls.py:899-922 回退块（度量/判据/λ 计划）、
  `LSNewtonSystem.residual`（newton_ls.py:205-234，试态评估）、
  Γ 一致性点（newton_ls.py:923-924 步后更新、:569-596 远场刷新）。
  `wake/multivalued.py` 无需动。
- 三大结构风险：① merit 单位混排（质量行 vs |q|² Kutta 行；
  conforming 用 `kutta_sigma` 归一，但那是 eliminated-Γ 结构，LS 没有
  对应物）；② 两表脱同步——验收/freeze/fail-fast 全读 inf-norm
  （newton_ls.py:600,662,696,784），L2 merit 可能接受 inf-norm 上升的
  态，与手调的 freeze 语义打架；③ 试步内 Γ/远场冻结（vortex 远场下
  试步 merit 非实态 merit；committed neumann 跑无此问题）。
- 成本注记：LS 每次试步付一整次 `newton_side_data` + `assemble_matrix`
  （newton_ls.py:215-226），merit 搜索比 conforming 贵。

## 4. 复检钩（GB31.3 之后填写）

GB31.3 的片终止 re-spec 若引入以下任一即复检，否则维持关闭：
trial 全非有限 / λ 打满地板成常态 / merit-inf 两表系统性脱同步。
