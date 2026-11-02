# 20261103-0420 — **预注册：第 7 轮 = IBL 内层地板独立归因（S5 前置研究）**

**性质**：研究轮预注册；本文件必须先于研究代码与正式测量提交。**不改生产库、
不改默认值、不改 D08/D13 判据、不刷新任何 committed gate 证据，不启动 S5 的 2×2。**
**由来**：使用者裁决单独研究 IBL；S5 第 6 轮判定 §6 要求先查 IBL 内层不收敛/既有地板。
**分支**：`kimi/phase7-ibl-floor`（stack 在尚未合入的 `kimi/phase7-xtr030` 顶端；
该顶端包含最新 `origin/main`，没有丢弃或复制 S5 的 11 个提交）。

---

## 1. 一个问题与不得预设的答案

问：**S5 当前 xtr030 的 B1/B16 末态，其 IBL 内层不收敛是否复现已登记的
TE 带 `(B,δ)` 公式层地板与伪时间控制器到底签名；在同一末态、同一边数据下原样重启动，
它究竟收敛、撞地板，还是走出另一种轨迹？**

已有事实只有：四条相关原生/对照腿的各外迭代 `ibl_converged` 全为 0；旧 V5 研究在
别的状态上证明过 medium IBL 地板。不得把“残差约为 1e-6”单独当成同一机制的判据，
也不得预设重启动仍不收敛。

本轮是**识别与归因**，不是闭包修复、控制器调参或 S5 重开。以下结果都合法：

- 两条 xtr030 主腿均复现旧地板签名；
- 两条主腿均在原配置重启动后达到 IBL 自身收敛判据；
- 两腿混合、签名不一致或出现新轨迹，判 UNDEFINED；
- 重建不能复现原生末态残差，判仪器 FAIL 并停止。

## 2. 固定出处、状态与仪器边界

### 2.1 四个已有状态，零条新耦合腿

| 角色 | S5 腿 | 输入末态 | 首 IBL 冻结边界 |
|---|---|---|---|
| 主腿 | `native_xtr030_b1_r1` | `_work/state_native_xtr030_b1_r1.npz` | `_work/first_native_xtr030_b1_r1.npz` |
| 主腿 | `native_xtr030_b16_r1` | `_work/state_native_xtr030_b16_r1.npz` | `_work/first_native_xtr030_b16_r1.npz` |
| 对照 | `control_xtr005_b1` | `_work/state_control_xtr005_b1.npz` | `_work/first_control_xtr005_b1.npz` |
| 对照 | `control_xtr005_b16` | `_work/state_control_xtr005_b16.npz` | `_work/first_control_xtr005_b16.npz` |

- 文件均位于已提交的 `bench/studies/xtr030_attractor/summary.csv` 所指本地 `_work/`。
  `_work` **只是仪器输入，不是证据**；每个输入文件记录 SHA-256，`U/phi/delta_star`
  数组 SHA 必须逐项等于 S5 committed summary，medium 网格直接文件 SHA 也必须相同。
- B1/B16 是状态的**历史出处标签**。本轮所有诊断 worker 都固定实际 BLAS=1、Numba=16，
  从库读 `blas_thread_count()` 并断言；不再研究 BLAS 响应。
- 不清 `__pycache__`，使用暖缓存；每腿在 fresh Python 进程中运行。缺少输入或指纹不符
  就响亮停止，**不在本轮重算 S5 耦合腿来补缓存**。
- 直接使用末态保存的 `ue_surf/U` 与该腿首 IBL 的 `inflow_mask/inflow_state`；
  `rho_e/mach_e` 由保存的 `ue_surf` 按生产公式重算。不得从最终 FP 重解后的 `phi`
  反推并冒充最后一次 IBL 的边数据。
- 旧 `phases/p2/bench/studies/v5_ibl_floor/results/` 只用于事前校准签名和作历史对照，
  不重跑旧 S1/S2/S3，不重算 dense SVD。

### 2.2 原生重建守卫（先于任何重启动）

用生产 `IBL3Solver` 在四个末态重建纯定常残差 `F=solver.residual(U)`：

- `maxabs(F)` 必须与相应 committed `history_*.csv` 末行 `ibl_final_residual`
  相符；CSV 只有 12 位有效数字，允许误差仅为 `5e-12 * recorded_residual`；
- 重启动的首个纯定常残差也必须通过同一复核；
- 实际 BLAS、数组/网格指纹或输入形状任一不符 ⇒ **FAIL（仪器无效）**，
  保存诊断，停止派生腿，不改输入去凑一致。

## 3. 测量与事前判据

### C1 — 末态残差解剖与闭包地板活动集（4 腿）

每腿输出全部 `(node,equation)` 残差、`x/c`、上下翼面、层流/湍流/入流钉扎标记，
以及 `DELTA_MIN/RE_D_MIN` 活动节点。区域事前固定为 `le_inflow(x<=0.02)`、
`le_band(0.02<x<=0.05)`、`mid(0.05<x<0.90)`、`te_band(x>=0.90)`；
`x>=0.96` 的近尾缘能量份额只作辅助读数。

“**旧地板残差签名**”必须同时满足：

1. 全局 `maxabs(F)` 所在节点位于 `te_band`；
2. 六方程的 L2 范数份额中，最大两项恰为 `B` 与 `delta`，且二者合并
   `norm(F[:,[B,delta]])/norm(F) >= 0.90`；
3. `DELTA_MIN` 与 `RE_D_MIN` 活动集均为空。

阈值由旧 S1/S2 committed CSV 事前校准：二者合并份额分别为 0.9605/0.9666，
全局最大均在 `x>=0.90`，活动集均为空。TE 总能量、近 TE 能量、逐节点排名完整记录，
但不另造阈值。签名不满足不等于“没有地板”，只表示不能归入这个已定义签名。

### C2 — 同状态、同边数据原样重启动（4 腿）

每腿从末态 `U` 出发，使用同一个重建 solver、原 `tol=1e-9/max_iter=100`、
默认 `cfl0=1.0/growth=2.0/cfl_max=1e8`，调用生产 `IBL3Solver.solve` 一次。
不改容差、不加迭代、不改阻尼、不换种子。保存逐迭代纯残差、verbose 控制器
CFL/拒绝记录、最终 CFL、终止迭代数和最终残差。

逐腿判别：

- **RESTART_CONVERGED**：生产返回 `converged=True`，且最终纯残差满足同一生产判据；
- **CONTROLLER_FLOOR**：生产返回不收敛，最终 CFL 到达 `cfl_min=1e-3`，
  出现步拒绝，且最后至多 10 个记录残差的
  `(max-min)/max(abs) <= 1e-3`；
- **OTHER_NONCONVERGED**：其余任何轨迹，包括仍在明显下降但撞迭代帽、周期/非周期
  振荡、未到底的停滞。完整记录，不硬塞进地板类。

C2 只诊断这个固定末态的 IBL 重启动。`RESTART_CONVERGED` 不证明松耦合外环能收敛，
`CONTROLLER_FLOOR` 单独也不证明公式层机制；后者必须与 C1 签名联合。

### C3 — 机制判定与 G-DOMAIN

| 总判定 | 必需证据 | 相反结果与不得越界 |
|---|---|---|
| **PASS_KNOWN_FORMULATION_FLOOR** | 两条 xtr030 主腿均同时满足 C1 旧签名 + C2 `CONTROLLER_FLOOR`，且至少一条 xtr005 当前对照也满足同一组合 | 只能结论“这些末态复现已知公式层地板签名”；不证明所有 IBL 状态都有地板，不证明全局无其他不动点 |
| **PASS_RESTART_CONVERGENCE** | 两条 xtr030 主腿均为 `RESTART_CONVERGED`，残差解剖与对照如实报告 | 只结论“两个固定末态都可在原配置重启动收敛”；不声称外层+IBL 已收敛，不自动启动 S5 |
| **UNDEFINED** | 主腿一收敛一不收敛、C1/C2 组合混杂、对照与主腿关系矛盾、出现 `OTHER_NONCONVERGED`，或任何结果不满足上两条 | 不挑选较好看的一腿，不把混合结果报成方向 |
| **FAIL** | §2.2 仪器守卫、实际 BLAS、指纹、verbose 解析或自测失败 | 不作物理归因，不通过修改输入或跳过错腿继续 |

**G-DOMAIN**：本研究定义域只有上述四个 committed-provenance 末态及一次固定重启动。
旧地板幅值不是跨状态等值判据；xtr005 是当前同机对照，不是“必须零残差”的负对照。
任何网格、转捩位置、容差、种子、迭代帽或闭包变化都在定义域外。

## 4. 仪器牙齿、时间盒与 kill criterion

- `--self-test` 在正式读取状态前验证：数组/文件指纹不符会失败；native 残差复核会拒绝
  错值；verbose 的 iter/reject/final 记录能解析；`RESTART_CONVERGED`、
  `CONTROLLER_FLOOR`、`OTHER_NONCONVERGED` 以及 C3 四条判定树合成案例均返回预期。
  任一失败 ⇒ 不启动正式诊断。
- 最多 **4 次末态残差重建 + 4 次 IBL 重启动**，每次仍最多 100 内迭代；
  **0 条 fresh FP/松耦合腿，0 次 dense SVD**。正式诊断总墙钟上限 **20 min**，到时不开新腿。
- 不扫 BLAS、omega、`tol_ds`、IBL 容差、迭代帽、网格、种子、人工扩散或闭包参数；
  不放宽判据，不把地板残差改名为收敛。
- 任何需要修改生产库、D08/D13 值钉或刷新 committed gate 证据才能继续的情况 ⇒ 停止，
  登记为新裁决项，不在本轮顺手实施。

## 5. 交付与验证

新增一次性研究 `bench/studies/ibl_floor_restart/`：`run.py`、README、`.gitignore`、
environment/run_state、输入出处与仪器自测 CSV、逐节点残差 CSV、逐迭代重启动 CSV、
残差/轨迹 PNG、机器可读 criteria/verdict metrics。登记 `bench/studies/README.md`；
phase 7 收口时归档，无测试加载。

收尾运行自测、语法检查和 `tests/F/test_F07_repo_invariants.py`。生产库和测试面不改，
不为本轮重跑贵重全套或 D08/D13；历史套件基线只能注明承接，不能写成当前机全套绿。
判定文件逐条回答 C1–C3，台账新增第 7 轮。若结果为已知地板，S5 继续保持阻塞；
若两腿重启动收敛，也只能另立生产外环复核预注册，不能直接续跑原 2×2。
