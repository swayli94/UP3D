# 20261103-0454 — **预注册：第 8 轮 = xtr030 `A/Psi` 残差的 Jacobian 投影与接受步分解**

**性质**：IBL 独立归因研究；本文件必须先于研究代码和正式测量提交。承接第 7 轮判定 §6
及使用者“下一步仍是针对 A/Psi 残差做 Jacobian 投影和接受步分解”的裁决。
**分支**：`kimi/phase7-ibl-ap-projection`（stack 在 `kimi/phase7-ibl-floor` 顶端；
`origin/main` 已 fetch 且仍是其祖先，不丢弃第 5–7 轮输入出处和证据）。

不改生产库、默认值、闭包、控制器、D08/D13 判据或 committed gate 证据；不启动 S5 2×2；
不新跑 FP/松耦合腿，也不新跑 SA/SST/CFL3D 对照。若以后确需补外部对照，先复用已有
committed SA/SST 数据；必须新跑时沿用 `tools/cfl3d_seq`，或先验证外部构建一致性并记录
二进制 SHA-256、构建参数和运行环境。

---

## 1. 一个问题与两个不得混同的读数

问：**第 7 轮 xtr030 相对 xtr005 增大的 `A/Psi` 残差，是否主要投影到固定边数据
`J_BL,BL` 的候选近零子空间？生产伪时间求解器实际接受的步，又是在消减还是放大各方程的
纯定常残差？**

两个测量轴必须分开：

1. **Jacobian 投影**回答残差相对线性算子range/候选近零左子空间的位置；不直接证明非线性
   不收敛由该子空间造成。
2. **接受步分解**回答生产 Armijo 伪时间 merit 下降时，纯定常六方程实际怎样变化；不接受
   “merit 下降 ⇒ 所有纯残差都下降”的预设。

以下结果均合法：投影强、投影弱、两腿不一致；接受步消减、放大或基本保持 `A/Psi`；某些腿
没有接受步。不得为得到机制故事而追加谱截断、改缩放或换状态。

## 2. 固定状态、算子和出处

### 2.1 四个既有状态与上一轮重启动

沿用第 7 轮完全相同的四个状态：

| 角色 | 腿 | 原生状态 | 接受步复核对象 |
|---|---|---|---|
| 主腿 | `native_xtr030_b1_r1` | S5 保存末态 `state_*.npz` + 首 IBL 冻结边界 `first_*.npz` | `ibl_floor_restart/_work/restart_*.npz` 与 committed restart history |
| 主腿 | `native_xtr030_b16_r1` | 同上 | 同上 |
| 对照 | `control_xtr005_b1` | 同上 | 同上 |
| 对照 | `control_xtr005_b16` | 同上 | 同上 |

- 当前输入文件 SHA-256、`U/phi/delta_star` 数组 SHA、medium 网格 SHA 必须同时符合第 7 轮
  committed `summary.csv/provenance.csv` 和其指向的 S5 committed summary；缺一停止。
- B1/B16 仅为历史出处标签。所有诊断 worker 实际固定 BLAS=1、Numba=16，并从库读
  `blas_thread_count()` 断言；不清缓存，每腿 fresh Python 进程。
- 边数据仍直接使用末态 `ue_surf` 和首 IBL `inflow_mask/inflow_state`，不从最终 FP 的 `phi`
  反推。重建 `F=solver.residual(U)` 的 `maxabs(F)` 必须满足第 7 轮同一
  `5e-12 * recorded_residual` native guard。
- 接受步正式记录结束后，最终 `U` 和最终逐节点纯残差必须与第 7 轮本地保存的
  `restart_*.npz` **逐位一致**；逐迭代纯残差与 committed CSV 的 12 位输出在
  `5e-12 * max(recorded,1)` 内一致，CFL 与 verbose 的 3 位输出按 `5e-4` 相对容差复核。
  任一不符均为仪器 FAIL，不作物理归因。

### 2.2 唯一 Jacobian 与均衡化定义

每腿在原生末态 `U`、固定上述边数据下调用生产 `solver.residual_jacobian(U)`（无
`dt_inv/U_old`），所得稀疏矩阵即旧研究及 `tight_driver` 定义的 `J_BL,BL`。不加入 FP 块、
edge-data 块或伪时间项。

沿用旧 V5 Q3 的一次 row-then-column 2-norm 均衡化，不迭代优化缩放：

```text
r_i = ||J[i,:]||2,   J1 = diag(1/max(r_i,1e-300)) J
c_j = ||J1[:,j]||2,  Js = J1 diag(1/max(c_j,1e-300))
f  = diag(1/max(r_i,1e-300)) F.ravel()
Js = L diag(s) V^T
```

对 `Js` 做一次完整 dense SVD（410 节点、2460 阶；`full_matrices=False`），奇异值按 NumPy
降序。候选子空间 `J20` 事前固定为**最小 20 个奇异值对应的右/左奇异向量**，与旧 V5
Q2/Q5 的 `N_NULL=20` 对齐；不得看完投影后改 K。完整奇异值谱和逐模态读数仍提交，
`sigma[n-20]/sigma[0]`、条件数及 `<1e-6/1e-8/1e-10*sigma_max` 计数用于判断它是否真有
谱软化；“bottom 20”本身不自动称为已证明的数学零空间。

`A/Psi` 行/列集合固定为每节点方程/未知量 `(1,3)`，历史 `(B,delta)` 对照集合固定为
`(2,0)`。对任意目标残差块 `g`（完整 `f`、仅保留 `A/Psi` 行的 `f_AP`、仅保留
`B/delta` 行的 `f_BD`），投影能量比例定义为：

```text
P20(g) = L[:,n-20:n] L[:,n-20:n]^T g
capture20(g) = ||P20(g)||2^2 / ||g||2^2
```

`g=0` 时响亮拒绝，不返回伪造的 0 或 NaN。每模态右侧 `A/Psi` 能量为
`sum(V[i,AP_columns]^2)`；bottom-20 `right_AP_share20` 是这 20 个单位右向量的总
`A/Psi` 能量除以 20。另记录 `right_AP>=0.5` 的模态数及其投影，但不据此改选主判据子空间。

## 3. 测量与事前判据

### C1 — 四状态 Jacobian 投影与主假设

每腿提交：完整奇异值谱、每模态 `sigma_rel`、右六变量能量、全残差/`AP`/`BD` 左投影系数，
以及 `capture20(full/AP/BD)`、`right_AP_share20`、均衡化前后条件数和缩放范围。

逐腿的 **DOMINANT_AP_PROJECTION** 必须同时满足：

1. `right_AP_share20 >= 0.50`：bottom-20 右候选子空间至少一半能量在 `A/Psi`；
2. `capture20(f_AP) >= 0.50`：`A/Psi` 残差至少一半能量落在对应 bottom-20 左子空间；
3. `sigma[n-20]/sigma[0] <= 1e-3`：候选子空间确实位于至少三个相对量级的软端。

0.50 是“主导”的预注册定义，不是收敛容差；`20/2460` 的维数基线和
`capture20/(20/2460)` 仅作读数，不另造显著性声明。`BD` 与 xtr005 是解释对照，不参与
主腿真假投票，不能拿较好看的对照替代主腿条件。

C1 判定：

- **PASS_AP_PROJECTION**：两条 xtr030 均满足全部三条；
- **FAIL_AP_PROJECTION**：两条 xtr030 均至少违反一条；
- **UNDEFINED**：一条满足而另一条不满足；
- **FAIL（仪器）**：SVD 重构/投影划分、native/provenance 或自测失败。与物理假设的 FAIL
  在 criteria 中分行，不混写。

若两腿均 PASS 且 xtr005 也均满足，只能结论“该签名在这些状态共有”，不得写成 xtr030
特异机制；若投影弱，只否定上述**主导投影**定义，不声称排除所有低维/非线性作用。

### C2 — 生产接受步的逐方程分解（4 腿）

用方法包装记录生产 `IBL3Solver.solve` 的原始调用，不重写非线性求解算法：在每次
`residual_jacobian` 保存当前点和伪时间残差，在每次 line-search `residual` 调用保存 trial
伪时间残差并额外计算同点纯定常残差。`omega` 由同一 trial 序列从 1 开始逐次减半；接受
条件仍由生产代码执行，记录器只用同一 Armijo 表达式标记事件。额外残差计算不得改变步长、
接受决策或返回状态，§2.1 的逐位重启动复核是硬守卫。

从原生末态按第 7 轮原配置各重启动一次（`tol=1e-9/max_iter=100`、其余默认）。每个接受步
提交：迭代号、CFL、omega、六方程 `F_before/F_after/dF`、六变量 `dU`、伪时间
`R_before/R_after`，以及完整/`AP`/`BD` 的 L2/inf 前后范数。拒绝迭代保留次数、trial 数和
CFL，不输出 12 次 trial 的全部大数组。

逐接受步分类仅描述实际纯残差：

- `AP_REDUCED`：`||F_AP_after||2 < ||F_AP_before||2 * (1-1e-12)`；
- `AP_AMPLIFIED`：`||F_AP_after||2 > ||F_AP_before||2 * (1+1e-12)`；
- `AP_FLAT`：其余（仅数值不动）；没有接受步的腿标 `NO_ACCEPTED_STEPS`，不投票方向。

每腿再汇总所有接受步首尾的 `sum ||F_AP||2^2` 变化及逐事件计数；这些是 C2 的 RECORDED
读数，不反改 C1，也不把接受伪时间 merit 的步改称“纯残差下降步”。至少一条主腿无接受步
不是仪器失败；两腿都无接受步时接受步机制总判定为 UNDEFINED，投影仍可独立判定。

### C3 — 合成结论与 G-DOMAIN

| 总判定 | 条件 | 允许结论与边界 |
|---|---|---|
| **PASS** | C1=`PASS_AP_PROJECTION`，四腿出处/native/replay 全过，C2 完整记录 | 仅“两个固定 xtr030 末态满足主导 AP 投影签名”；结合 C2 描述接受步，不自动宣称已证明非线性根因 |
| **FAIL** | C1=`FAIL_AP_PROJECTION`，仪器全过 | 仅“两个固定末态均不满足主导 AP 投影定义”；不排除弱投影、其他 K 或其他机制 |
| **UNDEFINED** | C1 主腿不一致，或 C2 因无接受步而无法回答且总叙述必须依赖接受步 | 报告每腿原始读数，不合并成方向 |
| **FAIL（仪器）** | 指纹、native、SVD、生产 replay、自测任一失败 | 停止派生解释，保留失败证据，不改输入/阈值换绿 |

**G-DOMAIN**：定义域只有这四个固定末态、上述固定边数据、一次均衡化、K=20 的完整 SVD
和一次原配置重启动。结论不外推到其他网格、转捩、种子、BLAS 吸引子、缩放、K、外层耦合或
所有 IBL 地板。C1 的反向结果（右子空间非 AP、AP 投影落在其余 2440 维、谱不软）均有明确
落点；C2 的反向结果（AP 实际下降、无接受步）也完整保留。

## 4. 仪器牙齿、时间盒与 kill criterion

`--self-test` 必须在正式读取状态前通过并提交 CSV：

- 非对角、非单位行列缩放的合成矩阵经定义均衡化后，完整 SVD 重构误差 `<=1e-10` 相对；
  投影能量与补空间能量之和为 1，AP/BD 行列 mask 不被展平顺序混淆；零目标响亮失败；
- 用 `IBL3Solver.solve` **原函数**驱动一个 6 自由度线性假求解器：一个正常 Newton 步必须
  被记录为接受且纯残差到零；一个恒定坏 Jacobian 的反向案例必须全部拒绝且零接受事件；
  记录器返回状态/轨迹与未包装生产调用逐位一致；
- Armijo 标记、`AP_REDUCED/AP_AMPLIFIED/AP_FLAT/NO_ACCEPTED_STEPS`、C1/C3 判定树均有
  正反案例；错误指纹、错误 replay 和错误 SVD 划分必须触发 FAIL。

正式预算最多 **4 次完整 dense SVD + 4 次原配置重启动**，0 条新耦合腿、0 次迭代特征值/
奇异值求解、0 组新参数。每腿 worker 上限 10 min，正式累计墙钟上限 30 min；每腿顺序执行，
SVD 数组不跨腿驻留。超出预算、需要改生产代码/换 K/换缩放/追加谱计算才能继续 ⇒ 照章停止，
标 UNDEFINED 或 RECORDED 阻塞，不调参换绿。

## 5. 交付与验证

新增一次性研究 `bench/studies/ibl_ap_projection/`：run.py、README、`.gitignore`、
environment/run_state、仪器自测、出处、投影谱/模态/摘要、逐接受步和逐迭代 CSV、投影与
接受步 PNG、机器可读 criteria/verdict metrics；登记 `bench/studies/README.md`，phase 7
收口归档，无测试加载。大数组仅进 gitignored `_work/`，提交 CSV/PNG 是结论证据。

收尾运行自测、语法检查和 `tests/F/test_F07_repo_invariants.py`。生产库与测试面不改，不为
本轮重跑贵重全套、D08/D13 或外部 SA/SST；历史套件基线只注明承接。判定文件逐条回答
C1–C3，台账新增第 8 轮，并明确是否值得另立“非线性根因/修复”研究；本轮本身不做修复。
