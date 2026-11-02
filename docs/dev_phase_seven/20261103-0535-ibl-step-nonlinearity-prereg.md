# 20261103-0535 — **预注册：第 9 轮 = xtr030 强回溯接受步的非线性余项与左补空间逐方程分解**

**性质**：IBL 独立归因研究；本文件必须先于研究代码和正式测量提交。承接第 8 轮判定
[20261103-0516-ibl-ap-projection-verdict.md](20261103-0516-ibl-ap-projection-verdict.md) §6。
**分支**：`kimi/phase7-ibl-ap-projection`（stack 在第 8 轮 `ed35cca` 顶端；`origin/main`
已 fetch 且仍是其祖先，不丢弃既有输入出处和证据）。

不改生产库、默认值、闭包、控制器、D08/D13 判据或 committed gate 证据；不启动 S5 2×2；
不新跑 FP/松耦合腿，不新跑 SA/SST/CFL3D 对照，不重跑第 7/8 轮四腿重启动，也不修改其
CSV/PNG。除 §2.2 所述一次不可替代的投影基重建外，不重复第 8 轮四腿 SVD 研究，不换 K、
缩放或状态寻找阳性结果。

---

## 1. 一个问题与两个不能混同的分母

问：**在第 8 轮唯一存在接受事件的 `native_xtr030_b1_r1` 上，实际强回溯接受步的纯定常
残差响应，在一阶模型 `F+JdU` 与非线性实际 `F(U+dU)` 之间有多大偏差；该偏差在旧
bottom-20 左补空间中由哪些方程承载？**

对第 `k` 个接受事件定义：

```text
P_k = J(U_k) dU_k                         # 一阶预测响应
D_k = F(U_k+dU_k) - F(U_k)               # 实际响应
Fpred_k = F(U_k) + P_k                   # 一阶预测的末残差
N_k = F(U_k+dU_k) - Fpred_k = D_k - P_k  # 非线性余项 = 两模型末残差之差
```

必须同时报告两种相对读数，不得只报好看的一个：

1. **末残差层面**：`||N|| / (||Fpred|| + ||F1||)`；它可能被较大的 `F0` 稀释，不能单独
   用来宣称一阶响应准确。
2. **步响应层面**：`||P||`、`||D||`、`||N||` 及其比值；这是本轮主判据。小步长不自动
   等于非线性余项小，伪时间 Armijo 被接受也不自动等于纯定常一阶模型方向正确。

以下结果均合法：一阶响应主导、非线性余项主导、混合/相消、逐方程结论不同、预测与实际
残差方向不一致。物理分类不是仪器 PASS/FAIL；不得为得到机制故事而改事件集合、投影顺序、
阈值或缩放。

## 2. 固定事件、状态重放与左补空间

### 2.1 唯一事件集合与离线接受路径重放

只使用第 8 轮 committed `accepted_event_metrics_native_xtr030_b1_r1.csv` 的 **18 个实际接受
事件**。其迭代号事前固定为：

```text
(0, 1, 2, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 21, 23, 25, 27)
```

“强回溯”事前定义为 `omega <= 2^-8 = 0.00390625`；上述 18 步均满足，逐步使用 committed
`cfl/omega`，不重新选择步长。其余三腿第 8 轮均为 `NO_ACCEPTED_STEPS`，本轮不为制造事件
而重启动或改控制器；它们不投票，结论不得外推到四腿。拒绝 trial 的大数组未保存，本轮不
补跑、不重建拒绝 trial，也不把它们混入接受事件。

从 S5 原生 `U0` 开始，按 committed 事件顺序做**离线接受路径重放**，不调用
`IBL3Solver.solve`，不运行迭代控制器或线搜索：

```text
U_old = 上一个接受事件的 U_before（首步 = U0）
dt_inv = q / (cfl_k * h)
Rpt0, Jpt = solver.residual_jacobian(U_before, dt_inv=dt_inv, U_old=U_old)
delta = spsolve(Jpt, -Rpt0.ravel())
scaled_dU = omega_k * delta.reshape(n,6)
U_after = U_before + scaled_dU
dU = U_after - U_before                   # 浮点状态的实际增量
Rpt1 = solver.residual(U_after, dt_inv=dt_inv, U_old=U_old)
F0 = solver.residual(U_before)
F1 = solver.residual(U_after)
_, J = solver.residual_jacobian(U_before)   # 纯定常；无 dt_inv/U_old
```

**首次正式重放的仪器勘误（未进入 SVD/物理测量）**：初版代码把 `omega*delta` 当作 committed
`dU`，但第 8 轮实际导出的是浮点 `U_after-U_before`。首次重放的 `F0/F1/Rpt0/Rpt1/dF`
回写 mismatch 全为 0、末态 `U/F` 均逐位一致，仅 `dU` 因加法舍入出现每步 1772–2178 个末位
字符串差异（绝对差最大 `2.64e-14`）。现将 `scaled_dU` 仅用于形成 `U_after`，实际 `dU`
和后续 `JdU` 一律使用 `U_after-U_before`；判据、容差、事件和预算不变，失败证据以
`instrument_failure_du_rounding_*` 留档，不用完整重启动或低精度状态替代。

拒绝迭代不改变 `U/U_old`，故不在离线路径中重算；其省略是否影响数值只能由下面的硬守卫
判定，不能靠口头假设。每一步的 `F0/F1/Rpt0/Rpt1/dU` 必须按第 8 轮相同 `f12` 格式与
committed `accepted_steps_native_xtr030_b1_r1.csv` **逐单元格字符串一致**，`dF=F1-F0`
也必须回写一致。末态 `U` 和最终逐节点纯残差必须与第 7 轮 `restart_*.npz` 逐位一致。

任一不符均为 **FAIL_INSTRUMENT**；不得退回低精度 CSV 累加状态、不得改容差、不得加跑一次
完整重启动来替换失败重放。原生输入文件、数组、medium 网格 SHA 与 native residual 继续满足
第 7/8 守卫；本轮另记录全部输入 CSV/NPZ 的 SHA-256。

正式预算因此只有 **18 次已知 CFL/omega 的稀疏 Newton 方向重算及相应残差/Jacobian 调用**，
0 次生产 `solve` 重启动、0 条新耦合/外部对照腿、0 个新接受或拒绝事件。

### 2.2 固定左投影基及一次重建预算

第 8 轮只提交了谱和投影系数 CSV，没有保存左奇异向量；要从新差值的普通方程空间进入其
**其余 2440 维左补空间**，投影基是不可替代输入，不能用 `B/delta` 等普通方程补集冒充。

因此只在 `native_xtr030_b1_r1` 原生 `U0` 重建一次与第 8 轮完全相同的
`J0 = residual_jacobian(U0)`、一次 row-then-column 2-范数均衡化和完整 dense SVD：

```text
Js0 = diag(row_scale0) J0 diag(col_scale0) = L0 diag(s0) V0^T
L20 = L0[:, 2460-20:2460]
FULL0(g) = row_scale0 * g.ravel()
C0(g) = (I - L20 L20^T) FULL0(g)
B0(g) = L20 L20^T FULL0(g)
```

这次计算**只重建固定投影工具**，不重新投票第 8 轮 C1，不重跑其余三腿 SVD，不在 18 个
接受状态上追加 SVD，也不改变 K=20 或缩放。完整谱须与 committed
`projection_modes_native_xtr030_b1_r1.csv` 的 `sigma` 按 `f12` 回写一致；SVD 重构、投影
划分误差均 `<=1e-10`，原生 `full/AP/BD` 的 `capture20` 与第 8 轮摘要按相同 `f12`
回写一致。若 LAPACK 符号选择不同，只比较符号不变的投影/奇异值，不要求单个左向量逐位同号。
任何不符为 FAIL_INSTRUMENT，禁止换 eigensolver 或迭代谱方法继续。

`row_scale0/L20` 随后对所有接受步冻结。对任一二阶数组 `X` 和源方程 `e`，先取源方程
`X_e = X*M_e`，再缩放和投影；本轮主分解是 **source-equation-first**：

```text
E_region,e(X) = ||region0(X_e)||2^2,  region ∈ {FULL, BOTTOM20, COMPLEMENT}
```

`M_e` 是每节点六方程展平顺序中的第 `e` 个方程 mask。先投影再取输出方程是另一个量，不得
与上述定义混写。每事件每源方程验证 `FULL = BOTTOM20 + COMPLEMENT`（相对误差 `<=1e-10`），
并验证六源方程在同一区域的**投影向量之和**等于直接投影 full 向量。零目标可以记录零能量，
但不得产生伪造的 0/NaN 比值或方向票。

**正式测量前勘误（合成自测发现，尚无正式数值）**：原文误写“六源方程能量之和在三个区域都
等于 full 能量”。投影是线性的，故 `sum_e region0(X_e)=region0(X)`；但不同源方程投影后的
向量一般不正交，`sum_e||region0(X_e)||²` 在 BOTTOM20/COMPLEMENT 通常不等于
`||region0(X)||²`，交叉项不能丢弃。自测已用已知正交基使错误能量守卫变红。现将硬守卫改为
上述**向量重构**恒等式；逐源方程能量仍按定义提交，但不得当作可加份额。每事件/区域的
`P/D/N` 另提交 `source_energy_sum=sum_e E_region,e`、直接 full 能量及
`cross_energy=full_energy-source_energy_sum`，event summary 的分类使用直接 full 向量，
equation summary 的分类使用各源方程自身跨步能量，两者不混写。

## 3. 测量、分类与 G-DOMAIN

### C1 — 逐步、逐方程、逐区域的非线性余项

对 18 步 × 6 方程 × 3 区域（FULL/BOTTOM20/COMPLEMENT）提交：

- `||P||, ||D||, ||N||, ||F0||, ||Fpred||, ||F1||` 及能量；
- `P·D`、余弦（两范数均非零时）、`||N||/(||Fpred||+||F1||)`；
- 响应尺度 `S=||P||+||D||` 和 `||N||/S`（`S>0` 时），以及 `||N||/||P||`
  （`||P||>0` 时）。

每事件/方程/区域按下列互斥规则分类：

1. `EXACT_ZERO`：`S==0`，此时必须同时有 `||N||==0`，否则为仪器失败；不投票方向。
2. `LINEAR_DOMINATED`：`S>0` 且 `||N|| <= 0.10 ||P||`；0.10 是事前固定的“次级
   修正”定义，不是收敛容差。
3. `NONLINEAR_DOMINATED`：不满足上条，且 `||N|| >= ||P||` 与 `||N|| >= ||D||`
   同时成立；即余项不弱于一阶和实际响应二者。
4. `MIXED`：其余全部情形，包括中间幅值及 `P/N` 相消使实际响应很小的情况。

跨 18 步汇总使用每方程/区域的 `sum||P||², sum||D||², sum||N||²`，以同一规则给出
`POOLED_*` 分类，同时提交 18 步各类计数、比值中位数/最小值/最大值和最大单步能量份额。
汇总不得抹掉逐事件反向案例；单步贡献和逐步类别是正式读数，不得只挑支持某一故事的步。

### C2 — `F+JdU` 与实际末残差的方向是否一致

在每事件/源方程/区域分别比较一阶预测末残差和实际末残差相对 `F0` 的范数变化，使用与第
8 轮一致的 `1e-12` 相对死区：

```text
REDUCED:  after < before * (1-1e-12)
AMPLIFIED: after > before * (1+1e-12)
FLAT:      其余
```

组合为 `BOTH_REDUCED / BOTH_AMPLIFIED / BOTH_FLAT / DIRECTION_MISMATCH /
PRED_FLAT_ONLY / ACTUAL_FLAT_ONLY`。`before=0` 时禁止构造归一化下降率；若 after 非零，
按 AMPLIFIED 记录，若二者均零为 BOTH_FLAT。逐步导出预测/实际分类及余弦，不根据总能量
方向反推单步方向。

C2 是方向测量，不把 BOTH_REDUCED 写成收敛，也不把 DIRECTION_MISMATCH 自动写成非线性
根因；幅值解释必须回到 C1 的 `P/D/N`。

### C3 — 判定树与允许结论

- **I0 PASS**：自测、输入/native、18 步离线重放、末态逐位、一次 SVD 回写、投影/能量划分
  全部通过；否则 **FAIL_INSTRUMENT**。
- **C1/C2 RECORDED**：I0 通过且 18×6×3 数据完整，按上述标签报告线性主导、非线性主导、
  混合和方向一致/不一致。物理反向结果不记仪器 FAIL。
- **总判定 RECORDED**：I0/C1/C2 完整；本轮是预注册诊断，不设“必须非线性主导才算 PASS”
  的假靶子。
- **总判定 UNDEFINED**：固定事件不存在或全部 `EXACT_ZERO`，无法形成响应比较；零接受步
  的其他腿不补跑、不投票。
- **总判定 FAIL（仪器）**：任一硬守卫失败或需要超出预算才能获得答案；保留失败证据，停止
  派生解释。

**G-DOMAIN**：定义域只有 xtr030 B1 原生末态出发的 18 个已接受强回溯事件、固定边数据、
原生点冻结的 bottom-20 左补空间、上述 source-equation-first 投影和纯定常 Jacobian。结论不
外推到拒绝 trial、其他三腿、其他 K/缩放/网格/转捩/BLAS 吸引子、外层耦合或所有 IBL 地板。
一阶主导、非线性主导和 MIXED 都有明确落点；即使非线性余项主导，也只说明这些步的一阶模型
误差结构，不单独证明它就是残差地板根因；即使一阶主导，也不证明继续迭代会收敛。

## 4. 仪器牙齿、时间盒与 kill criterion

`--self-test` 必须在读取正式状态前通过并提交 CSV，至少覆盖：

- 独立二次残差 `F(U)=J(U-Uref)+0.5Q(U)` 的解析 `P/D/N` 与数值计算一致；线性案例
  `Q=0` 得零余项，构造的非线性主导、MIXED、EXACT_ZERO 和零分母案例均分类正确；
- 已知正交 bottom 基上 FULL=BOTTOM+COMPLEMENT，六源方程投影向量求和恒等，并构造
  非零交叉项使旧的错误能量求和守卫变红；mask 展平顺序和“先 mask 后投影”不被换成另一种
  定义；错误投影或错误余项必须触发 FAIL；
- 用生产 `IBL3Solver.solve` 驱动合成求解器生成短轨迹后，离线重放可复现接受事件；改一个
  committed CFL/omega、末态或 CSV 单元格必须让重放/回写守卫变红，不得静默换路径；
- `BOTH_REDUCED/BOTH_AMPLIFIED/BOTH_FLAT/DIRECTION_MISMATCH` 及零 before 正反案例，
  C3 的 RECORDED/UNDEFINED/FAIL_INSTRUMENT 判定树均有牙齿。

正式累计上限 **1 次原生点 dense SVD + 18 次离线稀疏方向重算**，0 次生产 solve 重启动、
0 次其余三腿 SVD、0 条新耦合/外部对照腿、0 组新参数。单 worker 上限 10 min，累计上限
15 min。超出预算、重放不逐位、谱不能回写或需要改生产代码/改事件/改 K/改缩放才能继续，
照章停止并标 FAIL_INSTRUMENT 或 UNDEFINED；不以完整重启动、CSV 低精度状态或另一套投影
定义作为未注册退路。

## 5. 交付与验证

新增一次性研究 `bench/studies/ibl_step_nonlinearity/`：run.py、README、`.gitignore`、
environment/run_state、仪器自测、出处/重放/投影守卫、逐节点六方程 `P/D/N/Fpred/F1`、
逐事件×方程×区域 CSV、跨步摘要、机器可读 criteria/verdict metrics 及非线性分解 PNG；
登记 `bench/studies/README.md`，phase 7 收口归档，无测试加载。大数组仅进 gitignored
`_work/`，提交 CSV/PNG 是结论证据，第 7/8 轮已有文件不改。

收尾运行自测、语法检查和 `tests/F/test_F07_repo_invariants.py`。生产库与测试面不改，不为
本轮重跑贵重全套、D08/D13 或外部 SA/SST；历史套件基线只注明承接。判定文件逐条回答
C1–C3，台账新增第 9 轮，并明确非线性余项结果是否足以支持另立根因/修复研究；本轮本身
不做修复。
