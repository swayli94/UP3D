# 20261103-0442 — **判定（第 7 轮，IBL 独立研究）：四状态均控制器到底且无重启动收敛；xtr030 是否属于已知公式层地板仍 UNDEFINED**

轮次文件名：`20261103-0442-ibl-floor-verdict.md`
阶段 / 门：phase 7 / 无门（研究）——预注册：
[20261103-0420-ibl-floor-prereg.md](20261103-0420-ibl-floor-prereg.md)
测量代码：`5e9bea1`（分支 `kimi/phase7-ibl-floor`；本判定与证据随后提交）
证据：`bench/studies/ibl_floor_restart/`（run.py + README + 自测/出处/摘要/逐节点/逐迭代 CSV + 2 PNG）

---

## 1. 目的

回答 S5 当前 xtr030 B1/B16 末态的 IBL 内层不收敛，是否复现已登记的 TE 带 `(B,δ)`
公式层地板与伪时间控制器到底签名；在同一末态、同一边数据下原样重启动，它究竟收敛、
撞地板，还是走出另一种轨迹。预注册要求残差位置、六方程份额和闭包地板活动集共同成立，
不得只凭残差约为 1e-6 或 CFL 到底就命名为同一公式层地板。

本轮只做识别与归因，不修闭包、不调控制器、不改生产库/默认值/D08/D13，也不启动 S5 2×2。

## 2. 改了什么

- 新增一次性研究 `bench/studies/ibl_floor_restart/`：四个 S5 既有末态的出处核验、
  生产 `IBL3Solver` 纯定常残差重建、逐节点/逐方程分解、同状态同边数据重启动、
  verbose CFL/拒绝轨迹和预注册判定树。
- 新增 20 项合成仪器自测，登记 `bench/studies/README.md`。
- **生产库零改动、测试零改动、默认值零改动、D08/D13 判据零改动、committed gate
  证据零刷新。零条 fresh FP/松耦合腿，零次 dense SVD，S5 2×2 未启动。**

## 3. 怎么测的

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
python bench/studies/ibl_floor_restart/run.py --self-test
python bench/studies/ibl_floor_restart/run.py
python bench/studies/ibl_floor_restart/run.py --analyze-only
python -m pytest tests/F/test_F07_repo_invariants.py -q
```

- 主腿为 xtr030 B1/B16 的 S5 `r1` 末态，xtr005 B1/B16 为当前同机对照；
  B1/B16 只表示输入历史出处，所有诊断 worker 实际 BLAS 均为 1、Numba=16。
- 直接使用末态保存的 `ue_surf/U` 和各腿首 IBL 冻结边界；没有从最终 FP 重解后的 `phi`
  反推边数据。`U/phi/delta_star` 数组 SHA 和网格 SHA 均与 committed S5 summary 相同。
- 每腿先复算原生末态 `maxabs(F)`，再按原 `ibl_tol=1e-9/max_iter=100` 从该末态重启动一次。
  不增加迭代、不改容差、不换种子。判据和阈值均在正式测量前的 `2ba8658` 预注册中冻结。

## 4. 结果

### 4.1 I0：出处和原生重建全部通过

- 自测 **20/20**，包括指纹错误、残差位置错误、活动闭包地板、控制器分类、假收敛和
  完整判定树的反向案例，见 `instrument_self_test.csv`。
- 四腿 `provenance_guard/native_guard` 全部通过。重建残差相对 committed history 的
  最大相对差为 **1.1504e-13**，只相当于 12 位 CSV 的舍入；见 `provenance.csv`。
- 四腿 `DELTA_MIN/RE_D_MIN` 活动集均为空，最小 `delta=1.78162e-6`，
  最小 `re_d=2.21028`，不是这两个闭包数值地板直接夹死。

### 4.2 C1：xtr005 满足旧残差签名；两条 xtr030 主腿均不满足合并份额阈值

| 腿 | 原生 `maxabs(F)` | 全局最大 `x/c` | TE 能量份额 | `B+delta` 合并份额 | 旧签名 |
|---|---:|---:|---:|---:|---:|
| xtr030 B1 | 1.562489668912e-6 | 0.980207 | 0.803024 | **0.882496** | **0** |
| xtr030 B16 | 6.572426677031e-7 | 0.970307 | 0.717006 | **0.861668** | **0** |
| xtr005 B1 | 2.091980536021e-6 | 0.990105 | 0.674717 | **0.986194** | 1 |
| xtr005 B16 | 1.918127438964e-6 | 0.990105 | 0.654716 | **0.979630** | 1 |

四腿全局最大都在上翼面 TE 带，前两大方程均为 `B/delta`，但预注册要求合并份额
`>=0.90`。两条 xtr030 主腿分别为 0.882496/0.861668，**不满足旧签名**；
不得事后把阈值降到 0.86 来换 PASS。逐节点证据见 `residual_anatomy_*.csv`，
汇总见 `residual_summary.csv` 和 `residual_map.png`。

完整 L2 份额进一步显示 xtr030 的残差更混杂：

| 腿 | `delta` | `A` | `B` | `Psi` |
|---|---:|---:|---:|---:|
| xtr030 B1 | 0.508432 | **0.402899** | 0.721316 | **0.242444** |
| xtr030 B16 | 0.510775 | **0.380425** | 0.693961 | **0.335745** |
| xtr005 B1 | 0.546905 | 0.113943 | 0.820654 | 0.119847 |
| xtr005 B16 | 0.536554 | 0.150806 | 0.819624 | 0.132424 |

这些数由 committed `residual_anatomy_*.csv` 的六方程列直接重算。它们只支持
“xtr030 的 `A/Psi` 份额比当前 xtr005 对照更大、残差结构更混杂”。本轮没有测
Jacobian 谱、均衡化算子或残差投影，**不得进一步写成已证明旧 `(A,Psi)` 近零空间是原因**。

### 4.3 C2：四腿原样重启动均不收敛，全部满足 `CONTROLLER_FLOOR`

| 腿 | n_iter | 初始 → 最终纯残差 | 最终 CFL | 步拒绝 | 分类 |
|---|---:|---:|---:|---:|---|
| xtr030 B1 | 39 | 1.562489668912e-6 → 1.560038984365e-6 | 1.000e-3 | 20 | `CONTROLLER_FLOOR` |
| xtr030 B16 | 11 | 6.572426677031e-7 → 6.572426677031e-7 | 1.000e-3 | 10 | `CONTROLLER_FLOOR` |
| xtr005 B1 | 11 | 2.091980536021e-6 → 2.091980536021e-6 | 1.000e-3 | 10 | `CONTROLLER_FLOOR` |
| xtr005 B16 | 11 | 1.918127438964e-6 → 1.918127438964e-6 | 1.000e-3 | 10 | `CONTROLLER_FLOOR` |

末 10 个记录残差的相对范围均为 0；所有腿都没有达到生产 IBL 收敛判据。
这闭合的是**四个固定末态的原样重启动轨迹**：控制器到底且纯残差平台化，
不是“仍在下降、只是 100 步不够”。逐迭代残差/CFL/拒绝记录见 `restart_history_*.csv`，
摘要和图见 `restart_summary.csv`、`restart_trajectories.png`。

### 4.4 C3：两条预注册 PASS 路径都不成立

- `PASS_RESTART_CONVERGENCE` 不成立：两条 xtr030 主腿均未收敛。
- `PASS_KNOWN_FORMULATION_FLOOR` 不成立：虽然两主腿都是 `CONTROLLER_FLOOR`，
  但二者的 C1 旧残差签名均为 0。
- 因此按预注册总判定为 **UNDEFINED**。xtr005 两腿的“旧签名 + 控制器到底”组合均成立，
  只能作为当前对照证据，不能替 xtr030 补足缺失的签名。

机器可读结果见 `criteria.csv`、`verdict_metrics.csv`。四次正式重启动累计
**14.1588 s**，低于 20 min 时间盒；没有追加任何参数腿。

## 5. 判定

**UNDEFINED（xtr030 是否属于同一已知公式层地板不可判；四状态原样重启动均不收敛、
控制器到底且残差平台化已闭合）。**

一句话：**当前 xtr005 对照复现旧 TE `(B,δ)` 残差签名；xtr030 虽同样尾缘主导、
CFL 到底且不收敛，但 `B+delta` 份额未达到事前阈值且 `A/Psi` 份额更大，故不能
把两种状态合并成同一个已闭合机制，更不能靠放宽阈值宣称 S5 的前置问题已解决。**

对 S5 的直接后果：仍没有外层与 IBL 内层都收敛的两条腿，**S5 维持 UNDEFINED，
原 2×2 继续不启动**。本轮也没有证明“公式层地板在 xtr030 上永远不可越过”。

## 6. 下一步

- 若继续 IBL 归因，下一轮应另写预注册，针对 xtr030 相对 xtr005 增大的 `A/Psi`
  残差份额设计可证伪实验，例如均衡化 `J_BL,BL` 的残差投影与接受步前后逐方程分解；
  先证明仪器和定义域，再决定是否值得做新的谱计算。不得把本轮阈值改松。
- 若不再投入 IBL，按 phase 7 范围转入 S6 小尾巴或阶段收口；能力表述应写成
  “xtr005 旧地板签名复现，xtr030 混合残差结构的机制未判”，不是“xtr030 已证明同一地板”。
- D08/D13 的 BLAS=1 值钉和 committed 证据保持不动。

## 基线行

研究轮不改生产库和测试面，**不重跑贵重全套或 D08/D13**。自测 **20/20**；
F07 **7 passed**。历史开发机套件基线承 **839 / 0 / 35 / 2**，不是当前迁移机全套实测；
未运行 `bench/baseline_*.csv` 再基线，也不声称任何基线数值移动。
