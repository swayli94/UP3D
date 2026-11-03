# 20261103-0516 — **判定（第 8 轮，IBL 独立研究）：`A/Psi` 右软子空间存在，但残差几乎全在其左补空间；主导投影假设 FAIL**

轮次文件名：`20261103-0516-ibl-ap-projection-verdict.md`
阶段 / 门：phase 7 / 无门（研究）——预注册：
[20261103-0454-ibl-ap-projection-prereg.md](20261103-0454-ibl-ap-projection-prereg.md)
（`883b071`；raw 条件数预算勘误 `108a7a1`）
测量代码：`e94b07c`（分支 `kimi/phase7-ibl-ap-projection`；环境指纹见 environment.json；
后续提交只调整投影图布局并写判定，不重算数值）
证据：`bench/studies/ibl_ap_projection/`（run.py + README + 自测/出处/谱/模态/接受步 CSV + 2 PNG）

---

## 1. 目的

承接第 7 轮“xtr030 的 `A/Psi` 原始残差份额更大，但未做 Jacobian 投影”的开放问题，
在四个固定末态上检验：均衡化 `J_BL,BL` 的 bottom-20 候选近零子空间是否由 `A/Psi`
主导，以及 `A/Psi` 残差是否主要落入对应左子空间。并分解生产伪时间求解器实际接受步，
检查伪时间 merit 下降时六方程纯定常残差到底如何变化。

本轮不改生产库、默认值、闭包、控制器、D08/D13 或 committed gate 证据；不启动 S5 2×2，
不新跑 FP/松耦合或 SA/SST/CFL3D 对照，不换 K、缩放或状态寻找阳性结果。

## 2. 改了什么

- 新增一次性研究 `bench/studies/ibl_ap_projection/`：四状态出处/native 守卫、一次
  row-then-column 均衡化、每腿一次 2460 阶完整 dense SVD、六变量右能量和 full/AP/BD
  左投影、生产 `IBL3Solver.solve` 接受步包装记录及判定树。
- 每个接受步导出 `(iteration,node,equation)` 的 `F_before/F_after/dF`、伪时间
  `R_before/R_after` 和 `dU`；每腿另存全部迭代/trial/CFL 记录。xtr030 B1 的逐节点接受步
  表为 44,280 行，其余三腿只有表头（确实没有接受事件）。
- 新增 26 项合成仪器自测，登记 `bench/studies/README.md`。
- **生产库零改动、测试零改动、默认值零改动、D08/D13 和 committed gate 证据零改动。**
  正式测量为 4 次 SVD + 4 次固定重启动；0 条 fresh 耦合/外部对照腿，0 组新参数。

### 2.1 启动和仪器勘误（不进入物理判定）

- 首次误用 anaconda `base`，worker 的 BLAS 出处为 `None`，在 SVD/重启动前停止；
  `startup_failure_base_*` 留档。改用真正的 `up3d` 后，在实际暖调用前即读回 BLAS=1。
- 第二次首腿 SVD 和求解完成后，接受步导出函数遗漏 `numpy` 导入而停止；
  `instrument_failure_name_error_*` 留档。该腿未通过 replay、未进入物理判定。
  随后补 `trace_export_normal` 自测（25→26 项）并从头重跑四腿。
- 预注册代码前勘误：raw `J` 条件数会偷加第二次 SVD，故正式只做均衡化 `Js` 的完整谱；
  raw 谱/条件数没有计算或报告。

## 3. 怎么测的

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
python bench/studies/ibl_ap_projection/run.py --self-test
python bench/studies/ibl_ap_projection/run.py
python bench/studies/ibl_ap_projection/run.py --analyze-only
python -m pytest tests/F/test_F07_repo_invariants.py -q
```

- 四个输入与第 7 轮完全相同；当前文件/数组/网格 SHA 同时绑定第 7 轮和 S5 committed CSV。
  B1/B16 只表示历史出处，所有正式 worker 实际 BLAS=1、Numba=16。
- 每腿以生产 `solver.residual_jacobian(U)`（无伪时间项）构造 `J_BL,BL`，做一次旧 V5 定义的
  行列 2-范数均衡化。K=20、0.50 主导阈值及 `sigma20_rel<=1e-3` 均在测量前冻结。
- 接受步通过包装生产 `solve` 的方法记录，不重写求解算法。最终 `U` 和逐节点残差必须与第 7 轮
  本地缓存逐位一致；committed 历史只按已登记的输出舍入容差复核。
- 图布局在正式测量后由两面板改为“谱 / 右能量 / 对数残差投影”三面板；只运行
  `--analyze-only`，所有数值 CSV 和判定不变，未重跑 SVD 或重启动。

## 4. 结果

### 4.1 I0：仪器、出处和 replay 全部通过

- 自测 **26/26**，覆盖 SVD 重构、投影补空间、AP/BD 展平 mask、已知 bottom 投影、
  正常接受、全部拒绝、零接受步、CSV 导出及判定树反向案例。
- 四腿 `provenance_guard/native_guard/projection_guard/replay_guard` 全为 1。
- SVD 重构相对误差最大 **4.3963e-15**，投影划分误差最大 **6.6613e-16**，右变量能量
  划分误差最大 **3.5527e-15**。
- 四腿最终 `U` 和最终逐节点残差均与第 7 轮**逐位一致**；历史最大归一化差
  **4.8027e-19**，CFL 最大相对差 **3.2011e-4**，在 verbose 三位输出的预注册容差内。

### 4.2 C1：右软子空间确由 `A/Psi` 主导，但 `A/Psi` 残差不落入其左子空间

| 腿 | bottom-20 右 `A/Psi` 份额 | `sigma20/sigma0` | 均衡化 cond | `capture20(AP)` | capture/(20/2460) | 主导投影 |
|---|---:|---:|---:|---:|---:|---:|
| xtr030 B1 | 0.616273 | 4.18621e-4 | 5.24392e4 | **1.040857e-4** | 0.0128025 | **0** |
| xtr030 B16 | 0.602863 | 2.77925e-4 | 4.82378e4 | **3.982458e-5** | 0.00489842 | **0** |
| xtr005 B1 | 0.602937 | 2.31428e-4 | 3.99018e5 | 1.783390e-4 | 0.0219357 | 0 |
| xtr005 B16 | 0.600556 | 1.79202e-4 | 2.03751e5 | 1.712782e-4 | 0.0210672 | 0 |

四腿的前两条签名都成立：右 `A/Psi` 份额均超过 0.50，`sigma20_rel` 均小于 1e-3；
主腿分别只有 2/20、2/20 个 bottom 模态不是右 `A/Psi` 主导（对照各 1/20）。但第三条的
`capture20(AP)>=0.50` **四腿全部失败**。两条主腿的 `A/Psi` 残差能量分别有
**99.98959% / 99.99602%** 位于其余 2440 维左补空间。

维数基线为 `20/2460=0.00813008`。四腿投影不仅低于主导阈值，还低于该基线；主腿 B1/B16
分别约为基线的 1.28% / 0.49%。完整残差 `capture20(full)` 也只有
`3.58037e-5…1.59235e-4`，并非 `A/Psi` 单独选择造成的视觉现象。四腿在均衡化后
`<1e-6*sigma_max` 的奇异值计数均为 0。

因此，第 7 轮观察到的原始 `A/Psi` 方程份额增大，**不能**解释为残差主要占据旧
`(A,Psi)` 近零左空间。右软方向在四状态共有，但残差相对左软子空间的位置不支持该机制；
xtr005 也同样不满足，不是 xtr030 特异阳性。逐模态证据见 `projection_modes_*.csv`，
摘要和图见 `projection_summary.csv`、`projection.png`。

### 4.3 C2：唯一有接受事件的主腿逐步消减 `A/Psi`，但运动量很小；其余三腿零接受步

| 腿 | 接受步 | 拒绝迭代（含末次） | 接受 omega 范围 | `A/Psi` 逐步分类 | 接受步 `sum||F_AP||²` 相对变化 |
|---|---:|---:|---:|---|---:|
| xtr030 B1 | **18** | 21 | 4.8828125e-4 … 3.90625e-3 | **18 REDUCED / 0 AMPLIFIED / 0 FLAT** | **-5.384154e-4** |
| xtr030 B16 | 0 | 11 | — | `NO_ACCEPTED_STEPS` | — |
| xtr005 B1 | 0 | 11 | — | `NO_ACCEPTED_STEPS` | — |
| xtr005 B16 | 0 | 11 | — | `NO_ACCEPTED_STEPS` | — |

xtr030 B1 的 18 个接受事件满足伪时间 `R_inf` 全部下降，但纯定常残差并非一致下降：
`pure_inf` 有 2/18 步上升、完整 `pure_l2` 有 11/18 步上升，`B/delta` L2 有 12/18 步上升。
这正说明 Armijo 的伪时间 merit 不能冒充六方程纯残差判据。

与此同时，`A/Psi` L2 在 **18/18** 步都下降，单步比值范围 0.998653…0.999998；
合计能量从 `2.034692492499e-10` 到 `2.033596982721e-10`，只下降 **0.0538415%**。
`B/delta` 合计能量变化为 **-2.954740e-4**（由 committed accepted_event_summary 重算）。
所以本状态没有“接受步放大 `A/Psi`”的签名；实际是非常强回溯的小步带来极慢消减，
最后仍停在 `1.560038984365e-6`。逐事件见 `accepted_event_summary.csv`，逐节点六方程见
`accepted_steps_native_xtr030_b1_r1.csv`，图见 `accepted_steps.png`。

其余三腿 11 个迭代均全拒绝、零接受事件；`NO_ACCEPTED_STEPS` 不是消减或放大证据。
上一轮 verbose 的 10 条拒绝不含最后触发 `n_fail>10` 的拒绝，本表“拒绝迭代”包含末次，
两者口径不同，不是 replay 不一致。

### 4.4 C3：主假设被明确否定，但非线性根因没有自动闭合

两条 xtr030 均满足右子空间和谱软化条件，却都违反最关键的 `capture20(AP)>=0.50`，
故 C1=`FAIL_AP_PROJECTION`、总判定 **FAIL**。这是物理假设的 FAIL，不是仪器 FAIL。
C2 为 **RECORDED**；它进一步排除了“有接受步时 `A/Psi` 被放大”的描述，但只有一条主腿
存在接受事件，不能对四条腿作统一接受步方向结论。

## 5. 判定

**FAIL（两个固定 xtr030 末态均不满足预注册的 DOMINANT_AP_PROJECTION；仪器与 replay 全过）。**

一句话：**均衡化 Jacobian 的 bottom-20 右软子空间确实由 `A/Psi` 主导，但超过 99.989%
的 `A/Psi` 残差能量落在其左补空间；唯一有接受事件的 xtr030 B1 还在 18/18 步小幅消减
`A/Psi`，故“原始 `A/Psi` 份额增大 ⇒ 残差受近零左空间支配/接受步放大”的解释不成立。**

本轮只否定上述主导机制，不排除其他 K、其他缩放、其他非线性耦合或更广泛范围的作用；
但按 G-DOMAIN 不得为寻找阳性而事后换 K/缩放。S5 多吸引子根因仍为 **UNDEFINED**，
仍没有外层与 IBL 内层都收敛的两条腿，原 2×2 继续不启动。

## 6. 下一步

- 不再追加同类 dense SVD 或调 K 换答案；本轮已给出预注册投影判据的反向结果。
- 若继续 IBL 归因，另写预注册，在**补空间/实际强回溯步**上比较线性预测
  `F+JdU` 与非线性实际 `F(U+dU)` 的逐方程偏差，区分线性模型误差与非线性响应；
  本轮没有测该差值，不能提前命名为根因。
- 若不再投入 IBL，按 phase 7 范围转 S6 或阶段收口。能力表述保留为“xtr030 混合残差结构
  机制未判”，不得写成“已证明 `A/Psi` 近零空间导致地板”。
- 不需要补充 SA/SST 对照。若以后需要，先复用 committed 数据；必须新跑时沿用
  `tools/cfl3d_seq` 或先验证外部构建一致性并记录二进制指纹。D08/D13 保持不动。

## 基线行

研究轮不改生产库和测试面，**不重跑贵重全套或 D08/D13**。自测 **26/26**；
`up3d` 下 F07 **7 passed**。正式四腿累计 **38.0836 s**（4 SVD + 4 重启动），
低于 30 min 时间盒；0 新耦合/外部对照腿。历史开发机套件基线承 **839 / 0 / 35 / 2**，
不是当前迁移机全套实测；未运行 `bench/baseline_*.csv` 再基线，也不声称任何基线数值移动。
