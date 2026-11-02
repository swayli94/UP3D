# 20261103-0601 — **判定（第 9 轮，IBL 独立研究）：强回溯接受步在左补空间由一阶响应主导；`A/Psi` 全步线性，`delta/B` 混合但非线性不主导**

轮次文件名：`20261103-0601-ibl-step-nonlinearity-verdict.md`
阶段 / 门：phase 7 / 无门（研究）——预注册：
[20261103-0535-ibl-step-nonlinearity-prereg.md](20261103-0535-ibl-step-nonlinearity-prereg.md)
（`24232fc`；投影能量可加性勘误 `141daf2`；`dU` 浮点状态差仪器勘误 `94a4194`）
测量代码：`94a4194`（分支 `kimi/phase7-ibl-ap-projection`；环境指纹见 environment.json）
证据：`bench/studies/ibl_step_nonlinearity/`（run.py + README + 自测/出处/重放/投影/
逐节点与逐事件 CSV + 1 PNG）

---

## 1. 目的

承接第 8 轮“`A/Psi` 残差几乎全在 bottom-20 左补空间、唯一主腿只有强回溯小接受步”的
结果，在这 18 个实际接受步上比较纯定常一阶预测 `F+JdU` 与非线性实际 `F(U+dU)`，
逐源方程分解实际响应 `D`、一阶响应 `P=JdU` 和非线性余项 `N=D-P`。

本轮不改生产库、默认值、闭包、控制器、D08/D13 或 committed gate 证据；不启动 S5 2×2，
不新跑 FP/松耦合、SA/SST/CFL3D 或完整 IBL 重启动，不改第 7/8 轮已有证据。

## 2. 改了什么

- 新增一次性研究 `bench/studies/ibl_step_nonlinearity/`：18 个 committed CFL/omega 的离线
  接受路径重放、纯定常 `P/D/N`、原生点冻结 left-bottom20/complement 投影、逐源方程和
  直接 full 向量分类、方向一致性及交叉项导出。
- 只在原生 B1 状态重建一次与第 8 轮相同的 dense SVD，原因是第 8 轮未保存左奇异向量；
  完整谱和 `capture20` 必须回写旧证据。其余三腿 0 SVD，18 个接受状态 0 追加 SVD。
- 新增 27 项合成仪器自测，登记 `bench/studies/README.md`。
- **生产库零改动、测试零改动、默认值零改动、D08/D13 和 committed gate 证据零改动。**
  正式成功测量为 1 次投影基 SVD + 18 次离线稀疏方向重算，0 次生产 `solve` 重启动、
  0 条 fresh 耦合/外部对照腿、0 组新参数。

### 2.1 两项事前/仪器勘误（均不改变物理判据）

- **投影能量不可直接相加**：正式测量前的已知基自测证明，不同源方程投影后的向量一般不
  正交，原文“BOTTOM20/COMPLEMENT 的源方程能量之和等于 full 能量”错误。已改成投影向量
  求和守卫，并在 `event_summary.csv` 显式提交源能量和、直接 full 能量和交叉项；逐源方程
  能量不当作可加份额。勘误提交 `141daf2` 先于正式测量。
- **浮点 `dU` 定义**：首次正式重放将 `omega*delta` 误作 committed `dU`，而第 8 轮导出的是
  `U_after-U_before`。该次所有 `F/Rpt/dF` 回写已通过、末态逐位一致，仅 `dU` 末位字符串
  不同（绝对差最大 `2.64e-14`），在 SVD/非线性测量前停止。失败证据保留为
  `instrument_failure_du_rounding_*`；改用实际状态差并新增舍入回归自测后从头执行同一测量，
  未改容差、事件或预算，未使用完整重启动或低精度状态替代。

## 3. 怎么测的

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
python bench/studies/ibl_step_nonlinearity/run.py --self-test
python bench/studies/ibl_step_nonlinearity/run.py
python bench/studies/ibl_step_nonlinearity/run.py --analyze-only
python -m pytest tests/F/test_F07_repo_invariants.py -q
```

- 只使用第 8 轮 xtr030 B1 的 18 个接受事件，迭代号、CFL、omega 和输入 SHA 全部绑定 committed
  CSV；其他三腿零接受事件，不补跑、不投票，拒绝 trial 不在本轮定义域。
- 离线重放不调用生产 `solve`，每步以 committed CFL/omega 和生产稀疏 Newton 方向形成相同
  `U_after`；`F0/F1/Rpt0/Rpt1/dU/dF` 全部按旧 `f12` 字符串回写，末态 `U/F` 与第 7 轮逐位一致。
- 每步纯定常 `J` 不含伪时间项，`dU=U_after-U_before`。left bottom-20、row scaling 和
  source-equation-first 投影顺序在原生点冻结；直接 full 向量与逐源方程结果分别报告。

## 4. 结果

### 4.1 I0：自测、重放、谱回写和投影划分全部通过

- 自测 **27/27**，覆盖独立二次残差解析余项、线性/非线性/MIXED/零分母、已知投影、
  非零交叉项、方向正反例、离线重放及其破坏案例、完整导出和判定树。
- 18 步 `F0/F1/Rpt0/Rpt1/dU/dF` 回写 mismatch **全为 0**；最终 `U` 和逐节点残差均与第 7 轮
  **逐位一致**，`replay_roundtrip_mismatches=0`。
- 一次 SVD 的 2460 个 `sigma` 和三个 `capture20` 均与第 8 轮回写一致；SVD 重构相对误差
  **4.2351e-15**，投影划分误差 **2.2204e-16**。
- 响应 FULL=BOTTOM20+COMPLEMENT 最大相对误差 **1.1024e-15**，六源方程投影向量重构
  最大相对误差 **2.6545e-14**，均远低于 `1e-10` 守卫。

### 4.2 C1：补空间的 `A/Psi` 全步一阶主导；`delta/B` 混合，但没有任何非线性主导事件

下表均为固定左 **COMPLEMENT**、逐源方程的 18 步结果。`N/P` 是
`sqrt(sum||N||²/sum||P||²)`，不是逐步中位数；逐步分类次序为线性 / 混合 / 非线性主导。

| 源方程 | pooled 分类 | pooled `N/P` | 逐步线性 / 混合 / 非线性 | 方向不一致步数 |
|---|---|---:|---:|---:|
| `delta` | MIXED | 0.139405 | 13 / 5 / **0** | 0 |
| `A` | **LINEAR_DOMINATED** | **0.0395649** | **18 / 0 / 0** | 0 |
| `B` | MIXED | 0.129661 | 14 / 4 / **0** | **1** |
| `Psi` | **LINEAR_DOMINATED** | **0.0207072** | **18 / 0 / 0** | 0 |
| `Ct1` | LINEAR_DOMINATED | 0.0458598 | 18 / 0 / 0 | 0 |
| `Ct2` | LINEAR_DOMINATED | 0.0453376 | 17 / 1 / **0** | 0 |

`A` 的逐步 `N/P` 最大 **0.0753234**，`Psi` 最大 **0.0264923**，二者 18/18 均低于预注册
0.10 线性界。`delta/B` 的 pooled 比值略过 0.10 而进入 MIXED，但其非线性余项没有一步
同时不弱于 `P` 和 `D`；**全部 324 个事件×源方程×区域条目中，`NONLINEAR_DOMINATED`
计数为 0**。直接 full 向量与 complement 向量也都是 **18/18 LINEAR_DOMINATED**，逐步
`N/P` 范围分别约 `1.385e-3…9.964e-2` 和 `1.385e-3…9.964e-2`。

因此，强回溯并没有在这组接受步上制造“非线性余项压倒一阶响应”的签名；尤其不能用
非线性模型误差解释 `A/Psi` 的补空间响应。`delta/B` 的非线性修正大于 `A/Psi`，但仍只是
混合修正，不得命名为非线性主导根因。

### 4.3 C2：整体预测/实际同向下降；唯一源方程方向翻转不等于非线性主导

- 直接 full 和 complement 向量均为 **18/18 BOTH_REDUCED**；`A/Psi` 源方程在 complement
  也均为 **18/18 BOTH_REDUCED**，与第 8 轮原始 AP L2 逐步下降一致。
- 唯一 `DIRECTION_MISMATCH` 在迭代 **10、CFL=0.25、omega=1/512 的 `B` 源方程**：
  complement 预测下降而实际上升，但该步仍 `LINEAR_DOMINATED`，`N/P=0.0765282`。
  其预测/实际末范数都极贴近 before，说明近平坦处的小余项足以翻转符号，不能据此反称
  非线性主导。bottom20 同一源方程仍 BOTH_REDUCED，翻转位于 complement。
- `Ct2` 在迭代 3 的 complement 为 BOTH_AMPLIFIED、`N/P=0.112444`、MIXED；预测与实际
  同向放大，同样不是模型方向不一致。`delta` 的 12/18、`B` 的 10/18 complement 步为
  BOTH_AMPLIFIED，再次说明纯定常单方程方向不能被伪时间 Armijo 代替。
- bottom20 直接向量有 3/18 BOTH_AMPLIFIED 和 3/18 MIXED，而 complement 两者分别为
  0/18 和 0/18；这说明两个区域必须分开，不能用 bottom20 的局部方向替代全残差。

### 4.4 C3：回答了接受步的一阶/非线性分解，但没有闭合整体地板根因

本轮定义域只有 xtr030 B1 的 18 个**已接受**强回溯步。结果支持：这些步的纯定常响应，
特别是左补空间的 `A/Psi`，由一阶模型主导；`delta/B` 有较大的非线性修正和局部方向差异，
但没有非线性主导事件。

结果**不支持**把“实际强回溯接受步的一阶预测失效”归因为非线性余项主导，也不支持新立
以修正非线性余项为目标的求解器修复。它不排除拒绝 trial 的非线性、伪时间方向/ merit 与
纯定常目标不一致、步长/控制器或其他外层机制；整体 IBL 地板和 S5 多吸引子根因仍为
**UNDEFINED**。

## 5. 判定

**RECORDED（I0/C1/C2 完整；18 个接受步和一次投影基重建全部通过硬守卫）。**

一句话：**xtr030 B1 的 18 个强回溯接受步在 full 和左补空间均由一阶响应主导，`A/Psi`
逐源方程 18/18 线性且同向下降；`delta/B` 虽有混合余项和一次 `B` 方向翻转，但 324 条
事件×方程×区域分类无一非线性主导，故本轮没有“非线性余项压垮 `F+JdU`”的证据。**

## 6. 下一步

- 不追加同义的非线性余项测量，不把 MIXED 或近平坦符号翻转包装成非线性根因，也不据本轮
  修改生产控制器。
- 若继续 IBL 归因，另写预注册研究**为什么准确的一阶接受步仍只产生微小净消减**：固定同一
  轨迹，分解伪时间 Newton 方向、纯定常一阶响应与步长/接受 merit 的预算；拒绝 trial 未在
  本轮保存，若必须研究拒绝机制，须单独预注册新的记录重放，不能外推本轮接受步结论。
- 若不再投入 IBL，按 phase 7 范围转 S6 或阶段收口。S5 多吸引子根因仍 UNDEFINED，2×2 不启动；
  不需要补充 SA/SST 对照，D08/D13 保持不动。

## 基线行

研究轮不改生产库和测试面，**不重跑贵重全套或 D08/D13**。自测 **27/27**；
`up3d` 下 F07 **7 passed**。成功正式测量累计 **8.8334 s**（worker 8.6620 s；1 SVD +
18 离线方向重算），首次 `dU` 仪器失败 2.6008 s、0 SVD，留档但不进入物理判定；0 新耦合/
外部对照腿。历史开发机套件基线承 **839 / 0 / 35 / 2**，不是当前迁移机全套实测；
未运行 `bench/baseline_*.csv` 再基线，也不声称任何基线数值移动。
