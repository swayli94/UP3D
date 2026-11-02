# ibl_ap_projection —— `A/Psi` 残差 Jacobian 投影与接受步分解（phase 7）

预注册：[`20261103-0454-ibl-ap-projection-prereg.md`](../../../docs/dev_phase_seven/20261103-0454-ibl-ap-projection-prereg.md)。

**性质：一次性研究/测量；不改生产库、默认值、D08/D13 判据或 committed gate 证据。**

## 问题与范围

在第 7 轮相同的 xtr030 B1/B16 主腿和 xtr005 B1/B16 对照末态上：

1. 对固定边数据的生产 `J_BL,BL` 做一次 row-then-column 均衡化和完整 dense SVD，
   测量 bottom-20 右奇异子空间的 `A/Psi` 能量及对应左子空间对残差的投影；
2. 包装生产 `IBL3Solver.solve` 的原配置重启动，记录每个接受步前后六方程纯定常残差、
   伪时间残差和状态增量，并与第 7 轮末态及轨迹复核。

研究矩阵最多四次 SVD 和四次 IBL 重启动，**零条新 FP/松耦合腿、零条新 SA/SST/CFL3D
对照、零组新参数**。不启动 S5 的 2×2，不修改 K、缩放或控制器来寻找更好看的结论。

## 复算

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg

python bench/studies/ibl_ap_projection/run.py --self-test
python bench/studies/ibl_ap_projection/run.py
python bench/studies/ibl_ap_projection/run.py --analyze-only
```

- 输入是第 7 轮 `ibl_floor_restart/` 和 S5 `xtr030_attractor/` 已存在的 `_work/` 数组。
  `_work/` 不跟踪，只作仪器输入；当前输入 SHA、数组 SHA、网格 SHA 和原生残差必须同时
  符合第 7 轮与 S5 的 committed CSV。
- 接受步重启动的最终 `U` 和逐节点残差必须与第 7 轮本地缓存逐位一致；committed 逐迭代
  历史按预注册的输出舍入容差复核。缺输入或守卫失败会停止，不自动重算耦合腿。
- 所有诊断 worker 实际 BLAS=1、Numba=16；B1/B16 只表示输入历史出处。每腿顺序执行，
  单腿上限 10 min，累计上限 30 min，不清缓存。

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 出处守卫、自测、SVD 投影、生产接受步记录、判定树与证据导出 |
| `environment.json` / `run_state.json` | 代码/包/网格指纹、线程、预算、腿数和成本 |
| `instrument_self_test.csv` | 合成 SVD/投影、生产 solve 包装正反向案例与判定树牙齿 |
| `provenance.csv` / `summary.csv` | 每腿出处、native/replay 守卫及全部关键读数 |
| `projection_modes_*.csv` | 每腿完整谱、右六变量能量及 full/AP/BD 左投影系数 |
| `projection_summary.csv` / `projection.png` | C1 的机器可读摘要和谱/投影图 |
| `accepted_steps_*.csv` | 每个接受步、每节点、每方程的前后残差和状态增量 |
| `accepted_event_metrics_*.csv` / `iteration_trace_*.csv` | 接受事件范数与全部迭代/trial 记录 |
| `accepted_step_summary.csv` / `accepted_steps.png` | C2 的逐腿摘要和逐步残差比 |
| `criteria.csv` / `verdict_metrics.csv` | 预注册 C1–C3 与总判定 |

本研究不被测试加载；自测在正式复算前显式运行，收尾另跑 F07 仓库不变量。
phase 7 收口时随一次性研究归档。

## 启动与仪器勘误记录

- `startup_failure_base_*`：第一次误从 anaconda `base` 启动，worker 的 BLAS 出处为
  `None`，在任何 SVD/重启动前停止；改用 `conda activate up3d` 后实际 BLAS=1。
- `instrument_failure_name_error_*`：第二次首腿 SVD 和求解完成后，接受步导出函数遗漏
  `numpy` 导入而停止，未通过 replay、未纳入任何物理判定。已补 `trace_export_normal`
  自测，正式自测增至 **26/26**，随后从头重跑四腿。

两次失败均不作为物理证据；只有正式 `summary.csv` 中全部出处/native/SVD/replay 守卫通过的
腿才进入判定。

## 判定

第 8 轮判定：[`20261103-0516-ibl-ap-projection-verdict.md`](../../../docs/dev_phase_seven/20261103-0516-ibl-ap-projection-verdict.md) —— **FAIL（主导 `A/Psi` 投影假设；仪器全部通过）**。

- 四腿 bottom-20 右子空间的 `A/Psi` 份额均为 0.6006–0.6163，谱软化条件也全过；
  但 `capture20(AP)` 仅为 **3.982458e-5 … 1.783390e-4**，全部远低于 0.50，
  也低于 20/2460=0.00813008 的维数基线。两条主腿超过 99.989% 的 `A/Psi` 残差能量
  位于其余 2440 维左补空间。
- xtr030 B1 的 18 个接受步全部小幅消减 `A/Psi`，合计能量变化 **-5.384154e-4**；
  其余三腿 `NO_ACCEPTED_STEPS`。接受步伪时间 merit 全降不等于纯六方程残差全降。
- 四腿最终 `U` 和逐节点残差与第 7 轮逐位一致；SVD/native/replay/出处守卫全过。
- S5 根因仍 UNDEFINED，2×2 不启动；不追加换 K/缩放的谱计算，不新跑外部对照。

所有结论以 `criteria.csv`、`verdict_metrics.csv`、`provenance.csv`、`summary.csv`、
逐模态/逐接受步 CSV 和 PNG 为准；生产库、测试判据和 committed gate 证据未改。
