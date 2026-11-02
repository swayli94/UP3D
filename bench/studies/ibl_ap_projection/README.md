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

## 判定

正式测量尚未完成；最终判定、数值和证据指针以随后提交的 verdict 与本目录 CSV/PNG 为准。
