# ibl_floor_restart —— 固定末态 IBL 地板与重启动归因（phase 7）

预注册：[`20261103-0420-ibl-floor-prereg.md`](../../../docs/dev_phase_seven/20261103-0420-ibl-floor-prereg.md)。

**性质：一次性研究/测量；不改生产库、默认值、D08/D13 判据或 committed gate 证据。**

## 问题与范围

检查 S5 的 xtr030 B1/B16 末态是否复现已知 TE 带 `(B,δ)` IBL 公式层地板，
并在同一末态、同一边数据下用生产 `IBL3Solver.solve` 原样重启动。xtr005 B1/B16
是当前同机对照。研究矩阵只有四次残差重建和四次 IBL 重启动，**零条新 FP/松耦合腿、
零次 dense SVD**，不扫参数、不启动 S5 的 2×2。

## 复算

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg

python bench/studies/ibl_floor_restart/run.py --self-test
python bench/studies/ibl_floor_restart/run.py
python bench/studies/ibl_floor_restart/run.py --analyze-only
```

- 输入是 `bench/studies/xtr030_attractor/_work/` 中四个已存在的 S5 末态及其首 IBL
  冻结边界。`_work/` 不跟踪，只作仪器输入；脚本逐腿核对 committed S5 summary 的
  `U/phi/delta_star` 数组 SHA、网格 SHA 和末行原生 IBL 残差。
- 输入缓存缺失、指纹不符或原生残差重建失败会响亮停止；本研究不会自动重算 S5 耦合腿。
- 所有诊断 worker 都固定实际 BLAS=1、Numba=16；B1/B16 只表示输入状态的历史出处。
- 正式时间盒为四次重启动、累计 20 min。完成腿自动跳过，`_work/` 保存 worker 输出和
  本地数组，提交 CSV/PNG 才是结论证据。

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 出处守卫、自测、生产 IBL 残差重建与重启动、判定树、证据导出 |
| `environment.json` / `run_state.json` | 代码/包/网格指纹、线程与缓存策略、腿数和成本 |
| `instrument_self_test.csv` | 指纹、残差签名、verbose 解析、重启动分类和完整判定树牙齿 |
| `provenance.csv` / `summary.csv` | 每腿输入出处、原生残差复核、收敛和全部关键读数 |
| `residual_anatomy_*.csv` | 全节点、全六方程残差与闭包地板活动集 |
| `restart_history_*.csv` | 生产 verbose 的逐迭代残差、CFL 和步拒绝记录 |
| `residual_summary.csv` / `restart_summary.csv` | C1/C2 的机器可读摘要 |
| `criteria.csv` / `verdict_metrics.csv` | 预注册 C1–C3 判定与总判定 |
| `residual_map.png` / `restart_trajectories.png` | 残差空间分布和重启动轨迹 |

本研究不被测试加载；自测在每次正式复算前显式运行，收尾另跑 F07 仓库不变量。

## 判定

第 7 轮判定：[`20261103-0442-ibl-floor-verdict.md`](../../../docs/dev_phase_seven/20261103-0442-ibl-floor-verdict.md) —— **UNDEFINED（xtr030 是否属于同一已知公式层地板不可判）**。

- 四腿输入指纹和原生残差重建全部通过；四次原样重启动均不收敛，CFL 均到 `1e-3`，
  末段纯残差平台化，分类均为 `CONTROLLER_FLOOR`。
- xtr005 两腿复现旧 TE `(B,δ)` 签名；xtr030 两腿的 `B+delta` 合并份额为
  0.882496/0.861668，低于预注册 0.90，且 `A/Psi` 份额更大，不能事后降阈值合并机制。
- S5 仍没有外层与 IBL 内层都收敛的两条腿，2×2 继续不启动。

所有结论以 `criteria.csv`、`verdict_metrics.csv`、`provenance.csv`、`summary.csv`、
逐节点/逐迭代 CSV 和 PNG 为准；生产库、测试判据和 committed gate 证据未改。
