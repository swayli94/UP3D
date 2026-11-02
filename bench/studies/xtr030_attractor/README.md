# xtr030_attractor —— xtr030 松耦合跨 BLAS 落点归因（phase 7 S5）

预注册：[`20261103-0245-xtr030-attractor-prereg.md`](../../../docs/dev_phase_seven/20261103-0245-xtr030-attractor-prereg.md)。

**性质：研究/测量，不改生产库、默认值、D08/D13 判据或 committed gate 证据。**

## 复算

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MPLBACKEND=Agg

python bench/studies/xtr030_attractor/run.py --self-test
python bench/studies/xtr030_attractor/run.py --gate-check
python bench/studies/xtr030_attractor/run.py --baseline-only
python bench/studies/xtr030_attractor/run.py
python bench/studies/xtr030_attractor/run.py --analyze-only
```

- 工况与 D08 相同：NACA0012、M0.50/α2.0/Re3e6、committed medium 网格；
  主腿 xtr=0.30，对照 xtr=0.05。
- B1/B16 两臂的实际线程数都从库读取；B16 只在独立研究 worker 中使用
  `PYFP3D_ALLOW_BLAS_THREADS=1`，生产默认值不变。
- 原生生产容差为 `tol_ds=1e-3`。只有原生跨档分离成立时，才启动 B1 下的
  `tol_ds=1e-8` 种子 × 冻结入流 2×2；同输入 IBL 重放用于定位首差。
- 所有耦合臂调用生产 `run_loose_coupling`。研究 wrapper 只做种子/边界注入和记录，
  `_work/` 保存本地数组与 worker 输出，不提交、不充当证据。
- 时间盒：最多 10 条 fresh 松耦合腿 + 2 条 IBL 重放，正式耦合累计 75 min；
  D08/D13 两遍门检查另限 30 min。已完成腿自动跳过。

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 自测、门检查、worker、原生/对照/紧容差腿、同输入重放、比较与绘图 |
| `environment.json` | 版本、代码/网格指纹、线程与缓存策略、时间盒 |
| `instrument_self_test.csv` | 同一性判别、分离阈值、BLAS 守卫、种子与冻结边界注入牙齿 |
| `gate_checks.csv` | D08/D13 两遍定向证据锁结果；不重刷门证据 |
| `summary.csv` | 每腿出处、收敛、cl、残差、数组 SHA 与成本 |
| `history_*.csv` / `delta_*.csv` | 逐外迭代历史与最终表面 delta* |
| `replay.csv` / `comparisons.csv` | 同输入 IBL 重放、重复性、跨档与紧容差距离 |
| `verdict_metrics.csv` | 由预注册判据直接计算的机制判定及中间量 |
| `run_state.json` | fresh 腿数、累计成本与时间盒 |
| `trajectories.png` / `delta_profiles.png` | 轨迹/停止历史与最终 delta* 对照 |

判定尚待正式测量；完成后在此链接第 6 轮 verdict。
