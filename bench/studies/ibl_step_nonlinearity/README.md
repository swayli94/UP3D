# ibl_step_nonlinearity —— 强回溯接受步非线性余项与左补空间分解（phase 7）

预注册：[`20261103-0535-ibl-step-nonlinearity-prereg.md`](../../../docs/dev_phase_seven/20261103-0535-ibl-step-nonlinearity-prereg.md)
（含正式测量前的投影能量可加性勘误）。

**性质：一次性研究/测量；不改生产库、默认值、D08/D13 判据、committed gate 证据或第 7/8 轮产物。**

## 问题与范围

只在第 8 轮 `native_xtr030_b1_r1` 的 18 个已接受强回溯事件上比较：

```text
P = J dU
D = F(U+dU) - F(U)
N = D - P = F(U+dU) - (F + J dU)
```

逐步、逐源方程导出 FULL/BOTTOM20/COMPLEMENT 三个区域的 `P/D/N` 和预测/实际末残差。
左投影基在原生点冻结，补空间定义为第 8 轮 bottom-20 左奇异子空间的其余 2440 维。
其他三腿没有接受事件，不补跑、不投票；不研究拒绝 trial。

不同源方程投影后的向量一般不正交，因此 BOTTOM20/COMPLEMENT 的逐源方程能量**不是可加份额**。
仪器守卫检查投影向量之和，`event_summary.csv` 同时给出直接 full 能量、源方程能量和及交叉项；
两种聚合分类不混写。

## 复算

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg

python bench/studies/ibl_step_nonlinearity/run.py --self-test
python bench/studies/ibl_step_nonlinearity/run.py
python bench/studies/ibl_step_nonlinearity/run.py --analyze-only
python -m pytest tests/F/test_F07_repo_invariants.py -q
```

- 离线重放只使用 committed CFL/omega 重算 18 个接受方向，不调用生产 `IBL3Solver.solve`，
  不运行控制器或线搜索。所有相关单元格按第 8 轮 `f12` 格式回写一致，末态与第 7 轮缓存逐位一致，
  否则停止，不退回低精度 CSV 状态或完整重启动。
- 第 8 轮没有保存左奇异向量，故只在原生 B1 状态重建一次同一 dense SVD 投影基；完整谱和
  `capture20` 必须回写旧证据。不做其余三腿 SVD，不在接受状态追加 SVD。
- 正式预算：1 次原生点 dense SVD + 18 次离线稀疏方向重算，0 次生产 solve 重启动、
  0 条新 FP/松耦合/SA/SST/CFL3D 腿、0 组新参数。

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 自测、出处、离线接受路径、固定投影、非线性分解、判定与导出 |
| `instrument_self_test.csv` | 二次残差解析 oracle、投影/交叉项、重放/方向/判定树及完整导出牙齿 |
| `environment.json` / `run_state.json` | 代码/包/网格指纹、线程、预算与成本 |
| `provenance.csv` / `summary.csv` | 原生输入和全部上游 CSV/NPZ 指纹及总守卫 |
| `replay_guard.csv` / `projection_guard.csv` / `energy_guard.csv` | 接受路径、谱回写、投影向量重构守卫 |
| `accepted_nonlinearity_*.csv` | 每接受步、节点、方程的 `dU/F/JdU/Fpred/F1/D/N` |
| `event_equation_metrics.csv` | 18 步 × 6 源方程 × 3 区域的幅值、余项和方向分类 |
| `event_summary.csv` | 直接 full 向量分类及 `P/D/N` 源能量和、交叉项 |
| `equation_summary.csv` | 各方程跨 18 步能量、分类计数、比值分布和单步集中度 |
| `criteria.csv` / `verdict_metrics.csv` / `nonlinearity.png` | 机器可读判定与证据图 |

本研究不被测试加载；自测在正式复算前显式运行，phase 7 收口时随一次性研究归档。

## 仪器勘误记录

- `instrument_failure_du_rounding_*`：首次正式离线重放把 `omega*delta` 误作 committed `dU`，
  而第 8 轮导出的是浮点 `U_after-U_before`。所有 `F/Rpt/dF` 回写已通过、末态也逐位一致，
  仅实际状态差的末位字符串守卫触发；在 SVD 和非线性测量前停止。失败证据保留，修正后自测
  增至 **27/27**，不以完整重启动或低精度状态替代。
