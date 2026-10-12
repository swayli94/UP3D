# hex_cl_attribution —— hex 族 NACA0012 M0.80/α1.25 的 cl −24 % 归因（研究轮 Q2）

预注册：[`docs/dev_phase_six/20261012-1000-research-amg-hex-prereg.md`](../../../docs/dev_phase_six/20261012-1000-research-amg-hex-prereg.md) §Q2。
嫌疑登记：[`20261012-1000-w25-w26-verdict.md`](../../../docs/dev_phase_six/20261012-1000-w25-w26-verdict.md) §4.2
（hex medium 全收敛 cl 0.2582 vs tet committed 0.3413，x_shock 差 0.010 c，几何无辜）。
**性质：研究/测量，不建门、不改判据、不改库。**

## 复算（C3：命令 + 线程环境）

```
conda activate up3d
PYTHONNOUSERSITE=1 PYTHONPATH=<repo> NUMBA_NUM_THREADS=16 \
    python phases/p6/bench/studies/hex_cl_attribution/run.py            # 4 腿求解 + 分析 + 图（~2.5 min）
python phases/p6/bench/studies/hex_cl_attribution/run.py --analysis-only   # 不求解，从 CSV 重出分析
```

- 工况 M0.80/α1.25；生产配方与 D05 逐字一致（`upwind_c=1.5, m_crit=0.95,
  freeze_tol=1e-6, freeze_refresh_max=8, precond=direct, direct_refactor_every=4,
  n_newton_max=400, n_picard_seed=5`，见 `tests/D/test_D05_euler_naca0012.py`）。
- BLAS 由 pyfp3d 自身钉 1（脚本断言 `blas_thread_count()==1`；`ladder.csv`
  的 `blas_threads` 列全为 1）。
- tet 的 xcoarse/coarse/medium 三级值**抄自 committed**
  `cases/gates/D05_euler_naca0012/summary.csv`（证据纪律：不重算已提交 artifact）；
  tet medium 重算一次只为 Cp 曲线（committed 证据没有 Cp 数据列），其
  cl = 0.3413401232 与 committed **逐位一致**。
- Cp 分段边界事前写死（看数之前）：LE_peak [0,0.05]、front [0.05,0.45]、
  shock [0.45,0.75]（两族 x_shock 均在 ~0.60）、TE_unload [0.75,1.0]。

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 求解 + 分析 + 出图（单脚本，腿结果逐腿落盘） |
| `ladder.csv` | 4 腿全状态：cl/cd/x_shock/converged/residual/n_newton/线程出处 |
| `cp_curves.csv` | 各腿截面 Cp（family, level, side, x, cp） |
| `cp_segment_cl.csv` / `cp_segment_delta.csv` | 分段 cl 贡献 ∫(cp_l−cp_u)dx 与 hex−tet 逐段差 |
| `bubble_metrics.csv` | 上翼面 [0.2,0.55] max Cp（气泡计；正 = 中部减速气泡） |
| `grid_convergence_hex.csv` | 两族三级形态（W2.1a 仪器：参照侧 implied_order + 三分词汇） |
| `residual_history_<family>_<level>.csv` | 逐 Newton 步残差（kill criterion 记录） |
| `cp_tet_vs_hex_medium.png` | (a) 两族 medium Cp 对照（分段着色，CFL3D L4 参照） |
| `cp_hex_ladder_vs_tet.png` | (c) hex 三级 Cp 形态 vs tet medium |
| `cl_ladder.png` | (b)+(c) cl vs h ∝ N^(−1/2)，两族 + CFL3D L4 线 |

## 关键数字（2026-10-12 实测，BLAS=1，16 线程）

cl 三级序列（hex fresh / tet committed；CFL3D Euler L4 = 0.348881）：

| level | hex cl | tet cl | hex x_shock | tet x_shock |
|---|---|---|---|---|
| 1 | 0.296546 (coarse) | 0.296972 (xcoarse) | 0.5102 | 0.5329 |
| 2 | 0.258223 (medium) | 0.407108 (coarse) | 0.6094 | 0.6142 |
| 3 | 0.256059 (fine) | 0.341340 (medium) | 0.6295 | 0.5991 |

- (b) 差**不收缩**：medium −0.0831 → fine −0.0852（相对 tet medium）；hex 梯子
  在 ~0.256 走平（|Δ| 0.0383 → 0.0022），tet 梯子围绕 CFL3D 振荡。
  hex fine 的 h（0.00670）已小于 tet medium（0.00692），等 h 下仍差 −25 %。
- (a) 逐段 Δcl（hex medium − tet medium，总 −0.0827）：LE_peak −0.0003、
  **front −0.0341、shock −0.0509**、TE_unload +0.0026 —— 亏损 100 % 在
  前中段 + 激波区；LE 峰与尾缘无辜。
- 指纹（`bubble_metrics.csv`）：上翼面中部减速气泡（cp > 0）**随加细涌现并前移**
  —— hex coarse 无（max cp −0.297 @0.541），medium +0.199 @0.457，
  fine +0.383 @0.376；tet medium 无（−0.891）。下翼面两族皆干净。
- **kill criterion 触发**：hex fine `converged=False`（400 步帽，|R| 滞停在
  2.2–2.5e-7 振荡 ±5 % ≥200 步，非单调降、非发散 —— 激波结构极限环，
  残差历史见 `residual_history_hex_fine.csv`）。按预注册照记，未换网格。
  cl 对此停滞水平的稳健性有 D05 线程实验标定（|R|~1e-6 封顶态 vs 收敛态
  cl 差 0.15–0.36 %），「不收缩」结论不受影响。

## 结论

**② 系统性缺陷**（预注册判别式：加细不收缩 ⇒ ②），指纹 = hex 族收敛到
一个 tet 族从未出现的上翼面中部气泡/双激波结构，cl 亏损全部位于该结构与
被它移位的激波所占的两段。①（分辨率分布）被「等 h 仍差 25 %」排除；
③（tet 参考有偏）被 tet medium −2.16 % vs CFL3D L4 且 tet 梯子任何一档
都不近 0.26 排除。
**残余未分辨（诚实记录）**：② 内部 —— 缺陷在 hex 生成器（拓扑/分布/拼缝）
还是在该拓扑上被触发的求解器熵冻结（hex 各腿 `frozen_in_transient` 的
|dsigma| 0.21/0.19 vs tet 0.067，tail n_shock 搅动 127–410 vs 140–164），
(a)(b)(c) 三件仪器分不开；判别需动配方的探针轮，超出本轮范围。
