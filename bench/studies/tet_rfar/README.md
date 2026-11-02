# tet_rfar —— NACA0012 tet medium r_far 20/30 不收敛归因（phase 7 S4）

预注册：[`20261102-2030-tet-rfar-prereg.md`](../../../docs/dev_phase_seven/20261102-2030-tet-rfar-prereg.md)。
原始观察：`phases/p6/docs/dev_phase_six/20261012-1000-w25-w26-verdict.md` §4.2/§6c。

**性质：研究/测量，不改库、不改生产网格、不改判据、不刷新 committed gate 证据。**

## 复算

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MPLBACKEND=Agg

python bench/studies/tet_rfar/run.py --self-test      # 分类器牙齿
python bench/studies/tet_rfar/run.py --generate-only  # 6 变体 × 2 次生成 + SHA 校验
python bench/studies/tet_rfar/run.py --baseline-only  # coarse/medium × r15/20/30
python bench/studies/tet_rfar/run.py                  # 续跑未收敛 r20/30 medium 的诊断臂
python bench/studies/tet_rfar/run.py --analyze-only   # 从已写 CSV 重出图，不求解
```

- 流场工况 M0.80/α1.25，生产配方与 `tests/D/test_D05_euler_naca0012.py::RECIPE`
  逐字一致；BLAS 由 `pyfp3d` 钉 1，脚本断言实际线程数。
- 只改变 `r_far`，其他网格参数来自生产 `LEVELS`。r=15 重生成文件必须与 Git
  跟踪的 coarse/medium `.msh` SHA-256 相同；每个变体重复生成必须逐位一致。
- `_work/` 是 gitignored 的本地网格工作区，不是证据。脚本、参数、指纹、逐步历史
  和图才是提交证据。
- 时间盒：最多 10 条 fresh 顶层腿、累计求解 60 min；不扫参数、不跑 fine。
  已完成腿由 `summary.csv` 自动跳过，诊断臂只作用于未收敛的 r=20/30 medium。

## 文件

| 文件 | 内容 |
|---|---|
| `run.py` | 生成、求解、分类、自测、逐腿落盘与绘图 |
| `environment.json` | 版本、OS/CPU、线程、缓存策略、配方和时间盒 |
| `classifier_self_test.csv` | 收敛/clamp/周期/非周期/不动点五类合成历史自测 |
| `mesh_manifest.csv` | 6 变体的节点/单元/边界组计数、双生成与 committed SHA 对照 |
| `summary.csv` | 每腿收敛分类、力系数、激波、冻结/线性求解诊断与成本 |
| `history_*.csv` | 逐步残差、F、Gamma、clamp、步长、线性求解与熵修正量 |
| `cp_curves.csv` | 各腿上下翼面截面 Cp |
| `probe_levels.csv` | continuation 各 Mach level 的收敛与成本（若运行） |
| `run_state.json` | 顶层腿数、累计求解墙钟与时间盒状态 |
| `residual_histories.png` | coarse/medium 全腿残差历史 |
| `cp_medium.png` | medium r=15/20/30 的基线与诊断臂 Cp 对照 |

## 判定

第 5 轮判定：[`20261103-0050-tet-rfar-verdict.md`](../../../docs/dev_phase_seven/20261103-0050-tet-rfar-verdict.md) —— **UNDEFINED（根因不可判）**。

- r=20 medium：周期 5 非线性极限环，末 200 点残差/Gamma/线搜索步长同周期重复。
- r=30 medium：本次收敛，旧“撞 400 步帽”观察未复现，记 RECORDED；不归因给机器或吸引子。
- r=20 同网格 continuation 与 entropy-off 分别 clamp 和有界停滞，均未收敛，
  触发预注册 kill criterion；不扫额外参数换答案。

所有数字以 `summary.csv`、`probe_levels.csv`、`history_*.csv` 和判定文件为准。
continuation 的力系数/Cp 按实际 `m_final=0.7875` 后处理；非收敛曲线不作物理解。
