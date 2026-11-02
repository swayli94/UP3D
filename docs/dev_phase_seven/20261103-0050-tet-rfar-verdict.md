# 20261103-0050 — **判定（第 5 轮，研究轮 S4）：r=20 medium 周期 5 极限环实锤；两根因对照均未收敛，机制 UNDEFINED**

轮次文件名：`20261103-0050-tet-rfar-verdict.md`
阶段 / 门：phase 7 / 无门（研究）——预注册：
[20261102-2030-tet-rfar-prereg.md](20261102-2030-tet-rfar-prereg.md)
提交：分支 `kimi/phase7-charter`（本文件随判定提交）
证据：`bench/studies/tet_rfar/`（run.py + README + 指纹/摘要/逐腿历史 CSV + 2 PNG）

---

## 1. 目的

回答 NACA0012 tet medium 在 M0.80/α1.25、D05 生产配方下，r_far=20/30
重生成变体撞 400 步帽的机制：远场物理、网格实现、初值路径，还是熵修正参与。
预注册明确规定：同网格 continuation 与 entropy-off 均不改变不收敛状态时，
不得凭 r 扫描混合作出根因结论，判 **UNDEFINED**。

## 2. 改了什么

- 新增 `bench/studies/tet_rfar/`：可复算生成/求解/分类脚本、环境出处、网格 SHA、
  6 条基线 + r=20 两条同网格诊断臂的完整历史、Cp 与成本。
- `bench/studies/README.md` 登记一次性研究，phase 7 收口随阶段归档。
- **库零改动、生产网格零改动、测试/判据/默认值零改动、committed gate 证据零刷新。**

## 3. 怎么测的

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MPLBACKEND=Agg
python bench/studies/tet_rfar/run.py --self-test
python bench/studies/tet_rfar/run.py --generate-only
python bench/studies/tet_rfar/run.py --baseline-only
python bench/studies/tet_rfar/run.py
python bench/studies/tet_rfar/run.py --analyze-only
python -m pytest tests/F/test_F07_repo_invariants.py
```

- 工况与配方：M0.80/α1.25；D05 `RECIPE` 原样（direct、400 步帽、Picard seed 5）。
- G-ATTRACTOR 出处：16 Numba 线程、BLAS 实测 1、同一暖缓存策略；版本与 OS 入
  `environment.json`。每条腿的网格 SHA 同时入 `summary.csv`。
- 6 变体各生成两次，SHA 全部逐位一致；r=15 coarse/medium 与 Git 跟踪生产网格
  SHA 完全相同。Gmsh 4.15.2。完整哈希见 `mesh_manifest.csv`。
- 参考口径：r=15 是 committed 网格/历史防漂移锚，不是外部真值；r=20/30 是新鲜
  变体。非收敛腿的 cl/Cp 只描述末态，不当作可用物理解。

## 4. 结果

### 4.1 C1：网格出处通过；r=20 原失败复现，r=30 未复现

| level | r_far | n_nodes / n_tets | cl | x_shock | 最终 \|R\| / n | 分类 | wall_s |
|---|---:|---:|---:|---:|---:|---|---:|
| coarse | 15 | 5610 / 16386 | 0.407108 | 0.614243 | 8.67e-13 / 20 | converged | 3.69 |
| coarse | 20 | 6650 / 19476 | 0.416436 | 0.630041 | 6.01e-13 / 18 | converged | 4.57 |
| coarse | 30 | 9002 / 26466 | 0.433243 | 0.638450 | 9.91e-13 / 24 | converged | 6.46 |
| medium | 15 | 20888 / 61788 | 0.341340 | 0.599095 | 1.97e-13 / 46 | converged | 16.78 |
| medium | 20 | 25062 / 74250 | 0.502131* | 0.709232* | 9.05e-08 / 400 | **periodic_limit_cycle(p=5)** | 60.48 |
| medium | 30 | 33954 / 100800 | 0.766835 | 0.833535 | 6.44e-13 / 89 | converged | 38.03 |

`*` 非收敛末态，不是可用物理解。数据源：`summary.csv`；历史：`history_*.csv`。

- r=15 medium 新算 cl=0.3413404093583，与 committed D05 的 0.3413401232 只差
  2.86e-7（绝对）；本机 D05 锁的已知 Cp 微差不因本轮改变，也未刷新。
- **r=20**：末 200 点残差在 9.045e-8…2.826e-7 内以 5 步严格重复，折叠
  RMSE/极差=1.03e-9；Gamma 与线搜索步长也以同一周期重复。零 clamp、零 GMRES
  stall，三次冻结回退后末态未冻结。这是真实的非线性迭代极限环测量，
  不是“撞帽”别名，也不是线性求解失败的误命名。
- **r=30**：本次收敛到不同 Cp/激波形态，原“r=30 也撞帽”未复现，记 **RECORDED**。
  原始变体没有网格 SHA/逐腿历史，不能把差异归给机器、Gmsh 或吸引子；
  本次 r=30 收敛也不证明它与 r=15 是同一物理解。

### 4.2 C2：分类器有牙齿

`classifier_self_test.csv` 五类合成历史全过：converged、clamped、周期 4 极限环、
非周期有界停滞、不动点停滞。周期并列舍入在正式求解前修正为取最小基频；
period4 实测返回 p=4，而非 p=8。真实腿的原始逐步量全部保留，未用 cl 倒推分类。

### 4.3 C3：r=20 两个同网格诊断臂都没有收敛

| 臂 | M_final | 分类 | 最终 \|R\| | clamp / GMRES stall | cl（只作末态记录） | wall_s |
|---|---:|---|---:|---|---:|---:|
| 冷启动 entropy-on 基线 | 0.8000 | periodic_limit_cycle(p=5) | 9.05e-08 | 0 / 0 | 0.502131 | 60.48 |
| 同网格 Mach continuation | **0.7875** | **clamped** | 1.10e-02 | 5276+1162 / 0 | 0.844472 | 69.66 |
| 同网格冷启动 entropy-off | 0.8000 | **bounded_stall** | 3.80e-08 | 0 / 0 | 1.126030 | 64.11 |

continuation 的 level 历史（`probe_levels.csv`）：

| M | 结果 | n_newton | 最终 \|R\| |
|---:|---|---:|---:|
| 0.700 | converged | 6 | 4.35e-11 |
| 0.750 | converged | 26 | 3.68e-13 |
| 0.800 | failed | 17 | 9.48e-03 |
| 0.775 | converged | 17 | 1.89e-11 |
| 0.800 | failed | 13 | 4.43e-02 |
| 0.7875 | failed, clamped | 48 | 1.10e-02 |

⇒ 没有“continuation 收敛”证据，也没有“entropy-off 收敛”证据。
entropy-off 把周期环变成非周期有界停滞（周期分数 0.347 > 0.25 阈值），
只说明失败形态改变，**不满足预注册的“熵修正参与使收敛状态改变”判据**。

### 4.4 C4/C5 与时间盒

- C4 的两个同网格臂均未改变不收敛状态，触发预注册 **UNDEFINED** 分支。
  可以确认 medium-specific 数值失败及 r=20 的周期 5 形状，但不能闭合到远场物理、
  网格拓扑、初值或熵修正中的某一个根因，更不能声称多吸引子已实锤。
- C5 的摘要、网格、逐腿、level、Cp、成本与 PNG 均已落盘。8 条有效顶层腿求解
  263.70 s；一次仪器勘误复算 69.66 s；合计 333.36 s。上限 10 腿/3600 s 未触发，
  没有新增 r 值、seed/配方扫描或 fine 计算。

## 5. 仪器勘误

1. **预注册实施前（尚未流场测量）**：误假设 coarse/medium 有 manifest；实测该族
   coarse/medium 是 Git 跟踪文件，仅 fine 有 sidecar。改为直接比较 committed 文件
   SHA-256，强度不变；勘误已先于测量提交（`1c8ad5b`）。
2. **整理阶段**：continuation 第一版后处理误用目标 M=0.80，实际末态是 M=0.7875。
   收敛分类、残差、level 历史和成本不受影响；修正为统一按 `m_final` 计算力系数/Cp/
   激波后，仅重跑该臂。复算的残差与收敛路径逐位复现，表中 cl=0.844472 是修正后的
   M=0.7875 非收敛末态。错误的 M=0.80 后处理行未纳入提交证据。
3. Cp 图例显式标注分类与 M_final，避免把 clamp/stall 曲线读成收敛物理解。

## 6. 判定

**UNDEFINED（根因不可判；预注册 kill criterion 照章执行）。**

一句话：**r=20 medium 的 400 步失败已实锤为周期 5 非线性极限环，但同网格
continuation 与 entropy-off 分别 clamp 和有界停滞，均未收敛；r=30 旧失败本次未复现，
故本轮只交付可复算失败形态，不把远场半径、网格实现、初值或熵修正中的任何一个
写成已闭合根因。**

这不是“再扫一个参数就会绿”的授权；继续 S4 需要新证据、新预注册与使用者裁决。

## 7. 下一步

- 按阶段序列进入 S5（xtr030 多吸引子）前，把 r=20 的周期 5 历史作为候选仪器案例；
  不得把其非收敛 cl/Cp 当物理锚。
- r=30 的原观察不可复现只保留 RECORDED；若旧机器/旧变体原始 SHA 与历史以后补齐，
  可在固定出处下重开 A/B。
- 本轮暴露但不顺手修的问题：迁移机 11 项既有测试红、D05 证据锁的环境微差，
  均不属于 S4 根因结论。

## 基线行

研究轮不改库和测试面，**不重跑贵重全套**。分类自测 5/5、`py_compile` 通过，
F07 **7 passed / 0 failed @2.17 s**。历史开发机套件基线仍承 **839 / 0 / 35 / 2**；
这不是当前迁移机全套实测。当前迁移机已知的 11 项定向失败保持原样，未刷新证据、
未改阈值、未宣称全套转绿。
