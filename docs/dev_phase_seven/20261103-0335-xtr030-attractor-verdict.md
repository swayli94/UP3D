# 20261103-0335 — **判定（第 6 轮，研究轮 S5）：FP 首差已闭合；同输入首 IBL 无 BLAS 响应，多吸引子根因 UNDEFINED**

轮次文件名：`20261103-0335-xtr030-attractor-verdict.md`
阶段 / 门：phase 7 / 无门（研究）——预注册：
[20261103-0245-xtr030-attractor-prereg.md](20261103-0245-xtr030-attractor-prereg.md)
提交：分支 `kimi/phase7-xtr030`（本文件随判定提交）
证据：`bench/studies/xtr030_attractor/`（run.py + README + 自测/门检查/摘要/历史/比较 CSV + 2 PNG）

---

## 1. 目的

回答 D08 xtr030 medium 跨 BLAS 配置 cl 与外迭代数分家的机制，并检验历史“多吸引子”
是否真的是同一离散方程的两个紧收敛不动点。预注册明确要求先排除 FP 初始扰动、
IBL 直接敏感、提前停止与首个外迭代冻结入流的路径依赖；不得把单次 cl 差直接称为
多吸引子。历史 B24 数值不是本轮追跑目标，线程上限内只做 B1/B16。

## 2. 改了什么

- 新增 `bench/studies/xtr030_attractor/`：独立 worker、同档重复、xtr005 对照、
  同输入生产路径重放、种子/冻结边界包装、判定树与证据导出。
- 新增一次性研究登记 `bench/studies/README.md`；预注册含两项实施前/失败后仪器勘误。
- **生产库零改动、测试零改动、默认值零改动、D08/D13 判据零改动、committed gate
  证据零刷新。**紧容差 2×2 因原生不收敛而未启动。

## 3. 怎么测的

```bash
conda activate up3d
export PYTHONNOUSERSITE=1 NUMBA_NUM_THREADS=16
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MPLBACKEND=Agg
python bench/studies/xtr030_attractor/run.py --self-test
python bench/studies/xtr030_attractor/run.py --gate-check
python bench/studies/xtr030_attractor/run.py
python bench/studies/xtr030_attractor/run.py --analyze-only
```

- 工况与 D08 相同：NACA0012、M0.50/α2.0/Re3e6、committed medium、xtr030；
  xtr005 是对照，不预设零敏感。B16 只在独立研究 worker 中开启既有 BLAS 逃生舱。
- 每腿从库读取并断言实际 BLAS；Numba=16；同档重复为独立 fresh Python 进程。
  网格 SHA-256 全腿相同，见 `summary.csv`。
- C1 为 B1/B16 各两次原生腿，C2 为 xtr005 两档，C3 为生产首外迭代内的同输入
  B1/B16 重放。C4 仅在原生两腿均有资格判不动点时启动。
- 初次 standalone 重放触发仪器守卫后，保留 `replay_standalone_invalid.csv`，
  追加 1 条 `instrument_one_outer_b1` 诊断腿，再按预注册勘误改为生产路径重放。
  其余诊断比较见 `instrument_diagnostics.csv`；这些诊断不进入机制腿结论。

## 4. 结果

### 4.1 仪器与当前机门边界

- 自测 **15/15**：同一/不同/灰区/不收敛判别、BLAS 守卫、种子与全部 IBL 输入覆盖、
  四条完整判定树及原生不收敛护栏均通过，见 `instrument_self_test.csv`。
- D08/D13 定向证据锁两遍均各 **2 failed**，失败节点完全相同、墙钟 397.29/396.97 s，
  见 `gate_checks.csv`。它们登记当前机散布边界，不作为本轮吸引子机制证据；
  未改阈值或重刷证据。F07 7/7、F09 4/4 通过。

### 4.2 C1：同档逐位复现，跨档分离成立，但 B1 原生未收敛

| 腿 | 重复 | 外层 converged | n_outer | cl | 末 ds_change_rel |
|---|---:|---:|---:|---:|---:|
| B1 xtr030 | r1 / r2 完全相同 | **0** | 10 | 0.2653245956164 | 6.93725e-3 |
| B16 xtr030 | r1 / r2 完全相同 | 1 | 8 | 0.2543591890987 | 3.92781e-4 |

- 同档 `E_cl=E_ds=0`，跨档 `Delta cl=0.01096540652`、`max|Delta delta*|=8.301585e-4`，
  分离判据成立。k=0 的 `phi/gamma` 已不同。
- **但 B1 没有外层收敛，历史“B1 4 步收敛到 0.282068”未在本机复现。**
  B16 也不是历史 B24 的逐位复刻，预注册不要求追跑该值。当前两臂不能组成“两个
  收敛吸引子”的证据对，C4 因此不启动。

### 4.3 C2：xtr005 对照敏感较小，但不是零敏感

| 腿 | 外层 converged | n_outer | cl | 末 ds_change_rel |
|---|---:|---:|---:|---:|
| B1 xtr005 | 1 | 3 | 0.2826417959736 | 6.24932e-4 |
| B16 xtr005 | 1 | 4 | 0.2821403971504 | 7.90701e-4 |

跨档 `Delta cl=5.013988232e-4`、`max|Delta delta*|=1.348920e-4`。两档 k=0 FP
同样已分家；对照只说明 xtr030 的放大更强，不证明 xtr005 是唯一或紧不动点。

### 4.4 C3：首差在 FP；同输入首 IBL 对 BLAS 无直接响应

- k=0 B1/B16：`max|Delta phi|=4.681161084e-9`，
  `max|Delta gamma|=6.892261761e-11`，见 `comparisons.csv`。
- 修正后的生产路径重放中，B1 和 B16 都逐位复现 B1 原生首 IBL：
  `native_max_abs=0`、两腿 U 的 SHA-256 完全相同，见 `replay.csv`。
  ⇒ **这份首 IBL 输入下，IBL 求解没有直接 BLAS 响应；初始扰动先于 IBL，位于 FP。**
  该结论仅限首外迭代，不外推为所有后续 IBL 调用都无关。
- 重放本身 **100 次 IBL 内迭代未收敛**，残差 2.109157992e-6；逐腿历史也显示
  原生/对照各外迭代的 `ibl_converged` 均为 0。故 B16 的 converged=1 仅表示外层
  `tol_ds=1e-3` 判据通过，**不是 FP+IBL 全系统紧收敛**。

### 4.5 C4：未启动；多吸引子根因不可判

预注册要求排除提前停止和冻结边界差异后才能判多个紧不动点。当前 B1 原生未收敛，
B16 也只有外层收敛且 IBL 内求解未达自身判据 ⇒ 没有两个合格的紧收敛解可比较。
因此未运行种子 × 冻结入流 2×2，**提前停止、冻结入流路径依赖和真正多不动点均未判，
不能把未运行写成其中任何一项的负结果。**

### 4.6 仪器勘误与成本

1. 初版 standalone 重放与原生首 IBL 最大差 29.3718992945，守卫停止且证据保留。
   生产 1 外迭代诊断逐位复现原生首解；修正为生产路径内覆盖全部首 IBL 输入后，
   B1/B16 重放均逐位复现。详细检查值见 `instrument_diagnostics.csv`。
2. `run_state.json`：7 条 fresh 松耦合腿（6 原生/对照 + 1 诊断）累计 873.40 s；
   重放累计 107.94 s（含失败 standalone）；D08/D13 两遍门检查 794.25 s。
   均低于预注册时间盒；C4 未消耗，未扩扫参数。

## 5. 判定

**UNDEFINED（多吸引子根因不可判；FP 首差与首 IBL 直接响应已闭合）。**

一句话：**当前机上同档结果逐位可复现、跨档扰动首先出现在 FP，首 IBL 在同输入下
不受 BLAS 直接影响；但 B1 原生松耦合不收敛、B16 仅外层收敛且 IBL 内求解未收敛，
所以历史“两个多吸引子”在本轮不可判，不能宣称提前停止、冻结入流或真正多不动点。**

逐条机器可读判定见 `criteria.csv`；总判定与关键量见 `verdict_metrics.csv`。

## 6. 下一步

- 重开 S5 的前置证据应是：固定出处下得到两条**外层与 IBL 内层都收敛**的腿，再按
  原 2×2 预注册判别；先查 IBL 内层不收敛/既有地板，不扫 BLAS、seed 或松弛参数换答案。
- D08/D13 的 BLAS=1 值钉保持不动；当前机两遍同值红只作环境边界登记。
- 按 phase 7 范围可转入 S6 小尾巴或阶段收口；若要把本轮发现写入新能力边界，
  应明确写成“FP 扰动 + 当前不可判”，不得沿用“已证明多吸引子”的强句。

## 基线行

研究轮不改库和测试面，**不重跑贵重全套**。自测 15/15；F07 **7 passed**；
F09 **4 passed**。历史开发机套件基线承 **839 / 0 / 35 / 2**，不是当前迁移机全套
实测；当前机 D08/D13 定向证据锁两遍各 2 failed 已单列，未宣称全套转绿。
