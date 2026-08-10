# PRE_REGISTRATION — B32 ② weld 符号冻结修复 + ① conforming taper 生产采用

> **日期**：2026-07-22（跑前写就）
> **分支**：`kimi/b30-transonic-ceiling-attribution`（同 B29/B31 same-branch
> 先例：采用工作延续在裁决分支内）
> **触发**：B31 VERDICT §8 出口路由，用户裁决 2026-07-22 = **同意 ①
> （CONF taper 生产采用 + demo 锚刷新）+ ②（weld 符号冻结修复，① 的
> 配套前提）**。
> **勘误 2026-07-22（GB32.1 ✗，按 §3 判定树执行回滚）**：② 的逐步
> 刷新语义在 medium 复验第一腿即破坏定态系统适定性——G-re（0.82，
> σ_flips=0 健康种子，② 前 8 步 0 钳制收敛）在刷新下 80 步发散
> （23 lim + 13 flr，res 4.25e-5，weld_updates=16，证据
> `results/g1.log`）。机理：逐步刷新把 F(x; s) 的 s 变成状态相关的
> 不连续切换函数（s(x)=sign diag D(x)），Newton 在符号切换流形附近
> 振荡——冻结版解的是**固定**系统（适定），刷新版解的是**切换**
> 系统（病态）。F 腿（0.84 链式种子）的病根是**种子质量**
> （wrong-side 种子），不是 pin 语义——F2 健康种子已证同配方同库
> 代码收敛。⇒ **回滚至 B31 冻结语义**（pyfp3d/tests 逐位恢复
> 9822b60，29/29 回归绿）；0.84 fresh-workspace 高 M 种子隐患以
> **F2 健康种子模式**作操作解（demo 生产 ramp 的 level-0 冻结在
> 0.70 亚临界健康态，不触发该隐患）；② 剩余价值（切换系统的
> 适定化刷新）属另案研究，本 phase 不再提名。① 按原计划继续在
> B31 冻结语义上执行。

## 1. 范围

- **② weld 符号冻结修复（先）**：B31 F2 探针确诊——0.84 链式种子
  失败是 `kutta_weld_sign` 首次残差冻结到瞬态翻号站所致（可修复，
  非估计器结构缺陷）。修复语义跑前锁定见 GB32.1。
- **① conforming taper 生产采用（后）**：b18 demo 生产配方
  （`run_demo.py` CONF legs，`kutta_estimator="pressure"`）接入
  `tip_taper`（vanish_smooth，r_c=0.05·b_semi，全 M 统一）；conforming
  腿重解 + demo 锚刷新 + 天花板爬升记录。
- **明确不做**：③ B10 roll-up / 显式尖涡 rescope（用户已问清，模型级
  另案）；`freeze_max_clamped` 重定规格（用户保留）；LS 侧任何改动；
  pytest 库默认变更（`tip_taper=None` 默认不动，default-off 逐位一致
  纪律不变）；demo LS 腿（A/C）不重解（taper 只影响 conforming 路径）。

## 2. 既有事实锁（证据在案，本 phase 不重测）

| # | 事实 | 出处 |
|---|---|---|
| G1 | taper 是因：2a 因子实验 0.83 死→治配对（同种子同配方 ±taper），治疗腿 0 钳制、mmax 2.85→1.82 | b31 `results/g2_conf_taper.csv` |
| G2 | 生产 pressure+taper 0.83 濒死级治愈（strict、0 钳制、11 步 90s）；0.84 健康种子（2a-D 态）strict 收敛（13 步 86s、0 钳制、cl_p 0.2761、corrM 1.180） | b31 `results/g2c_medium_pressure_taper.csv` 腿 E/F2 |
| G3 | 0.84 链式种子（E 的 0.83 态）失败 = 首次残差冻结 3 个翻号站的 weld 符号，尖 Γ 被反向 pin（−1.5e-4），97 lim+63 flr 极限环 | 同上腿 F；`newton.py:355-376` |
| G4 | `kutta_blocks` 每个 Newton 步重建精确 D（`newton.py:459` `cvs.newton_rows`）→ **逐步刷新符号零额外装配成本**；σ 幅值是 merit 权重（消除步中抵消），保持冻结不动 | `newton.py:430-472` |
| G5 | 采用代价（已标注）：medium 0.82 cl_p −3.00%（翼身，超 F3 单翼带 −1.1..−1.6%）；coarse M0.5 −1.05%；尖 Γ 卸载至 ~2%、尖峰 M 0.994→0.803（coarse） | b31 VERDICT §3；`g2b_coarse_pressure_taper.csv` |
| G6 | demo conforming 锚：medium 0.50/0.65/0.75/0.79 strict（cl_p 0.2173/0.2321/0.2483/0.2579，"M0.80+ stalls; recorded"）；coarse ramp 0.84 reached；cross-model M0.65 gate ≤5%（现 1.1%）、M0.75 RECORDED（1.1%）；cl_fus/cl_p = 16% @0.79 | b18 `results/checks.csv` |
| G7 | B31 库改动全部 default-off 逐位一致；`tip_taper=None` 时生产行逐位不变（pytest 断言持有） | b31 VERDICT §6；`tests/test_b31_pressure_taper.py` |

## 3. Gates（判定阈值跑前锁定）

### GB32.1 — ② weld 符号逐步刷新（库改动 + 单测 + 0.84 链式种子复验）

**修复语义（跑前锁定）**：`kutta_weld_sign` 由"首次残差冻结"改为
**每个 Newton 线性化步从 `kutta_blocks` 已重建的精确 diag(D) 刷新**
（G4：零额外装配成本）；σ 幅值维持首次冻结（merit 权重语义不变）；
**滞回地板**：`|d_j| < 0.1·median|d|`（与 σ 地板同源）时保持上一步
符号，避免瞬态过零站抖振；翻号计数改为记录相对初始冻结态的逐步
翻号历史（诊断量）。`tip_taper=None` 与 taper=1 区域逐位不变。

**单测**（`tests/test_b31_pressure_taper.py` 扩充或新增
`test_b32_weld_sign.py`）：合成翻号态符号跟踪断言；滞回地板断言；
t=1 / None 逐位不变断言；② 前后 G/E 等价态逐位断言（σ_flips=0
路径不受刷新影响）。

**medium 复验**（`run_g1.py`，缓存自 b31 种子）：

- G/E 腿（0.82/0.83，σ_flips=0）：② 后路径**逐位一致**（断言
  residual_history 逐位）；
- F 腿（0.84 链 E 种子，② 前 97+63 极限环）：目标 = **strict 收敛、
  0 钳制**，复现 F2 健康种子结果（cl_p ≈ 0.276、corrM ≤ 1.3）。

**判定**：
- **✓**：单测全绿 + G/E 逐位 + F 收敛 0 钳制 + corrM ≤1.3 → 进 GB32.2。
- **◐**：F 收敛但 G/E 路径漂移（刷新语义有副作用）→ 代价表交用户。
- **✗**：F 仍不收敛 → 回滚修复（保留 B31 状态 + F2 种子方案记录），
  报告用户再裁。

### GB32.2 — ① taper 生产采用 + demo 锚刷新（config 变更 + 重解）

**config 变更**：`cases/demo/b18_wingbody_transonic/run_demo.py` 的
CONF `newton_kw` 增加 `tip_taper=tip_taper_factors(station_z, B_SEMI,
"vanish_smooth", 0.05·B_SEMI)`——全 M 腿统一（G5：M0.5 代价 ~−1%
已知情接受）。库默认 `None` 不动。

**重解与锚刷新**（后台，全缓存）：

1. conforming medium 链 0.50/0.65/0.75/0.79 strict 重解 + coarse ramp
   0.84 重解；
2. **天花板爬升记录腿**：0.79 起以 dm=0.0125 爬升（生产 ramp 语义），
   记录新 reached 天花板（G2 预期 ≥0.83，**以实测为准**，不预支
   声称）；生产 ramp 的 workspace 冻结发生在 level 0（0.70 亚临界，
   符号健康），B31 的 0.84 fresh-workspace 链式种子隐患**不适用于
   本路径**（勘误 §头；若 0.84 级仍停滞，按实测记录 m_reached，
   m_top ≥ 0.80 为 PASS 下限）；
3. `checks.csv` 锚刷新：GB18.1 各腿 cl_p/mmax、GB18.3 cross-model gap
   （conf cl 降 → gap 预注册预期 ~2%，gate ≤5% 不变）、GB18.4 峰值
   记录行、GB18.5 cl_fus 比例、天花板记录语（"M0.80+ stalls" 按实测
   改写）；
4. `demo_report.md` / `docs/roadmap/track_b.md` / `docs/agent-rules.md`
   相关数字同步。

**护栏**：每腿 corrM ≤1.3；钳制计数全报（honesty，非零不记"干净"）；
每腿 cl_p 代价落盘 `results/g2_adoption_cost.csv`（预期 M0.5 ≈ −1%、
0.79 ≈ −3%，超 −4% 即 ◐）；pytest 全树扫——若有锚受 ②/① 影响，
逐处重钉并在 commit message 记录 Δ。

**判定**：
- **✓**：demo checks 全绿 + 代价表落盘 + 回归全树绿 → 采用完成，
  台账勾 GB32.x。
- **◐**：代价超 −4% / 钳制非零 / cross-model gap 超 5% → 呈用户再裁
  （gap gate 语义变更只能用户批）。
- **✗**：任一 demo 腿不收敛 → 回滚 config，保持 B31 状态，报告。

### GB32.3 — 收尾（证据性）

`VERDICT.md`（本目录）+ track_b 条目/台账 + `cases/analysis/README.md`
登记 + 合并回归（b1/b2/m2/v0/p8/p13/p14/b31/b32，线程 16）+
demo checks.csv 8/8 复核。

## 4. 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|---|
| V1 | 逐步符号刷新引入新极限环（符号随迭代抖动） | 低-中 | 中 | 滞回地板（预注册语义）+ F 腿门禁 + G/E 逐位断言 |
| V2 | demo 锚 churn 范围超预期（pytest 树 pin 受影响） | 低 | 低 | 全树扫；逐处 Δ 记录进 commit；库默认 None 不动 |
| V3 | 天花板爬升不及 0.83（生产 ramp 语义 ≠ 分析 strict 单级） | 中 | 低 | 以实测记录，不预支声称；天花板行是 RECORDED 不是 gate |
| V4 | cross-model gap 超 5%（conf cl 降 ~3% → gap 1.1%→~4% 上沿） | 低 | 中 | 预注册预期 ~2%；实测超 5% → ◐ 呈用户（gate 变更用户批） |
| V5 | 重解成本失控 | 低 | 低 | 全缓存后台；预算 §6 |

## 5. 预注册的后继路由（不在本 phase 范围）

- **③ B10 roll-up / 显式尖涡 rescope**：LS 侧 C 类关闭后的余路，
  模型级根治，另立 phase（用户已问清内涵，未裁决启动）。
- **`freeze_max_clamped` 重定规格**：用户保留项，本 phase 不动。
- taper 半径/形式扫描（代价分解）：仅当用户认为 −3.00% 代价需优化
  时另立。

## 6. 成本预算（T2 类）

GB32.1 = 库改动 + ~4 单测 + F 腿 1 个 medium 解（~8 min；G/E 走缓存
逐位断言）；GB32.2 = conforming demo 腿重解 ~15–30 min（后台）+
爬升 ~4 级（~10 min）；GB32.3 = 回归 ~2 min。求解 wall ≤ ~1 h。

## 7. 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- 库改动 default-off 逐位一致（`tip_taper=None` 路径 pytest 断言）；
  ② 只改 taper 激活路径的符号来源，不改行结构。
- 证据纪律：结论以 committed CSV/PNG 为准；npz 小缓存按 b23–b31
  分析缓存惯例随 phase 提交；demo 侧种子缓存不提交（2026-07-19
  用户指令）。
- 不重算贵重已提交 artifact；b31 缓存（g2a/g2c 种子）直接复用。
