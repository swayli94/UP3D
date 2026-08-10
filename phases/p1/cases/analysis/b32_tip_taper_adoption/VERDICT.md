# VERDICT — B32 ② weld 符号冻结修复 + ① conforming taper 生产采用

> **日期**：2026-07-22
> **分支**：`kimi/b30-transonic-ceiling-attribution`（B29/B31 same-branch
> 先例：采用工作延续在裁决分支内）
> **预注册**：`PRE_REGISTRATION.md`（跑前写就，含 GB32.1–32.3 判定树与
> GB32.1 勘误头部；本文件按其执行）
> **判定**：**GB32.1（②）✗ 回滚**——逐步刷新 weld 符号把固定系统变成
> 状态相关**切换系统**（s(x)=sign diag D(x)），病态：健康 0.82 种子
> （σ_flips=0，修复前 8 步 0 钳制）在刷新下 80 步发散（23 lim+13 flr，
> 证据 `results/g1.log`）；按判定树 ✗ 分支回滚至 B31 冻结语义
> （pyfp3d/tests 逐位恢复 9822b60，回归绿）。F 腿（0.84 链式种子）病根
> 重定性为**种子质量**，0.84 fresh-workspace 隐患以 **F2 健康种子模式**
> 作操作解。**GB32.2（①）✓ 采用落地**——b18 生产配方接入 tip_taper
> （vanish_smooth，r_c=0.05·b_semi，全 M 统一），conforming medium 链
> 0.50–0.79 strict 重解全腿 0 钳制，**天花板爬升 0.79 → 0.84 REACHED**
> （cl_p 0.2738、Mmax 2.14、0 钳制、res 1.9e-14），"M0.80+ stall" 记录
> 退役；demo checks **8/8 PASS**；采用代价 ≈ −1.3%（F3 带内，7/7 腿
> 0 钳制、corrM ≤1.180）；跨模型 gap 改善（M0.65 1.1→0.3%、M0.75
> 1.1→0.2%）。**GB32.3 = 本文件**（合并回归 170 passed, 4 skipped）。
> **证据**：`results/g1.log`（GB32.1 失败腿；`*.log` 按仓库惯例不入库，
> 由 `run_g1.py` 重跑再生成，失败终态缓存
> `results/g1_m0_82_ptaper_b32.npz` 已提交）、
> `results/g2_adoption_cost.csv`+`.png`（由 `run_g2.py` 重跑自 demo
> 缓存可再生成）、demo 侧 `cases/demo/b18_wingbody_transonic/results/`
> （checks.csv / cl_vs_mach.csv / cross_model.csv / PNG /
> `b32_demo_rerun.log`）。

---

## 1. GB32.1 — ② weld 符号逐步刷新：**✗ 回滚（按判定树）**

**修复语义（预注册锁定）**：`kutta_weld_sign` 由首次残差冻结改为每个
Newton 线性化步从 `kutta_blocks` 已重建的精确 diag(D) 刷新（零额外装配
成本）；σ 幅值维持首次冻结；滞回地板 `|d_j| < 0.1·median|d|` 时保持
上一步符号。单测（合成翻号跟踪 / 滞回地板 / t=1 逐位）曾写并通过。

**medium 复验第一腿即失败**：G-re（0.82，σ_flips=0 健康种子，修复前
8 步 0 钳制收敛）在刷新语义下 **80 步发散**——23 lim + 13 flr，
res 4.25e-5，weld_updates=16（`results/g1.log`）。预注册风险 V1
（逐步刷新引入新极限环）实发。

**机理**：逐步刷新把 F(x; s) 的 s 变成状态相关的不连续切换函数
s(x)=sign diag D(x)——冻结版解的是**固定**系统（适定），刷新版解的是
**切换**系统（病态），Newton 在符号切换流形附近振荡。G/E 腿逐位断言
的前提（σ_flips=0 路径不受刷新影响）不成立：符号源一旦状态相关，
即使未翻号，迭代映射本身已变。

**F 腿病根重定性**：0.84 链式种子失败的根因是**种子质量**
（wrong-side 种子在首次残差时把 3 个翻号站的符号冻进 weld），不是
pin 语义缺陷——B31 F2 已证同配方同库代码换健康种子即 strict 收敛。

**回滚与操作解**：

- pyfp3d/ 与 tests/ **逐位恢复 9822b60**（`git diff --stat 9822b60` 对
  这些路径为空），回滚后 29/29 回归绿；`tests/test_b32_weld_sign.py`
  随库改动一并移除。
- 0.84 fresh-workspace 高 M 种子隐患 → **F2 健康种子模式**作操作解；
  生产 ramp 的 level-0 workspace 冻结发生在 0.70 亚临界健康态，**不
  触发**该隐患（GB32.2 爬升腿实测确认，§2）。
- ② 剩余价值（切换系统的适定化刷新）属另案研究，**本 phase 不再
  提名**。

## 2. GB32.2 — ① taper 生产采用 + demo 锚刷新：**✓（8/8 PASS）**

**config 变更**（`cases/demo/b18_wingbody_transonic/run_demo.py`）：CONF
腿 `newton_kw` 增加 `tip_taper=tip_taper_factors(station_z, B_SEMI,
"vanish_smooth", 0.05·B_SEMI)`，全 M 统一；库默认 `tip_taper=None` 不动
（default-off 逐位一致纪律不变）；LS 腿（A/C）不碰、走缓存。

**重解结果**（日志 `cases/demo/b18_wingbody_transonic/results/b32_demo_rerun.log`）：

| 腿 | cl_p | 钳制 | 备注 |
|---|---|---|---|
| medium 0.50 | 0.2143 | 0+0 | strict；锚重钉（原 B9/CL_M05 历史锚 0.2173） |
| medium 0.65 | 0.2290 | 0+0 | strict |
| medium 0.75 | 0.2450 | 0+0 | strict |
| medium 0.79 | 0.2545 | 0+0 | strict |
| **medium 爬升** | **0.2738 @ 0.84** | **0+0** | **m_reached=0.84**，Mmax 2.14，res 1.9e-14 |
| coarse 0.60 | 0.2155 | 0+0 | strict |
| coarse 0.84 | 0.2590 | 0+0 | ramp reached（proof-of-concept） |

- **天花板 0.79 → 0.84 REACHED**：生产 ramp 链（level-0 冻结在 0.70
  亚临界健康态，GB32.1 隐患不适用）从 0.79 爬至 0.84，超预注册下限
  m_top ≥ 0.80 与 G2 预期 ≥0.83——**以实测为准，未预支声称**。
  "M0.80+ stalls; recorded" 记录语就此退役。
- demo checks **8/8 PASS**（GB18.1/18.2 PASS 门；GB18.3 M0.65 PASS
  ≤5%；M0.75/18.4/18.5 RECORDED 语同步刷新）。

## 3. 采用代价与护栏（`results/g2_adoption_cost.csv`+`.png`）

| 腿 | cl_p 新 | 参考 | 代价 | corrM | 钳制 |
|---|---|---|---|---|---|
| medium M0.50 | 0.2143 | 0.2173（B9/CL_M05 历史锚） | **−1.37%** | 0.617 | 0+0 |
| medium M0.65 | 0.2290 | 0.232114（旧 conf 缓存） | **−1.35%** | 0.837 | 0+0 |
| medium M0.75 | 0.2450 | 0.248301（旧 conf 缓存） | **−1.32%** | 1.009 | 0+0 |
| medium M0.79 | 0.2545 | 0.257883（旧 conf 缓存） | **−1.31%** | 1.074 | 0+0 |
| medium 爬升 0.84 | 0.2738 | —（新纪录，无参考） | — | 1.180 | 0+0 |
| coarse M0.60 | 0.2155 | 0.217796（旧 conf 缓存） | **−1.04%** | 0.732 | 0+0 |
| coarse M0.84 | 0.2590 | 0.261749（旧 conf 缓存） | **−1.03%** | 1.076 | 0+0 |

- **判定 ✓**：无腿触 −4% ◐ 线；corrM 全部 ≤1.3（最大 1.180，与 B31 F2
  健康种子 0.84 的 1.180 一致）；钳制计数全 0（honesty 全报）。
- **两种代价语义并陈**：本表 = **生产 ramp 链**语义（level-0 冻结在
  0.70、温链式种子）→ medium ≈ −1.3%，落 F3 单翼带（−1.1..−1.6%）内；
  B31 分析网格 G 腿 = **单级 strict**语义（b30 0.81 种子直接解 0.82）
  → −3.00%。两者均为实测记录；−3.00% 不是本采用配置的生产读数，
  机理差异（温种子链 vs 单级）本 phase 不分解、不预支。
- coarse 两腿 −1.03/−1.04% 与 G5 既有事实锁（coarse M0.5 −1.05%）
  一致。

## 4. 跨模型与锚重钉

- **M0.65**：conf 0.2290 vs LS flat-frag 0.2296 → gap **0.3%**（采用前
  1.1%），gate ≤5% PASS——conf cl 降向 LS 靠拢，gap 改善而非恶化
  （预注册 V4 未触发）。
- **M0.75**：0.2450 vs 0.2456 → **0.2%**（采用前 1.1%）RECORDED。
- **M0.5 锚重钉**：0.2173 → **0.2143**（B9/CL_M05 历史锚按预注册
  GB32.2 item 3 让位）；对 LS 0.2184 的 gap = **1.9%** RECORDED
  （无 gate）。
- **GB18.5**：cl_fus 0.0423 = 17% of 0.2545 @M0.79（RECORDED 语保留）。

## 5. 成本与护栏

- 求解成本：GB32.1 库改动 + 单测 + G-re medium 复验（~8 min）；
  GB32.2 demo 重解（后台，~40 min 量级，一次性，不重跑）+ 爬升；
  GB32.3 合并回归 110 s。在预注册 §6 预算内。
- 护栏：corrM ≤1.3 ✓（§3）；honesty ✓（GB32.1 失败腿连同 23+13
  钳制与 weld_updates=16 全报；爬升 0.84 为 strict 收敛态才计
  天花板）；pytest 全树扫——锚 churn 限于 demo results 的
  CSV/PNG 与三处文档现状段（track_b / agent-rules / demo_report +
  overview 勘误句），历史 phase 记录（B9/B17/B18/B26/B27 等 closed
  条目内数字）一律不改。
- **库零净改动**：pyfp3d/ 与 tests/ 对 9822b60 逐位一致（GB32.1
  回滚后）；B31 的 default-off taper 移植原样持有
  （`tests/test_b31_pressure_taper.py` 54+15+35 绿）；本 phase 唯一
  生产变更 = `run_demo.py` 的 CONF 腿 config。
- 合并回归（v0/b1/b2/m2/p8/p13/p14/b31，线程 16）：**170 passed,
  4 skipped**。

## 6. 预注册判定树终态

- GB32.1 **✗**（F 腿目标未达且健康腿被破坏）→ 按 ✗ 分支：回滚修复、
  保留 B31 状态 + F2 种子方案记录、勘误入预注册头部 →
- GB32.2 **✓**（demo checks 8/8 PASS + 代价表落盘 7/7 绿 + 回归绿）
  → 采用完成，台账勾 GB32.x →
- GB32.3 = 本文件。

## 7. 出口路由（交用户）

1. **③ B10 roll-up / 显式尖涡 rescope**：LS 侧 C 类关闭后的余路
   （B31 §8-③，用户已问清内涵），模型级根治，**未启动**——如需抬
   LS 天花板另立 phase。
2. **`freeze_max_clamped` 重定规格**：用户保留项，本 phase 未动
   （LS 0.80@c=2.0 一单元之差的窗语义问题原样保留）。
3. **weld 符号切换系统的适定化刷新**：GB32.1 证明朴素逐步刷新
   病态；若要消除 fresh-workspace 高 M 种子隐患需切换系统层面的
   适定化设计，另案研究，**不提名**。
4. **taper 半径/形式扫描（代价分解）**：仅当用户认为 −1.3% 生产
   代价（或 −3.00% 单级读数）需优化时另立。

**demo 锚刷新：已执行**（① 采纳条款，本 phase 内完成）；LS 侧零改动；
pyfp3d/ 零净改动已核验（对 9822b60 diff 为空）。
