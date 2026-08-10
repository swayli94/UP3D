# PRE-REGISTRATION — B27 B18 demo 刷新（LS 腿复活：inboard clip）+ Track B 文档收尾

> **日期**：2026-07-20
> **分支**：`kimi/b25-inboard-fragment-clip`（叠在 B25/B26 上，未 push 链）
> **状态**：预注册，待用户裁决后实施
> **前置**：B26 VERDICT（B26-A 成立：袋愈后 LS 天花板 coarse 0.84
> reached / medium 0.50→0.7625，死于 (b) 类翼尖-激波 Newton 停滞；
> T1 独立发现 = A 侧锚点背离系 B21/B22 freeze-capture 效应）、
> B25 VERDICT（`inboard_clip` 已入库）、2026-07-19 全库 inspection
> （D1：B15–B18 回顾性章节普遍未勘误的系统性 stale 债）

---

## 1. 物理推导链（为什么是这一条）

B26 已在 analysis 层面判定：**袋是 LS 翼身跨声速天花板的限制器，袋愈后
LS 天花板与 conforming 同址**。但仓库的门面证据仍讲旧故事：

- `cases/demo/b18_wingbody_transonic/` committed checks.csv / cl_vs_mach.csv /
  cross_model.csv / PNG 全部仍是"LS junction-limited（medium 死 0.50、
  coarse 死 0.55）"的 B18 叙事（post-B20、pre-B21 基线）。B22 刷新了
  B15/B14 demo，**B18 是唯一未在 B21/B22 态重跑过的 demo**。
- `run_demo.py:12-19` docstring 数字是初版遗留，与 committed checks.csv
  本就不一致（B26 §1 已记）；现在连 checks.csv 也过时了。
- 文档五处落点（demo_report / design_track_b §18 / agent-rules 相位表 /
  overview B-track 表 / roadmap B-track 行）均写"LS 交界受限
  （closed-negative）"，且 B26 的 T1 发现（A 侧 0.50 级带袋收敛）尚未在
  任何文档留痕。

证据纪律：结论以 committed artifact 为准——B26 的 CSV/PNG 已提交，但
**demo 是交付门面**，demo 不刷新则 B26-A 不落账。本 phase = 预注册 §3
B26-A 出口的执行：**B18 demo 刷新（LS 腿复活）+ GB9.4/GB20.5 重设 +
Track B 文档收尾**。

**与 B26 的关系**：B26 是判定实验（analysis case，A/C 同码对照）；B27 是
交付刷新（demo 门面重述 + 文档勘误），不含新物理判定。唯一的**新测量**
是跨模型升档：B26 之前 medium 无共同跨声速 Mach（LS 离不开 0.5），现在
LS+clip 可到 0.75——GB18.3 的"跨模型只剩 M0.5"可以升级。

**不在本 phase 范围**（边界承诺）：

- (b) 类新天花板归因（C medium 0.775 vs conforming 0.80+ 是否同机制）——
  下一 phase 候选，本 phase 只记录现象。
- cl_fus out-band ×2 的治理（P11/曲面壁元路线的输入，本 phase 只记录）。
- B9/B17 的 M0.5 锚（不动）；conforming 配方（不动）；库代码（零改动，
  `make_inboard_clip` B25 已入库）。

## 2. 实验设计

### 2.1 改动面（预注册承诺）

| 文件 | 改动 |
|---|---|
| `cases/demo/b18_wingbody_transonic/run_demo.py` | docstring 重写（新叙事）+ `ls_mesh` 增加 clip 侧 builder（`inboard_clip=make_inboard_clip(FUS)`）+ LS 腿改 A/C 对 + 新增跨模型腿 + 图/CSV 刷新 |
| `cases/demo/b18_wingbody_transonic/results/*` | 全量重生成（checks.csv、cl_vs_mach.csv、cross_model.csv、3 PNG） |
| `docs/demo_report.md` | B18 行刷新 |
| `docs/design_track_b.md` | 新 §22（B26/B27）+ §18 erratum 指针（§19–21 体例） |
| `docs/agent-rules.md` | B26 ✓ CLOSED（B26-A）+ B27 ✓ CLOSED 条目；next-phase 优先级更新 |
| `docs/overview.md` | B-track 表 GB20.5 故事勘误（袋归因：刻面 G1.6 → 自由边奇点已愈；G1.6 退居 cl_fus 嫌疑） |
| `docs/roadmap.md` | B-track 行 B18 描述刷新 |
| `cases/analysis/b27_b18_demo_refresh/` | 本文件 + `results/g27_consistency.csv`（demo 重跑 vs B26 committed 的逐项 diff） |

### 2.2 腿（全部重型、显式启动、后台）

本 worktree **无任何 demo npz 缓存**（gitignored，主侧手工同步政策），
故 demo 全腿重解；conforming 腿顺带完成 B21/B22 态的首次复测
（B22 只刷了 B15/B14）。

| # | 腿 | 说明 | 预估 |
|---|---|---|---|
| L1 | conf coarse 0.84（0.70 起）+ coarse 0.60 cross | GB18.1 复测 | ~10–25 min |
| L2 | conf medium 0.65（0.60）+ 0.79（0.70） | GB18.1 复测 | ~20–50 min |
| L3 | **conf medium 0.75（0.70）新点** | 跨模型第二横坐标 | ~10–25 min |
| L4 | LS A/C ceiling probe coarse 0.50→0.84 | = B26 A/C coarse（bit 复测） | ~4 min |
| L5 | LS A/C ceiling probe medium 0.50→0.84 | = B26 A/C medium（bit 复测） | ~66 min |
| L6 | LS+clip medium 定向 0.50→0.65、0.50→0.75 | 跨模型 LS 侧（末级 strict） | ~10 min |

LS ceiling probe 配方 = B26 冻结配方逐字（LS_RAMP_KW、dm=0.05、
dm_min=0.01、α=3.06）；conforming 配方 = B18 原样（CONF_SEED_KW /
CONF_RAMP_NK，freeze_tol=1e-5）。合计 ~2–2.5 h（T2 后台）。

### 2.3 度量

1. **bit 一致性**（主度量）：L1/L2/L4/L5 对 committed 锚逐项 diff →
   `g27_consistency.csv`（demo 重跑 vs B18 committed checks.csv 的
   conforming 行、vs B26 committed g1_summary/g1_levels 的 LS 行）。
2. 跨模型：0.65/0.75 两侧 cl_p（LS+clip 用定向 ramp 的**收敛末态**
   phi_ext 积 cl_p；conforming 用 conf_ramp 收敛态），|gap|%。
3. 袋表征（GB18.4 重答）：A 侧濒死级峰位（交界条带）vs C 侧濒死级峰位
   （翼尖）+ 走廊 corrM（引用 B26 committed g1_peaks，demo 内不重复
   topk 机器，直接 hardcode 引用并注出处）。
4. cl_fus（GB18.5 刷新）：C 侧新天花板态的 cl_fus 与 out-band 分解
   （引用 B26 committed g1_summary：0.078 / out 0.057）。
5. cl(M) 全表：conforming 0.50/0.65/0.75/0.79(medium)+0.84(coarse)、
   LS+clip 0.50/0.65/0.75/0.7625(medium)+0.84(coarse)。

### 2.4 判定（预注册 gates）

- **GB27.1（PASS 判据）**：conforming 腿复现 committed 锚（cl_vs_mach
  4 位小数一致）。B21/B22 只动 LS 路径，期望 bit 一致；漂移 →
  RECORDED 并作 B21/B22 独立发现披露（不判 B27 失败）。
- **GB27.2（PASS 判据）**：LS ceiling A/C 复现 B26 committed——
  C coarse reached 0.84；C medium m_last=0.7625、死于 0.775 (b) 类；
  A medium m_last=0.50、死于 0.5125 (a) 类；A coarse m_last=0.82。
  同码同线程期望 bit 一致；漂移 → 披露于 consistency CSV。
- **GB27.3（PASS 判据 / RECORDED 兜底）**：新跨模型——0.65：
  |cl_p gap| ≤ 5% 判 PASS（B9/B17 M0.5 的 2.6% 为参照，跨声速放宽）；
  0.75：RECORDED（激波区敏感，不设阈值，只披露）。
- **GB27.4（RECORDED）**：GB18.4 重答（袋 = 自由边奇点，C 侧已愈，
  残余限制器 = 翼尖 P13 类 + 高 M Newton）+ GB18.5 刷新（cl_fus 新顶
  值 + out-band ×2 → P11 输入）。
- **GB27.5（RECORDED）**：T1 勘误入 demo（A 侧复测 vs B18 committed
  锚 = B21/B22 freeze-capture 效应；A medium 0.50 级带袋 Mmax 5.22
  现收敛，袋真实杀伤线 0.55 Mmax 13.1 > freeze_max_clamped=8）。

demo 叙事（docstring + checks 注）预注册为：**conforming 到 0.84/0.79；
LS 无 clip 仍袋限（medium (a) 类死 0.5125、coarse (b) 类死 0.84）；
LS+inboard clip 到 0.84 coarse / 0.7625 medium——两侧天花板同址，
交界袋曾是限制器（B25/B26），残余限制器与 conforming 同类（高 M
Newton/激波，峰在翼尖）；跨模型升档为 M0.5(2.6%)+0.65+0.75**。

## 3. 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|---|
| T1 | conforming 腿漂移（B21/B22 本应惰性） | 低 | 中 | GB27.1 条款：漂移记独立发现，不判失败 |
| T2 | 成本 ~2–2.5 h（无缓存全腿重解） | 高 | 低 | 后台跑；LS 腿成本已由 B26 实测（69 min） |
| T3 | L6 定向 ramp 0.75 末级 strict 行为与 ceiling probe 中 loose 通过不同 | 低 | 中 | B26 0.7625 strict 已收敛作先验；若死，跨模型第二点退 0.70（另补 conf 0.70） |
| T4 | 文档范围蠕变（5 处落点外的连锁引用） | 中 | 低 | 只改 2.1 表列文件；其余引用留 inspection 债记录 |
| T5 | demo 重跑与 B26 committed 数值漂移（线程/numba 非确定性） | 低 | 中 | 同 16 线程环境；GB27.2 条款披露于 consistency CSV |

## 4. 工程纪律

- 库代码零改动；默认路径 bit-identical（A 侧 = 现码原样）。
- 开跑前置：b1/b2/m2/v0 复核一次（本 phase 仍零库改动，跑前跑一次即可）。
- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- demo 重型腿维持 `PYFP3D_TRANSONIC_GATES=1` 门控惯例。
- 结论以 committed CSV/PNG 为准；npz/msh gitignored。
- 不 push；push/PR 先获用户确认。

## 5. 时间与产出

- 实现（demo 改造 + consistency diff 机器）：~0.5 天
- 求解：~2–2.5 h（后台）
- 文档五处 + VERDICT：~0.5 天
- 产出：刷新的 B18 demo（代码 + committed CSV/PNG）、
  `results/g27_consistency.csv`、五处文档勘误、本 phase VERDICT.md

## 6. 与其他路线的关系（备忘）

- B27 落地后：LS 翼身全包线评估真正开题（B25 §7 落地）；conforming
  退居交叉验证；Track V 片拓扑前置全部就位。
- 下一 phase 候选（本 phase 不做）：(b) 类天花板归因——LS+clip 与
  conforming 的 medium 天花板若同机制，后续投入转向 Newton 鲁棒性
  （两路径共享），否则 LS 侧仍有独立问题。
- cl_fus out-band ×2 → P11（曲面壁元）输入监视项。
- P11 与本 phase 判据正交（P11 修几何根因，本 phase 交付操作后果）。
