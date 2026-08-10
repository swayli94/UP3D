# PRE_REGISTRATION — B31 C 类翼尖根治（片终止重定规格）+ LS step 语义配套评估

> **日期**：2026-07-21（跑前写就）
> **分支**：`kimi/b30-transonic-ceiling-attribution`（用户指令 2026-07-21：
> B31 在 B30 分支内开展——同 B28→B29 same-branch 先例）
> **触发**：B30 VERDICT §7 出口路由，用户裁决 = **① C 类翼尖根治为主 +
> ④ LS line search 配套评估**。B30 定论：两路径 (b) 类死法同机制
> （B30-SAME，钳制 100% 翼尖、O(10)），耗散杠杆 LS ◐ / CONF ✗——
> 耗散非约束，翼尖片终止奇点是根。

## 1. 范围

- **① C 类根治（主）**：LS 片终止函数空间 re-spec（B8 遗留，B30 提名）
  + conforming taper×片终止交互（G13.2 药效已证、生产未采用）。
- **④ LS step 语义评估（配套）**：探查证实 LS Newton **已有**
  best-of-tried 回退线搜索（newton_ls.py:899-922）——"line search
  移植"框架退役，降级为证据性评估门（GB31.4）。
- **明确不做**：roll-up / 显式尖涡（B10 rescope，模型级，另案）；
  `freeze_max_clamped` 重定规格（honesty gate，用户保留）；active-set
  滞回仪表（B30 已按条款不提名）。

## 2. 既有事实锁（证据在案，本 phase 不重测）

| # | 事实 | 出处 |
|---|---|---|
| F1 | LS 片终止 = Heaviside：末切环带有限跳量 \|δ\|≈0.026（h/taper 无关），跨一个单值单元归零；honest 尖指数 p=+0.62（去 sliver +0.367），与 conforming +0.52 同物同量 | `pyfp3d/wake/multivalued.py:94-100`, `:819-829`；track_b B8 |
| F2 | G13.2 离散机制：求解器只见最外 TE 站，Γ_last 在末单元上作集中涡脱落，诱导 ~Γ_last/h；判据 q≥1 | track_p.md:1317-1326；design.md:413-425 |
| F3 | conforming taper（Γ_eff=F·Γ_Kutta，vanish_smooth，r_c=0.05·b_semi）根治尖边：指数 +0.592→+0.009，fine M0.84 成真解；代价 cl −1.1…−1.6%、局部（η<0.95 不变） | track_p.md:1332-1343 |
| F4 | taper **与 `kutta_estimator="pressure"` 不兼容**（NotImplementedError），且生产翼身配方（probe/pressure 均含）从未携带 taper | `pyfp3d/solve/newton.py:146-151`；b18 run_demo.py:121-128 |
| F5 | B8 死路（约束侧，实测）：TE row blend 无效（病灶在**下游终止环**不在 TE 节点）；span_blend 焊死环跳量但 −20% 全局升力（唯一全局 Γ 模态被 re-level，加密更差） | track_b.md:463-501；run_b8_span_blend.py |
| F6 | B30：LS 0.7875 死 = freeze 窗语义（残差 1e-13、钳制 5+1 卡死；c=2.0 时 0.80 死于 9>8 一单元）；CONF 0.83 死 = 真 Newton 停滞（0 limited、残差 8e-6 锯齿、加耗散更差）；两腿钳制 z≈1.197–1.198（≥ B_SEMI=1.1963，骑尖边/圆帽带） | b30 VERDICT §2–3 |
| F7 | LS Newton 已有 safety backtrack（inf-norm、strict 减、λ≥0.05、best-of-tried）；conforming 有 merit 线搜索仍停滞 | newton_ls.py:899-922；newton.py:938-972 |
| F8 | 翼身网格默认圆帽（墙到 z=1.2185，片止于 1.1963）；G13.3：圆帽不造奇点但**放大**设计尖锐 TE 奇点 | wingbody.py:41,246；track_p.md:2205-2213 |

## 3. Gates（判定阈值跑前锁定）

### GB31.1 — 翼尖终止图谱（cache-only，零求解零库改动）

输入 = B30 已提交缓存（g2_ls_m0_7875 死态、g2_conf_m0_82 收敛/0_83 死态、
g3 杠杆收敛级）。交付（`results/g1_tip_atlas.csv` + `.png`）：

1. **钳制归属**：每个死级钳制单元归类 = on-cap-wall（贴圆帽壁）/
   wake-adjacent（贴片）/ beyond-tip straddler（B20 mixed-plain 类）；
   距片尖边与距尖 TE 角距离。——回答"死细胞住在帽壁还是片尖边"。
2. **LS 终止环 δ(q) 剖面**：末 15% 展向逐环跳量（cut 单元 main−aux），
   honest-metric 纪律（F1 教训：峰值只读 main 场，with/without sliver
   并报）；straddler 侧读/main 读 M² 对照（B20 类污染普查）。
3. **CONF Γ(z) 尖部剖面**：末 10 站 Γ + Γ_last（脱落集中涡强度，
   F2 判据量）；叠加 vanish_smooth r_c=0.05·b_semi 的 taper 预测
   （GB31.2 的 cure leverage 预估）。

**PASS** = 图谱落盘且每个死级钳制归属明确；**FAIL** = 归属不明 →
gate 重设计再动手。

### GB31.2 — conforming taper×生产配方（先因子实验，后采用工程）

- **2a（零库改动，因子设计）**：0.83 strict，种子 = G2 链 0.82 态，
  probe 估计器 **无 taper（对照腿）** vs probe+`tip_taper`
  （vanish_smooth，r_c=0.05·b_semi，F3  proven 参数，**治疗腿**）。
  对照腿隔离估计器切换本身的影响（生产 = pressure，F4）。
  - **✓**：对照死且治疗收敛（生产接受语义）且尖部钳制对 0+2 基线
    缩小 → taper 是因 → 以 dm=0.0125 爬升（≤0.84）记新天花板，
    触发 2b。
  - **◐**：治疗收敛但天花板不动 / 对照也收敛（估计器效应）→
    代价-收益表交用户。
  - **✗**：治疗不收敛或尖部钳制不减 → CONF 侧 C 类关闭，余路
    = B10 roll-up（§5）。
- **2b（仅 2a✓ 触发，库改动，default-off bit-identical）**：taper 兼容
  `kutta_estimator="pressure"`（语义：F_j·pressure-Kutta 行 +
  (1−F_j)·Γ-weld，F=1 处与现状逐位一致；参照 F5 死路——这是
  conforming 侧**有逐站 Γ 目标**的 blend，与 LS 侧不同）。单测 +
  生产配方（pressure+taper）重验 0.83/爬升；**采用裁决交用户**。

### GB31.3 — LS 片终止函数空间 re-spec（库改动，default-off bit-identical）

候选阶梯（预注册顺序；F5 两条约束侧死路**明确排除**）：

1. **C1 尖外渐隐 fringe**：spanwise clip 外扩 w 环（尖**外**），fringe
   内 aux 耦合渐隐 1→0（B8 命名方向：终止处函数空间/亚单元终止）。
   机制论证（写进实现注记）：aux 行是上游对流，尖**外** weld 不回灌
   升力区——与 F5 span_blend（焊尖**内**升力区）机械区分。
2. **C3 片过尖延伸**：polyline 过尖延伸（新增 outboard extend，对比
   `inboard_clip` 的 B25 语义"流体中不留自由边"）；远场
   `farfield_aux="pin_gamma"` 分支切割语义核查 + P5 branch-ray 回归门。

门禁（每候选）：LS **0.7875 strict，生产 c=1.5**，种子 demo 0.775：

- **✓**：收敛且钳制 < 基线（5+1）且以 dm=0.0125 爬升 ≥1 级且护栏全绿
  （inboard q<0.95 span Γ/cl 漂移 ≤1%；corrM ≤1.3；honest 尖指数不
  恶化；P5 回归干净）→ 采用候选，交用户。
- **◐**：收敛但天花板不动 / 钳制不降 → 代价表交用户。
- **✗**：无效应或更差 → 下一候选；全 ✗ → LS 侧 C 类关闭（阴性定论），
  余路 = B10 roll-up rescope。

### GB31.4 — LS step 语义评估（④ 配套，证据性）

交付（`results/g4_step_semantics.md`）：

1. 事实陈述落盘：LS 已有 backtrack（F7）——"移植"框架正式退役；
2. B30/G3 traces 复核：LS 濒死级残差到 1e-13 仍死（step 接受非约束）、
   CONF 带 merit 线搜索仍停滞（线搜索非约束）；
3. 仅当 GB31.3 改变 LS step 动力学（新病理：trial 全非有限 / λ 打满
   地板 / merit-inf 两表脱同步）时复检。

**判定**：无新病理 = **关闭 ④（证据性）**；有 = 提名 merit 语义升级
（独立 phase，库改动，自带测试）。

### GB31.5 — 裁决与出口

VERDICT.md + track_b 条目/台账 + 采用决定（生产配方/默认变更一律交
用户）+ demo 锚刷新（仅当任一 ✓ 被用户采纳后）。

## 4. 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|---|
| U1 | 2a 估计器切换混淆归因（probe≠pressure） | 中 | 中 | 对照腿同种子同配方仅差 taper；probe-vs-pressure 在 0.82 收敛级的 cl_p/Γ 差 RECORDED |
| U2 | C1 fringe 重演 span_blend 全局 re-level（F5） | 低-中 | 高 | 渐隐全部在尖**外**（q>span_length）+ 上游对流论证 + inboard Γ/cl ≤1% 硬护栏 |
| U3 | C3 延伸触发 P5 远场 branch-ray 伪影 | 中 | 中 | P5 回归门（远场跳跃审计）；pin_gamma 语义先核查再跑 |
| U4 | honest-metric 事故重演（B8 ×5 伪影） | 中 | 中 | 峰值只读 main 场（B20 语义）；with/without sliver 并报；钳制归属含 straddler 类 |
| U5 | medium 解成本失控（~10–30 min/级） | 中 | 低 | 全缓存；后台；预算 §6 |
| U6 | 库改动破既有门禁 | 中 | 高 | default-off bit-identical（关 = 逐位一致，测试断言）；每步库改动前后跑 b1/b2/m2/v0 + 受影响测试文件 |

## 5. 预注册的后继路由（不在本 phase 范围）

- **B10 roll-up / 显式尖涡 rescope**（模型级根治，track_p 早已记录在
  案）：C 类双腿全 ✗ 时提名，用户裁决。
- **freeze_max_clamped 重定规格**：仍为用户保留项，本 phase 不动。
- **merit 语义升级**：仅 GB31.4 检出 LS step 新病理时提名。

## 6. 成本预算（T2 类）

GB31.1 ≈ 0（缓存）；GB31.2a = 2–4 个 medium CONF 解；2b = 库改动 +
~4 解；GB31.3 = 库改动 + 每候选 2–6 个 LS medium 解；GB31.4 ≈ 0。
求解 wall ≤ ~5 h（后台）；实现 ~1.5 天。

## 7. 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- 库改动一律 **default-off 且关态逐位一致**（B20/B25 先例），开跑前
  b1/b2/m2/v0 前置门禁全绿；改动后受影响测试文件全绿。
- 结论以 committed CSV/PNG 为准；求解缓存 `results/*.npz`（~1 MB/个，
  本分支新算，按 b23/b25/b28/b30 分析缓存惯例随 phase 提交）；msh/
  日志本地；demo 侧缓存不提交（2026-07-19 用户指令）。
- 不 push；push/PR 先获用户确认（AGENTS.md）。

## 8. 证据

`results/g1_tip_atlas.csv` + `.png`（GB31.1）、
`results/g2_conf_taper.csv` + `.png`（GB31.2）、
`results/g3_ls_termination.csv` + `.png`（GB31.3）、
`results/g4_step_semantics.md`（GB31.4）、
VERDICT.md 按 §3 判定树裁决。求解缓存 `results/*.npz`。
