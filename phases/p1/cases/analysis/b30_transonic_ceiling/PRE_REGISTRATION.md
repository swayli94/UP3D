# PRE-REGISTRATION — B30 翼身跨声速 (b) 类天花板归因 + López 耗散杠杆

> **日期**：2026-07-21（跑前写就）
> **分支**：`kimi/b30-transonic-ceiling-attribution`（基于 origin/main `ff02d10`，
> 含 B23–B29 已合入链）
> **执行依据**：B26 VERDICT §6.2（"(b) 类高 M Newton 停滞与 conforming medium
> 0.80+ 停滞是否同一机制" = 命名的下一 phase 候选）+ 2026-07-20 审查报告
> §1.4/§2.2（(b) 类归因未裁决）——**用户 2026-07-21 裁决开工**（"我接受你的
> 建议"：现码基线 → 归因 → 便宜杠杆 → 再决定 C 类）。
> **前置锚（全部 committed，只读不重算）**：B29 刷新的
> `cases/demo/b18_wingbody_transonic/results/checks.csv` 8/8 及其 npz 缓存。

## 1. 问题复述

B29 之后翼身跨声速天花板两路径同址但不同死状（committed checks.csv）：

| 腿 | m_last | 死于 | 死状（demo 分类） | 峰位 |
|---|---|---|---|---|
| conforming medium | 0.79（strict，res 2.2e-14） | 0.80+ stall | (b) Newton 停滞 | B18 叙事 = "更锐激波/交界相互作用"（**未普查**） |
| LS flat-frag C medium | 0.775 | 0.7875 | (a+dm)：res < 1e-6 但 clamp > freeze_max_clamped=8 | 翼尖 M3.98 @ z=1.20（P13 类，B29 live 实测） |
| 两侧 coarse | **0.84 reached** | — | — | — |

审查报告 §2.2 的假设：两路径死于**同一机制** = 翼尖 P13 自由边奇点 +
高 M 下 Newton 的 lim/flr 振荡（"(b) 类"）。但 conforming 侧的 lim/flr
落点从未被普查（B18 只记了叙事），LS 侧的峰位实测在翼尖——**同机制
目前只是疑似**。本 phase 的第一问 = 用同一套普查机器把两路径濒死级
的钳制单元落点钉死；第二问 = López 的跨声速配方（高耗散爬升 + 到顶
降耗散，Table 4.13）能否不动数值代码地抬升天花板。

**与 B26 的区隔**：B26 测"袋愈后天花板动没动"（操作后果）；B30 测
"(b) 类死因的落点归因 + 耗散杠杆的操作后果"，均不裁决奇点的几何根
治（C 类 = 片终止方式重定规格 / conforming taper 交互，留给用户裁决，
见 §5）。

## 2. 锚（committed，只读）

- `cases/demo/b18_wingbody_transonic/results/checks.csv`：GB18.1/18.2
  各行（上表）+ GB18.4（C 侧 dying peak M3.98 @ tip z=1.20，走廊干净）。
- `cases/demo/b18_wingbody_transonic/results/ls_flat_medium_084.npz`：
  LS C 侧 `phi_ext` @ m_last=0.775 + levels_json（含 dying 级 res/nlim/
  nflr）；`conf_medium_079.npz`：conforming `phi` @ 0.79。
- B26 `results/g1_summary.csv` / `g1_peaks.csv`：斜片时代的 (b) 类死状
  对照（C medium 死 0.775 res ~2e-6 8/7 振荡，峰 M4.18 @ z=1.20）。
- 配方（逐字冻结，同 B18/B26/B29）：
  - LS：`LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma",
    freeze_tol=1e-4, freeze_max_clamped=8, intermediate_tol=1e-3, n_seed=30,
    direct_refactor_every=1000, n_newton_max=80)`；构建 = B29 生产侧
    `WakeLevelSet(te_polyline(FUS), direction=(cosα,sinα,0),
    sheet_direction=(1,0,0))` + `CutElementMap(inboard_clip=
    make_inboard_clip(FUS))`。
  - conforming：`CONF_SEED_KW` + `CONF_RAMP_NK`（precond="direct",
    freeze_refresh_max=8, n_newton_max=80, farfield_spanwise_gamma=True）
    + `freeze_tol=1e-5, intermediate_tol=1e-4, kutta_estimator="pressure"`；
    构建 = `cut_wake(read_mesh(CONF_DIR/medium.msh))`。
  - 公共：α=3.06，upwind 默认 c=1.5 / m_crit=0.95 / m_cap=3.0 /
    rho_floor=0.05。

## 3. 实验设计

### GB30.1 — 基线锚定（纯缓存，无求解）

从 §2 npz 缓存重建两路径锚点字段（m_last / cl_p / Mmax / dying 级
res/nlim/nflr），与 checks.csv 引用值逐位核对；同时验证缓存可作
GB30.2 的 warm-start 种子（形状与当前算子 DOF 数一致）。
**PASS** = 全部引用值在引用精度内一致且种子可加载；**FAIL** = 缓存
与 committed CSV 背离（= 独立发现，走 T1 条款，phase 继续但锚重记）。

### GB30.2 — (b) 类归因普查（2 个濒死级 strict 解 + 普查）

每路径一条腿，single-level strict 解 warm-started 自 §2 缓存（零库改动）：

| 腿 | 构建 | 级 | 种子 | 配方 |
|---|---|---|---|---|
| LS | B29 生产侧（flat+clip） | **0.7875** strict | ls_flat_medium_084.npz `phi_ext` | LS_RAMP_KW（去 ramp-only 键）+ n_seed=0 |
| CONF | conforming | **0.80** strict | conf_medium_079.npz `phi`（Γ 从零重解） | CONF_RAMP_NK + freeze_tol=1e-5 + kutta_estimator="pressure" + n_picard_seed=0 |

**续爬条款**：若某腿单级 strict 收敛（= 与级联死法路径相关，本身记为
独立发现），以 dm=0.01 续爬至**首个失败级**（≤ 0.84），普查在首个
失败级执行；若一路爬到 0.84 reached，同样记为独立发现并按 reached
处理。

**普查机器**（`wb30.py`，脚本侧复刻求解器自己的 monitor 数学，
`physics.isentropic.limit_q2_field`（cap = q²(M=3)）+ `rho_tilde`
floor 语义；LS 侧逐侧复刻 `multivalued.py:472-486`，conforming 侧复刻
workspace 钳制计数）：

1. **钳制掩码**：最终迭代态的 limited（局部 M > 3）/ floored
   （ρ̃ = 0.05）单元集合。**自校验**：掩码计数必须等于该次求解记录的
   n_limited / n_floored（不等 = 普查代码错，gate FAIL，先修再判读）。
   **2026-07-21 语义勘误（探针验证）**：conforming 的 `clamp_history`
   比返回 φ 旧一步（newton.py:672 vs :992），死级末态簿记有 ≤2 漂移；
   LS 求解器在末态重评估（newton_ls.py:928-932），保持精确匹配。
   修正后自校验 = (i) 链上最后**收敛**态掩码必须 0/0（conforming 收敛
   门禁拒收钳制态，0.82 缓存实测 0/0）且 (ii) 死级掩码与记录差 ≤2；
   LS 腿维持精确等值。
2. **落点分区**（单元质心，区域定义照 B26/B23 机器）：
   tip box = z > 0.95·b_semi（≈1.136；on/off-body 分列）；
   交界走廊 = z < 0.5 且 x > 0.8；近机身 = |dist_fus_surface| < 0.10；
   其余 = field。每钳制类 × 每区域报 计数 + 峰 M + 代表坐标。
3. **振荡结构**：step_records / clamp_history 的 n_lim/n_flr 时间序列
   （末 20 步 min/max/mean、是否单调、toggling 计数幅度）。
   **预注册限制**：逐迭代的单元级集合不可得（需库改动），本 gate 只交
   付计数级振荡 + 最终态落点；若后续需要逐迭代集合，走 §5 候选。
4. 对照量：全场峰 M 及位点、走廊 corrM、cl_p/γ（对锚）。

**判定（阈值跑前锁定）**：

- **B30-SAME**：两路径钳制集各自 **≥ 2/3 落在 tip box**（on/off-body
  合计）且总数 O(10)（< 50，非片级）→ (b) 类同机制成立 → GB30.3
  两腿同做。
- **B30-SPLIT**：conforming 钳制集主落点在交界/激波带而 LS 在 tip box
  → 归因分裂，两路径分别立项（LS = P13 类奇点；conforming = 交界/
  激波 Newton），GB30.3 只做 LS 腿，conforming 侧另报。
- **B30-AMBIG**：介于其间 → 落点直方图全量报告，用户裁决。

### GB30.3 — López 耗散杠杆（零库改动）

机制假设：López M6 以 **μc=2.0** 爬升（UP3D 默认 1.5），人工粘性加
大 → 翼尖过冲区缩小 → clamp 数落入 freeze 捕获窗 / Newton 停滞消除；
到顶后 `upwind_c_post=[1.8, 1.6]` 逐级收紧（两求解器该槽位均在位，
生产配方**从未用过**，capability_review 记录零外部调用）。

协议（对 GB30.2 的首个失败级及其后续）：

1. **L1 爬升腿**：以 `upwind_c=2.0` 重跑 GB30.2 首个失败级（同种子——
   即 G2 链中该失败级实际所用种子，T2 实测补记见 §4；同其余配方）。接受语义 = 生产语义（strict：res < tol 且 0 clamp，
   或 freeze 捕获 ≤ freeze_max_clamped）。
2. 若收敛：以 dm=0.0125 续爬（≤ 0.84），记录新 m_last。
3. **L1 降阶腿**：在新 m_last（或原天花板，若爬升未动但本级收敛）以
   `upwind_c_post=[1.8, 1.6]` 逐级 strict 重解。
4. 判定树分支 SAME → LS + CONF 两腿；SPLIT → 仅 LS 腿（conforming
   侧按 SPLIT 归因另报，若其钳制在激波带则同一协议同样适用并记录）。

**度量**：新 m_last；接受级 clamp 数；降阶存活性；**耗散代价**（对
committed 锚的 cl_p/γ/Mmax 漂移 + Cp 截面对照，RECORDED）；护栏：
走廊 corrM ≤ 1.3、n_te=76/150 断言、sliver 普查、 honesty 字段
（不普查任何未 reached 的状态为"收敛"）。

**判定（跑前锁定）**：

- **B30-L1-✓**：天花板抬升 ≥ 1 级（LS medium > 0.775、CONF medium
  > 0.79，按生产接受语义）**且** 1.6 降阶在新顶级收敛 → López 配方
  在翼身成立 → 刷新 demo 锚 + 文档收尾（耗散代价 RECORDED）。
- **B30-L1-◐**：爬升动了但降阶失败 → 天花板只能靠高耗散维持 →
  代价定量报告，用户裁决是否接受为新生产配方。
- **B30-L1-✗**：天花板不动 → 耗散非约束 → C 类（片终止重定规格 /
  奇点根治）成为下一命名候选，§5 路由。

### GB30.4 — 裁决与出口

VERDICT.md + track_b 条目/台账 + demo 锚刷新（若 L1-✓）。出口按
§5 路由，最终裁决 = 用户。

## 4. 风险登记

| # | 风险 | 概率 | 影响 | 缓解 |
|---|---|---|---|---|
| T1 | 缓存与 committed CSV 背离（B29 后代码漂移） | 低 | 中 | GB30.1 先验；背离记独立发现，锚重记 |
| T2 | 单级 warm-start 与级联死法路径相关（续爬条款触发） | 中 | 低 | **2026-07-21 实测触发**：CONF strict 暖启动链 0.80→0.81→0.82 全收敛、0.83 才死（lim=0，flr 1–2 振荡，峰 z=1.197 翼尖）——conforming "M0.80+ stall" 是级联路径相关，非硬天花板。GB30.3 CONF 腿种子 = G2 链上 0.82 态（φ+Γ，`g2_best_converged` 实现）；LS 腿无收敛链级，回退 demo 0.775 种子（= G2 该级原种子，"同种子"语义不变） |
| T3 | 普查掩码与求解器计数不等（复刻误差） | 中 | 中 | GB30.2 自校验条款；不等即 FAIL 先修 |
| T4 | medium 解成本失控（每 strict 级 ~10–40 min） | 中 | 低 | 后台跑；n_newton_max=80 不变；预算 §6 |
| T5 | upwind_c=2.0 改变已收敛物理（激波抹宽） | 高 | 低 | 正是降阶腿存在的理由；代价 RECORDED |
| T6 | Γ 从零重解使 CONF 腿轨迹与生产级联不同 | 中 | 低 | ~~φ warm 携带激波/场~~ **2026-07-21 实测证伪**：Γ=0 种子同时杀掉尾迹跳量与远场涡 BC，0.80 级 iter 0 即 lim 22k、iter 1 \|R\|=1e10 发散——改从缓存切网格 φ 用 `kutta_targets(φ, wc)` 重建 Γ 种子（与生产 ramp 的 Γ 传递同机制，探针估计、Newton 抛光） |

## 5. 预注册的后继路由（不在本 phase 范围）

- **C 类奇点根治**（LS 片终止方式重定规格 = B8 遗留 re-spec；
  conforming taper 交互）——L1-✗ 或 SPLIT-LS 出口时提名，用户裁决。
- **limiter/floor active-set 滞回 / 逐迭代集合仪表**（库改动）——
  仅当 GB30.2 计数级振荡证据明确要求时提名，自带测试与 bit-identical
  默认值，独立 phase。
- **freeze_max_clamped / m_cap 重定规格**（honesty gate 变更）——
  只能由用户裁决，本 phase 不动。
- **LS Newton line search 移植**——L1-✗ 后排序时与 C 类一并评估。

## 6. 成本预算（T2 类）

GB30.1 ≈ 0；GB30.2 = 2 个 medium strict 解（LS ~10–40 min、CONF
~10–30 min）+ 普查分钟级；GB30.3 = 每腿 1–6 个 medium 解
（≤ ~2 h/腿）。合计 wall ≤ ~4 h（后台）。实现 ~0.5 天。

## 7. 工程纪律

- 线程：NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16。
- **零 `pyfp3d/` 改动**（普查 = 脚本侧复刻 + 自校验；杠杆 = 既有
  kwargs）；开跑前 b1/b2/m2/v0 前置门禁全绿。
- 结论以 committed CSV/PNG 为准；msh gitignored；日志本地。求解缓存
  `results/*.npz` 为小规模（~1 MB/个）本分支新算分析缓存，按
  b23/b25/b28 惯例随 phase 提交（支撑"缓存重跑再生成 CSV"）；demo 侧
  种子缓存（`cases/demo/*/results/*.npz`）仍按 2026-07-19 用户指令不
  提交。
- 不 push；push/PR 先获用户确认（AGENTS.md）。

## 8. 证据

`results/g1_anchor.csv`（GB30.1）、`results/g2_census.csv` +
`results/g2_regions.csv` + `results/g2_oscillation.png`（GB30.2）、
`results/g3_levers.csv` + `results/g3_levers.png`（GB30.3）、
VERDICT.md 按 §3 判定树裁决。求解缓存 `results/*.npz`（~1 MB/个，
本分支新算，按 b23/b25/b28 分析缓存惯例随 phase 提交；demo 种子缓
存不提交，政策见 §7）。
