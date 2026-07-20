# PRE-REGISTRATION — B28 cl_fus 带外 flat-vs-tilted 解耦 + GB9.4 重设

> **日期**：2026-07-20（跑前写就）
> **分支**：`kimi/b28-cl-fus-gate-respec`（基于 origin/main 0f27900，含 B23–B27 已合入链）
> **执行依据**：B23 VERDICT §(c)（GB9.4 阈值重设建议，未执行）+ B25 注 2
> （flat-vs-tilted 模型差归因，"哪个模型更物理超出本 phase 判据范围"）
> + B26 §4 监视项（C 侧 cl_fus ≈ A 侧 ×2，增量主要在 out-band）
> **用户指令**：执行 B23 §(c) 的 gate 重设 + 补 flat-fragment 解耦实验腿

## 1. 问题复述

GB9.4 的 ≤5% 阈值已被 B23 判定误设（物理 carryover ≈ 10% of cl_p 两路径
共有）。B23 §(c) 建议重设为"袋带外（out-band）机身升力两路径一致"，但
B25 的 conforming oracle 分解（`w2_conf.csv`）实测该前提**本身不成立**：

| 配置（medium M0.5 α=3.06，committed） | cl_fus_out |
|---|---|
| conforming oracle（平片，缝贴 y=0 水线） | **0.03514** |
| A（LS q≥0 截断 + 斜片，带袋） | 0.02138 |
| C（LS fragment clip + 斜片，袋愈） | **0.05035** |

A/C 近乎等距夹住 oracle（B25 注 2：A 总量契合是 band 虚高与 out 虚低
**巧合相消**）。C-vs-oracle 带外差 0.0152（43%），B25 归因为
**flat-vs-tilted 片位置模型差**（口头归因，未判别）；B26 在跨声速重测
同型（C 侧 out 0.057–0.068 vs A 0.022–0.033）。

**本 phase 的单一问题**：C-vs-oracle 的带外差是否由片位置（平 vs 斜）
唯一决定？判别腿 = **F 侧：fragment clip（拓扑同 C）+ 平片
`sheet_direction=(1,0,0)`（位置同 conforming）**，与 C 的唯一变量是片
位置，与 conforming 的唯一变量是离散路径/拓扑。

## 2. 锚（committed，只读，不重算）

- `cases/analysis/b25_inboard_fragment_clip/results/f1_summary.csv`：
  A/C 两腿 medium α∈{0,2,3.06} + coarse α∈{0,2,3.06} 全套度量。
- `cases/analysis/b25_inboard_fragment_clip/results/w2_conf.csv`：
  conforming oracle 分解（medium α=3.06：all 0.03561，band 0.00048，
  out 0.03514，poles 0.00089）。
- `cases/demo/b9_wingbody/results/checks.csv`：GB9.4 XFAIL 现状
  （conf 0.164 / LS 0.187，legacy 配方）。

## 3. 实验设计

### 3.1 新腿（全部 M0.5，LS Picard，SOLVE_KW = freestream + pin_gamma，逐字同 B25）

| 腿 | level | α | 角色 |
|---|---|---|---|
| F | medium | 3.06 | **决定性腿** |
| F | medium | 2.0 | 趋势腿 |
| F | coarse | 3.06 | level 端点 |

α=0 **不设新腿**：`direction=(cos0,sin0,0)=(1,0,0)`，平≡斜由构造恒等，
B25 的 A/C α=0 惰性证据直接继承。

F 侧构建：`WakeLevelSet(te_polyline(FUS), direction=(cos α, sin α, 0),
sheet_direction=(1,0,0))` + `CutElementMap(inboard_clip=make_inboard_clip(FUS))`。
`sheet_direction` 是 B28 新引入的 `WakeLevelSet` knob（默认 `None` =
等于 `direction`，bit-identical 由 test_b1 新增同款测试锁定）：片的
**几何**（s 场、法向、cut 分类）由 `sheet_direction` 决定，跳量
**对流** `u_hat`（`wake_ls_coo`）仍取 `direction`（来流）——相对 C
侧单变量 = 片位置。`inboard_clip` 谓词（`fuselage.py:185-193`）只消费
穿越点坐标，与片坐标系正交，无需改。

### 3.2 度量（逐字复用 B25 `measure_f1`/`wb24.measure_e1` 核心）

- **主量**：`cl_fus_out`（带外，bw=0.06、x>1.0，同 B23/B25 定义）。
- **副量**：`cl_fus_band`（澄清带内语义：C 0.02019 vs oracle 0.00048）、
  `cl_fus` 总量、`cl_fus_poles`。
- **护栏**（对照 C 同 α 腿）：|Δcl_p| ≤ 2%；|Δγ| ≤ 5%；root te_jump
  失真 ≤ 5%；n_outer ≤ 1.5×C；res（B25 注 1 语义）；nlim/nflr = 0；
  sliver min 二面角不劣于 C − 5°；strip aux jump ≤ 1.5×γ；
  n_te_nodes = 150/76 逐腿断言；走廊 corrM ≤ 1.3 且 n_sup_corr = 0
  （袋必须保持治愈）；tip Mmax 记录（P13 监视）；n_aux_farfield census。

### 3.3 判定树（阈值先验锚定，跑前锁定）

容差 TOL = **15%**，先验锚：带外 cl_fus 自身的 coarse→medium 实现噪声
A 侧 −8.1%（0.02325→0.02138）、C 侧 −9.7%（0.05575→0.05035）；跨模型
一致性阈值不能紧于该量自身网格噪声，取 15%（≈1.5× 最大观测噪声）。

- **F1（位置唯一因子）**：|F_out − oracle| ≤ 0.15·|oracle|（=0.0053）
  **且** |F_out − C_out| > 0.15·|C_out|（=0.0076）。
  → C-vs-oracle 带外差全部归因于片位置；平-fragment 是跨模型一致配置。
- **F2（位置排除）**：|F_out − C_out| ≤ 0.15·|C_out| **且**
  |F_out − oracle| > 0.15·|oracle|。
  → 差异在拓扑/根涡强度；下一嫌疑 = M2 台账遗留（交界最内 TE 节点
  CV fan 含机身面）。
- **F3（混合/模糊）**：其余情形。报告位置因子占比
  `r = |C_out − F_out| / |C_out − oracle|`（C-oracle 差 = 0.0152），
  按 r 定性。

α=2 与 coarse 腿用于一致性复核（同向性），不独立触发分支。

副判定（不进主树）：F_band vs C_band（0.02019）vs oracle（0.00048）——
若 F_band 落到 oracle 附近，C 的带内高载同属位置效应（斜片骑上子午面
y>0 侧）；若 F_band 仍高，属 fragment 拓扑的合法近片载荷。

### 3.4 各分支的 gate 重设动作（B23 §(c) 落地方式，跑前锁定）

- **F1 落地**：B9 demo 的 LS 腿改采 flat-fragment 配置（平片 +
  inboard_clip；legacy farfield 配方不动，保持 B9 committed 对照链），
  GB9.4 重设为两款：
  (i) **硬 gate**：medium `|cl_fus_out(conf) − cl_fus_out(LS)| ≤
  0.15·|cl_fus_out(conf)|`（coarse RECORDED）；
  (ii) **RECORDED**：两路径 band/out/poles/总量分解全量报告。
  旧 ≤5% 阈值以 erratum 退役（物理 carryover 不再当误差，B23 §(c)）。
  斜片配置的带外偏差在 gate note + track_b erratum 中文档化为
  **片位置模型敏感性**（"哪个更物理"的裁决留给用户，同 B25 注 2 原意）。
  demo LS 腿需重解（coarse ~2 min + medium ~15 min）。
  附预注册条款：若重解后 GB9.5（cl_p 跨模型 <1%）因配置更换翻出界，
  记为 re-spec erratum 的一部分重锚（B27 刷新 B18 demo 同例），不算
  回归。
- **F2/F3 降级落地**：demo 配置不动；GB9.4 改写为 RECORDED 分解报告
  （≤5% 前提正式退役，带外跨模型差按实测值文档化为未决模型差），
  并在 VERDICT 提出下一归因 phase（F2 → M2 接线审计；F3 → 按 r 定）。
  B9 demo LS 腿是否顺带 demo 化 inboard_clip（B25 遗留"demo 化"项）
  留用户裁决，不在本 phase 执行。

### 3.5 风险登记

- **R1 平片非顺流**：α>0 时平片不沿来流（斜切 3.06°）——这正是
  conforming 近似本身，对照即目的；护栏验证翼面解不动。若平片腿收敛
  退化（n_outer > 1.5×C 或不收敛），先在 α=2 弱状态判读，收敛性问题
  与物理判定分开记录。
- **R2 远场相互作用**：片出口位置移到 y=0；pin_gamma 锚定出口跳量；
  n_aux_farfield census 记录。
- **R3 片擦机身水线**：与 conforming 同型（片沿 y=0 水线贴体）；
  sliver 护栏覆盖。
- **R4 knob 回归**：`sheet_direction` 默认 None ⇒ bit-identical，由
  test_b1 新增同款测试 + b1/m2/v0 前置门禁锁定。
- **R5 对流模型敏感性**：F 腿保持 û=来流（斜）而片平，û 与片不平行
  本身是一个模型选择。应急腿 **F′**（script-only，无需 knob：
  `direction=(1,0,0)` 连对流一起压平，= conforming 对流语义）仅在
  F1/F2 判读模糊（两检验比都贴近 15% 阈值，或 α=2/coarse 与决定性腿
  不同向）时启用，预算 +15 min。

### 3.6 成本预算（T2 类）

medium α=3.06 ≈ 15 min + medium α=2 ≈ 13 min + coarse ≈ 0.5 min +
CutElementMap 构建 ≈ 4 min + knob 测试 ≈ 1 min → **≤ 40 min wall**
（F1 落地追加 demo 重解 ≈ 20 min）。

### 3.7 证据

`results/f2_summary.csv`（全度量）、`results/f2_decomposition.png`
（A/C/F/oracle 分解对照）；求解缓存 `results/*.npz`（gitignored 同
B25 政策）；VERDICT.md 按 §3.3 树裁决。
