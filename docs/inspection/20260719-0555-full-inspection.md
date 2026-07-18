# 全面检查报告（Kimi 第二轮独立审查）

- **Date**: 2026-07-19 05:55 CST
- **Auditor**: Kimi Code CLI（独立于主 authoring agent）
- **Basis**: branch `kimi/inspection-b20-capability-review` @ `origin/main` = `3db08b9`
  （B20 re-baseline）。覆盖 A3 响应落实、B16–B20 全部新工作、文档系统、
  规划系统、能力矩阵。
- **Method**: (a) 独立复跑默认全套件；(b) 四路并行只读审查——B15–B20 文档↔提交
  CSV 逐数字核对 / B19–B20 代码改动公式级复推导 / 能力矩阵与开放项普查 / A3 响应
  逐条落实核实；(c) 所有头条发现由主审在代码与 CSV 上亲自二次核实。
  **未重跑任何贵重 demo**（P4 heavy、P5 medium、B15 ramp 等），以 committed
  CSV/PNG 为准；除本报告外未修改任何文件。

---

## 0. 头条结论

1. **代码：A3 的 12 项修复逐条核实全部落实；B19/B20 的数值修复本身经公式级
   复推导为正确，测试锁定真实。** 本轮新发现 **1 个真实 bug 并已实测确认**：
   `freeze_side_state` 漏打 B20 补丁（N1，3-D only；M6 coarse M0.70 下
   shipped 捕获与 patched 捕获的上游选择 83 个元素不同、分支 9 个不同，全部
   落在 aux-touching mixed-plain 类上），它是开放项 **GB20.7 的头号机制候选**
   ——一行可修，判别实验便宜。
2. **基线独立复跑吻合**：`465 passed + 22 skipped + 2 xfailed`（1179.17 s
   @16 线程），与文档逐字一致。门控套件（`PYFP3D_TRANSONIC_GATES=1`）复跑
   结果见 §7 补记。
3. **文档：B19/B20 新写的内容与提交 CSV 精确一致（抽查约 40 个数字全部对
   上）；但 B20 re-baseline 之后，B15–B18 的回顾性章节普遍未勘误**——新旧两套
   数字并存于同一文档集，GB20.7 这个开放项在 7 个文档面中只披露了 2 个。
   另有两处 pre-existing doc error（cond1 数值、design.md §12 陈旧 spec）和
   一处 A3 响应文档的夸大声明（M2 ledger 实际未修）。
4. **计划：A3→B19→B20 的决策路径健康、纪律有效**（alias bug 靠"测量异常而非
   解释异常"抓到，没有被记成物理发现）。当前最关键开放项是 GB20.7；建议
   下一步顺序：N1 修复 + GB20.7 判别 → 3-D LS 数字测试锁 → 文档勘误一波 →
   LS fine route / Track V / P11（P11 的受益面在 B20 后变大，见 §4）。

---

## 1. 代码正确性

### 1.1 A3 修复逐条核实：12/12 落实（file:line 证据）

| 项 | 结论 | 证据 |
|---|---|---|
| C2（conforming fail-fast 跨 epoch） | ✓ 按描述修复 | `solve/newton.py:726/772/833` 三处 `r_level_best = np.inf` reset |
| C3（freeze_max_reverts 回移） | ✓ | `newton.py:454`（=3）、`:647`、`:729-730` 永久 disarm |
| C4（reader 丢未命名面组） | ✓ + 测试 | `mesh/reader.py:146-162` placeholder；`tests/test_mesh_reader_roundtrip.py:135-157,174-180` |
| C5（reader 体组崩溃） | ✓ | `reader.py:194-197` `volume_<i>` padding；测试 `:159-172` |
| C6（`dy == 0.0` 精确相等） | ✓ | `constraints/dirichlet.py:93` 纯 mask 成员判定，与审计建议逐字一致 |
| C7a（TE 无 aux DOF 无守卫） | ✓ scoped | `wake/cut_elements.py:233-236`（`n_ext_dofs > 0` 限定；首版过宽已诚实记录） |
| C7b（TE 节点强制 shift） | ✓ | `cut_elements.py:146` `shift_mask |= te_mask` |
| T1（B3 测试 2% vs gate 0.3%） | ✓ | `tests/test_b3_lifting.py:124-131`（0.003，引 gate + 实测 0.1441%） |
| T2（validate_coloring 返回丢失） | ✓ | `tests/test_mesh_coloring.py:59` |
| P1（section_cp 丢 gamma） | ✓ | `post/unified.py:110/122` → `section_cut.py:376/405/229/268`，默认 1.4 |
| F0（A1 计时 flake） | ✓ | `tests/test_a1_instrumentation.py:33-34`（REPO_ROOT）、`:76`（20%） |
| P2（physics-bounds 矛盾） | 确认未修 = 符合 BACKLOG 声明 | `isentropic.py:370-371` 原样；但 backlog 状态**只记在 response doc**——`track_a.md` 的 A3 条目未列，代码处无注释指针，建议补 |

**A3 响应文档的一处夸大**：response 称 finding 2（M2 权威矛盾）"resolved in
track_m's title, checkbox **and ledger**"——title/checkbox 确已修，但
**ledger 行 `docs/roadmap/track_m.md:247` 仍 `| M2 | ◐ |`**，`docs/roadmap.md:13`
与 `docs/overview.md:131` 也仍写 "M2 ◐"。权威级矛盾以缩小的形式残留至今。

### 1.2 B19/B20 修复正确性（本轮重点复核，公式级）

复核方式：对残差 `R_a = ρ̃(|g_read|²)·V·(g_scatter·B_a)` 逐元素类重推 Terms 2/3
的行/列因子与行列映射。

- **J = dR/dφ 现在在全部四类元素上成立**：cut（readvec≡dofvec，
  `wake/multivalued.py:590-592` 有断言）、te_lower（两映射一致地走 aux，
  正确地不补丁——那是真 aux 读）、single-side plain（side 场 ≡ main 场）、
  mixed-side plain（B20 后两个映射都回到 main）。**不再存在任何
  side 梯度×main DOF 的交叉路径**——组合后坍缩为单场 Newton。
- **knob 移除干净**：`pyfp3d/`、`tests/` 内无 `plain_density` 残留；
  `_apply_main_density` 无条件、只作用于 mp mask 内；历史 A/B 脚本以
  `SystemExit` + 复现指针禁用（如 `cases/analysis/c1_ls_jacobian_fd/run_legb_apply.py:24-30`）。
- **alias 修复完整**：18 个 `PicardOperator.velocities` 消费者逐一检查，全部
  "下一次调用前消费"或已 `.copy()`；其它共享缓冲模式（`_R`、`_rho_tilde` 等）
  同样安全。
- **测试锁定真实、非空转**：`tests/test_b19_jacobian_3d.py` 有 premise asserts
  （`len(msp)>0`、`touched/control>0`、`n_nu_active>50`）；erratum 处理正确——
  断言的是"129 个 mixed 元素中 0 个触 cut 节点"这个真不变量，而非旧的假命题。
- **GB20.3 机制成立**：quasi-2D 上 mp 元素不触 cut 节点 ⇒ 补丁写入逐位相同值；
  mask = sign-change candidates − cut − te_lower，实测不过包不欠包。

### 1.3 新发现

**N1 [bug，已实测确认，3-D only] `freeze_side_state` 漏打 B20 补丁——GB20.7 头号嫌疑**

- `pyfp3d/wake/multivalued.py:739-745`：`freeze_side_state` 在**未补丁的 side
  场**上算 `q2l/rho` 并捕获 `(up_sel, branch)`；而 `newton_side_data`
  （`:667-668`）的 live 和 frozen 两个分支都先过 `_apply_main_density`。
- ⇒ docstring 承诺的不变式"the state is captured at the SAME q2l/rho
  `newton_side_data` computes, so the frozen sweep reproduces the live density
  **bitwise** at the freeze point"（`:730-732`）在 mixed-side plain 触-cut 元素
  上破裂。
- **本轮已用探针实测确认**（`docs/inspection/20260719-n1-freeze-probe.py`，
  M6 coarse、M0.70 seeded 跨声速态，与 GB19.6 同一工况）：shipped 捕获与
  B20-patched 捕获的选择差异——upper 侧 **upstream 83 个元素不同、branch 9 个
  不同**，且差异**全部**落在 aux-touching mixed-plain 类上（83/83、7/7；
  该类总数 252，与 GB19.6 一致）；lower 侧无差异。⇒ 不变式破坏不是理论可能，
  是在该工况下实际发生：shipped `freeze_side_state` 捕获的上游/分支选择，
  有 83+7 个元素不是 live 系统（B20 后）在该状态会做的选择。
- 后果有界：冻结相自身的 J 与自身 R 一致（不是 Jacobian 错误），是**冻结选择
  错配**——3-D 冻结 Newton 收尾（B15 recipe，M6 ramp 每级 armed）从冻结点起
  就在一个与 live 系统不一致的选择上迭代。
- 现有测试看不见：`tests/test_b15_ls_newton_freeze.py:126-143` 的"冻结≡live
  逐位"测试只跑 2.5-D 网格（`M0_DIR = naca0012_2.5d`，本文件 `:48`），那里
  side ≡ main，对该缺陷空转。
- **与 GB20.7 的关系**：stall 案例的记录是"freeze armed every level with zero
  reverts"——N1 给出了一个具体、已实测存在的机制候选（为什么 freeze 没救场）。
  **与 GB20.7 的因果仍未证明**；建议的判别实验：给 `freeze_side_state` 补一行
  `_apply_main_density`（与 `newton_side_data` 对齐），重跑 M6 medium ramp——
  恢复到 M0.84 则 N1 即机制；仍 stall 则回到 Claude 的 recipe 重标定路线
  （freeze_tol / dm）。修复后还应把探针的断言固化进测试（3-D 网格上
  patched≡shipped 选择一致）。

**N2 [测试缺口]** B20 的两处补丁（`readvec[mp]=el[mp]`、grad/q2 patch）没有
ungated 结构锁，唯一覆盖在 gated FD 测试里。建议加 ungated 断言（便宜）。

**N3 [流程缺口，B20 自己记录]** 3-D LS 数字不受任何测试锁定：B15 的测试全在
2.5-D，M6 medium γ = 0.088338 是 demo 数字——re-baseline 移动它们时套件不报警。
建议给核心 3-D LS 数字加锁（哪怕 gated）：M6 coarse/medium ramp 的 m_final、γ、
M_max、clamp 计数。

**N4 [陈旧注释 5 处]** `kernels/cut_assembly.py:417-421` 与 `:444-445`、
`wake/multivalued.py:561-571`（pre-B20 叙事当现状写）、
`tests/test_b19_jacobian_3d.py:5-9`（pre-erratum 措辞，与文件内第二个测试
自相矛盾）。另见 §2 D9/D10 的 docs 面同类问题。

---

## 2. 文档正确性

### 2.1 核对为一致的数字链（抽查到 CSV 行，约 40 个数字）

- **B17 三角**：coarse 0.2089/0.1853/0.2086/0.2087、medium
  0.2173/0.2165/0.1691/0.2117、GB17.3 双求解器 0.1691/0.1690、GB17.6 vortex
  +2.5% — `triangle_coarse.csv`、`triangle_medium.csv`、`lift_ab_medium.csv` ✓
- **B18 conforming**：cl_p 0.2617@M0.84 coarse、0.2579@M0.79 medium、cl(M) 升
  0.2173/0.2321/0.2579 — `cl_vs_mach.csv`、`checks.csv` ✓（B20 不触 conforming）
- **B19 全套**：targeted 1.145684e-01→1.333699e-08、control 6.33e-10、eps 翻转
  spread 1.00→131.5、partial 1.4697e-02、GB19.4 A/B（γ 0.07212068、M_max
  1.134235、+3.6% 墙钟）、GB19.6（252 元素、Δρ 0.4474/45.3%、q² 3.2229 vs
  1.3379）、blind-spot erratum 129/0 — `c1_fd_probes*.csv`、`b19_*.csv` ✓
- **GB20.1–20.6 与 re-baseline 表全部**：B7 M_max 1.453→1.392、B16 legacy
  limited 3690→11、pin floored 3→0、ramp 0.7875→M0.84 converged、M6 medium γ
  0.088338→0.071909、B15 17/20（3 FAIL 行已亲自核对：m_final=0.6625、2/5
  levels、M_max 1.5822 vs 2.4549）、B14 5/7、B9 头条未动一位 ✓
- **套件计数**：465+22+2 由本轮独立复跑确认 ✓

**结论：B19/B20 新写的内容在数字层面全部准确。** 问题集中在旧相位的回顾章节。

### 2.2 发现的问题（按严重度）

**D1 [系统性 stale] B18 旧故事在 6 个面未勘误。** re-baseline 重新生成的
`b18.../checks.csv` 已带更正注记（GB18.4 行："ATTRIBUTION CORRECTED by
B20/GB20.5: it is NOT the mixed_plain aux artifact"，已亲自核对），但
`track_b.md:1205-1212`、`agent-rules.md:129-153`、`overview.md:92-102`、
`demo_report.md:68`、`demo_report/track_b.md:1760-1784`、`design_track_b.md:1485`
仍断言旧归因（"G1.6/GB9.4/**B8 mixed-plain**"、Mmax 伪影 3.96、nlim 43/nflr 40）。
两个未披露的事实移动：**coarse LS 腿回退**（0.575/Mmax 1.44 → 0.55/Mmax 1.31，
`m_last_conv=nan`）和 **M0.6 cross-model 点已存在**（`cross_model.csv` 第 3 行
0.2178 vs 0.2174 = 0.2%，已亲自核对）而所有文档仍说 "SKIPPED"。
注意：**B18 的结论本身仍然成立**（LS 过不了 M0.5、medium 无共同跨声速 Mach）——
要勘误的是机制归因和细节数字。

**D2 [stale] B15 GB15.4 成功故事 5 处未勘误；GB20.7 披露面不足。**
`track_b.md:878-895`、`:8` header、`:1015-1018`/`:1716-1717`（B14/GB16.2 锚点
仍引 γ 0.088338 / M_max 2.4938）、`demo_report/track_b.md:1305-1326`、
`demo_report.md:61`（"19/19" → 现状 17/20）仍是成功叙事。GB20.7 只在
`track_b.md` + `agent-rules.md` 披露；**overview.md、PROJECT_STRUCTURE.md、
demo_report.md、demo_report/track_b.md、design_track_b.md 五面缺席**。

**D3 [doc error] ledger 自相矛盾（track_b.md 内部）。**
`:1684` track-status 行："B20 ✓ CLOSED … **NOT adopted, user's call**"、"B19 …
NOT adopted"——与同文件 `:1480-1487` / `:1701` 的 "PERMANENTLY, with no
switch (user-arbitrated)" 直接矛盾。`:1245`/`:1460` B19/B20 section 头仍
"◐ OPEN"；`:8` 顶部 banner 止于 B18；GB15.3/GB19.2/GB19.3/GB19.4 checkbox
未勾（GB19.4 头还写 "A/B in progress"，其正文 `:1360` 已写 "A/B COMPLETE"）；
`:1569-1572` GB20.6 正文仍是 pre-permanence 措辞（"the knob is inert until
someone asks for it"）。

**D4 [doc error，A3 遗留] M2 ledger 仍 ◐** — 见 §1.1 末段。

**D5 [doc error，pre-existing] B16 cond1**：文档（`track_b.md:984-985,1010`、
overview、agent-rules）写 **6.36e18**，提交 CSV（两个 epoch）为 **9.1e18**。

**D6 [doc error] re-baseline 表错标**：`track_b.md:1640` "B16 medium pin clamps
42/40 churn → 0/0"——B16 自己的 pre-B20 committed CSV 显示其 medium pin 为
0 lim/0 flr；42/40 属 B17 pin_gamma medium / B18 side 腿。`agent-rules.md:42`
同错。

**D7 [stale] B17 0.2115 → 0.2114**（re-baseline 后）六处未更新
（`track_b.md:1063,1139,1163,1719,:8`、`agent-rules.md:111`、`overview.md:132`、
`demo_report.md:67`、`demo_report/track_b.md:1681,1698`）；且 B17 的"junction
churn 只限残差地板"框架已过时（post-B20 同一轨迹收敛到 res ~1e-13）。
披露点 `track_b.md:1630` 还把 0.2114 列在 "*Unchanged:*" 条目下，自相矛盾。

**D8 [stale] B16 多个数字位移未披露**：legacy churn res 7.95→5.25、pin res
5.88e-14→1.49e-13、outer jumps 5.3e-15→1.8e-15、max|R| 84.457→84.8（文档说
"reproduced to the digit"）、medium pin res 7.03e-6→1.79e-13、neumann
unbounded 2.6e8→5.8e5。B16 章节（`track_b.md:982-1032` 等 6 面）全部保留旧值；
artifact 内部也有 stale（GB16.4 行仍 "OPEN/XFAIL" 已被 B17  supersede；GB16.3
medium 行值注自相矛盾）。

**D9 [stale] 索引与规范文件缺 B19/B20**：`demo_report.md` 无 B19/B20 行（违反
其自身规则）；`demo_report/track_b.md:8` header 缺 B20；**`design_track_b.md`
止于 §18（B18）——B20 永久改变了离散（mixed-plain 密度来源），数值规范文件
不再描述 shipped 方案**；`PROJECT_STRUCTURE.md:961` "A3 ◐"、`:493`
"test_b1..test_b18"、`cases/analysis/` 不在目录树；`track_a.md:8` "no A3
scoped" 与同文件 "A3 ✓ CLOSED" 并存；`roadmap.md:14` B 行漏 B16/B17；
`overview.md:14-16`（knob 存在）vs `:37`（knob 移除）同文件矛盾、GB20.7 缺席、
track 表 B 行止于 B18、回归基线仍是 B19 归因（计数仍对）；`agent-rules.md:5`
knob 与 `:30` "knob is REMOVED" 同文件矛盾。

**D10 [pre-existing] 三处更老的 doc 债**：`design.md` §12 risk 1 仍写"TE nodes
not duplicated"（P2 早已重定规格为 duplicated，单值 TE 是"定量错误"）；P13
动机段仍说 tip 修复"交给 Track B / B10 rescope"（G13.2 已记 superseded，B10
条目无回指）；`docs/analysis/capability_review_2026-07-15.md` 无 superseded
banner（其 LS medium M0.84"已演示"恰是 GB20.7 现状 FAIL 的案例）。

**D11 [trivial]** 时序漂移：GB15.3 文档 41.9/7.5/5.6×/6.5 vs 现 CSV
41.7/7.7/5.4×/6.6（物理数字逐位不变）；GB15.2 残差 8.5e-13 vs 1.16e-12；
GB17.2 "~22%" vs CSV 20.1%/19.4%（参照基不同）。

**流程建议**：close-out 五面清单管得住"新章节"，管不住"re-baseline 后旧章节的
勘误"。建议加一条：**任何 re-baseline 提交必须附"引用旧数字的章节清单"，逐一
勘误或标注**——本轮 D1/D2/D7/D8 全部是这一条缺失的产物。

---

## 3. 计划合理性评估

- **决策路径健康**。A3（先验证后修）→ B19（分两条腿；只修一半时如实记
  PARTIAL 而非四舍五入，逼出第二个缺项）→ B20（A/B 测量后用户裁决永久化，
  而非留一个"明知有缺陷的选项"）→ re-baseline（范围实测而非假设）。alias bug
  靠"测量异常而非解释异常"被抓出，没有被记成物理发现——纪律在发挥作用。
- **GB20.7 是当前唯一开放的 gate 级技术问题**，两个机制候选：
  (a) Claude 的 recipe-mismatch 假设（B15 recipe 针对污染场标定）；
  (b) 本轮 N1（冻结选择仍捕获于污染场——更具体，且直接解释了"freeze armed
  every level with zero reverts"为什么没救场）。判别实验便宜，见 §1.3。
- **后续排序建议**（最终裁决是用户的）：
  1. **N1 修复 + GB20.7 判别**（小、解锁 re-baseline 收尾）；
  2. **N3 测试锁**（小，防下一次无声 re-baseline）；
  3. **文档勘误一波**（D1–D9，小，但现在是讨论前的最佳时机）；
  4. **LS fine route**（中，GB20.7 之后才有意义）；
  5. **Track V vs P11**：上一轮我给的是"V 看物理收益、P11 看 gate 卫生"。
     **B20 改变了权衡**：交界 pocket 被归因为 G1.6 类刻面几何误差 ⇒ P11
     （曲面壁元）的潜在受益面从"sphere 一个 gate"扩大到 **sphere + 翼身
     fuselage Cp（GB9.4）+ LS 翼身跨声速**三个已知伤口。若翼身跨声速 LS
     是目标，P11 的优先级应上调。

---

## 4. 能力矩阵（截至 3db08b9，证据均为 committed artifacts）

图例：**V** = 对参考数据验证过；**D** = 已演示（committed 自检证据，无外部参考）；
**L** = 受限；**F** = 失败/阴性；**N** = 未尝试。

### 4.1 翼型 NACA0012（2.5-D，M0 wake-embedded + M3 wake-free 族）

| 域 | conforming | level-set |
|---|---|---|
| 不可压 | **V**：G2.3 cl_p 0.47858 vs Hess–Smith 0.482556 = −0.82%；Kutta 2 次更新收敛 | **V**：B3/B4 隐式 pressure-Kutta；same-mesh Γ 差 0.1–0.7%；wake-free 复现 0.3% |
| 亚声速 M0.5 | **V**：cl 0.284372 = −0.33%（PG/KT bracket 内）；**网格收敛** 2.71→0.33→0.03%（G9.2） | **V**：cl_KJ 0.2828 bracket 内；GA3.5 重锁 0.1441% |
| 跨声速 M0.74–0.82 | **D**：真离散解（G8.1/G4.3）；**medium M0.80 无孤立解**（FP fold，模型极限） | **D/L**：coarse gate MET（M3 wake-free Γ +0.9%）；medium fold Picard stall −18.8%；same-mesh γ 差 **−7.4% 未归因** |

### 4.2 机翼 ONERA M6（3-D，M1/M4/M5 族）

| 域 | conforming | level-set |
|---|---|---|
| 亚声速 M0.5 | **V**：repo 唯一 3-D Richardson（round-tip 自相似梯，p = 2.31，cl_KJ 外推 0.2050） | **D**：4 法 spread 0.322%（A1）；workflow cl 0.2129 |
| 跨声速 M0.84 | **D**：P14 pressure-Kutta 0.2776/0.2823，**跨模型 vs LS 0.17%/0.36%**；tapered fine cl_KJ 0.2866 = 参考的 99.5%（"strongly indicated, NOT earned"）；AGARD Cp 定性 overlay | **D/L**：coarse 到 M0.84 收敛（GB20.4，post-B20）；**medium 回退到 M0.6625（GB20.7 开放）** |
| fine（~2.5M tets） | **D**：M0.84 真离散解（M_max 2.818，0 clamps，44.6 min） | **N**：路线设计了未建（AMG O(n) + 薄带 LU） |

### 4.3 翼身组合体（3-D，M2 族）

| 域 | conforming | level-set |
|---|---|---|
| 亚声速 M0.5 | **D**（B9 新能力，0 lim/flr） | **D**：post-B20 medium res 1.1e-13、6 clamps |
| — 跨模型 M0.5 | colspan：**V** — GB9.5：0.2173/0.2188 vs 0.2165/0.2175 = **0.4%/0.6%**（B20 未动一位）；pin_gamma 后 LS 低 2.6%（远场截断） | |
| 跨声速 | **D**：coarse **M0.84** cl_p 0.2617；medium **M0.79** 0.2579 strict；M0.80+ stall（非 sliver，记录不追） | **F（closed-negative）**：交界 pocket 真实存在——M0.5 收敛（res 1.1e-13）但带 **Mmax 5.22** 尖峰；coarse ceiling 0.55；归因 = G1.6 类刻面几何，随加密恶化 |

### 4.4 验证形状与其它

- **Sphere**：不可压 Cp **F**（G1.6 strict xfail 11.6%，变分罪）；压缩性 **V**
  （M0.3 Cp peak 0.32%）；非升力 Newton 入口 G10.1 开放。
- **Cylinder 2.5-D**：**D**（Cp 误差 9.1→4.5→2.2%，slope 1.02）。
- **NLR7301 双元件**：**N**（B9 多元件腿 2026-07-17 superseded，未排期）。
- **Track V（粘性）**：设计完整（Drela IBL3 + transpiration BC），零实现。

### 4.5 开放问题清单（精确位置）

1. **GB20.7**（Track B 唯一开放 checkbox）：M6 medium M0.84 LS ramp 只到
   M0.6625；recipe mismatch vs 能力损失未裁决；本轮新增机制候选 N1。
   （`track_b.md` B20；committed FAIL 在 `b15.../checks.csv` ×3、`b14.../checks.csv` ×2）
2. **G1.6** sphere Cp（strict xfail 11.6%）：唯一在案路线 = Option C gate
   重定义 + P11 曲面壁元。受益面已扩到 GB9.4/交界 pocket。
3. **G13.3 跨声速 Richardson**：阴性开放（round fine 死于 M0.75，尖 tip TE）；
   "0.019 gap = 分辨率" 仍 *strongly indicated, NOT earned*。
4. **G10.1** 非升力 Newton 入口开放；**P11** 条件性未开；**P12** backlog；
   **B10** shelved（直尾迹 O(θ²) 记录）。
5. **翼身跨声速不对称**（B18 closed-negative，B20 修正归因后结论不变）。
6. **LS fine-mesh 路线缺口**（B14 只关 medium 瓶颈）。
7. **GB9.4** 机身带 16–20% 机翼升力且 LS 值随加密增长（XFAIL → G1.6）。
8. **B17 残差 2.6%**：pin_gamma medium 三角单调闭合但低于 conforming（远场
   截断）；vortex 从 +2.5% bracket，不闭合。
9. **B18 conforming medium M0.80+ stall**（res 2–7e-6，0 clamp，记录不追）。
10. **LS-vs-conforming 离散差**：NACA coarse M0.80 γ 低 −7.4%（已量化未归因，
    用户裁决不追）。
11. **M2 附带验证项**：交界最内 TE 节点 fan 含机身面——上下 CV 是否只取
    机翼侧元素（`track_m.md` M2）。
12. **N1**（本轮新增，已实测确认）：`freeze_side_state` 未打 B20 补丁；
    M6 coarse M0.70 下选择差异 83 upstream + 9 branch（探针：
    `docs/inspection/20260719-n1-freeze-probe.py`）。
13. **流程缺口**：3-D LS 数字无测试锁（N3）；re-baseline 旧章节勘误无流程
    （§2 流程建议）。

---

## 5. 给 Claude 的讨论要点（按优先级）

1. **N1**：`multivalued.py:739-745` 补 `_apply_main_density`（与
  `newton_side_data:667-668` 对齐），然后把 GB20.7 判别实验跑掉（M6 medium
  ramp）。若恢复 M0.84，GB20.7 关闭；若仍 stall，走 freeze_tol/dm 重标定。
2. **文档勘误一波**：D1（B18 六面）、D2（B15 五面 + GB20.7 披露）、D3
  （ledger 矛盾）、D4（M2 ledger）、D5（cond1）、D6（42/40 错标）、D7
  （0.2114）、D8（B16 位移）、D9（demo_report 行 / design_track_b §19-20 /
  PROJECT_STRUCTURE / overview / agent-rules 自相矛盾）。
3. **N3**：给 3-D LS 核心数字加测试锁。
4. **P11 优先级重估**（受益面扩大，见 §3 第 5 条）。
5. **close-out ritual 加 re-baseline 勘误条款**（§2 流程建议）。

---

## 6. 本轮未做

- 未重跑任何贵重 demo（P4 heavy ~40 min、P5 medium 45–75 min、B15 ramp
  ~22 min 等），以 committed CSV 为准；N1 与 GB20.7 的因果未实测（给出的是
  判别实验，不是结论）。
- 门控套件（`PYFP3D_TRANSONIC_GATES=1`）复跑结果见 §7。
- 未修改 `pyfp3d/`、`tests/` 或任何 docs（本报告除外）。

## 7. 补记（门控套件复跑）

（`PYFP3D_TRANSONIC_GATES=1` 全套件复跑在本报告撰写时进行中；结果出来后补记于此。）
