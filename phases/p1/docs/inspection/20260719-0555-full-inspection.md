# 全面检查报告（Kimi 第二轮独立审查）

- **Date**: 2026-07-19 05:55 CST（**2026-07-19 12:30 更新**：纳入 Claude 的 GB20.7
  结论提交 `37a0799`，分析见 §8）
- **Auditor**: Kimi Code CLI（独立于主 authoring agent）
- **Basis**: branch `kimi/inspection-b20-capability-review` @ `origin/main` =
  `37a0799`（GB20.7 结论；初稿基于 `3db08b9` B20 re-baseline）。覆盖 A3 响应
  落实、B16–B20 全部新工作、文档系统、规划系统、能力矩阵。
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
   落在 aux-touching mixed-plain 类上）。**（更新）Claude 的 GB20.7 扫掠
   （`37a0799`）只测试了 freeze 的"武装时机"轴（freeze_tol），未测试 N1 的
   "捕获内容"轴——N1 仍是开放且独立的机制候选，且无论 GB20.7 结论如何都
   是确定 correctness 修复**（§8）。
2. **基线独立复跑吻合**：`465 passed + 22 skipped + 2 xfailed`（1179.17 s
   @16 线程），与文档逐字一致。门控测试（`PYFP3D_TRANSONIC_GATES=1`，12 个
   门控文件）独立复跑 **85 passed 全绿**（§7）——同时实测确认 3-D LS 数字
   无测试锁的 N3 缺口在门控全开下依然存在。
3. **文档：B19/B20 新写的内容与提交 CSV 精确一致（抽查约 40 个数字全部对
   上）；但 B20 re-baseline 之后，B15–B18 的回顾性章节普遍未勘误**——新旧两套
   数字并存于同一文档集，GB20.7 这个开放项在 7 个文档面中只披露了 2 个。
   另有两处 pre-existing doc error（cond1 数值、design.md §12 陈旧 spec）和
   一处 A3 响应文档的夸大声明（M2 ledger 实际未修）。
4. **计划：A3→B19→B20 的决策路径健康、纪律有效**（alias bug 靠"测量异常而非
   解释异常"抓到，没有被记成物理发现）。**（更新）GB20.7 已由 Claude 给出
   第一轴答案："不是 freeze_tol-recipe 失配，是真实能力损失"（附诚实边界：
   只扫了 freeze_tol），并留下 GB15.4/GB14.4 重定规格给用户裁决；N1 修复
   是该结论收口前最后一个未做的便宜判别实验。**建议下一步顺序：N1 修复 +
   M6 medium ramp 复验 → GB15.4/GB14.4 重定规格（用户裁决）→ 3-D LS 数字
   测试锁 → 文档勘误一波 → LS fine route / Track V / P11（P11 的受益面在
   B20 后变大，见 §4）。

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
  恢复到 M0.84 则 N1 即机制；仍 stall 则"真实能力损失"结论排除最后一个
  已知代码侧嫌疑、正式收口。修复后还应把探针的断言固化进测试（3-D 网格上
  patched≡shipped 选择一致）。**（2026-07-19 更新）Claude 已在 freeze_tol 轴
  定谳 GB20.7（§8），但该扫掠在原理上测不到 N1——它仍是开放的独立机制
  候选。**

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
**（2026-07-19 更新）`37a0799` 把 GB20.7 条目改为已定谳并声明 GB15.4 条款
转 NEGATIVE / GB14.4 superseded，但 B15 章节正文的成功叙事仍未勘误**——
现在同一文件里"GB15.4 ✓ 达到 M0.84"（:878-895）与"GB15.4 的该条款现在是
NEGATIVE"（GB20.7 条目）直接并存，矛盾从"披露不足"升级为"同文件自相矛盾"。

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
- **（更新）GB20.7 第一轴已有答案**（`37a0799`，2026-07-19）：Claude 扫
  freeze_tol 1e-3→1e-6，ceiling 只动 0.6625→0.6750 ⇒ "**不是
  freeze_tol-recipe 失配，是真实能力损失**"，附诚实边界（只扫了 freeze_tol；
  dm/n_newton_max/m_start 未动），并把 GB15.4/GB14.4 的重定规格留给用户
  裁决——方法论核对通过（忠实复现 B15 committed 调用、每变体新建 operator、
  不踩 demo-cache；CSV 与文档数字一致，§8）。**但第二个机制轴仍未测试：
  本轮 N1（freeze 捕获内容本身在未补丁场上计算）——它解释同一组症状，
  且 freeze_tol 扫掠在原理上看不见它**（扫的是"何时武装"，N1 是"武装时
  捕获了什么"）。判别实验便宜，见 §1.3 与 §8。
- **后续排序建议**（最终裁决是用户的）：
  1. **N1 修复 + M6 medium ramp 复验**（小；无论 GB20.7 结论如何都是确定性
     correctness 修复，且是"真实能力损失"结论收口前最后一个未做的便宜
     判别实验）；
  2. **GB15.4/GB14.4 重定规格**（用户裁决；Claude 已按 G14.7 先例提议
     re-lock 到新包线，留待用户）；
  3. **N3 测试锁**（小，防下一次无声 re-baseline）；
  4. **文档勘误一波**（D1–D9，小，但现在是讨论前的最佳时机）；
  5. **LS fine route**（中，GB20.7 收口之后才有意义）；
  6. **Track V vs P11**：上一轮我给的是"V 看物理收益、P11 看 gate 卫生"。
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
| 跨声速 M0.84 | **D**：P14 pressure-Kutta 0.2776/0.2823，**跨模型 vs LS 0.17%/0.36%**；tapered fine cl_KJ 0.2866 = 参考的 99.5%（"strongly indicated, NOT earned"）；AGARD Cp 定性 overlay | **D/L**：coarse 到 M0.84 收敛（GB20.4，post-B20）；**medium 包线回退到 ≈M0.675**（committed recipe 0.6625；GB20.7 已定谳"真实能力损失（freeze_tol 轴）"，GB15.4/GB14.4 重定规格待用户裁决，N1 轴未测） |
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

1. **GB20.7（更新：第一轴已定谳，第二轴开放）**：M6 medium M0.84 LS ramp
   只到 M0.6625（committed recipe；freeze_tol 扫掠最佳 0.675）。Claude
   `37a0799` 定谳"**真实能力损失，非 freeze_tol-recipe 失配**"（附边界：
   只扫了 freeze_tol），GB15.4 "reaches M0.84" 转 NEGATIVE、GB14.4 条款
   superseded，**重定规格留待用户裁决**。仍未测试的第二机制轴 = N1
   （§8）。
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

1. **N1（仍是第一优先级，GB20.7 定谳后尤甚）**：`multivalued.py:739-745` 补
  `_apply_main_density`（与 `newton_side_data:667-668` 对齐）。两点理由：
  (a) 不变式破坏已实测，是确定性 correctness 修复，与 GB20.7 结论无关；
  (b) Claude 的 freeze_tol 扫掠在原理上测不到它（§8.3）——补完后跑一次
  committed M6-medium 调用，若恢复 M0.84 则"真实能力损失"结论翻案，若仍
  stall 则该结论排除最后一个已知的代码侧嫌疑、正式收口。
2. **文档勘误一波**：D1（B18 六面）、D2（B15 五面 + GB20.7 披露）、D3
  （ledger 矛盾）、D4（M2 ledger）、D5（cond1）、D6（42/40 错标）、D7
  （0.2114）、D8（B16 位移）、D9（demo_report 行 / design_track_b §19-20 /
  PROJECT_STRUCTURE / overview / agent-rules 自相矛盾）。
3. **N3**：给 3-D LS 核心数字加测试锁。
4. **P11 优先级重估**（受益面扩大，见 §3 第 6 条）。
5. **close-out ritual 加 re-baseline 勘误条款**（§2 流程建议）。

---

## 6. 本轮未做

- 未重跑任何贵重 demo（P4 heavy ~40 min、P5 medium 45–75 min、B15 ramp
  ~22 min 等），以 committed CSV 为准；N1 与 GB20.7 的因果未实测（给出的
  是判别实验，不是结论）。**（更新）Claude 的 GB20.7 扫掠证据已由我核对
  （§8），其结论同样基于 committed CSV，我未复跑那 4 条 ramp（每条
  ~22–31 min）。**
- 门控套件（`PYFP3D_TRANSONIC_GATES=1`）复跑结果见 §7。
- 未修改 `pyfp3d/`、`tests/` 或任何 docs（本报告除外）。

## 7. 补记（门控套件复跑，2026-07-19 09:40 CST）

`PYFP3D_TRANSONIC_GATES=1` 下，对**全部 12 个含门控测试的文件**独立复跑：
`test_p4_transonic / p5_onera_m6 / p8_newton / p8_jacobian / b6_newton /
b6_transonic / b11_linear_ls / b14_schur_ls / b16_farfield_aux /
b17_farfield_pin_gamma / b18_wingbody_transonic / b19_jacobian_3d` ——
**85 passed in 4910.55 s（1:21:50），全绿**。这独立证实了 B20 re-baseline
提交声明的"no test lock breaks at all"：包括 b14 的 medium A/B 门控在内，
没有任何测试锁定被 re-baseline 打破。

同时这也实测确认了 **N3 缺口本身**：b14 门控 A/B 比较的是同一代码下的两条
求解路径（都停在 M0.6625，彼此一致），而 B15/B14 demo 的 17/20、5/7 FAIL
只存在于 committed CSV——**套件在门控全开的最大配置下依然看不见这次 3-D LS
数字回退**，3-D LS 数字无测试锁的缺口是真实且当前存在的。

（附注：首次全套件门控复跑在 2 小时上限被杀——与 Claude 在 UP3D 的自身
测试存在 CPU 争用；改为只跑门控文件后完成。Claude 声称的"67/67 gated 3-D
LS 测试绿"与本结果相容：本运行的 85 含各文件 ungated 测试。）

---

## 8. GB20.7 更新分析（2026-07-19 补，commit `37a0799`）

Claude 在我报告初稿之后提交了 GB20.7 的第一轴结论。本节是对该提交的
独立核对与评估。

### 8.1 核对：证据与文档一致，方法论严谨

- **数字一致**：`cases/analysis/c1_ls_jacobian_fd/results/gb207_recipe_sweep.csv`
  四行（freeze_tol 1e-3/1e-4 → m_final 0.6625、2/5、assignment_cycle；
  1e-5/1e-6 → 0.675、3/6、+tol；wall 1348–1886 s）与 commit message、
  `track_b.md` GB20.7 条目、`agent-rules.md` 的引用**逐位一致**。
- **方法严谨**（`run_gb207_recipe.py`）：忠实复现 B15 committed M6-medium
  调用（wakefree medium、`farfield="neumann"`、`n_seed=40`、`n_newton_max=80`、
  `tol_residual=1e-10` + `B_NEWTON_M6_DEFAULTS`，只动 freeze_tol）；每变体
  新建 operator（避免 state-dependent cache 串扰）；直接调求解器，不经过
  demo 的 npz cache——正确避开了它自己记录的 demo-cache 陷阱；判定规则
  预先写死（任一档到 M0.84 ⇒ mismatch；全部不到 ⇒ real loss，并显式标注
  "NOT proof: only freeze_tol was varied"）。

### 8.2 结论中成立的部分

- "**freeze_tol 是 contributor 不是 cause**"——由数据支持：降 freeze_tol
  动了 ceiling（0.6625→0.675）且出现了按 `tol` 收敛的级别，但不到 M0.84。
- "**真实能力损失**（在 freeze_tol 轴上）"的表述附了诚实边界（dm、
  n_newton_max、m_start 未扫）——准确。
- "**contamination 是无意的稳定器**"（与 GB20.5 翼身同型：收敛更干净、
  爬升更近）——与 committed 证据一致。
- 旧 M0.84 态的 M_max 2.45–2.49 只对同样被污染的 LS Picard 验证过
  （common-mode），conforming 参考 1.995——与 `3db08b9` re-baseline 记录
  一致，"trade"的表述（旧的走得更远但状态可疑；新的停得早但干净）成立。
- GB15.4/GB14.4 的重定规格留给用户而非单方面 re-spec——符合项目规则
  （重新定义 committed capability claim 是用户裁决），G14.7 先例引用恰当。

### 8.3 结论未覆盖的部分：N1 轴仍然开放（关键）

- 扫掠变的是"freeze **何时**武装"（freeze_tol）；N1 是"武装时**捕获了
  什么**"：只要 freeze 武装，capture 就在未补丁场上计算
  （`multivalued.py:739-745`），**与 freeze_tol 取值无关**。⇒ freeze_tol
  扫掠在原理上不能判别 N1。
- 值得注意的呼应：Claude 的假设措辞是"arming the freeze **before the
  selection settles** and locking a bad assignment"；N1 说的是"**即使
  selection 已 settled，捕获到的也不是 live 系统的那一个**"（实测 83+9
  个元素）。两者是同一不变式的两个独立破坏轴——时机 vs 内容。时机轴
  已测（contributor）；内容轴未测。
- 因此 GB20.7 的准确当前状态 = "**非 freeze_tol 失配；真实能力损失**"
  成立，但带两个显式 qualifier：(a) 其它 recipe 旋钮未扫（Claude 自己
  声明）；(b) **N1 capture bug 未修未测**（本轮发现）。N1 修复是收口前
  最后一个便宜的判别实验，且无论结果如何都是必须做的 correctness 修复。
- 附带：N1 修复后，B15 recipe 的校准前提（churn floor 在污染场上测得）
  整体需要重测——freeze 捕获修正 + freeze_tol 重标定是同一个 follow-up
  phase 的自然内容，不是障碍。

### 8.4 文档债的连带变化

- GB20.7 答案仍只披露于 `track_b.md` + `agent-rules.md`（2/7 面）；
  overview、PROJECT_STRUCTURE、demo_report×2、design_track_b 依然缺席
  （D2/D9 延续）。
- B15 章节正文的成功叙事未动 ⇒ `track_b.md` 同文件内"GB15.4 ✓ 达到
  M0.84"（:878-895）与"GB15.4 该条款现在是 NEGATIVE"（GB20.7 条目）
  **直接并存**——D2 从"披露不足"升级为"同文件自相矛盾"（已在 §2.2 D2
  更新）。
- 新增用户待裁决项：**GB15.4/GB14.4 按新包线重定规格**（已加入 §4.5
  与 §5）。
