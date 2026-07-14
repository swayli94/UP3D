# pyFP3D 总览（快照 + 文档地图）

> **快照日期：2026-07-15。** 本文件是给人读的高层总览，**不是**权威进度源：
> 阶段/gate 状态以 [roadmap.md](roadmap.md)（track 索引）+ [roadmap/](roadmap/)
> 各 track 文件（含各自的进度台账）为准；当前阶段以 [agent-rules.md](agent-rules.md)
> 为准；证据在 [demo_report.md](demo_report.md) 索引 + [demo_report/](demo_report/)。
> 若本文件与它们冲突，以它们为准。gate 关闭时按 CLAUDE.md 工作流更新台账后，顺手刷新本文件。

## 一句话状态

四条 track：求解器主线（P）与网格线（M）基本关闭，level-set 尾迹线（B）是当前
工作面——**B15 已关（2026-07-15，LS Newton 跨声速 ramp + N5 冻结选择，M6 medium
M0.84 的 Picard plateau 消除、38.4→11.0 min）**，**下一步 = B9（翼身组合体 LS
求解，M∞0.5；M2 网格已交付）**。粘性耦合线（V）设计完毕、零实现。

## Track 状态表

| Track | 状态 | 开放项 / 下一步 |
|-------|------|----------------|
| **P — 求解器**（[roadmap/track_p.md](roadmap/track_p.md)） | P0–P9 ✓（P1 仅 G1.6 以 strict xfail 挂起）；P10 ◐（G10.2/G10.3 ✓）；P13 ◐（G13.1 ✓、G13.2 conforming ✓、G13.3 亚声速 Richardson ✓ p=2.31） | G10.1（非升力 Newton 入口，无顺序约束）；G13.3 **跨声速阴性开放**（圆帽 fine 的 ramp 死于 M=0.75，site=尖 tip TE）；P11 条件性未开（仅剩 G1.6 理由）；P12 backlog |
| **M — 网格**（[roadmap/track_m.md](roadmap/track_m.md)） | M0、M1(+M1b 自相似阶梯)、M3、M4、M5（圆顶翼尖盖）✓ | M2 ◐：翼身网格 ✓（2026-07-13），求解腿 = B9 |
| **B — level-set 尾迹**（[roadmap/track_b.md](roadmap/track_b.md)） | B1–B5、B7、B8（characterized-not-cured）、B11–B13、B15 ✓；B6 ◐（coarse gate ✓；medium 定量项由 GB15.4 补上） | **B9 = NEXT**（翼身 LS，M∞0.5）；B14 designed-not-scheduled（Schur+AMG，fine 触发）；B10 搁置 |
| **V — 粘性耦合**（[roadmap/track_v.md](roadmap/track_v.md)） | 设计完整（Drela IBL3 + transpiration BC），零实现 | V1 依赖 P6（已满足），预算等同一个 Track-P 阶段 |

## 文档地图（每份文档的职能、权威范围、何时更新）

| 文档 | 职能 | 权威范围 | 何时更新 |
|------|------|----------|----------|
| [roadmap.md](roadmap.md) | track 索引 + working rules + gate 编号约定 | 各 track 一行状态 | track 状态行变化时 |
| [roadmap/track_{p,m,b,v}.md](roadmap/) | 各 track 的 phase 条目 + gate 清单 + 进度台账 | **阶段/gate 状态的唯一权威** | gate 开/关时 |
| [design.md](design.md) | 理论与数值参考（方程、离散、内核规则、求解策略、V0–V6 验证阶梯） | 数值方法（conforming 路径 + 共享理论） | 方法/勘误变化时 |
| [design_track_b.md](design_track_b.md) | Track B（level-set 尾迹）数值方案 + 逐阶段技术结论 | Track B 数值 | Track B 阶段推进时 |
| [demo_report.md](demo_report.md) + [demo_report/](demo_report/) | 已关阶段的证据档案（每阶段一个自检 demo + 提交的图/CSV） | 证据；**无 committed 工件的断言不是证据** | 阶段关闭时加节 |
| [agent-rules.md](agent-rules.md) | 每 session 注入的当前阶段 + 操作纪律（经 CLAUDE.md `@include`） | 当前阶段行 | 阶段变化时 |
| overview.md（本文件） | 人读总览 + 文档地图 | 无（快照） | 顺手刷新 |
| [analysis/](analysis/) | 分析/审查类报告（capability review 等），非规范文档 | 无（报告注明快照日期） | 新报告放这里 |
| [archive/](archive/) | 历史归档（勿作规范；rule 11 同样适用） | 无 | 只进不改 |
| `docs/references/` | 外部文献（López dissertation PDF 等） | — | — |

已删除：`docs/STATUS.md`（2026-07-15，本文件取代）、`docs/discussion_notes/`
（2026-07-14，commit 0e4895a；历史经 `git show 8aa4aee:docs/discussion_notes/<file>`）。

## 回归基线

**395 passed + 18 skipped + 2 xfailed @ B15 关闭（2026-07-15，实测 968.29 s
@16 线程）**；同日勘误提交（8ccd9b2）又加 1 条非门控测试
（`test_freeze_max_clamped_relaxes_the_convergence_semantics`）⇒ 现基线
**396 passed + 18 skipped + 2 xfailed**（2026-07-15 全套件复测实测 988.73 s
@16 线程）。
重 gate 走 `PYFP3D_TRANSONIC_GATES=1`；M6 `.msh` gitignored，16 条 M1 测试在
本地未生成网格时跳过（`cases/meshes/onera_m6/generate_onera_m6.py`，~30 s）。
内核/装配改动后先跑 `tests/test_v0_freestream.py`。

**基线演进**（近期；完整记录在各 track 台账）：182+8+2（P8）→ 184（P10 G10.2）→
218（B1）→ 229（B2）→ 276+17+2（B7，719 s）→ 291（P13/G13.2）→ 294（M1b）→
350（B8 关闭：+M5/B8/M2）→ 358（B8 backlog：M2 census 锁）→ 375+18+2（B11）→
384（B12+B13，此二者曾漏记，395=375+9+11 才对账）→ 395（B15）→ 396（B15 勘误）。

## 长期挂起项（勿反复重提）

- **G1.6 球面 Cp <2%**（strict xfail，11.6%）：根因 = 平坦面片壁上的自然边界条件
  （变分罪）；h 加密、恢复调参、Nitsche、边界数据修正均已用证据排除。唯一在案
  路线 = Option C gate 重定义 + P11 曲面壁元。
- **V6 <1%**：已定性为随加密消失的真 O(h) 地板（6.30%→3.29%→**1.41%**，P13 实测），
  fine 上已逼近目标；非 P11 问题。
- **G13.3 跨声速 / M0.84 Richardson**：两种几何皆未挣得（平帽序列非渐近；圆帽无
  M0.84 fine 态——尖 tip TE 被圆帽放大）。"0.019 缺口=分辨率"= strongly
  indicated、NOT earned（2026-07-14 用户裁决措辞）。
- **LS-vs-conforming 离散差**：NACA coarse M0.80 上 γ 低 −7.4%（strict-to-strict，
  B15 首次可测；B6 曾记 ~13%）——已量化、未归因。
- **G13.1/G9.1 原文的 LS 指数（1.34/×2.28）已判为 ×5 度量伪影**（诚实 +0.62）；
  引旧文以 2026-07-14 勘误为准。**B15 未让 G9.1 复活**（其 conforming limited
  单元前提未重测）。
- 后掠 TE Kutta 探针跨站共享（P5 记录）；`element_densities` mixed-plain junk
  权重修复（B8 backlog 记录未排期）；B14 触发条件 = medium 仍慢或真 fine 战役。
