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
| **P — 求解器**（[roadmap/track_p.md](roadmap/track_p.md)） | P0–P9 ✓（P1 仅 G1.6 以 strict xfail 挂起）；P10 ◐（G10.2/G10.3 ✓）；P13 ◐（G13.1 ✓、G13.2 conforming ✓、G13.3 亚声速 Richardson ✓ p=2.31） | G10.1（非升力 Newton 入口，无顺序约束）；G13.3 **跨声速阴性开放**（圆帽 fine 的 ramp 死于 M=0.75，site=尖 tip TE）；P11 条件性未开（仅剩 G1.6 理由）；P12 backlog；**P14 ◐ 2026-07-17 开启并落地（用户指示）**：壁面邻接 CV 压力相等 Kutta 估计器（A2 路由的修复）。**G14.1–G14.6 ✓**——S1/S2 一次换掉：M0.84 Γ(z) roughness 0.0970→**0.0043**（coarse）/ 0.0365→**0.0024**（medium，均达/优于 LS 带），全站 raw TE Cp gap 0.2206→**0.0040** / 0.1585→**0.0024**（**55×/67×**，走 G14.6 **主条款**，回退未动用；★ 首版误引 A2 的 *section 末点* 数 0.318/0.228 报成 80×/95×，当日据用户提问勘误——A2 有两套 TE gap 度量，须引用自己实际跑的那套）。**★ 跨模型对照（V14.6，`cross_model_medium_m084.csv`）= 本相位最强证据**：LS 路径**一直**用压力相等 Kutta（B4），故若升力位移真是 Kutta *形式*，压力路径必须落到 LS 的答案上——它做到了：medium M0.84 conforming-pressure **0.2776/0.2823** vs level-set **0.2772/0.2813**，**差 0.17%/0.36%**（探针路径当年比 LS **低 4.5%/4.3%**），而且是不同尾迹模型、不同 DOF 空间、**不同网格族**（`onera_m6_wakefree`）。⇒ 长期存在的 conforming-vs-LS 升力分歧**就是** Kutta 形式误差，已消除。保留告诫：跨模型非同网格 A/B；LS 态带 1 lim/2 flr（B15 caveat）而压力态 0/0；"两者一致"≠"两者都对"（刚性平面尾迹等共有模型误差对两者是共模）。**★ G14.7 XFAIL 且不改带**：medium M0.84 cl_KJ 0.2823（+4.85% vs G8.2 **探针路径**锁）——机理在一级已测并在二级跑前**预注册**（两闭合逐点只差探针自身 O(h) 读数偏差，cross-read medium 0.79%；Kutta 映射 b≈0.93 把它 1/(1−b)≈14× 放大进 Γ）；**方向（记录，非 gate）**：\|cl_KJ−0.288\| 0.0188→0.0057，**P9 那 0.019 gap 关掉 69%**（P9 看不见它：两套网格共用同一估计器，偏差对其 Richardson 是共模）。**非**网格收敛断言（M6 fine 不是离散解）、**非**"0.019 是分辨率"翻案（仍 *strongly indicated, NOT earned*）、**非**"压力升力就是对的"。**待用户裁决**：认作发现（G14.7 改锚压力路径锁）还是当缺陷追。★ **自我更正（V14.7 实测 2026-07-17）**：先前各版都断言 TE Cp **spike** 不会被 P14 触及（"两路共有的 P1 恢复伪影"）——那是从 A2 归因**推**的，没测。实测 medium M0.84 raw：探针 **0.1143** → 压力 **0.0533**（2.1×），且**低于 LS 的 0.0743**。A2 对的部分：确有共有残余（~0.05 = 真正的恢复地板）；需修正的部分：conforming 相对 LS 的**超出量**也是 Kutta 形式误差（Kutta 错⇒TE 流场真的错，末点偏离趋势有物理原因，共模度量分不开）。旁证：P6 平滑在压力路径上不再有效（0.0533→0.0660→0.0626；A2 在探针路径上测的是 0.147→0.081）。**教训：别把上一相位的归因当结论带进新测量——去测。** 诚实记录：判别器 D=7.33→**1.80**（落在 A2 的 inconclusive 区，非 O(1)）|
| **M — 网格**（[roadmap/track_m.md](roadmap/track_m.md)） | M0、M1(+M1b 自相似阶梯)、M3、M4、M5（圆顶翼尖盖）✓ | M2 ◐：翼身网格 ✓（2026-07-13；**机身+远场 2026-07-16 按用户指示重定规格并重生成**：5 倍翼根弦长、机翼居中、2 倍直径椭球机鼻、蒙皮 h_body=2h_wall + 两端按半径加密；**R_FAR 15→25 MAC**，h_far 与所有固定加密距离同比放大 ⇒ 2.78× 域几乎不要钱；★需 `Mesh.OptimizeNetgen` 治尾迹 corridor 在对称面压出的细带 sliver **抽签**），求解腿 = B9 |
| **B — level-set 尾迹**（[roadmap/track_b.md](roadmap/track_b.md)） | B1–B5、B7、B8（characterized-not-cured）、B11–B14、B15 ✓；B6 ◐（coarse gate ✓；medium 定量项由 GB15.4 补上） | **B9 = NEXT**（翼身 LS，M∞0.5）；**B14 ✓ 2026-07-17**（`precond="schur"` Schur 消元 aux + AMG(SPD 主块)：A1 的 precond 瓶颈消失，M6 medium M0.84 42.6%→2.6%、ramp 1.43×/亚声速 2.08×，γ=已提交 GB15.4；★ 小规模反而更慢，fine 内存受限路线仍是未建的设计用例）；B10 搁置 |
| **V — 粘性耦合**（[roadmap/track_v.md](roadmap/track_v.md)） | 设计完整（Drela IBL3 + transpiration BC），零实现 | V1 依赖 P6（已满足），预算等同一个 Track-P 阶段 |
| **A — 校验与分析**（[roadmap/track_a.md](roadmap/track_a.md)） | 2026-07-15 新建；**A1 ✓ 2026-07-16**（GA1.1–GA1.5：四求解器统一计时插桩 + conforming×level-set × Picard×Newton 耗时基准） | **A2 ✓ 2026-07-17 关闭**（TE/Kutta 保真度归因，GA2.1–GA2.5）：**S1 定谳**——conforming Γ(z) 逐站抖动是逐站探针差势跳 Kutta target **估计器**的测量伪影（fixed-Γ 判别量 D=7.33/25.70 coarse/medium，把抖动从光滑场里重新生出来；闭合残差 ≤0.6% 排除"未闭合"、抖动局域于 TE 邻层 0.02–0.07× 排除"流场"），**非流场内容**；**S2 分解**——TE Cp 突跳=势跳 Kutta 形式误差（conforming 独有,同估计器 34×/133× vs LS）+ P1 末点恢复伪影（两路共有）；2.5-D `a1_cp.png`（有 S2 无 S1）证明两者是不同机制。修复路由至 **P14**（无 `pyfp3d/` 改动）；工作目录 `cases/analysis/`（区别于 `cases/demo/`）。★A1 结论：3-D 下两条 Newton 路径都是 **precond（LU 分解）受限**（~40% 墙钟，lagged LU 已开），2.5-D 的"seed 是成本"**不外推**——引用主导相位必须带网格 |

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

现基线 **421 passed + 18 skipped + 2 xfailed**（2026-07-17 P14 一级+二级：
+15 = `tests/test_p14_te_pressure.py`；实测 1015.17 s @8 线程）。
重 gate 走 `PYFP3D_TRANSONIC_GATES=1`；M6 `.msh` gitignored，16 条 M1 测试在
本地未生成网格时跳过（`cases/meshes/onera_m6/generate_onera_m6.py`，~30 s）。
内核/装配改动后先跑 `tests/test_v0_freestream.py`。

**基线演进**（近期；完整记录在各 track 台账）：182+8+2（P8）→ 184（P10 G10.2）→
218（B1）→ 229（B2）→ 276+17+2（B7，719 s）→ 291（P13/G13.2）→ 294（M1b）→
350（B8 关闭：+M5/B8/M2）→ 358（B8 backlog：M2 census 锁）→ 375+18+2（B11）→
384（B12+B13，此二者曾漏记，395=375+9+11 才对账）→ 395（B15）→ 396（B15 勘误）
→ **399+18+2（M2 机身+远场重定规格 2026-07-16：+3 = 比例规则锁 + 已提交 census CSV 锁
+ 远场净空锁；实测 973.59 s @8 线程）**。★ 该次实测报的是 **406**：当时工作区里还带着
Track A **尚未提交**的 7 个 A1 测试（`tests/test_a1_instrumentation.py`，单独实测 7 passed），
406 − 7 = 399 才是 M2 这次的账；A1 落地后基线即为 406（已兑现）
→ **421+18+2（P14，2026-07-17：+15 = `tests/test_p14_te_pressure.py`；实测
1015.17 s @8 线程；406 + 15 = 421 逐项对账，零回归）**。

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
  权重修复（B8 backlog 记录未排期）；B14 ✓（2026-07-17 建成，`precond="schur"`）——
  剩余未建 = fine 内存受限路线（AMG O(n) + 薄带 LU，无全尺寸 splu）。
