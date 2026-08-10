# 全面审查报告：文档/代码一致性、翼身跨声速能力、Track V 可行性

> **日期**：2026-07-20（第三轮 Kimi 独立审查；前两轮 = 20260717 三份 + 20260719 full-inspection）
> **分支**：`kimi/b25-inboard-fragment-clip` @ `a537324`（B23–B27 链，领先 origin/main 18 个提交）
> **方法**：通读全部五个 track 文件 + design.md + design_track_b.md + overview/agent-rules/demo_report
> + 已删除的 DN2/DN6（git 历史恢复）；两个只读子代理分别核实翼身证据（committed CSV **及 npz 缓存
> 逐位核对**，b1+b2+m2+v0+b18 = 90 条测试实跑全绿）与 Track V 代码审计（全仓 grep + 装配层结构核对）。
> **性质**：横向审查报告，非 gate 文档，不改变任何 roadmap 状态；最终裁决 = 用户。

---

## 0. 头条结论

1. **翼身跨声速的全部关键声称经 npz/CSV 逐位核实成立**（conforming 0.84/0.79、
   LS+clip 0.84/0.7625、跨模型 2.4–2.6%、`inboard_clip` 默认逐位不变）。未发现数值正确性问题。
2. **主要问题是文档债**：B23–B27 的"五面收尾"缺四面——权威源 `docs/roadmap/track_b.md`
   完全没有 B23–B27 条目/台账行，且其 B18 条目仍陈述已退役叙事（无勘误指针），与已更新的
   索引面（roadmap.md / overview.md / agent-rules.md / design_track_b.md §22 /
   demo_report.md）**直接矛盾**。这是 B22 刚立的纪律 #10/#11 的同类失守；
   B27 VERDICT §5 自认"留待合入收尾"。**本报告落盘后即执行统一收尾。**
3. **Track V 可行**：设计绑定清晰、前置（P6/P8）全关闭、solver 侧约 70% 钩子已在位、
   参考数据齐备；最大工作量 = `viscous/` 包本体（Track-P 规模），最大物理风险 = IBL
   边缘速度输入的 P1 精度天花板与 TE/交界区已知污染。**V1（2.5-D 阶梯）现在就可开**；
   V3/V4（翼身 VII）建议排在 (b) 类天花板归因与 cl_fus 归宿之后。

---

## 1. 文档与代码的错误、矛盾、未解决问题

### 1.1 实质性矛盾（close-out 债，权威源 vs 索引面）

| 面 | 状态 | 问题 |
|---|---|---|
| `docs/roadmap/track_b.md`（**唯一权威**） | ✗ | 全文无 B23–B27 条目与台账行；B18 条目（`:1230`）仍写"LS junction-limited / G1.6 刻面 / closed-negative"，无 B26/B27 勘误指针——与 `roadmap.md:14`、`overview.md:190`、`agent-rules.md:549`、`design_track_b.md:1563`（§22）、`demo_report.md:69` 冲突 |
| `docs/demo_report/track_b.md` | ✗ | 头部（`:8`）停在 B21（连 B22 都没进头部）；B18 详章（`:1756-1819`）仍是旧叙事无加注；无 B26/B27 节 |
| `PROJECT_STRUCTURE.md` | ✗ | footer 写"Track B — B1–B9, B11–B22 ✓"；目录树无 b23–b27 分析目录 |
| `cases/analysis/README.md` | ✗ | 有 b23/b24/b25 行，缺 b26/b27 行 |
| `cases/demo/README.md` | ✗ | b18 行仍写 "1 PASS + 6 RECORDED"；B27 刷新后 = 8 gates 8/8 PASS |
| `docs/agent-rules.md` 基线行 | ✗ | 仍写 473+25+2（2026-07-19 P11）；B24/B25 向既有测试文件新增约 6–8 条测试（`test_b1_cut_elements.py:519-609` `TestInboardFragmentClip` 4 条；`test_m2_wingbody.py:152-162` 水线延伸；`test_b1_cut_elements.py:149-160` corner-theft 回归），真实基线已变 |
| `docs/overview.md` "现基线"行 | ✗ | 写 465+25+2（B22），比 P11 的 473 已旧一档；其台账末行才是 473 |

### 1.2 陈旧数字 / 失效路由

- **`docs/roadmap/track_v.md:62-67,84-85` scope guard**："M6 CL 0.245 vs FP-reference 0.288，
  归因于 sharp-TE/LE P1 壁面梯度地板 → P11，P9 先判别"——三处全过时：① 0.245 已被 P14
  （medium 0.2823）与 P13 tapered fine（0.2866 = 参考的 99.5%）超越；② P11 曲面壁元
  2026-07-19 实测**阴性**关闭；③ P14/G14.7 测得 0.019 gap 的 69% 是 Kutta 估计器偏差。
  定性结论（VII 把 CL 往下带、不把无粘缺口记到 Track V）仍成立且更强，但引用数字与路由失效。
- **DN6 §8.3** 声称 "kernels/residual.py 中壁面项的位置已经有占位（g=0 时为零）"——代码中
  不存在：`kernels/residual.py:43-76,144-189` 纯体积装配，无任何边界面循环。DN6 是历史文档，
  但 track_v.md 把 DN2/DN6 列为设计记录，开工时会误导（正确起点 = `solve/wall_correction.py`
  的 RHS 装配模板，见 §3.2）。
- `docs/inspection/20260719-0555-full-inspection.md` §4.5 开放清单中 #1（GB20.7）、#5（翼身
  跨声速不对称）、#12（N1）、#13（N3）已被 B21/B22/B26 解决（快照报告，不补注，读者注意）。

### 1.3 代码层小问题（非 bug，记录）

- `beyond_tip_elems` 语义拓宽：`inboard_clip` 生效时该集合同时收集**内侧被拒**单元
  （`wake/cut_elements.py:206-212`，仅注释说明），名字仍叫 beyond_tip；下游若按"翼尖外"
  解读会误读。现有测试用默认 clip，未受影响。
- B18 demo `checks.csv` GB18.4 行文压缩：`dist_fus=0.005` 属次级峰 M3.53，主峰 M6.17 的
  dist_fus=0.076（`g1_peaks.csv`）——不算错但易误读。
- `solve_multivalued_newton_transonic` 的 `target_reached` 在 `upwind_c_post` 之前计算
  （`newton_ls.py:1144-1146`），post 级失败不回收该标志；conforming `timings` 只覆盖末级
  （已标注 footgun），总量用 `timings_total`。
- pyfp3d/ 内无 TODO/FIXME；B25 默认逐位不变由 `test_default_none_is_bit_identical` 锁定。

### 1.4 未解决问题清单（精确现状）

| 问题 | 状态 |
|---|---|
| G1.6 球面 Cp 11.6% | strict xfail；P11 阴性 + 重归因（P1 场固有 max-norm 能力 @ h=0.08）；路线三岔口（Option C 重定规格 / isoparametric P2 壁层 / 接受）**未裁决** |
| **cl_fus 机身虚假升力**（GB9.4，16–20% wing cl_p；B26：C 侧 out-band ×2） | **当前无在案路线**——原指向的 P11/曲面壁元已阴性关闭；"G1.6 类"标签被 P11 摘除球锚；B23 给了判别但归宿悬空（B23 VERDICT §(c) 建议把 GB9.4 阈值重设为"袋带外跨模型一致"，未执行） |
| conforming medium M0.80+ stall 与 LS+clip medium 0.775 死 (b) 类 | B26 证据指向同机制（翼尖 P13 奇点 + 高 M Newton lim/flr 振荡）；**(b) 类归因 = 命名的下一 phase 候选，未裁决** |
| G13.3 跨声速 Richardson | 阴性开放；"0.019 gap = 分辨率" 仍 *strongly indicated, NOT earned* |
| 翼尖-TE 奇点（P13） | conforming 有 tip_taper（G13.2 ✓）；**LS 路径未愈**（B8 characterized-not-cured） |
| LS fine 路线（AMG O(n) + 薄带 LU） | B14 设计了未建 |
| B17 残差 2.6%（远场截断） | 开放；vortex 从 +2.5% bracket 不闭合 |
| LS-vs-conforming γ −7.4%（NACA coarse M0.80） | 已量化未归因，用户裁决不追 |
| M2 遗留：交界最内 TE 节点 CV fan 含机身面 | `track_m.md` M2 台账仍记 open |
| G10.1 非升力 Newton 入口 / B10 curved wake / P12 backlog | 开放 / shelved / backlog |
| 分支状态 | B23–B27 在未 push 链上（kimi/wingbody-junction-discriminator → b24 → b25 三分支叠放）；push/PR 需用户确认（AGENTS.md 纪律） |

---

## 2. 翼身跨声速能力评估

### 2.1 能力矩阵（全部经 npz 逐位核实）

| | conforming | level-set（默认 q≥0 裁切） | level-set + `inboard_clip`（B25） |
|---|---|---|---|
| 亚声速 M0.5 | ✓ medium 0.2173（B9 跨模型基准 0.4%/0.6%） | ✓ 0.2117（B17 pin_gamma） | 同左 |
| 跨声速 coarse | **M0.84 reached**，cl_p 0.2617，M_max 2.15，res 2.8e-12（proof-of-concept，欠解析，网格含 27 slivers） | 死 0.84（交界条带袋 M3.28 q≈0 卡 Newton，(b) 类） | **M0.84 reached**，cl_p 0.2542 |
| 跨声速 medium | **M0.79 严格**，cl_p 0.2579，res 2.2e-14；cl(M) 单调 0.2173/0.2321/0.2483/0.2579 | 死 0.5125（(a) 类：交界袋 Mmax 6.17，突破 freeze_max_clamped=8） | **M0.7625**（cl_p 0.2475），死 0.775 (b) 类：峰在**翼尖** M4.18（P13 类），交界走廊 corrM ≤1.10 干净 |
| 跨模型一致性 | colspan：M0.5 2.6%（B9/B17）· **M0.65 medium 2.4% PASS** · M0.75 2.5% RECORDED · coarse M0.6 2.1%（C 侧，欠解析）——gap 全 Mach 平在 ~2.5% = B17 cl_p↔cl_kj 口径差带 | | |

证据：`cases/demo/b18_wingbody_transonic/results/checks.csv` 8/8（含 3 条 RECORDED 语义性
PASS）、`cases/analysis/b26_ls_transonic_ceiling/results/g1_*.csv`、
`cases/analysis/b27_b18_demo_refresh/results/g27_consistency.csv`（336/336 逐位）。
配方：conforming 翼身 ramp 需 `freeze_tol=1e-5`（`run_demo.py:156`，coarse/medium 统一）；
LS 需 `solve_multivalued_newton_transonic` + `CutElementMap(inboard_clip=make_inboard_clip(FUS))`
（`wake/cut_elements.py:122,190-204` + `meshgen/fuselage.py:155-193`；默认 None = 逐位不变，
测试锁）。

### 2.2 结论性判断

1. **B18 的"能力不对称"叙事已被 B23–B27 推翻**：交界袋的根因不是 G1.6 刻面几何，而是尾迹片
   **inboard 自由边奇点**（默认 `q>=0` 裁切把片停在交界站位，自由边悬在流体内部，违反
   Helmholtz；B23 D5 定位：med |s|=0.0058 贴尾流面、med q=−0.003 正是 polyline 内端、
   折痕命中 0.00）。`inboard_clip` 把片裁成 conforming fragment 拓扑（贴机身水线→尾锥→
   对称面）后袋愈（走廊 corrM 14.66→0.63），**LS 天花板与 conforming 同址**：
   coarse 0.84=0.84；medium 0.7625≈0.79。
2. **残余限制器两路径同类**：(b) 类 = 翼尖 P13 自由边奇点 + 高 M Newton（lim/flr 振荡）。
   conforming medium M0.80+ stall 与 LS+clip medium 0.775 死疑似同机制——翼身跨声速真正的
   下一道墙是**两路径共享的 Newton 鲁棒性**，不再是某一尾迹模型的问题。
3. **共同缺陷**：翼尖-TE 奇点未根治（conforming taper 是模型修正非几何修复；LS 侧 B8 阴性）；
   无 fine 级翼身证据（G13.3 阴性 + LS 无 fine 逃逸）；无实验对照（翼身无 reference data，
   所有验证是跨模型/跨分辨率自洽）；cl_fus 无路线。
4. **工程成熟度**：conforming 是"默认即正确"路径（B18 配方 `kutta_estimator="pressure"`）；
   LS 翼身跨声速能力依赖**非默认的 `inboard_clip`**——默认裁切仍带袋。LS 翼身全包线评估
   开题时，第一个裁决点 = inboard_clip 是否提升为翼身默认。

---

## 3. Track V（粘性耦合）可行性评估

### 3.1 设计状态

设计完整且绑定清晰：Drela IBL3 六方程（δ, A, B, Ψ, C_τ1, C_τ2）**表面 Galerkin FE**
（壁面+尾迹片，不做流线积分）+ transpiration BC（无网格运动，RHS-only，δ*=0 逐位回无粘），
V1 松耦合 → V2 拟同时（可选）→ V3 增广 Newton → V4 尾迹片 IBL。注意历史设计笔记的内部张力
已由 track_v.md 裁决：DN6 §10.2 曾推荐 V1 = Green lag-entrainment 流线积分（其模块结构
streamline.py 等已过时），**绑定设计 = IBL3 表面 FE**——开工时以 track_v.md 为准重述 V1
交付物，勿从 DN6 模块图起步。

### 3.2 代码侧：约 70% 的 solver 钩子已存在（子代理逐点核实）

**已有可复用**：
- 壁面 Neumann RHS 装配模板 `solve/wall_correction.py:156`（P11 通量修正，结构同构：3 点积分、
  面积/法向/属主 tet 预计算）；Picard RHS 通道 `solve/picard.py:31 body_source_rhs`；LS Newton
  的 `b_base` 槽（`newton_ls.py:148`）——松耦合下主矩阵/AMG/消元全不动。
- 表面基础设施：壁面三角提取（`wall`/`fuselage` 组）、面积/法向、`wall_triangle_adjacency`
  （`post/surface.py:121`）、着色并行装配（`mesh/coloring.py` 接口通用）、三角积分规则。
- V3 先例齐全：conforming (φ,Γ) 增广 Newton（`newton.py:97`，Woodbury/AMG/ILU 三档 + EW
  forcing + 线搜索/PTC/freeze）；LS 扩展 DOF + 精确分块 Jacobian + `schur_ls.py` 块预处理
  （"AMG-φ + ILU-BL"的结构原型）。
- 网格资产齐全（M0/M3 quasi-2D 做 GV1.1，M1/M4/M5 机翼、M2 翼身两族）；**参考数据已备**：
  `cases/reference_data/rae2822_experiment/`（Case 7 = 经典粘性跨声速标模）、
  `naca0012_experiment/`、`onera_m6_experiment/`。
- 无新第三方依赖（IBL3 = numpy/scipy/numba）。

**缺失需新建**（按工作量排序）：
1. **`viscous/` 包本体 = V1 主体工作量**（六方程非线性表面 FE 装配 + Drela 闭包 + 数值通量/
   稳定化 + Newton 化 + 紧凑表面编号；roadmap 自评 "Track-P-sized effort"）。
2. transpiration 装配函数（照 wall_correction 模板，小）；**conforming Newton 无 RHS 通道**
   （`newton.py:267-347` 需加 `R_free -= (Tᵀb)[free]`，小但触及核心求解器）。
3. cut-mesh 语义决策（IBL 面网格建在 uncut 壁上时需 master 映射把 cut 解喂回——V1 数据布局
   第一设计点）；翼身 `wall`+`fuselage` 表面在交界曲线的缝合（现有邻接只单组工作）。
4. **LS 路径尾迹片源**：`pyfp3d/wake/` 无任何片面积分机制（零等值面多边形积分 = 全新几何代码；
   或体积带近似 = 需设计裁决，偏离"片 RHS"表述）。conforming 路径片源结构完全支持
   （显式 `wake_minus/plus` 面 + slave→master 折叠机制本来就是通量连续的弱形式）。

### 3.3 物理与验证风险（按严重度）

1. **边缘速度输入质量天花板**：P1 壁面梯度在 LE/TE 带 O(h) 偏差（V6 地板、G1.6 重归因的固有
   max-norm 能力）；`wall_tangential_gradient_quadratic` 在楔形角 <~6° 的锐 TE **直接抛
   ValueError**（`post/surface.py:534-570`）——M6 类翼型 TE 区 IBL 输入只能用线性+光滑路径，
   而 TE 正是 IBL 最敏感处（决定 δ*_TE 与脱腔效应）。翼身交界区壁面梯度还带着 cl_fus 未决
   污染；LS 路径翼尖奇点未愈。IBL 消费 u_e **和** du_e/ds——若不先定量输入误差带，
   GV1.x/GV3.x 对标会把"粘性模型误差"和"无粘输入误差"混在一起（next_phase_priorities §3
   论点 1，部分仍成立）。
2. **有效性包线**：attached/mildly-shocked，DN6 建议 M_shock<1.3；M6 M0.84 的激波正在包线
   边缘，翼尖奇点区（M_max 2.5–2.8）远在包线外（局部，但 IBL 种子点布署要避开）。松耦合
   近分离发散是 DN6 自评高概率风险 → GV1.2"5–10 次收敛"在跨声速可能要 under-relaxation 才达。
3. **gate 可验证性**：GV1.1 需要 XFOIL δ* 参考——仓库没有、XFOIL 不是依赖；需外部生成并
   commit（有 `generate_panel_reference.py` 先例，可行但要在 phase 计划里列出）。GV3.3 的
   "CL 降到 0.26–0.27"依赖的无粘基线已变成 0.2823/0.2866（track_v 的 0.245 陈旧，见 §1.2）
   ——**这反而强化了 Track V 的前提**：无粘离散已基本洗干净，剩下到实验的 ~0.02 差确实主要
   是粘性，VII 归因空间比设计书写成时干净得多。
4. **耦合机制交互**：transpiration 注入质量改变 Γ 与远场涡修正（DN6 开放问题 #4：J_φ,δ*
   是否含远场项未分析）；V3 的 BL 块 O(6×壁面节点) 远大于现有 aux 薄带，schur 精确消元
   未必划算，块预处理要新写（构件齐）。

### 3.4 时机判断

- 2026-07-19 的 next_phase 分析把 V 排在 P11 之后，其论据 1（"几何误差随加密恶化，
  transpiration 源直接坐在污染 Cp 上"）**已被 B25/B26 部分削弱**（交界袋愈）；P11 已执行。
  时机的剩余约束 = cl_fus 无归宿、(b) 类天花板未归因、LS 翼尖奇点未愈——三个都在**翼身**线上。
- **V1 的验证阶梯（2.5-D NACA + sphere/cylinder 冒烟）完全不接触翼身伤口**，与 (b) 类归因
  （A-track 小项）无资源冲突之外的耦合——现在开工 V1 可行。
- **V3/V4（翼身 VII、尾迹片 IBL）建议排在 (b) 类归因和 cl_fus 归宿之后**——否则 GV3.3 的
  M6/翼身对标会再次陷入"粘性模型 vs 已知无粘伤口"的归因混淆。
- **建议路径**：V1 先行（2.5-D IBL 阶梯 + transpiration RHS 通道 + GV1.3 的 δ*=0 逐位 gate），
  并行做 (b) 类天花板归因；V1 收尾时用实测 u_e 误差带定量决定 V3 是否进翼身。B25 的
  `inboard_clip` 已把 LS 片拓扑对齐到 conforming，V4 两路径统一处理的物理前置已就位
  （design_track_b §22.3 的自我定位一致）。

---

## 4. 工程纪律核对（本轮审查自身）

- 全程只读；子代理实跑 b1+b2+m2+v0+b18 = **90 passed**（34 s，本机网格在）；线程上限 16。
- 数字引用以 committed CSV/npz 为准；VERDICT 五份（B23/B24/B25/B26/B27）与 demo checks 交叉一致。
- 本报告落盘后执行文档债统一收尾（§1.1 各面 + §1.2 track_v 勘误 + 基线行重测），
  收尾提交与本报告同批；push/PR 待用户确认。

---

## 5. 文档债收尾记录（2026-07-20，本报告同批执行）

§1.1 各面全部落地，逐面结果：

- `docs/roadmap/track_b.md`：标题 B1–B21 → **B1–B27**；L8 状态行追加 B23–B27 一句话条目；
  B18 条目加 ★★★ 勘误块（"junction-limited" 退役）+ Consequence 段 B27 跨模型注记
  （M0.65 2.4% PASS / M0.75 2.5%，M0.5-only 退役；coarse 0.60 旧 0.2% 行退役 → C 侧 0.2133/2.1%）
  + GB18.2/18.3/18.4 superseded 标注；新增 B23–B27 五个 phase 条目（Progress ledger 之前）；
  台账状态行/"as of" 行/台账表格（B23–B27 五行，最新在前）同步。
- `docs/demo_report/track_b.md`：B18 详章顶部 erratum 指针块 + 章内 "1 PASS + 6 RECORDED"
  标注 superseded；文末新增 B23/B24/B25/B26/B27 五个紧凑节。
- `docs/demo_report.md`：索引表新增 B23–B27 五行（B18 行此前已是 B27 刷新版）。
- `docs/agent-rules.md`：Track-B 状态行补 B23/B24/B25 条目（B26/B27 已在）；B18 专条加
  ★★★ erratum 指针 + gates 计数标注；基线行更新（见下）。
- `docs/overview.md`：B-track 状态列 B11–B22 → B11–B27；行内补 B23/B24/B25 一句话；
  "现基线"行与演进链更新（见下）。
- `docs/roadmap.md`：Track-B 行状态列 B11–B21 → B11–B27，补 B23/B24/B25 一句话。
- `CLAUDE.md`：doc-map track_b B1–B21 → B1–B27 + B23–B27 收尾注；workflow 第 3 条基线更新。
- `PROJECT_STRUCTURE.md`：目录树 cases/analysis 段加 b23–b27 五行；`meshgen/fuselage.py`
  注释加 `make_inboard_clip`；`wake/cut_elements.py` 注释加 `inboard_clip`；footer 基线行更新。
- `docs/roadmap/track_v.md`：scope guard 2026-07-20 erratum 块（0.245 陈旧 → 0.2823/0.2866；
  P11 阴性；DN6 §8.3 residual.py 占位无代码对应）——本轮复核已在位。
- 基线重测（全套件，16 线程）：**479 passed + 25 skipped + 2 xfailed，1100.63 s**
  ——与 473+25+2（P11）对账 +6 passed = `TestInboardFragmentClip`（4）+ b1 foot-preference 锁（1）
  + m2 水线延伸锁（1）（`git diff origin/main...HEAD -- tests/` 逐项核对）。
- 历史章节（design_track_b §18、track_b B18 条目/台账行、agent-rules B18 专条）保留原文，
  均已有节首/行首勘误指针指向 B26/B27；2026-07-19 前 inspection 文档为日期快照，不回改。

