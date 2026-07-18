# pyFP3D 总览（快照 + 文档地图）

> **快照日期：2026-07-18。** 本文件是给人读的高层总览，**不是**权威进度源：
> 阶段/gate 状态以 [roadmap.md](roadmap.md)（track 索引）+ [roadmap/](roadmap/)
> 各 track 文件（含各自的进度台账）为准；当前阶段以 [agent-rules.md](agent-rules.md)
> 为准；证据在 [demo_report.md](demo_report.md) 索引 + [demo_report/](demo_report/)。
> 若本文件与它们冲突，以它们为准。gate 关闭时按 CLAUDE.md 工作流更新台账后，顺手刷新本文件。

## 一句话状态

四条 track：求解器主线（P）与网格线（M）基本关闭，level-set 尾迹线（B）是当前
工作面。

**最新：B19 已关闭（2026-07-18，用户指示）——LS Newton 的 Jacobian 在 3-D 上现已精确，
并测清了残差自身的不对称**（执行 A3/GA3.6 的 C1 发现，按用户裁决**分两条腿**）。

★★ **Leg A 是两处缺项，不是一处。**（1）**DOF 映射**：Terms 2/3 行列共用了质量守恒的
**散射**映射，而列必须跟随 `side_potentials` 的**逐节点读取**映射（在 cut 元素上两者恒等
——`readvec` 复现 `dofs_upper/lower`，已断言——只在 mixed-side plain 这个**仅存于 3-D**
的元素类上分歧）。（2）**梯度因子**，同一二元性下沉一层：残差是
`ρ̃(读取场梯度)·V·(散射场梯度·B_a)`，行因子须用 `grad_row`、列因子保持侧场梯度，而代码
两处都用了侧场。★ **只修（1）后仍是 1.4697e-02——8 倍改善但依旧与 ε 无关 ⇒ 记为 PARTIAL
而非四舍五入成通过，正是这个决定逼出了（2）。** ★ 块隔离（`|FD23−J23| ≡ |FD23−J2|` 精确
相等 ⇒ Term 3 在那里毫无贡献）找到了（2）——**我据行分类的第一次推断指向 Term 3，是错的**；
若照此动手，会在正确代码上制造新 bug。

**GB19.1 ✓** targeted 探针 **1.145684e-01 → 1.333699e-08**（control 6.33e-10 不变）；
★★ **ε 判别量翻转**：修复前三个 ε 恒为 1.532e-01（spread 1.00 = 真实缺项），修复后
1.6e-09/2.1e-08/2.2e-07（spread 131.5，按 ~1/ε = 纯 FD 舍入）。**GB19.2 ✓**
`max|ΔR| = 0.000e+00` 逐位不变（两次修复后各以 git stash A/B 验证）⇒ 已提交的 LS 结果
一个都不会动。★ **GB19.4 ✓ 阴性如实记录——没有任何收敛收益**：γ 0.07212068 到 8 位、
M_max 1.134235 到 6 位完全相同，步数不变，**+3.6% 墙钟**；那个平台是 **B15 的选择 churn
极限环**，精确导数治不了不连续选择 ⇒ **本阶段不得被记为收敛改善**，它买到的是**正确性**
（现在是 Newton 而非 quasi-Newton）。**GB19.5 ✓** 新增 `tests/test_b19_jacobian_3d.py`（+3）
堵住盲区；★★ **勘误：盲区此前被说错**（C1 原文与我自己的第一版测试都错）——quasi-2D 有
**129** 个 mixed-side plain 元素而非零，它为零的是**能读到 aux** 的那一类（129 个里 0 个
接触 cut 节点），这才是真正的不变量。

★★ **GB19.6（Leg B）——残差的不对称并不无害。** 真正读 aux 的只有 **252** 个元素
（占体积 **0.19%**），但那里 max|ρ_side − ρ_main| = **0.4474（45.3%）**，且决定性的是
**侧场读出 q² 高达 3.2229（M≈1.80，已顶到 M_cap 限制器）而 main 场只有 1.3379** ——
**虚假超声速状态**，人工密度开关随后就作用在这个虚构的 q² 上，即污染进入了**求解器的密度**
而不只是诊断量。★ 同一元素类第三次咬人（B8 的 ×5 度量伪影 → Jacobian → 现在是残差）。
★ **假设而非结论**：可能是 B18 翼身"随加密恶化"虚假袋（M_max 伪影 3.96）的贡献者之一；
但测的是**纯机翼** M6，未证因果，已记下具名验证方案。**不予采纳**——改密度来源即改 R、
会移动所有已收敛答案 ⇒ 单独立项。证据 `cases/analysis/c1_ls_jacobian_fd/`。

**A3 已关闭（2026-07-18，用户指示）——响应 2026-07-17 Kimi 独立审查**
（`docs/inspection/` 三份：文档一致性 17 项、代码审查 C1–C7/T1/T2/P1/P2、规划评估）。
所有发现先在当前树逐条复核，**全部仍然成立**。★★ **头条：C1 被审查自己要求的验证
证实了——LS Newton 的 Jacobian 在 mixed-side plain 元素上不是其残差的导数**（M6 coarse
M0.70：targeted 探针 ‖Jv−FD‖/‖FD‖ = **1.146e-01** vs control **6.33e-10**，差 8 个数量级；
且 **与 eps 无关**——eps 1e-6/1e-7/1e-8 恒为 1.532e-01，max/min **1.00**，跨三个数量级
⇒ 缺项而非 FD 噪声）。**后果有界：R 未被触及 ⇒ 所有已收敛 LS 状态、γ、cl、M_max 与
gate 数字全部成立；退化的是收敛速率——3-D 下 LS Newton 是 quasi-Newton。**
**记录，不修**（side-aware 列映射是已发布内核改动，会移动已提交的步数轨迹 ⇒ 单独立项、
用户裁决）。另：C2/C3 两个 B15 期 LS 修复回移到 conforming `newton.py`（悬置三个 phase）；
reader C4/C5 加固（C4 是**静默**的：只命名部分物理面组会丢弃其余三角形，链条终点是
**Γ(root) 被静默钉为 0** 且无任何报错）；C6/C7/P1/T1/T2/F0 修复；**零 pyfp3d/ 数值改动**。
close-out 流程扩为**五个面 + backport check**。响应报告
[inspection/20260718-response-to-kimi-inspection.md](inspection/20260718-response-to-kimi-inspection.md)。

**B18 已关闭（2026-07-18，用户指示，追加于 B17 之后；执行 GB16.6 债）：翼身
组合体跨声速 M0.84——conforming 到、level-set 交界受限**。翼身跨声速能力**不对称，这个
不对称本身就是结论**：conforming（Newton+压力 Kutta+Mach 续接）是翼身跨声速路——coarse 到
**M0.84（cl_p 0.2617）**、medium 到 **M0.79 严格（0.2579）**，干净 cl(M) 升 0.2173/0.2321/0.2579
@ M0.50/0.65/0.79（medium M0.80+ 停滞：非 sliver，更锐激波/交界相互作用，记录不追；★ 需把
freeze_tol 抬到翼身 churn 地板 1e-6→1e-5，B17 教训）。level-set（B15 ramp+B17 pin_gamma）**不能**
翼身跨声速：交界虚假超声速袋（G1.6/GB9.4，M0.5 就 M²≈1.27）**随加密恶化**——coarse ~M0.575、
medium 死于第一级 ~M0.5（Mmax 伪影 1.4→4.0）；GB9.4 的同向 closed-negative，纪律 #8 不追。⇒
medium 无共同跨声速 Mach，跨模型只留 M0.5（2.6%）。GB18.1 PASS + GB18.2–5 RECORDED；偿还 GB16.6
证据债；**零 pyfp3d/ 数值改动**。tests `test_b18_wingbody_transonic.py`(4)、demo
`cases/demo/b18_wingbody_transonic/`(7 gate)。
**B17 已关闭（2026-07-18，用户指示，追加于 B16 之后；解决 GB16.4）：远场 aux
钉扎必须携带 jump=γ，而非 0**。★★ **GB16.4 不是非收敛，是 B16 freestream pin 的边界条件
建模错误**：B16 把出流尾迹跳跃强制为 **0**，抹掉尾迹携带到边界的物理环量（medium cl_p
0.2165→0.1690，−22% 分辨率相关误差；coarse "吻合" 是巧合——jump=0 恰抵消 coarse legacy 的
外层 tet 垃圾）。**判别性证据**：给 Picard 驱动加同样 freestream pin（`solve_multivalued_lifting`
新旋钮），medium Picard-pin 干净收敛（res 7.5e-8）到 cl_p **0.1691**，与"停滞"的 Newton-pin
**0.1690** 吻合 0.1% ⇒ 两个独立求解器落到同一 BC 决定的态，**非** Newton 停滞。修复 =
`farfield_aux="pin_gamma"`（aux=host φ∞−side·γ，jump→γ，两个求解器新默认）：三角单调闭合向
conforming——coarse 0.2087、medium **0.2117（Picard）=0.2115（Newton）**，两求解器吻合 0.1%。
★ **B16 混淆了两个正交问题**：远场近奇异**条件数**（pin 治，jump 值无关）与出流**环量**（需
jump=γ）；翼身交界 churn（medium 仍 nlim 42/nflr 40）是第三个既有问题，只限残差地板不污染升力。
★ **后处理排查（用户疑点，GB17.2）**：cl_p（压力面积分）与 cl_KJ（环量积分）同步差 ~22% ⇒ 真实
流场态变化非伪影；"Cp 目测对齐却 cl_p 差 22%" 是 Cp 轴尺度错觉，绘制的展向 sectional cl 是 Γ 基
`2Γ/(u·c)`。**默认值（用户裁决）**：pin_gamma 两求解器新默认，仅作用 freestream、vortex/neumann
惰性（committed 逐位不变）；B9 coarse-12.8% erratum 为远场污染（medium 头条仍立）。vortex 从
+2.5% 另一侧 bracket，不闭合缺口。
**B16 ✓ 关闭 2026-07-18（churn 修复；GB16.4 由 B17 解决）**：翼身 LS-Newton churn 根因 = 近奇异
远场 aux 块（8 个远场 MAIN 行 max|R|=**84.457** 逐位复现，cond1 **6.36e18→8.70e6**）；pin 治
条件数（coarse res 5.88e-14、0 limited vs legacy 7.95/3690），但 freestream 的 jump=0 值对升力
是错的（见上）。**B9 ✓ 关闭 2026-07-17（重定规格）：翼身跨模型验证** LS（Picard）+
conforming（全新能力，Newton）在中网格 M0.5 升力一致到 cl_p 0.4% / cl_kj 0.6%。粘性耦合线
（V）设计完毕、零实现。

## Track 状态表

| Track | 状态 | 开放项 / 下一步 |
|-------|------|----------------|
| **P — 求解器**（[roadmap/track_p.md](roadmap/track_p.md)） | P0–P9 ✓（P1 仅 G1.6 以 strict xfail 挂起）；P10 ◐（G10.2/G10.3 ✓）；P13 ◐（G13.1 ✓、G13.2 conforming ✓、G13.3 亚声速 Richardson ✓ p=2.31） | G10.1（非升力 Newton 入口，无顺序约束）；G13.3 **跨声速阴性开放**（圆帽 fine 的 ramp 死于 M=0.75，site=尖 tip TE）；P11 条件性未开（仅剩 G1.6 理由）；P12 backlog；**P14 ✓ 2026-07-17 关闭（用户指示，当日开+关）**：壁面邻接 CV 压力相等 Kutta 估计器（A2 路由的修复），G14.1–G14.7 ✓、demo 28 PASS。**头条：修改后的 conforming 结果与 level-set 相同、计算结果合理**——S1/S2 一次换掉：M0.84 Γ(z) roughness 0.0970→**0.0043**（coarse）/ 0.0365→**0.0024**（medium，均达/优于 LS 带），全站 raw TE Cp gap 0.2206→**0.0040** / 0.1585→**0.0024**（**55×/67×**，走 G14.6 **主条款**，回退未动用；★ 首版误引 A2 的 *section 末点* 数 0.318/0.228 报成 80×/95×，当日据用户提问勘误——A2 有两套 TE gap 度量，须引用自己实际跑的那套）。**★ 跨模型对照（V14.6，`cross_model_medium_m084.csv`）= 本相位最强证据**：LS 路径**一直**用压力相等 Kutta（B4），故若升力位移真是 Kutta *形式*，压力路径必须落到 LS 的答案上——它做到了：medium M0.84 conforming-pressure **0.2776/0.2823** vs level-set **0.2772/0.2813**，**差 0.17%/0.36%**（探针路径当年比 LS **低 4.5%/4.3%**），而且是不同尾迹模型、不同 DOF 空间、**不同网格族**（`onera_m6_wakefree`）。⇒ 长期存在的 conforming-vs-LS 升力分歧**就是** Kutta 形式误差，已消除。保留告诫：跨模型非同网格 A/B；LS 态带 1 lim/2 flr（B15 caveat）而压力态 0/0；"两者一致"≠"两者都对"（刚性平面尾迹等共有模型误差对两者是共模）。**★ G14.7 ✓ 关闭——开相位时锚探针 G8.2 锁、XFAIL 且不改带，用户裁决后改锚 level-set oracle**：medium M0.84 压力路径 cl_p/cl_KJ 0.2776/0.2823 vs LS 0.2772/0.2813 = **0.15%/0.34%**（<1% 通过）。相对旧探针锁的 +4.85% 位移是**发现**：机理在一级已测并在二级跑前**预注册**（两闭合逐点只差探针自身 O(h) 读数偏差，cross-read medium 0.79%；Kutta 映射 b≈0.93 把它 1/(1−b)≈14× 放大进 Γ），\|cl_KJ−0.288\| 0.0188→0.0057，**P9 那 0.019 gap 有 69% 是 Kutta 估计器偏差**（P9 看不见它：两套网格共用同一估计器，对其 Richardson 是共模）。关闭 G14.7 断言的是**两路一致**，**非**网格收敛（M6 fine 不是离散解）、**非**"0.019 是分辨率"翻案（仍 *strongly indicated, NOT earned*）。★ **自我更正（V14.7 实测 2026-07-17）**：先前各版都断言 TE Cp **spike** 不会被 P14 触及（"两路共有的 P1 恢复伪影"）——那是从 A2 归因**推**的，没测。实测 medium M0.84 raw：探针 **0.1143** → 压力 **0.0533**（2.1×），且**低于 LS 的 0.0743**。A2 对的部分：确有共有残余（~0.05 = 真正的恢复地板）；需修正的部分：conforming 相对 LS 的**超出量**也是 Kutta 形式误差（Kutta 错⇒TE 流场真的错，末点偏离趋势有物理原因，共模度量分不开）。旁证：P6 平滑在压力路径上不再有效（0.0533→0.0660→0.0626；A2 在探针路径上测的是 0.147→0.081）。**教训：别把上一相位的归因当结论带进新测量——去测。** 诚实记录：判别器 D=7.33→**1.80**（落在 A2 的 inconclusive 区，非 O(1)）|
| **M — 网格**（[roadmap/track_m.md](roadmap/track_m.md)） | M0、M1(+M1b 自相似阶梯)、M3、M4、M5（圆顶翼尖盖）✓ | M2 ◐：翼身网格 ✓（2026-07-13；**机身+远场 2026-07-16 按用户指示重定规格并重生成**：5 倍翼根弦长、机翼居中、2 倍直径椭球机鼻、蒙皮 h_body=2h_wall + 两端按半径加密；**R_FAR 15→25 MAC**，h_far 与所有固定加密距离同比放大 ⇒ 2.78× 域几乎不要钱；★需 `Mesh.OptimizeNetgen` 治尾迹 corridor 在对称面压出的细带 sliver **抽签**），求解腿 = B9 |
| **B — level-set 尾迹**（[roadmap/track_b.md](roadmap/track_b.md)） | B1–B5、B7、B8（characterized-not-cured）、B9、B11–B16 ✓；B6 ◐（coarse gate ✓；medium 定量项由 GB15.4 补上） | **B16 ✓ 关闭 2026-07-17（用户指示，追加于 B15 之后；执行 B9 的 recorded follow-up）**：LS Newton 远场 BC 通用化——远场 aux DOF 钉扎。★ 翼身 LS-Newton churn 根因 = 近奇异远场 aux 块（尾迹片贯穿远场边界的 aux DOF 只受巨型外区单元的 wake-LS 行约束）：8 个远场 MAIN 行 max\|R\|=**84.457** 逐位复现，cond1 **6.36e18→8.70e6**。`farfield_aux="pin"`（默认，按模式适配）使 **coarse** 翼身 freestream Newton 达 **res 5.88e-14、0 limited**（升力≈conforming 0.1%），legacy churn 7.95/3690 limited；neumann 逐位不变。★★ **GB16.4 由 B17 解决（2026-07-18）**：不是非收敛，是 freestream pin 的 BC 建模错误（jump=0 杀出流环量）。⚠ 提案机理被自我修正（Picard 升力也用 wake_ls，差别在不动点吸收非 closure）。 · **B17 ✓ 关闭 2026-07-18（用户指示；解决 GB16.4）**：远场 pin 必须携带 jump=γ 而非 0——B16 的 jump=0 抹掉出流尾迹环量（medium −22%）；独立 Picard-pin 收敛到 Newton-pin "停滞" 的同一 0.169（两求解器同 BC 一致，非停滞）。`farfield_aux="pin_gamma"`（jump→γ，两求解器新默认）三角单调闭合向 conforming（coarse 0.2087、medium 0.2117 Picard/0.2115 Newton）；仅作用 freestream、vortex/neumann 惰性；vortex 从 +2.5% bracket；B9 coarse-12.8% erratum 为远场污染。GB17.1–17.4 ✓、17.5/17.6 RECORDED；demo `cases/demo/b17_farfield_pin_gamma/`、tests `test_b17_farfield_pin_gamma.py`(6)。**B18 ✓ 关闭 2026-07-18（用户指示；执行 GB16.6 债）**：翼身跨声速 M0.84——conforming 到（coarse 0.2617、medium M0.79 0.2579、cl(M) 升 0.2173/0.2321/0.2579）、level-set 交界受限（coarse ~M0.575、medium 死 ~M0.5，交界袋随加密恶化 Mmax 1.4→4.0，closed-negative）；medium 无共同跨声速 Mach ⇒ 跨模型只留 M0.5（2.6%）；GB18.1 PASS+18.2–5 RECORDED；零 pyfp3d/ 改动。**B9 ✓ 关闭 2026-07-17（重定规格）**：翼身跨模型 LS+conforming 一致 0.4%/0.6%；GB9.4 XFAIL⇒G1.6。B10 搁置 |
| **V — 粘性耦合**（[roadmap/track_v.md](roadmap/track_v.md)） | 设计完整（Drela IBL3 + transpiration BC），零实现 | V1 依赖 P6（已满足），预算等同一个 Track-P 阶段 |
| **A — 校验与分析**（[roadmap/track_a.md](roadmap/track_a.md)） | 2026-07-15 新建；**A1 ✓ 2026-07-16**（GA1.1–GA1.5：四求解器统一计时插桩 + conforming×level-set × Picard×Newton 耗时基准） | **A2 ✓ 2026-07-17 关闭**（TE/Kutta 保真度归因，GA2.1–GA2.5）：**S1 定谳**——conforming Γ(z) 逐站抖动是逐站探针差势跳 Kutta target **估计器**的测量伪影（fixed-Γ 判别量 D=7.33/25.70 coarse/medium，把抖动从光滑场里重新生出来；闭合残差 ≤0.6% 排除"未闭合"、抖动局域于 TE 邻层 0.02–0.07× 排除"流场"），**非流场内容**；**S2 分解**——TE Cp 突跳=势跳 Kutta 形式误差（conforming 独有,同估计器 34×/133× vs LS）+ P1 末点恢复伪影（两路共有）；2.5-D `a1_cp.png`（有 S2 无 S1）证明两者是不同机制。修复路由至 **P14**（无 `pyfp3d/` 改动）；工作目录 `cases/analysis/`（区别于 `cases/demo/`）。★A1 结论：3-D 下两条 Newton 路径都是 **precond（LU 分解）受限**（~40% 墙钟，lagged LU 已开），2.5-D 的"seed 是成本"**不外推**——引用主导相位必须带网格。**A3 ✓ 2026-07-18 关闭**（GA3.1–GA3.6：响应 Kimi 独立审查）：★★ **C1 证实**——LS Newton Jacobian 在 mixed-side plain 元素上≠dR/dφ（targeted 1.146e-01 vs control 6.33e-10，且 eps 无关 max/min=1.00）⇒ 3-D 下是 quasi-Newton；**R 未动 ⇒ 所有已收敛状态与 gate 数字成立**；记录不修，修复单独立项。C2/C3 回移 conforming Newton、reader C4/C5（静默丢面组 ⇒ Γ(root) 静默为 0）、C6/C7/P1/T1/T2/F0；文档 17 项全部处置；close-out 扩为五面 + backport check；零 pyfp3d/ 数值改动 |

## 文档地图（每份文档的职能、权威范围、何时更新）

| 文档 | 职能 | 权威范围 | 何时更新 |
|------|------|----------|----------|
| [roadmap.md](roadmap.md) | track 索引 + working rules + gate 编号约定 | 各 track 一行状态 | track 状态行变化时 |
| [roadmap/track_{p,m,b,v,a}.md](roadmap/) | 各 track 的 phase 条目 + gate 清单 + 进度台账 | **阶段/gate 状态的唯一权威** | gate 开/关时 |
| [design.md](design.md) | 理论与数值参考（方程、离散、内核规则、求解策略、V0–V6 验证阶梯） | 数值方法（conforming 路径 + 共享理论） | 方法/勘误变化时 |
| [design_track_b.md](design_track_b.md) | Track B（level-set 尾迹）数值方案 + 逐阶段技术结论 | Track B 数值 | Track B 阶段推进时 |
| [demo_report.md](demo_report.md) + [demo_report/](demo_report/) | 已关阶段的证据档案（每阶段一个自检 demo + 提交的图/CSV） | 证据；**无 committed 工件的断言不是证据** | 阶段关闭时加节 |
| [agent-rules.md](agent-rules.md) | 每 session 注入的当前阶段 + 操作纪律（经 CLAUDE.md `@include`） | 当前阶段行 | 阶段变化时 |
| overview.md（本文件） | 人读总览 + 文档地图 | 无（快照） | 顺手刷新 |
| [analysis/](analysis/) | 分析/审查类报告（capability review 等），非规范文档 | 无（报告注明快照日期） | 新报告放这里 |
| [archive/](archive/) | 历史归档（勿作规范；rule 11 同样适用） | 无 | 只进不改 |
| `docs/references/` | 外部文献（López dissertation PDF 等）；**gitignored、未纳入版本库** ⇒ 新 clone 不含此目录，而 design_track_b/track_b 大量按章节引用（"López eq. 3.33–3.34"、"López p.57"），需自备 PDF | — | — |

已删除：`docs/STATUS.md`（2026-07-15，本文件取代）、`docs/discussion_notes/`
（2026-07-14，commit 0e4895a；历史经 `git show 8aa4aee:docs/discussion_notes/<file>`）。

## 回归基线

现基线 **465 passed + 22 skipped + 2 xfailed**（2026-07-18 B19 LS-Newton Jacobian 精确化：
+2 passed / +1 skipped = `tests/test_b19_jacobian_3d.py`；实测 1101.50 s @16 线程）。
上一档 463+21+2（2026-07-18 A3 审查响应，+3 = `tests/test_mesh_reader_roundtrip.py`）。上一档 460+21+2（2026-07-18 B18 翼身跨声速，+4 =
`tests/test_b18_wingbody_transonic.py`，均 ungated）。
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
1015.17 s @8 线程；406 + 15 = 421 逐项对账，零回归）**
→ 429+19+2（B14 Schur+AMG：+8/+1）→ 442+20+2（B9 翼身跨模型：+13/+1，1084.20 s
@16 线程）→ 450+21+2（B16 远场 aux 钉扎，2026-07-17：+8 passed / +1 skipped =
`tests/test_b16_farfield_aux.py`；8 ungated（其中 2 条依赖翼身网格，CI 无网格时跳过）
+ 1 门控 GB16.3 跳过）→ 456+21+2（B17 远场 pin_gamma，2026-07-18：+6 =
`tests/test_b17_farfield_pin_gamma.py`，6 ungated；1097.11 s @16 线程；解决 GB16.4）
→ **460+21+2（B18 翼身跨声速，2026-07-18：+4 = `tests/test_b18_wingbody_transonic.py`，
4 ungated；翼身跨声速 ramp 在门控 demo；执行 GB16.6 债）**。

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


> **V14.6 的两组数字（A3 2026-07-18 溯源注）**：跨模型一致性同时以 **0.17%/0.36%**（对 A1 缓存的 LS 全精度值，V14.6 行）和 **0.15%/0.34%**（对 `run_demo.py` 中四位小数的 `LS_REF` 常量，G14.7 行）出现，两者都在已提交的 `checks.csv` 里。**这是同一次测量**，~0.02pp 之差是参照值舍入，对 "< 1%" 结论无影响；引用时以全精度的 **0.17%/0.36%** 为准。
