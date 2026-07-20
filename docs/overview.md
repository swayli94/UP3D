# pyFP3D 总览（快照 + 文档地图）

> **快照日期：2026-07-19。** 本文件是给人读的高层总览，**不是**权威进度源：
> 阶段/gate 状态以 [roadmap.md](roadmap.md)（track 索引）+ [roadmap/](roadmap/)
> 各 track 文件（含各自的进度台账）为准；当前阶段以 [agent-rules.md](agent-rules.md)
> 为准；证据在 [demo_report.md](demo_report.md) 索引 + [demo_report/](demo_report/)。
> 若本文件与它们冲突，以它们为准。gate 关闭时按 CLAUDE.md 工作流更新台账后，顺手刷新本文件。

## 一句话状态

四条 track：求解器主线（P）与网格线（M）基本关闭，level-set 尾迹线（B）是当前
工作面。

**最新：P11 已关闭（2026-07-19，用户指示当日开+关；sphere 腿）——曲面壁元路线
实测阴性，且 G1.6 被重归因。** ★★ 头条：经验证的曲面壁邻层（tet10 几何 + mapped-P1
场 + ΔA 增量装配，`pyfp3d/solve/curved_wall.py`，默认逐位不变）只把 medium 球
Cp 误差移动 **11.56%→11.33%**（= G1.4 边界数据 oracle 天花板）——预注册风险触发：
mapped-P1 在二次几何上线性重现丢失 **O(h)**（实测 0.138），与被移除的法向误差同阶。
★★ 对照实验重归因 G1.6：同样平坦面片的结构化 icosphere 壳 φ 阶 **1.67/1.98**（二阶，
h≈0.036 时 Cp 2.14%）；h_min 扫掠的阶坍塌（0.88/0.56/0.42，脚本化精确复现）=
**固定中远场网格的污染地板**（固定 h_min=0.03 只加密远场 ⇒ φ_wall 降 3.17×、
argmax 从 r=1.53 回壁面、阶恢复 1.89）⇒ **medium 的 11.6% ≈ P1 场在 h=0.08 的
固有 max-norm 能力**，几何罪份额 ≈0.2pp。G1.6 xfail 保持；路线三岔口
（Option C 重定规格实测可通过形式 / isoparametric P2 壁层 / 接受为长期边界）
= **用户裁决**；翼身"G1.6 类"交叉归因失去球锚、需自己的判别实验。
demo `cases/demo/p11_curved_walls/`（14 PASS+2 XFAIL）、tests
`test_p11_curved_walls.py`（8）。

**B22 已关闭（2026-07-19；执行 B21 的 recorded follow-up + Kimi N3/§2/§5）——
B21 态证据全面刷新 + 3-D LS 数字首次上测试锁。** B15 demo **20/20**（B20 时代
17/20；缓存删净、零 `cached` 行）、B14 demo **7/7**（曾 5/7；medium schur
1.47×、precond →1.8%）；★ coarse ramp 也被 B21 移动（γ 0.0848→**0.084931**，
如实披露）。★★ **N3 缺口关闭**：`tests/test_b22_ls_3d_anchors.py`（+2 gated）
重解 committed coarse/medium ramp 并绝对断言 m_final/γ/M_max/钳位——两天内两次
re-baseline 套件全绿的"无报警"状态到此为止。★ re-baseline 勘误清单入流程
（CLAUDE.md 步骤 5 + 纪律 #11）。★ 下一相位分析
（`analysis/next_phase_priorities_2026-07-19.md`）：**推荐 P11 优先**（G1.6
现在拥有三个随加密恶化的伤口、翼身线上数值嫌疑已清空），LS fine 第二，
Track V 置于 P11 后（V1 阶梯可并行）；M_max 2.4818-vs-1.995 是跨网格族比较、
不得引作缺陷。**下一相位 = 用户裁决。**（★ 后记：P11 已于同日执行——其"三个
伤口同根因"表的**球锚已被 P11 重归因**，翼身两伤口的"G1.6 类"标签待自己的
判别实验，见最新块。）

**B21 已关闭（2026-07-19；执行 Kimi 第二轮审查的 N1 发现）——
`freeze_side_state` 漏打 B20 补丁（在未打补丁的 side 场上捕获冻结选择，3-D 下
与 live 系统不一致，实测 83+9 个元素选择不同），补齐后 **M6 medium M0.84 ramp
能力恢复**：committed recipe 重新到 M0.84（γ 0.088343、res 9.0e-14、clamps 0/1、
515 s，比 pre-B20 的 657 s/3 clamps 更快更干净）。⇒ **GB20.7 的"真实能力损失"
定谳被推翻——损失是 B20 补丁打得不完整（N1），不是修复的内在代价**；
"contamination 是无意稳定器"的 synthesis 同时证伪。GB15.4 能力条款恢复，
数字轻微重基线（γ 0.088338→0.088343）。**

**B20 已关闭（2026-07-18，用户指示）——mixed-side plain 元素的密度改读 main 场，
且 B19 Leg B 的假设有了答案：一个分裂的答案。** 为测 A/B 临时引入的
`plain_density` 开关已随采纳**移除**（见下方 ★★ 永久采纳段）。★ **报告层早已
如此裁决**——`element_mach2` 自 2026-07-14 起默认 `mixed_plain="main"`——B20
让**装配**与**诊断**一致。

★★ **一个缓冲区别名 bug 是靠"测量异常而非解释异常"抓到的**：`PicardOperator.velocities`
返回的是**共享缓冲区的视图**，我在密度路径里重算 main 梯度就把调用方的 side 值原地覆盖了——
quasi-2D "移动"0.77、亚声速 Γ **翻三倍**、Jacobian 退化，全是**同一个 bug**（2940 个元素被
污染，而掩码只有 129），`.copy()` 修复。*若当时接受"Γ 翻三倍 ⇒ Leg B 影响很大"这个看似合理的
叙述，一个纯别名 bug 就会被记成物理发现。*

**GB20.1 ✓** quasi-2D 逐位不变、M6 只动 164/~12k 行；**GB20.2 ✓** main 下 Jacobian 仍**精确**
（8.07e-09／6.29e-10）⇒ Leg A ∘ Leg B 正确复合；**GB20.3 ✓** 2.5-D 亚声速 Γ **+0.0000%**
（该元素类只存在于 3-D ⇒ 所有已提交 quasi-2D 锁一个不动）；**GB20.4 ✓** M6 coarse ramp→M0.84：
side **停 m0.7875 不收敛** vs main **到 M0.84 收敛**。

★★ **GB20.5 —— 假设分裂，这才是发现**。B18 medium 翼身 @M0.5：side res 6.8e-5／**82 个 cell
被钳**／Mmax 3.920（**被钳位的、未收敛的**数）vs main **res 1.1e-13／仅 6 个被钳**／Mmax 5.220
（**真正收敛的解**）。⇒ **收敛病**（churn／钳位）**大部分确实是这个污染**；但**交界袋是真实的**，
B19 的字面假设**被证伪**——去掉污染反而**解除了钳位**，露出亚声速来流下真实的 M≈5.2 尖峰，
即 **G1.6/GB9.4 的刻面几何误差**，不是 mixed-plain 密度。main 仍过不了 M0.5。
**"袋是 mixed-plain 污染"这句话不要再重复——已实测为假。**

★★ **已永久采纳，开关移除**（用户裁决 2026-07-18）。开关当初只为把 A/B 测出来而存在；
既然 side 读法已被确立为**内部不自洽**（一个方程用两个速度场，而且这个单元根本没有尾迹跳跃
穿过），把它留作可选项就等于明知有缺陷还当选项发布。**接受的代价**：3-D 已提交 LS 数字重基线
（2.5-D 不动，实测 +0.0000%）；交界袋不在此列（那是 G1.6 几何问题）。旧值经由 git 与
before/after CSV 保持可追溯。

★★ **重基线（2026-07-19）与 GB20.7→B21 的翻案链。** 重基线实测：套件 465+22+2
不变、gated 3-D LS 67/67 绿（★ 但 3-D 数字本就无测试锁——开放流程缺口 N3）；
所有移动的数字向 B20 预测方向走（B7 M_max 1.453→1.392、B16 legacy limited
3690→11、B18 翼身 res 6.8e-5→1.1e-13、B17 medium Newton 0.2115→0.2114、M6
coarse ramp 0.7875-不收敛→M0.84 收敛），**唯一回退 = M6 medium ramp 停在
M0.6625**。GB20.7 扫 freeze_tol 1e-3→1e-6 只把天花板挪到 0.675，遂定谳"真实
能力损失"——**次日被 B21 推翻**：真机制是 `freeze_side_state` 漏打 B20 补丁
（Kimi/N1），修复后 committed recipe 重回 M0.84。教训入档：freeze_tol 扫的是
"何时武装"，N1 是"武装时捕获了什么"——扫掠在原理上看不见它。

**B19 已关闭（2026-07-18，用户指示）——LS Newton 的 Jacobian 在 3-D 上现已精确，
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
翼身跨声速：交界虚假超声速袋（G1.6/GB9.4，M0.5 就 M²≈1.27）**随加密恶化**——close-out 时
coarse ~M0.575、medium 死于第一级 ~M0.5（Mmax 伪影 1.4→4.0）；★ B20 重基线勘误
（2026-07-19）：coarse 天花板 **~M0.55（Mmax 1.31）**、medium ~M0.5 但 Mmax **5.22 为真实
未钳值**（旧 3.96 是钳位伪影；GB20.5 证实袋=G1.6 几何、非 mixed-plain 污染）；GB9.4 的
同向 closed-negative，纪律 #8 不追。⇒
medium 无共同跨声速 Mach，跨模型只留 M0.5（2.6%）；★ coarse M0.6 跨模型点在重基线工件中
已存在（0.2178 vs 0.2174 = 0.2%，原记 skipped）。GB18.1 PASS + GB18.2–5 RECORDED；偿还 GB16.6
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
conforming——coarse 0.2087、medium **0.2117（Picard）=0.2114（Newton；pre-B20 0.2115）**，
两求解器吻合 0.1%。
★ **B16 混淆了两个正交问题**：远场近奇异**条件数**（pin 治，jump 值无关）与出流**环量**（需
jump=γ）；翼身交界 churn（close-out 时 medium nlim 42/nflr 40）是第三个问题，只限残差地板
不污染升力——★ B20 勘误：该 churn 即 mixed-plain 污染，重基线后同一轨迹收敛到 |R|~1e-13、
钳位归零（GB20.5）。
★ **后处理排查（用户疑点，GB17.2）**：cl_p（压力面积分）与 cl_KJ（环量积分）同步差 ~22% ⇒ 真实
流场态变化非伪影；"Cp 目测对齐却 cl_p 差 22%" 是 Cp 轴尺度错觉，绘制的展向 sectional cl 是 Γ 基
`2Γ/(u·c)`。**默认值（用户裁决）**：pin_gamma 两求解器新默认，仅作用 freestream、vortex/neumann
惰性（committed 逐位不变）；B9 coarse-12.8% erratum 为远场污染（medium 头条仍立）。vortex 从
+2.5% 另一侧 bracket，不闭合缺口。
**B16 ✓ 关闭 2026-07-18（churn 修复；GB16.4 由 B17 解决）**：翼身 LS-Newton churn 根因 = 近奇异
远场 aux 块（8 个远场 MAIN 行 max|R|=**84.457** 逐位复现，cond1 **9.1e18→8.70e6**——勘误
2026-07-19：旧值 6.36e18 是 CSV 前预跑值，committed CSV 为准）；pin 治
条件数（coarse res 5.88e-14、0 limited vs legacy 7.95/3690），但 freestream 的 jump=0 值对升力
是错的（见上）。**B9 ✓ 关闭 2026-07-17（重定规格）：翼身跨模型验证** LS（Picard）+
conforming（全新能力，Newton）在中网格 M0.5 升力一致到 cl_p 0.4% / cl_kj 0.6%。粘性耦合线
（V）设计完毕、零实现。

## Track 状态表

- **P — 求解器**（[roadmap/track_p.md](roadmap/track_p.md)） — P0–P9 ✓（P1 仅 G1.6 以 strict xfail 挂起）；P10 ◐（G10.2/G10.3 ✓）；
  P13 ◐（G13.1 ✓、G13.2 conforming ✓、G13.3 亚声速 Richardson ✓ p=2.31） — G10.1（非升力 Newton 入口，无顺序约束）；G13.3 **跨声速阴性开放**（圆帽 fine 的 ramp 死于 M=0.75，
  site=尖 tip TE）；**P11 ✓ 2026-07-19 关闭（用户指示当日开+关；sphere 腿）**——G11.1 未达（曲面层 11.56%→11.33% = oracle 天花板；superparametric O(h) 风险触发）、G11.2 阴性+
  前提被驳（阶坍塌=固定中远场地板，E8 3.17×/1.89 阶；结构化壳平坦面片 ~2 阶）⇒ **G1.6 重归因为 P1 固有能力**，路线三岔口=用户裁决；P12 backlog；**P14 ✓ 2026-07-17 关闭（用户指示，当日开+关）**：
  壁面邻接 CV 压力相等 Kutta 估计器（A2 路由的修复），G14.1–G14.7 ✓、demo 28 PASS。**头条：修改后的 conforming 结果与 level-set 相同、计算结果合理**——
  S1/S2 一次换掉：M0.84 Γ(z) roughness 0.0970→**0.0043**（coarse）/ 0.0365→**0.0024**（medium，均达/优于 LS 带），全站 raw TE Cp gap 0.2206→**0.0040** /
  0.1585→**0.0024**（**55×/67×**，走 G14.6 **主条款**，回退未动用；★ 首版误引 A2 的 *section 末点* 数 0.318/0.228 报成 80×/95×，当日据用户提问勘误——
  A2 有两套 TE gap 度量，须引用自己实际跑的那套）。**★ 跨模型对照（V14.6，`cross_model_medium_m084.csv`）= 本相位最强证据**：LS 路径**一直**用压力相等 Kutta（B4），故若升力位移真是 Kutta *形式*，
  压力路径必须落到 LS 的答案上——它做到了：medium M0.84 conforming-pressure **0.2776/0.2823** vs level-set **0.2772/0.2813**，**差 0.17%/0.36%**（探针路径当年比 LS **低
  4.5%/4.3%**），而且是不同尾迹模型、不同 DOF 空间、**不同网格族**（`onera_m6_wakefree`）。⇒ 长期存在的 conforming-vs-LS 升力分歧**就是** Kutta 形式误差，已消除。
  保留告诫：跨模型非同网格 A/B；LS 态带 1 lim/2 flr（B15 caveat）而压力态 0/0；"两者一致"≠"两者都对"（刚性平面尾迹等共有模型误差对两者是共模）。**★ G14.7 ✓ 关闭——开相位时锚探针 G8.2 锁、XFAIL 且不改带，
  用户裁决后改锚 level-set oracle**：medium M0.84 压力路径 cl_p/cl_KJ 0.2776/0.2823 vs LS 0.2772/0.2813 = **0.15%/0.34%**（<1% 通过）。
  相对旧探针锁的 +4.85% 位移是**发现**：机理在一级已测并在二级跑前**预注册**（两闭合逐点只差探针自身 O(h) 读数偏差，cross-read medium 0.79%；Kutta 映射 b≈0.93 把它 1/(1−b)≈14× 放大进 Γ），\|
  cl_KJ−0.288\| 0.0188→0.0057，**P9 那 0.019 gap 有 69% 是 Kutta 估计器偏差**（P9 看不见它：两套网格共用同一估计器，对其 Richardson 是共模）。关闭 G14.7 断言的是**两路一致**，**非**网格收敛
  （M6 fine 不是离散解）、**非**"0.019 是分辨率"翻案（仍 *strongly indicated, NOT earned*）。★ **自我更正（V14.7 实测 2026-07-17）**：先前各版都断言 TE Cp **spike** 不会被 P14 触及
  （"两路共有的 P1 恢复伪影"）——那是从 A2 归因**推**的，没测。实测 medium M0.84 raw：探针 **0.1143** → 压力 **0.0533**（2.1×），且**低于 LS 的 0.0743**。
  A2 对的部分：确有共有残余（~0.05 = 真正的恢复地板）；需修正的部分：conforming 相对 LS 的**超出量**也是 Kutta 形式误差（Kutta 错⇒TE 流场真的错，末点偏离趋势有物理原因，共模度量分不开）。
  旁证：P6 平滑在压力路径上不再有效（0.0533→0.0660→0.0626；A2 在探针路径上测的是 0.147→0.081）。**教训：别把上一相位的归因当结论带进新测量——去测。** 诚实记录：判别器 D=7.33→**1.80**（落在 A2 的
  inconclusive 区，非 O(1)）
- **M — 网格**（[roadmap/track_m.md](roadmap/track_m.md)） — M0、M1(+M1b 自相似阶梯)、M2、M3、M4、M5（圆顶翼尖盖）✓ — **M2 ✓（求解腿由 B9 于 2026-07-17 关闭；
  台账勘误 2026-07-19——A3 曾称已改而 ledger 行仍 ◐）**：翼身网格 2026-07-13 交付；**机身+远场 2026-07-16 按用户指示重定规格并重生成**（5 倍翼根弦长、机翼居中、2 倍直径椭球机鼻、蒙皮 h_body=2h_wall +
  两端按半径加密；**R_FAR 15→25 MAC**；★需 `Mesh.OptimizeNetgen` 治 sliver 抽签）；遗留验证项（交界最内 TE 节点 CV fan）在 track_m 记录
- **B — level-set 尾迹**（[roadmap/track_b.md](roadmap/track_b.md)） — B1–B5、B7、B8（characterized-not-cured）、B9、B11–B22 ✓；
  B6 ◐（coarse gate ✓；medium 定量项由 GB15.4 补上，B21 恢复、B22 上锁） — **B16 ✓ 关闭 2026-07-17（用户指示，追加于 B15 之后；执行 B9 的 recorded follow-up）**：
  LS Newton 远场 BC 通用化——远场 aux DOF 钉扎。★ 翼身 LS-Newton churn 根因 = 近奇异远场 aux 块（尾迹片贯穿远场边界的 aux DOF 只受巨型外区单元的 wake-LS 行约束）：
  8 个远场 MAIN 行 max\|R\|=**84.457** 逐位复现，cond1 **9.1e18→8.70e6**（勘误 2026-07-19：旧文 6.36e18 是 CSV 前预跑值）。`farfield_aux="pin"`（默认，按模式适配）使
  **coarse** 翼身 freestream Newton 达 **res 5.88e-14、0 limited**（升力≈conforming 0.1%），legacy churn 7.95/3690 limited；
  neumann 逐位不变。★★ **GB16.4 由 B17 解决（2026-07-18）**：不是非收敛，是 freestream pin 的 BC 建模错误（jump=0 杀出流环量）。⚠ 提案机理被自我修正（Picard 升力也用 wake_ls，差别在不动点吸收非
  closure）。 · **B17 ✓ 关闭 2026-07-18（用户指示；解决 GB16.4）**：远场 pin 必须携带 jump=γ 而非 0——B16 的 jump=0 抹掉出流尾迹环量（medium −22%）；
  独立 Picard-pin 收敛到 Newton-pin "停滞" 的同一 0.169（两求解器同 BC 一致，非停滞）。`farfield_aux="pin_gamma"`（jump→γ，两求解器新默认）三角单调闭合向 conforming（coarse 0.2087、
  medium 0.2117 Picard/0.2114 Newton——pre-B20 0.2115）；仅作用 freestream、vortex/neumann 惰性；vortex 从 +2.5% bracket；
  B9 coarse-12.8% erratum 为远场污染。GB17.1–17.4 ✓、17.5/17.6 RECORDED；demo `cases/demo/b17_farfield_pin_gamma/`、tests
  `test_b17_farfield_pin_gamma.py`(6)。**B18 ✓ 关闭 2026-07-18（用户指示；执行 GB16.6 债）**：翼身跨声速 M0.84——conforming 到（coarse 0.2617、medium M0.79 0.2579、
  cl(M) 升 0.2173/0.2321/0.2579）、level-set 交界受限（close-out coarse ~M0.575、medium 死 ~M0.5；B20 重基线：coarse ~M0.55/Mmax 1.31、medium Mmax 5.22
  真实未钳，袋=G1.6 几何而非 mixed-plain——GB20.5 勘误，closed-negative）；medium 无共同跨声速 Mach ⇒ 跨模型只留 M0.5（2.6%）+ coarse M0.6 点（0.2%，重基线工件）；
  GB18.1 PASS+18.2–5 RECORDED；零 pyfp3d/ 改动。**B19 ✓ 关闭 2026-07-18**：LS-Newton Jacobian 3-D 精确化（两处缺项；探针 1.146e-01→1.33e-08、ε 判别翻转；
  R 逐位不变、无收敛收益——GB19.4 阴性）；Leg B 测得 mixed-plain 侧场密度污染（虚假超声速 q² 3.22 vs 1.34）→ 路由 B20。**B20 ✓ 关闭 2026-07-18、永久采纳（用户裁决，开关移除）、重基线
  2026-07-19**：mixed-plain 密度改读 main 场；所有移动数字向好，唯一回退 = M6 medium ramp（GB20.7 曾定谳"真实能力损失"）。**B21 ✓ 关闭 2026-07-19（执行 Kimi/N1）**：
  `freeze_side_state` 漏打 B20 补丁——修复后 **M6 medium M0.84 恢复**（γ 0.088343、res 9e-14、515 s），GB20.7 翻案；3-D freeze 捕获一致性测试锁落地。
  **B22 ✓ 关闭 2026-07-19（执行 B21 follow-up + Kimi N3/§2/§5）**：B15 demo 20/20 + B14 7/7 刷新（coarse ramp γ 0.084931 披露）；
  **N3 关闭**——gated 绝对锚锁 `test_b22_ls_3d_anchors.py`；re-baseline 勘误清单入流程；下一相位分析推荐 P11（用户裁决）。**B9 ✓ 关闭 2026-07-17（重定规格）**：
  翼身跨模型 LS+conforming 一致 0.4%/0.6%；GB9.4 XFAIL⇒G1.6。B10 搁置
- **V — 粘性耦合**（[roadmap/track_v.md](roadmap/track_v.md)） — 设计完整（Drela IBL3 + transpiration BC），零实现 — V1 依赖 P6（已满足），预算等同一个 Track-P 阶段
- **A — 校验与分析**（[roadmap/track_a.md](roadmap/track_a.md)） — 2026-07-15 新建；**A1 ✓ 2026-07-16**（GA1.1–GA1.5：四求解器统一计时插桩 + conforming×level-set
  × Picard×Newton 耗时基准） — **A2 ✓ 2026-07-17 关闭**（TE/Kutta 保真度归因，GA2.1–GA2.5）：**S1 定谳**——conforming Γ(z) 逐站抖动是逐站探针差势跳 Kutta target **估计器**
  的测量伪影（fixed-Γ 判别量 D=7.33/25.70 coarse/medium，把抖动从光滑场里重新生出来；闭合残差 ≤0.6% 排除"未闭合"、抖动局域于 TE 邻层 0.02–0.07× 排除"流场"），**非流场内容**；
  **S2 分解**——TE Cp 突跳=势跳 Kutta 形式误差（conforming 独有,同估计器 34×/133× vs LS）+ P1 末点恢复伪影（两路共有）；2.5-D `a1_cp.png`（有 S2 无 S1）证明两者是不同机制。
  修复路由至 **P14**（无 `pyfp3d/` 改动）；工作目录 `cases/analysis/`（区别于 `cases/demo/`）。★A1 结论：3-D 下两条 Newton 路径都是 **precond（LU 分解）受限**（~40% 墙钟，lagged LU
  已开），2.5-D 的"seed 是成本"**不外推**——引用主导相位必须带网格。**A3 ✓ 2026-07-18 关闭**（GA3.1–GA3.6：响应 Kimi 独立审查）：★★ **C1 证实**——LS Newton Jacobian 在 mixed-side
  plain 元素上≠dR/dφ（targeted 1.146e-01 vs control 6.33e-10，且 eps 无关 max/min=1.00）⇒ 3-D 下是 quasi-Newton；**R 未动 ⇒ 所有已收敛状态与 gate 数字成立**；
  记录不修，修复单独立项。C2/C3 回移 conforming Newton、reader C4/C5（静默丢面组 ⇒ Γ(root) 静默为 0）、C6/C7/P1/T1/T2/F0；文档 17 项全部处置；close-out 扩为五面 + backport check；
  零 pyfp3d/ 数值改动

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

现基线 **465 passed + 25 skipped + 2 xfailed**（2026-07-19 B22 3-D LS 锚锁：
+2 skipped = gated `tests/test_b22_ls_3d_anchors.py`；实测 1127.38 s @16 线程）。
上一档 465+23+2（2026-07-19 B21 freeze 捕获对齐：
+1 skipped = `tests/test_b15_ls_newton_freeze.py` 的 gated 3-D 捕获一致性锁；
实测 1105.87 s @16 线程）。
上一档 465+22+2（2026-07-18 B19 LS-Newton Jacobian 精确化：
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
4 ungated；翼身跨声速 ramp 在门控 demo；执行 GB16.6 债）**
→ 463+21+2（A3 审查响应，2026-07-18：+3 reader 锁）→ 465+22+2（B19，2026-07-18：
+2/+1 = `test_b19_jacobian_3d.py`；B20 重基线 2026-07-19 复测不变）
→ 465+23+2（B21 freeze 捕获对齐，2026-07-19：+1 skipped = gated 3-D 捕获锁；
1105.87 s @16 线程）
→ 465+25+2（B22 3-D LS 锚锁，2026-07-19：+2 skipped = gated
`test_b22_ls_3d_anchors.py`；实测 1127.38 s @16 线程）
→ **473+25+2（P11 曲面壁元，2026-07-19：+8 passed = ungated
`test_p11_curved_walls.py`；实测 1124.94 s @16 线程）**。

## 长期挂起项（勿反复重提）

- **G1.6 球面 Cp <2%**（strict xfail，11.6%）：★★ **根因已于 2026-07-19 由 P11
  重归因**——旧说"平坦面片壁的自然边界条件（变分罪）"被实测推翻（罪的份额
  ≈0.2pp；11.6% ≈ P1 场在 h=0.08 的固有能力；阶坍塌是混淆扫掠的中远场地板）。
  仍然成立：恢复非主导、Nitsche 死、边界数据修正无可修正、h 加密（同族）贵而不达。
  新增死路：mapped-P1（superparametric）曲面壁元。在案路线三岔口（P11 close-out，
  用户裁决）：Option C 重定规格（实测可通过形式：全尺度同步加密族阶 ≥1.8 +
  h_min 0.03 mean-Cp <1%）/ isoparametric P2 壁层（唯一能触及字面 2%-max@medium
  的路线）/ 接受为长期边界。
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
- **N3 流程缺口 ✓ 已由 B22 关闭（2026-07-19）**：3-D LS 核心数字（M6 coarse/medium
  ramp 的 m_final/γ/M_max/钳位）现由 `tests/test_b22_ls_3d_anchors.py` 的 gated
  绝对锚锁重解+断言——下一次无声 re-baseline 会使套件报警；B15/B14 demo 已按
  B21 数字刷新（20/20、7/7）。仍开放：M6 medium M_max 2.4818 只有 LS-vs-LS 共模
  验证（conforming 记录 1.995 是**跨网格族**比较，勿引作缺陷；同族核对 = 便宜
  A-track 小项，见 analysis/next_phase_priorities_2026-07-19.md §4）。


> **V14.6 的两组数字（A3 2026-07-18 溯源注）**：跨模型一致性同时以 **0.17%/0.36%**（对 A1 缓存的 LS 全精度值，V14.6 行）和 **0.15%/0.34%**（对 `run_demo.py` 中四位小数的 `LS_REF` 常量，G14.7 行）出现，两者都在已提交的 `checks.csv` 里。**这是同一次测量**，~0.02pp 之差是参照值舍入，对 "< 1%" 结论无影响；引用时以全精度的 **0.17%/0.36%** 为准。
