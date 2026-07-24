# pyFP3D 总览（快照 + 文档地图）

> **快照日期：2026-07-22（B28–B32 收尾）。** 本文件是给人读的高层总览，**不是**权威进度源：
> 阶段/gate 状态以 [roadmap.md](roadmap.md)（track 索引）+ [roadmap/](roadmap/)
> 各 track 文件（含各自的进度台账）为准；当前阶段以 [agent-rules.md](agent-rules.md)
> 为准；证据在 [demo_report.md](demo_report.md) 索引 + [demo_report/](demo_report/)。
> 若本文件与它们冲突，以它们为准。gate 关闭时按 CLAUDE.md 工作流更新台账后，顺手刷新本文件。

## 一句话状态

四条 track：求解器主线（P）与网格线（M）基本关闭，level-set 尾迹线（B）是当前
工作面。

**最新（2026-07-22）：B30/B31/B32 关闭——翼身跨声速能力实质推进。** conforming 翼身
跨声速 **medium 天花板 M0.79 → M0.84 达成**（B31 翼尖片终止奇点 cure + B32 生产采纳
tip_taper；cl_p 0.2738、0 钳制、代价 ≈ −1.3%）；(b) 类天花板归因（B30：两路径同机制 =
翼尖 P13 自由边奇点 + 高 M Newton，非某尾迹模型的袋）；LS 侧 C-class 实测阴性关闭
（B31 C1/C3，`outboard_fringe` 保留 default-inert）；B32 ② weld-sign per-step refresh
回滚（ill-posed switching）。B28（cl_fus 解耦 + GB9.4 重定规格，"机身虚假升力"标签退役）
+ B29（flat-fragment 升格翼身 LS 生产配置）于 2026-07-20 关闭。B31 新增的 Newton
Gamma-pin 行 blend 已 FD 验证（`test_blend_jacobian_fd_phi/_gamma`）。审计报告见
[inspection/20260722-0335-b28-b32-audit-pre-trackv.md](inspection/20260722-0335-b28-b32-audit-pre-trackv.md)。

**P11 已关闭（2026-07-19，用户指示当日开+关；sphere 腿）——曲面壁元路线
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
freeze_tol 抬到翼身 churn 地板 1e-6→1e-5，B17 教训）。**★ B32 勘误（2026-07-22）**：停滞归因 =
翼尖片终止奇点（B31），生产配方接入 tip_taper（vanish_smooth 0.05·b_semi）后 medium 爬升**达
M0.84**（cl_p 0.2738，0 钳制），代价 ≈ −1.3%，锚重钉 0.2143/0.2290/0.2450/0.2545 @ 0.50/0.65/0.75/0.79。level-set（B15 ramp+B17 pin_gamma）**不能**
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

- **P — 求解器**（[roadmap/track_p.md](roadmap/track_p.md)） — P0–P9 ✓（P1 仅 G1.6 以 strict xfail 挂起）；P10 ◐（G10.2/G10.3 ✓）；P13 ◐（G13.1 ✓、
  G13.2 conforming ✓、G13.3 亚声速 Richardson ✓ p=2.31） — G10.1（非升力 Newton 入口，无顺序约束）；G13.3 **跨声速阴性开放**（圆帽 fine 的 ramp 死于 M=0.75，
  site=尖 tip TE）；**P11 ✓ 2026-07-19 关闭（用户指示当日开+关；sphere 腿）**——G11.1 未达（曲面层 11.56%→11.33% = oracle 天花板；superparametric O(h) 风险触发）、
  G11.2 阴性+前提被驳（阶坍塌=固定中远场地板，E8 3.17×/1.89 阶；结构化壳平坦面片 ~2 阶）⇒ **G1.6 重归因为 P1 固有能力**，路线三岔口 **2026-07-22 裁决：Option C 重定规格 ADOPTED**（`TestG16Respec` PASS；字面 2%-max 保留 xfail）；P12 backlog；
  **P14 ✓ 2026-07-17 关闭（用户指示，当日开+关）**：壁面邻接 CV 压力相等 Kutta 估计器（A2 路由的修复），G14.1–G14.7 ✓、demo 28 PASS。**头条：
  修改后的 conforming 结果与 level-set 相同、计算结果合理**——S1/S2 一次换掉：M0.84 Γ(z) roughness 0.0970→**0.0043**（coarse）/ 0.0365→**0.0024**（medium，
  均达/优于 LS 带），全站 raw TE Cp gap 0.2206→**0.0040** / 0.1585→**0.0024**（**55×/67×**，走 G14.6 **主条款**，回退未动用；
  ★ 首版误引 A2 的 *section 末点* 数 0.318/0.228 报成 80×/95×，当日据用户提问勘误——A2 有两套 TE gap 度量，须引用自己实际跑的那套）。**★ 跨模型对照（V14.6，
  `cross_model_medium_m084.csv`）= 本相位最强证据**：LS 路径**一直**用压力相等 Kutta（B4），故若升力位移真是 Kutta *形式*，压力路径必须落到 LS 的答案上——它做到了：
  medium M0.84 conforming-pressure **0.2776/0.2823** vs level-set **0.2772/0.2813**，**差 0.17%/0.36%**（探针路径当年比 LS **低 4.5%/4.3%**），
  而且是不同尾迹模型、不同 DOF 空间、**不同网格族**（`onera_m6_wakefree`）。⇒ 长期存在的 conforming-vs-LS 升力分歧**就是** Kutta 形式误差，已消除。保留告诫：跨模型非同网格 A/B；
  LS 态带 1 lim/2 flr（B15 caveat）而压力态 0/0；"两者一致"≠"两者都对"（刚性平面尾迹等共有模型误差对两者是共模）。**★ G14.7 ✓ 关闭——开相位时锚探针 G8.2 锁、XFAIL 且不改带，
  用户裁决后改锚 level-set oracle**：medium M0.84 压力路径 cl_p/cl_KJ 0.2776/0.2823 vs LS 0.2772/0.2813 = **0.15%/0.34%**（<1% 通过）。
  相对旧探针锁的 +4.85% 位移是**发现**：机理在一级已测并在二级跑前**预注册**（两闭合逐点只差探针自身 O(h) 读数偏差，cross-read medium 0.79%；
  Kutta 映射 b≈0.93 把它 1/(1−b)≈14× 放大进 Γ），\|cl_KJ−0.288\| 0.0188→0.0057，**P9 那 0.019 gap 有 69% 是 Kutta 估计器偏差**（P9 看不见它：两套网格共用同一估计器，
  对其 Richardson 是共模）。关闭 G14.7 断言的是**两路一致**，**非**网格收敛（M6 fine 不是离散解）、**非**"0.019 是分辨率"翻案（仍 *strongly indicated, NOT earned*）。
  ★ **自我更正（V14.7 实测 2026-07-17）**：先前各版都断言 TE Cp **spike** 不会被 P14 触及（"两路共有的 P1 恢复伪影"）——那是从 A2 归因**推**的，没测。实测 medium M0.84 raw：
  探针 **0.1143** → 压力 **0.0533**（2.1×），且**低于 LS 的 0.0743**。A2 对的部分：确有共有残余（~0.05 = 真正的恢复地板）；需修正的部分：
  conforming 相对 LS 的**超出量**也是 Kutta 形式误差（Kutta 错⇒TE 流场真的错，末点偏离趋势有物理原因，共模度量分不开）。旁证：P6 平滑在压力路径上不再有效（0.0533→0.0660→0.0626；
  A2 在探针路径上测的是 0.147→0.081）。**教训：别把上一相位的归因当结论带进新测量——去测。** 诚实记录：判别器 D=7.33→**1.80**（落在 A2 的 inconclusive 区，非 O(1)）
- **M — 网格**（[roadmap/track_m.md](roadmap/track_m.md)） — M0、M1(+M1b 自相似阶梯)、M2、M3、M4、M5（圆顶翼尖盖）✓ — **M2 ✓（求解腿由 B9 于 2026-07-17 关闭；
  台账勘误 2026-07-19——A3 曾称已改而 ledger 行仍 ◐）**：翼身网格 2026-07-13 交付；**机身+远场 2026-07-16 按用户指示重定规格并重生成**（5 倍翼根弦长、机翼居中、2 倍直径椭球机鼻、
  蒙皮 h_body=2h_wall + 两端按半径加密；**R_FAR 15→25 MAC**；★需 `Mesh.OptimizeNetgen` 治 sliver 抽签）；遗留验证项（交界最内 TE 节点 CV fan）在 track_m 记录
- **B — level-set 尾迹**（[roadmap/track_b.md](roadmap/track_b.md)） — B1–B5、B7、B8（characterized-not-cured）、B9、B11–B32 ✓；
  B6 ◐（coarse gate ✓；medium 定量项由 GB15.4 补上，B21 恢复、B22 上锁） — **B16 ✓ 关闭 2026-07-17（用户指示，追加于 B15 之后；执行 B9 的 recorded follow-up）**：
  LS Newton 远场 BC 通用化——远场 aux DOF 钉扎。★ 翼身 LS-Newton churn 根因 = 近奇异远场 aux 块（尾迹片贯穿远场边界的 aux DOF 只受巨型外区单元的 wake-LS 行约束）：
  8 个远场 MAIN 行 max\|R\|=**84.457** 逐位复现，cond1 **9.1e18→8.70e6**（勘误 2026-07-19：旧文 6.36e18 是 CSV 前预跑值）。`farfield_aux="pin"`（默认，
  按模式适配）使 **coarse** 翼身 freestream Newton 达 **res 5.88e-14、0 limited**（升力≈conforming 0.1%），legacy churn 7.95/3690 limited；
  neumann 逐位不变。★★ **GB16.4 由 B17 解决（2026-07-18）**：不是非收敛，是 freestream pin 的 BC 建模错误（jump=0 杀出流环量）。⚠ 提案机理被自我修正（Picard 升力也用 wake_ls，
  差别在不动点吸收非 closure）。 · **B17 ✓ 关闭 2026-07-18（用户指示；解决 GB16.4）**：远场 pin 必须携带 jump=γ 而非 0——B16 的 jump=0 抹掉出流尾迹环量（medium −22%）；
  独立 Picard-pin 收敛到 Newton-pin "停滞" 的同一 0.169（两求解器同 BC 一致，非停滞）。`farfield_aux="pin_gamma"`（jump→γ，
  两求解器新默认）三角单调闭合向 conforming（coarse 0.2087、medium 0.2117 Picard/0.2114 Newton——pre-B20 0.2115）；仅作用 freestream、vortex/neumann 惰性；
  vortex 从 +2.5% bracket；B9 coarse-12.8% erratum 为远场污染。GB17.1–17.4 ✓、17.5/17.6 RECORDED；demo `cases/demo/b17_farfield_pin_gamma/`、
  tests `test_b17_farfield_pin_gamma.py`(6)。**B18 ✓ 关闭 2026-07-18（用户指示；执行 GB16.6 债）**：翼身跨声速 M0.84——conforming 到（coarse 0.2617、
  medium M0.79 0.2579、cl(M) 升 0.2173/0.2321/0.2579）、level-set 交界受限（close-out coarse ~M0.575、medium 死 ~M0.5；B20 重基线：
  coarse ~M0.55/Mmax 1.31、medium Mmax 5.22 真实未钳，袋=G1.6 几何而非 mixed-plain——GB20.5 勘误，closed-negative）；
  medium 无共同跨声速 Mach ⇒ 跨模型只留 M0.5（2.6%）+ coarse M0.6 点（0.2%，重基线工件）；GB18.1 PASS+18.2–5 RECORDED；零 pyfp3d/ 改动。
  ★★ **B26/B27 勘误（2026-07-20）**：上行"交界受限 closed-negative"叙事**已退役**——袋根因 = **B23 inboard 自由边奇点**（非 G1.6 刻面；G1.6 退居 cl_fus 嫌疑，
  C 侧 out-band ×2 → P11 输入），B25 `inboard_clip` 愈袋（corrM 14.66→0.63）；**B26-A：LS 天花板与 conforming 同址**（C 侧 = +clip：
  coarse 0.84 reached / medium 0.7625 死 0.775，(b) 类翼尖 P13 + 高 M Newton；A 侧复测爬过 B18 锚 = B21/B22 freeze-capture 效应，
  袋真实杀伤线 A medium 0.55 / Mmax 13.1）；**B27 demo 刷新**：conforming 逐位复现 B18 锚、LS A/C 逐位复现 B26 锚（g27_consistency.csv 336/336 bit），
  跨模型升档 **M0.5(2.6%) + M0.65(2.4% PASS ≤5%) + M0.75(2.5%)** ≈ 全 Mach 的 B17 口径差带，demo 8/8 PASS。**B29 ✓ 2026-07-20**：
  flat-fragment 升格为翼身 LS 生产配置（用户裁决 B28 §6）——B18 C 侧 = clip+平片（`sheet_direction=(1,0,0)`），M0.5 锚重钉
  0.2115/0.2184；medium 天花板 0.7625→**0.775**（垂死峰 M3.98 @ 翼尖实况）；跨模型 **0.5/1.1/1.1%**（M0.5/0.65/0.75，原
  2.6/2.4/2.5%）；GB18.5 实况平片分解 cl_fus **0.0382**（带 −0.0006/带外 0.0388/极 0.0007）@0.7875 vs conf 0.0423——B26
  斜片 ×2 读数退役（B28 位置敏感性），demo 8/8 PASS。**B19 ✓ 关闭 2026-07-18**：
  LS-Newton Jacobian 3-D 精确化（两处缺项；探针 1.146e-01→1.33e-08、ε 判别翻转；R 逐位不变、无收敛收益——GB19.4 阴性）；
  Leg B 测得 mixed-plain 侧场密度污染（虚假超声速 q² 3.22 vs 1.34）→ 路由 B20。**B20 ✓ 关闭 2026-07-18、永久采纳（用户裁决，开关移除）、重基线 2026-07-19**：
  mixed-plain 密度改读 main 场；所有移动数字向好，唯一回退 = M6 medium ramp（GB20.7 曾定谳"真实能力损失"）。**B21 ✓ 关闭 2026-07-19（执行 Kimi/N1）**：
  `freeze_side_state` 漏打 B20 补丁——修复后 **M6 medium M0.84 恢复**（γ 0.088343、res 9e-14、515 s），GB20.7 翻案；3-D freeze 捕获一致性测试锁落地。
  **B22 ✓ 关闭 2026-07-19（执行 B21 follow-up + Kimi N3/§2/§5）**：B15 demo 20/20 + B14 7/7 刷新（coarse ramp γ 0.084931 披露）；**N3 关闭**——
  gated 绝对锚锁 `test_b22_ls_3d_anchors.py`；re-baseline 勘误清单入流程；下一相位分析推荐 P11（用户裁决）。**B23 ✓ 关闭 2026-07-19（预注册判别）**：
  翼身交界袋 = **尾迹 inboard 自由边奇点**（升力/尾流耦合——α=0 两级干净、袋随 α 超线性增长、峰在机身后方 x=2.13；非 G1.6 刻面），P11 close-out 判别落地。
  **B24 ✓ 关闭 2026-07-19（阴性）**：袋跟自由边走（机制二次证实），但水线延伸两变体（B1 贴面/B3 离锥）都换出同等或更糟奇点（B1 medium corrM 78.56 非收敛）——判定树出口 3，(b)-1 路线关闭。
  **B25 ✓ 关闭 2026-07-19（治愈 C-A）**：`inboard_clip`（`meshgen/fuselage.py:make_inboard_clip` + `wake/cut_elements.py`）把片内侧边界移到机身面/对称面 
  = conforming fragment 拓扑——medium α=3.06 走廊 corrM **14.66→0.63**、n_sup 88→0、cl_p +0.38%（∈[A, oracle]），全物理护栏干净；默认 None 逐位不变；
  次级护栏带外 cl_fus +135% 经 oracle 归因 flat-vs-tilted 片模型差，记录不阻塞 → P11 监视。**B9 ✓ 关闭 2026-07-17（重定规格）**：
  翼身跨模型 LS+conforming 一致 0.4%/0.6%；GB9.4 XFAIL⇒G1.6（**B28 2026-07-20 更正**：cl_fus=尾流片位置敏感性，
  非 G1.6 误差；gate 重设为带外跨模型一致 ≤15%，medium 差 7.0% PASS，demo 8/8）。B10 搁置
- **V — 粘性耦合**（[roadmap/track_v.md](roadmap/track_v.md)） — 设计完整（Drela IBL3 + transpiration BC）；**V1 ✓ CLOSED
  2026-07-22**（GV1.1 9 PASS / 2 FAIL，(a)×2 = 闭包族不动点物理，recorded FAIL 接受）；**V2 ✓ CLOSED 2026-07-22**
  （GV2.1 23 PASS / 0 FAIL：cylinder Fourier blowing 对解析 relmax 严格降、阶 1.650/1.640，ṁ=0 五路驱动逐位一致，
  lagged ṁ 下 Newton Jacobian 逐位不变 + FD 6.6e-09–7.2e-08——transpiration 通道三路驱动全部打通）；**V3 ✓ CLOSED
  2026-07-22**（GV3.1/3.2 2 PASS / 4 FAIL / 23 RECORDED · GV3.3 0 PASS / 2 FAIL / 7 RECORDED：`viscous/coupling.py`
  松耦合 + committed XFOIL 参考；PASS Δcl 比 0.542∈[0.5,2.0]、松环 4–5 次外迭代 ω=1.0（跨声速 M0.72 记录点 4 次
  无调参）；honest FAIL 局域化：cf 仅转捩后首站 +44%（XFOIL e^N 斜坡 vs 瞬时切换）、δ* H 族偏移 ≤27.9%；GV3.3
  旋成体三轮调试稳定化（尾带钉扎 + 钉扎带 ṁ 掩蔽 + FP 护栏），中段轴对称优秀、尾锥 σ/μ 0.55/横流 0.26 FAIL、
  环不收敛 = 实测尾部失稳——V4 跳过判据按字面满足（GV3.2），GV3.3 反方证据入台账，**V4 ⊘ 跳过
  2026-07-22（用户定；重开触发 = V5 受挫或闭体粘性提前进范围）**）；**V5 ◐ OPEN 2026-07-23**
  （**GV5.0 ✓ EXECUTED 16 RECORDED / 0 FAIL**，`cases/analysis/v5_m6_bridge/`：M6 亚声速松耦合桥——
  桥答案 = 松环在 3-D 升力翼上不够用：coarse 根部上翼面 TE 分离斑块（H 4–5.5）δ*↔ṁ↔u_e 失稳
  ṁ_max ×12.4（GV3.3 尾部同类），medium 加密消除斑块但留有界 δ* 极限环（2–12 %/k）不达 tol 1e-3；
  ΔCL 双估计量下行（coarse −5.2 %/−4.8 %，medium −2.4 %/−2.1 % 输入受限）；横流首次活体 3-D 演习
  max|B|/|A| ≤ 0.072；翼尖掩蔽有效；新增 `viscous/coupling.py::build_wing_case` +
  `tests/test_v5_wing_case.py` (5)；δ*(z) CSV 喂 GV5.3 带预注册；medium 壁时被外部负载污染，引用须带旗标）
  （gate 按 B32/A4 现状重定规格，同日三分重排：V1 独立 IBL3 核心（GV1.1 解析/自相似对标）· V2 transpiration 通道
  （GV2.1 精确性 + ṁ=0 逐位 + FD）· V3 松耦合（GV3.1 NACA0012 对 committed XFOIL 引 A4 输入带 · GV3.2 松耦合 ≤10 次 →
  V4 跳过判据 · GV3.3 机身旋成体冒烟，唯一机身-alone 项）；V4 ⊘ 跳过 2026-07-22（原可选 quasi-simultaneous，判据满足）；V5 紧耦合 ◐ OPEN（入口 GV5.0 M6 亚声速
  松耦合桥 ✓ EXECUTED 2026-07-23 16R/0F——桥答案 = 松环在 3-D 升力翼上不够，紧耦合动机证据；**GV5.1 ✓ EXECUTED 9P/1F/36R**——精确增广 (φ, Γ, BL) Newton 交付
  且两级 FD 验证（worst 甜点 2.2e-8 coarse / 5.1e-9 medium），二次尾段 HONEST FAIL = IBL 稳态残差在 cond(J_BL,BL)~4e10 近零流形上的内禀地板（standalone 伪时间
  同地板，非紧耦合缺陷），N_total 14/13 vs 松环 4/5；已提交 GV3.1 medium 不动点不可复现（IBL 地板轨迹散布，诊断已提交，HEAD 重生成种子经用户裁决接受）；**IBL 地板 follow-up 诊断
  ✓ EXECUTED 2026-07-24**（14 RECORDED 无 band，`cases/analysis/v5_ibl_floor/`：近零簇在松环收敛态持续（S1 500/1236 cond 1.3e11、S2 1082/2460 cond 4.0e13、s1/s3 谱几乎重合），由湍流 (A,Ψ)
  变量承载、mid-chord→TE 分布；原始 cond 4e10–4e13 主要是缩放 artifact（行列均衡 → 2e4/7e5/1e7、亚 1e-6 计数 501/500/1082 → 0/0/2、无精确零方向），均衡后余 1e5–1e7 真 (A,Ψ)
  刚度 = GV5.1b/GV5.4 的真正靶子；F_BL 地板住 TE 带 (B,δ) 方程且几乎全在 J 值域内；闭包地板活动集全空（假说死）、eps_diff ×4 地板仅移 ≤6 %（非人工粘性截断）、伪时间控制器触底
  （cfl 钉 cfl_min、残差自第 0 迭代冻结）= 公式化地板经控制器表现，单靠全局化过不了地板）；**GV5.1b ✓ EXECUTED 2026-07-24**
  （2P/0F/7R 裁决后（执行时 1P/1F/7R 保留于 commit 1c55906），`cases/analysis/v5_1b_scaled_newton/`：scaled+damped 机构交付且精确——求解器内部行列均衡 + Levenberg 对角阻尼 +
  floor-reached 停止类，旗标默认关 = legacy 逐位，套件 28 green；medium 活体 e2 恒等式 1.96e-10 超自设 ≤1e-10 阈值
  = cond~1e10 下 SuperLU 主序舍入的机器地板（阈值非预注册，非代数错）——**2026-07-24 用户裁决 cond-aware 读 PASS**
  （tol = max(1e-10, 10·κ₁·eps)，κ₁~1e10 → ~1e-5 量级，4 个 decade 余量，VERDICT §3；run.py 阈值改由 κ₁ 一范数估计现算）；amended 种子自第 0 迭代即坐在
  10× 地板带内侧（F_BL = 1.00× 地板）⇒ 构造上无 above-band 收缩段，走预注册 fallback：medium floor_reached 第 5 迭代
  同 merit 收官（9.074e-11 ≈ 9.025e-11，取代 GV5.1 的 10 步 λ-collapse 爬行）、coarse 末 merit 2.044e-10 < GV5.1
  2.068e-10 仍在降、k=1 standalone F_BL −31 %/merit 2.3× 更深、μ 拒绝重试 0 次——缩放是活性成分，阻尼臂惰性；
  窗口问题被重构为需 above-band 种子的协议 → **GV5.1c ✓ EXECUTED 2026-07-24**（2P/1F/7R，
  `cases/analysis/v5_1c_above_band_window/`：标定 above-band δ 扰动种子（ε = 1e4 → 种子 F_BL ≈ 1e4× 地板带）真正读
  地板前 slope-2 窗——**地板之上处处无二次收缩**：干净下降段全是线搜索封顶折半（λ = 0.5 → p = 1.00 构造值），
  随后中程停滞（F_BL ~ 3e-2 → 1.3e-2/2.2e-2，10 迭代内从未进带，距地板 4262×/12867×），binding medium
  median p = 0.56 honest FAIL；紧 Newton 的障碍不止地板——其上 3–4 个 decade 还有一道中程下降屏障；
  近带种子是否有局部二次盆 = 后续问题 → **GV5.1d ✓ EXECUTED 2026-07-24**（2P/1F/7R，
  `cases/analysis/v5_1d_near_band_window/`：近带种子 T1 = [1e-4, 1e-3]（coarse 5.4×/medium 35× 带）
  ——**地板紧邻处也无二次盆**：coarse 一次封顶折半后爬行（λ → 6e-5）至 24× 地板仍未进带，medium 首个
  接受步把 F_BL 推离带（6.0e-4 → 9.8e-4，merit 靠块再平衡换得）后爬行至 493×；binding medium
  median p = 1.17 honest FAIL；μ 拒绝重试第三次为 0——平坦/锯齿 merit 邻域向下延伸到距地板
  ~1.5 个 decade 内，盆地搜寻穷尽（GV5.1b/1c/1d），**GV5.5 成为破地板唯一在册路线**）；μ 拒绝重试再为 0；band (a) PASS 两腿
  （cond-aware e2 容差本次预注册）；8 线程临时约束下执行（runner 默认 16，壁时标记）；medium 不动点在
  8 线程下再次散布（第 4 个不动点 cl 0.28245999，coarse 逐位一致））；破地板本身登记为独立项
  **GV5.5 TE 带 (B,δ) 公式层**（2026-07-24 用户定，未开工，排序待裁决）；GV5.2/5.3/5.4/5.5 排序待
  用户裁决；GV5.3 锚定 committed Cp——实验 CL 无 committed 来源）；
  V6 尾迹面片；翼身 VII 延后至 LS 侧翼尖 cure）— 依赖 P6+A4（均已满足），预算等同一个 Track-P 阶段，V4–V6 尚无实现。
  参考文献在手：Drela 2013 = AIAA 2013-2437（`docs/references/` 本地，gitignored）
- **A — 校验与分析**（[roadmap/track_a.md](roadmap/track_a.md)） — 2026-07-15 新建；**A1 ✓ 2026-07-16**（GA1.1–GA1.5：
  四求解器统一计时插桩 + conforming×level-set × Picard×Newton 耗时基准） — **A2 ✓ 2026-07-17 关闭**（TE/Kutta 保真度归因，GA2.1–GA2.5）：**S1 定谳**——
  conforming Γ(z) 逐站抖动是逐站探针差势跳 Kutta target **估计器**的测量伪影（fixed-Γ 判别量 D=7.33/25.70 coarse/medium，把抖动从光滑场里重新生出来；闭合残差 ≤0.6% 排除"未闭合"、
  抖动局域于 TE 邻层 0.02–0.07× 排除"流场"），**非流场内容**；**S2 分解**——TE Cp 突跳=势跳 Kutta 形式误差（conforming 独有,同估计器 34×/133× vs LS）+ P1 末点恢复伪影（两路共有）；
  2.5-D `a1_cp.png`（有 S2 无 S1）证明两者是不同机制。修复路由至 **P14**（无 `pyfp3d/` 改动）；工作目录 `cases/analysis/`（区别于 `cases/demo/`）。★A1 结论：
  3-D 下两条 Newton 路径都是 **precond（LU 分解）受限**（~40% 墙钟，lagged LU 已开），2.5-D 的"seed 是成本"**不外推**——引用主导相位必须带网格。
  **A3 ✓ 2026-07-18 关闭**（GA3.1–GA3.6：响应 Kimi 独立审查）：★★ **C1 证实**——LS Newton Jacobian 在 mixed-side plain 元素上≠dR/dφ（targeted 
  1.146e-01 vs control 6.33e-10，且 eps 无关 max/min=1.00）⇒ 3-D 下是 quasi-Newton；**R 未动 ⇒ 所有已收敛状态与 gate 数字成立**；记录不修，修复单独立项。
  C2/C3 回移 conforming Newton、reader C4/C5（静默丢面组 ⇒ Γ(root) 静默为 0）、C6/C7/P1/T1/T2/F0；文档 17 项全部处置；close-out 扩为五面 + backport check；
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

现基线 **627 passed + 25 skipped + 2 xfailed**（2026-07-24 Track V **V5 GV5.1d
执行**（近带种子读地板紧邻处二次盆：同样无盆——近带种子立即停滞，coarse 爬至
24× 地板未进带，medium 首步推离带；binding medium median p = 1.17 honest FAIL——
VERDICT `cases/analysis/v5_1d_near_band_window/VERDICT.md`）：全套件实测 627
@1340.77 s **@8 线程**（本 session 临时 8 核约束，用户定；与 16 线程账目不可直接
比；壁时明显低于同线程数 GV5.1c 账目的 3903 s——机器/缓存条件不同，标记引用）；
+7 vs 下档 620 = `tests/test_v5_near_band_seed.py`（7））。
上一档 620+25+2（2026-07-24 Track V **V5 GV5.1c
执行**（above-band 种子读地板前 slope-2 窗：地板之上无二次收缩——λ 封顶折半 +
中程停滞，binding medium median p = 0.56 honest FAIL——VERDICT
`cases/analysis/v5_1c_above_band_window/VERDICT.md`）：全套件实测 620 @3903.16 s
**@8 线程**（本 session 临时 8 核约束，用户定；与 16 线程账目不可直接比，
机器空载）；+9 vs 下档 611 = `tests/test_v5_above_band_seed.py`（9））。
上一档 611+25+2（2026-07-24 Track V **V5 GV5.1b
执行**（scaled+damped 增广 Newton：机构精确交付，band (b) 窗口问题重构——VERDICT
`cases/analysis/v5_1b_scaled_newton/VERDICT.md`）：全套件实测 611 @6556.77 s
@16 线程（wall 受合租负载 ~70–80 污染，标记引用；GV5.1 时空载为 1537 s）；
+8 vs 下档 603 = `tests/test_v5_tight_scaled.py`（8））。
上一档 603+25+2（2026-07-23 Track V **V5 GV5.1 执行**（增广紧耦合 (φ, Γ, U)
Newton：band (a) FD 精确性两级 PASS，band (b) 二次尾段被 IBL 地板挡住 HONEST
FAIL——VERDICT `cases/analysis/v5_tight_coupling/VERDICT.md`）：全套件实测
603 @1537.09 s @16 线程；+20 vs 下档 583 = `tests/test_v5_tight_jacobian.py`（8）+
`tests/test_v5_tight_edge.py`（7）+ `tests/test_v5_tight_system.py`（5））。
上一档 583+25+2（2026-07-23 Track V **V5 GV5.0 执行**（M6 亚声速松耦合桥，
RECORDED 入口检查）：全套件实测 583 @1218.05 s @16 线程；+5 vs 下档 578 =
`tests/test_v5_wing_case.py`（5））。
上一档 578+25+2（2026-07-22 Track V **V3 松耦合交付 + GV3.1/3.2/3.3 执行**：
全套件实测 578 @1637.39 s @16 线程；+7 vs 下档 571 =
`tests/test_v3_coupling.py`（7））。
上一档 571+25+2（2026-07-22 Track V **V2 transpiration 通道 + GV2.1**：实测
571 @1321.89 s；+17 vs 554 = `tests/test_v2_transpiration.py`（9）+
`tests/test_v2_newton_rhs_channel.py`（8）；NOJIT 路 17/17 绿）。
上一档 554+25+2（2026-07-22 Track V **V1 IBL3
core 交付 + GV1.1 执行**：全套件实测 554 @1462.64 s @16 线程；+35 vs 下档 519 =
`tests/test_v1_surface_mesh.py`（13）+ `tests/test_v1_closures.py`（17）
+ `tests/test_v1_ibl3.py`（5）；NOJIT 路 35/35 绿）。
上一档 519+25+2（2026-07-22 B28–B32 收尾 **+ G1.6
Option C 重定规格**：全套件实测 516 @1223.39 s，+ 3 条 `test_laplace_sphere.py::TestG16Respec`
断言（读 P11 已提交 sweep CSV、无交互）= 519。516 明细：B28–B32 收尾：
+37 passed vs B25 的 479 = B28 cut-from-fragment 锁（`test_b1_cut_elements.py`，+4）
+ B31 `test_b31_pressure_taper.py`（13）+ `test_b31_tip_fringe.py`（19）
+ `test_p14_te_pressure.py` 锁（1）；B29/B30/B32 未加测试；实测 1223.39 s @16 线程）。
上一档 479+25+2（2026-07-20 B25 inboard fragment clip：+6 passed =
`test_b1_cut_elements.py::TestInboardFragmentClip`（4）+ 同文件 foot-preference 锁（1）
+ `test_m2_wingbody.py` 水线延伸锁（1）；实测 1100.63 s @16 线程）。
上一档 473+25+2（2026-07-19 P11 曲面壁元：+8 passed = ungated
`tests/test_p11_curved_walls.py`；实测 1124.94 s @16 线程）。
上一档 465+25+2（2026-07-19 B22 3-D LS 锚锁：
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
`test_p11_curved_walls.py`；实测 1124.94 s @16 线程）**
→ **479+25+2（B25 inboard fragment clip，2026-07-20：+6 passed =
`test_b1_cut_elements.py::TestInboardFragmentClip`（4）+ foot-preference 锁（1）+
`test_m2_wingbody.py` 水线延伸锁（1）；实测 1100.63 s @16 线程）**。

## 长期挂起项（勿反复重提）

- **G1.6 球面 Cp <2%**（strict xfail，11.6%）：★★ **根因已于 2026-07-19 由 P11
  重归因**——旧说"平坦面片壁的自然边界条件（变分罪）"被实测推翻（罪的份额
  ≈0.2pp；11.6% ≈ P1 场在 h=0.08 的固有能力；阶坍塌是混淆扫掠的中远场地板）。
  仍然成立：恢复非主导、Nitsche 死、边界数据修正无可修正、h 加密（同族）贵而不达。
  新增死路：mapped-P1（superparametric）曲面壁元。★★ **路线三岔口已于 2026-07-22
  裁决（用户指示）：(a) Option C 重定规格 ADOPTED**——G1.6 活跃 gate 改为可达实测标准
  （全尺度加密族 φ_w 阶 ≥1.8 + h_min 0.03 mean-Cp <1%；E6/E8 实测 1.98/1.89 + 0.60%），
  由 `test_laplace_sphere.py::TestG16Respec` 读 P11 已提交 sweep 断言 **PASS**；字面
  2%-max@medium xfail **保留 = 记录的 P1 边界**（需 O(h²) 壁速 @h=0.08，超任何 P1 场方法）。
  未取路线：(b) isoparametric P2 壁层（唯一触及字面标准，且会收紧 Track V 的 u_e 输入带 A4）
  / (c) 接受为长期边界。
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
