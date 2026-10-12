# `phases/` —— 已结束阶段的归档

这里放的是 **phase 1–6** 的材料。当前阶段的工作面在仓库原位置
（`pyfp3d/`、`tests/`、`bench/`、`cases/`、`docs/`）。

**归档是目录重组，不是取消证据**：所有归档文件仍被 git 跟踪，已提交的 CSV/PNG 都还在原位
（归档里）。要一个**保证能跑**的历史树，用恢复点而不是从归档里跑脚本：

```bash
git worktree add ../up3d-prereorg d224223     # 2026-08-10 大整理之前（phase 1/2 全貌）
git worktree add ../up3d-prebench 1dc4647     # 2026-08-24 bench 大归档之前
git worktree add ../up3d-prep6    1544c8c     # 2026-10-12 phase 6 归档之前
```

**归档脚本的可运行性如实说**：移动会打断脚本找路的方式（深度表达式、拼出来的路径、
写死的绝对路径）。历史轮次修过能机械验证的三类，但**没有逐一验证到「能真的跑起来」**。
⇒ 证据是已提交的 CSV/PNG（纪律：没有提交产物的说法不算证据）；脚本是出处（provenance）。

---

## 各阶段一览

| 阶段 | 时间 | 一句话 | 归档 |
|---|---|---|---|
| p1 | 2026-07-06 → 07-25 | 从零建起求解器：五条 track（P 求解 / M 网格 / B level-set 尾迹 / V 粘性 / A 校验） | `p1/` |
| p2 | 2026-07-28 → 08-10 | 纲领重建：五条产品指标、S0–S6、激波台 S1、熵修正、第一版能力边界 | `p2/` |
| p3 | 2026-08-10 → 08-16 | 删 level-set（D5）；hex 剖分路线判死；M0.80 多解实测；GS4.0 审计修复 | `p3/` |
| p4 | 2026-08-18 → 08-20 | GS4.1：二维 strip IBL / 剪应力滞后（粘性重启的第一块） | `p4/` |
| p5 | 2026-08-21 → 08-23 | M1 二维无粘标定（网格收敛阶梯、产品指标第一次设门） | `p5/` |
| p6 | 2026-08-24 → 10-12 | 全面梳理 + 门重建：门号治理、CFL3D/XFOIL 对标、BLAS 确定论、能力边界换代 | `p6/` |

---

## phase 1（`p1/`）—— 建 solver

**做了什么**：从空仓库到可用的三维全势流求解器。Track P（求解器：Picard → Newton、
Mach 延拓、人工密度迎风、Kutta/尾迹、P14 壁压相等 Kutta 估计器）、Track M（网格：
挤出 2.5D、翼身、M6）、Track B（level-set 尾迹线 B1–B32，含翼身跨声速 M0.84）、
Track V（粘性耦合 V1–V6：IBL3 核心、transpiration 通道、松耦合、紧耦合增广 Newton、
wake sheet）、Track A（校验：A1 计时、A2 TE/Kutta 归因、A3 独立审查响应）。

**主要问题（已结论的）**：G1.6 球面 Cp 误差 = P1 场在 h=0.08 的固有能力（P11 重归因，
Option C 重定规格采纳）；conforming-vs-LS 升力分歧 = Kutta 形式误差（P14 消除）；
翼尖片终止奇点（B31 cure + B32 tip_taper）；松耦合回路 run-to-run 不可复现。

**归档内容**：`cases/`（24 demo + 23 analysis 证据链）、`docs/`（roadmap + 五 track
台账、demo_report 证据档案、design_track_b.md、inspection 审查报告、archive）、
**18 个纯 level-set 测试**（唯一能归档的测试——它们那条路后来被删除了）。

## phase 2（`p2/`）—— 纲领与激波台

**做了什么**：2026-07-28 奠基审计（无粘亚声速健康；跨声速激波位置由调参和网格决定）
⇒ 立第二阶段纲领：五条产品指标 M1–M5、八条原则、S0–S6 阶段划分。S0 地基（wheel 56 文件、
bench 基线 + bitcheck 逐位工具）、S1 激波判据台（Laval 喷管精确解）、GS1b.3 熵修正实现
（三级 cl 跳幅 +40.5% → +13.8%）、预条件子决策（DECISION-2026-08-02，EW forcing 1e-6
至今是库默认）、**第一版能力边界**（PHASE_TWO_CAPABILITY_BOUNDARY.md——仍是 phase 2 那
些声明的唯一原件）。

**主要问题**：M1 如实 FAIL（激波 0.6581 带外、两级 cl +101%）——可信包线实测 M_max ≲ 1.2，
M1 工况落在包线外；鲁棒性随加密变差；PTC 被否证并删除；σ 光滑化、σ 绑 freeze 两条路线
均为实测负结果（kill criterion 照章）。

**归档内容**：`bench/`（45 个 Track A/V studies 脚本 + gate_results/capability）、
`docs/dev_phase_two/`（92 个文件：progress.md、roadmap.md（裁决 D1–D5）、
PHASE_TWO_CAPABILITY_BOUNDARY.md、DECISION-2026-08-02-precond.md、**_TEMPLATE.md**
（轮次格式模板，至今在用）、README.md）、2026-07-28 奠基审计 + exp1–exp6。

## phase 3（`p3/`）—— 删路线、判死路线、多解

**做了什么**：执行裁决 D5 **删除 level-set 路线**（9 库文件 / −4624 行，账目逐位闭合）；
六面体→四面体剖分路线 (A) 两条判据 FAIL、kill criterion 成立（壁面精度由 O(h) 刻面几何
决定，与网格结构无关；熵修正的供体确定性问题**本来不存在**）；加密悖论调查答成了
**M0.80/α1.25 离散解不唯一**（固定网格只换种子，x_shock 撑开 0.168c）；GS4.0 =
2026-08-16 独立审计的六条活缺陷修复。

**归档内容**：`bench/`（run_task3_* 26 + run_hex_* 3 + run_gs40* 10）、
`docs/dev_phase_three/`（68 个轮次文件）、`docs/inspection/`（2026-08-16 独立审计，
触发 GS4.0）。

## phase 4（`p4/`）—— 二维 strip IBL

**做了什么**：GS4.1：二维 strip IBL / 剪应力滞后——粘性重启的第一块（S4 的进入前置
GS4.0 在 p3）。二维条带积分边界层求解器、剪应力滞后闭包、对 XFOIL 的标定腿。

**归档内容**：`bench/studies/gs41_*`（14 个 studies）、`docs/dev_phase_four/`（82 个轮次文件）。

## phase 5（`p5/`）—— M1 二维无粘标定

**做了什么**：M1 产品指标的二维无粘标定：网格收敛阶梯、三级 cl 一致性的第一次设门尝试。
27 轮（R11–R23、M1a/M1b/M1e、N1/N2）。

**读它之前必读**：`p5/docs/dev_phase_five/progress.md` 开头的 **§A 仍然成立 / §B 已撤回**
两张表——引用 phase 5 任何数字之前先读那两张表。

**归档内容**：`bench/studies/`（r11..r22、m1a/m1b/m1e、n1/n2 共 17 个）、
`docs/dev_phase_five/`（61 个文件，含 progress.md）。

## phase 6（`p6/`）—— 全面梳理 + 门重建

**做了什么**（使用者 2026-08-24 裁决立项：梳理代码、文档、测试、算例，为测试案例建合理的
门；50 轮 + 两次独立审计）：

- **门号治理**：AST 逐条读完 556 条断言做归类；测试全量重编号进 `tests/A..F`（A/B 未门控
  机制、C/E 快层、D 门控全集、F 仪器）；证据锁机器 `tests/_gate_evidence.py`（16/19 证据
  CSV 曾只写不读，现已新鲜复算对已提交 CSV）；删掉 gitignored 的 `artifacts/`（它作为证据
  形式自相矛盾），C/D 门证据落到被跟踪的 `cases/gates/<gate_id>/`。
- **CFL3D / XFOIL 对标**：三层参照按模型层级对齐（无粘↔CFL3D Euler 设门 D05/D06/D07、
  有粘↔RANS 设门 D08、对实验记录偏置不设门 D11/D12、IBL↔XFOIL D13）；M6 结构化 C-H 网格
  生成器交付；查出两条 M6 阶梯全程跑的是一阶（NITFO≥NCYC 缺陷）；加第四档翻掉 31% 判定
  ⇒ 判定改三分制（yes/no/unstable/undefined）。
- **能力边界换代**：`p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md`（2026-10-11）
  取代 phase 2 原件成为**现行能力声明**。三个实测真边界：α=0 可信上界 M 0.78–0.80（D14，
  Steinhoff–Jameson 非唯一）；有粘耦合在举力机翼上跑不起来（D09/D10）；生产配置在精确解
  上的已登记代价（C05）。
- **确定论**：H30 —— BLAS 由库运行时钉单线程（pytest 里设 env 不生效，OpenBLAS 只在
  numpy 首载时读）；numba 编译缓存出处态是第三个环境变量 ⇒ D08/D13 锁拆两层 + 暖态纪律
  （重编译后先白跑再测）。
- **两个归因研究**：C05 上 amg 不收敛 = AMG 层级建在不含激波耦合的 Picard 块上
  （激波腿不该用 amg 当默认）；hex 族 NACA medium 全收敛下 cl 比 tet 低 24% = 收敛到
  tet 从未出现的双激波/气泡结构（残余：生成器缺陷 vs 熵冻结交互未分开）。
- 基线：开工 595/12/2 → 收口 **831 passed / 33 skipped / 2 xfailed / 0 failed**（16 线程、
  BLAS=1）；门控全集 812/21/3/1 XPASS；快层 5/5。

**遗留（登记为「留到下一阶段」，见能力边界 §9 与收口判定 §5）**：M3b 速度指标锁、
W3.1 G-TEETH 定量化、W3.3 ZPG 湍流平板锚、W3.4 非线性 MMS、W3.6 E01/E03 升格评估、
hex 族能力边界门（红=好消息式）、C05 预条件子文书 + 「amg 能否修」研究、xtr030 多吸引子、
tet medium 在 r_far 20/30 重生成变体上不收敛。

**归档内容**：`docs/dev_phase_six/`（54 个文件：立项、承接台账、两次审计、50 轮
prereg/verdict、PHASE_SIX_CAPABILITY_BOUNDARY.md、progress.md 轮次总账）、
`bench/studies/`（4 个一次性研究：ibl_xfoil_recon、gate_audit_20260905、
c05_amg_investigation、hex_cl_attribution）、`bench/gate_results/`（m3a_smoothing_ab/
shock.csv，仅被 m3a 勘误文书引用）、`docs/`（**phase 6 结束时退役的四份纲领文档旧版**：
`agent-rules.md`、`design.md`——一直维护到 phase 6 的活文档的归档；另两份内容停在
phase 1 的（`overview.md`、`design_track_v.md`）按内容归属入 `p1/docs/`）。

---

## 归档口径（量出来的规则）

1. **相位归属 = 文件首次提交日期**（边界 2026-07-28 = phase two 纲领建立日）。
2. **留在原位的判据是可测的**：如果一个仍然存在的测试或脚本读它，它就留在原位。
   注意**传递闭包**（一跳规则曾漏掉 `v3_loose_coupling`/`v5_1b_scaled_newton`，全套复跑
   才抓到），且路径可以是**逐段拼出来的**——只认字面量的判据不完整。
3. **测试代码不归档**（phase 2 起）：每一个测试都在检验活的 `pyfp3d/` 代码；把一条活的
   收口门指进归档比没有这道门更糟（它会跑、会绿、并为一条没人维护的路线断言能力）。
   phase 1 的 18 个测试能归档，是因为它们那条路（level-set）被删除了。
4. **归档保持旧布局是故意的**：归档是历史快照，不追活的一侧的改名（活的在
   `bench/studies/`，历史的在 `phases/p1/cases/analysis/`）。
5. **搬迁必须修两类断裂**：根相对路径字面量（绑定文本行）与相对 markdown 链接
   （`](../…)` 按被引用文件自己的目录解析，搬迁后深度全变）。方法学顺序：**先取基线 →
   `git mv` → 先证明仪器有牙 → 改写 → 复验回基线**；替换前必须排除「已落在正确路径里」
   的出现（朴素子串替换曾造出 37 处 `phases/p1/phases/p1/…`）。
6. **归档脚本里的引用修复**只做能机械验证的，其余如实登记不装死——见顶部可运行性口径。

## 本次（p6）归档做了什么（2026-10-12）

- `git mv`：`docs/dev_phase_six/`（54 文件）→ `p6/docs/dev_phase_six/`；
  `bench/studies/` 的 4 个一次性研究 → `p6/bench/studies/`（`bench/studies/README.md`
  自己点名的「不被测试加载」二分；其余 7 个被活测试运行时加载，留下）；
  `bench/gate_results/m3a_smoothing_{ab,shock}.csv` → `p6/bench/gate_results/`；
  退役纲领文档：`agent-rules.md`、`design.md` → `p6/docs/`，`overview.md`、
  `design_track_v.md` → `p1/docs/`（内容归属）。
- 引用修复：活树 55 处 `docs/dev_phase_six/` 根相对引用重指；归档内 13 处 `](../../`
  链接深度 +4 个研究的证据指针；`pyfp3d/viscous/` 的 `design_track_v.md` § 引用重指到
  `p1/docs/` 归档原件；`bench/studies/README.md` 的一次性名单清零。
- **留在原位**（有活引用）：`tests/` 全部（含 F07–F10 治理锁——它们断言活仓库状态）、
  `cases/gates/`（活测试证据目录）、`cases/reference_data/cfl3d/`（D05–D10/F08 运行时读）、
  两个新网格族（rae2822_2.5d、cylinder_hex_2.5d）、`bench/` 根的 12+7 结构与
  `bench/gate_results/` 其余、`pyfp3d/` 全部。
- `tools/`（CFL3D/XFOIL 工具链）从不被 git 跟踪，无档可归。
