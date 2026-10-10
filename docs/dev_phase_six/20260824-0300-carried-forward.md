# 20260824-0300 — **承接台账**：`docs/` 清空之后，**仍然在管事的东西在这里**

**性质**：承接文书，**phase 6 唯一的"当前状态"落点**。
使用者裁决 2026-08-24：**phase 2–5 的文档全部移入 `phases/p*/`，`docs/inspection/` 删除；
phase 6 在自己的文档里指向或跟踪它们 —— 以前的分析和文档都只是参考。**

★★★ **这份文件为什么必须存在**：被归档的 6 个 phase-2 文件里，有 **4 个当时是"活的"**
（`phases/p2/docs/dev_phase_two/README.md` 那份 2026-08-10 检查单挡住了它们，理由有两条是**可测的事实**）。
裁决把它们变成参考，**那么它们承载的"当前状态"必须有一个新的落点**，否则就不是归档而是丢失。
**这里就是那个落点。**

★ 口径：**下面每一条要么是我实测的，要么明确标为"引用归档"**。
不从被归档的文档里抄数字 —— 那正是本项目记过多次的「把过期的数字改成新数字，下次还会过期」。

---

## 1. ★★★ 库的当前默认值 —— **实测，2026-08-24**

用 `inspect.signature` 直接读出来的，不是抄的：

| 参数 | 实测默认 | 它为什么在这里 |
|---|---|---|
| `ew_eta0` / `ew_eta_max` | **1e-06** / **1e-06** | ★★ 这一条是 `DECISION-2026-08-02-precond.md` 归档前唯一被挡住的**可测理由**：EW forcing 1e-6 **是当前默认**。决策记录（含**被排除的选项与推翻条件**）在 [`phases/p2/docs/dev_phase_two/DECISION-2026-08-02-precond.md`](../../phases/p2/docs/dev_phase_two/DECISION-2026-08-02-precond.md) |
| `upwind_c` | **1.5** | 人工密度强度 `C`。M1c 扫 `C ∈ {1.0, 1.5, 3.0}` |
| `m_crit` | **0.95** | 人工密度开关的 `M_c` |
| `m_cap` | **3.0** | 限幅的是**局部 Mach**（不是 M²）；判据是 `q2 >= q2_at_mach(m_cap, ...)` |
| `rho_floor` | **0.05** | 限幅的是**人工**密度 ρ̃，库不返回它 |
| `precond` | **`'amg'`** | M6 fine 的 `'direct'` 是 4h39m/26GB 的 splu 陷阱 |
| `entropy_correction` | **True** | ★★ **模型修正，不是稳定化技巧**。实测它让峰更**强**（M²max 6.426 ON vs 5.607 OFF） |
| `kutta_estimator` | **`'probe'`** | ★ 生产的翼身用 `'pressure'`；R19 的坑就是改了 `pressure` 那一行而 NACA 配方走 `probe` |
| `tip_taper` | **None** | ★ 生产的**保形翼身**显式传 `("vanish_smooth", 0.05·b_semi)`，代价 **−1.3 % cl_p**，是**模型偏置** |
| `gamma_target` | **None** | phase 5 R19 加的（规定环量）；`None` = 逐位沿用原表达式 |
| `n_picard_seed` | **0** | 冷启动回退是 `_SEED_FALLBACK = 5`，只在"无种子 ∧ 无热启动 ∧ 结果被钳制"三条同时成立时触发一次 |
| `wing3d` 的 `tip_cap` | **`"round"`** | ★ flat 必须**按名字显式索取** —— P13/G13.3 实测 flat 的凸边在加密下**发散**（p = +0.321） |

★★ **这张表是 phase 6 的责任面**：现在**没有任何东西**检查这些默认值有没有被静默改掉。
★★★ 但**不在这里立门** —— 使用者 2026-08-24 指示：门与测试案例的最终设置**还没有指示，且要大改**。
⇒ 本表的作用是**事实登记**，供那次设计当输入；**擅自把它变成一批锁会把旧结构锚死**。

## 2. 仍然在管事的裁决 D1–D5

原件：[`phases/p2/docs/dev_phase_two/roadmap.md`](../../phases/p2/docs/dev_phase_two/roadmap.md)。
★ 归档前它被挡住的理由是「**没有任何后续计划文件承接 D1–D5**」——**本节就是那次承接**：

| 裁决 | 内容 | 现在的状态（实测） |
|---|---|---|
| **D5** | **放弃 level-set 路线** | ✓ **已执行**：`pyfp3d/wake/` 不存在（9 文件 / 4624 行，PR #26）。`from pyfp3d.wake import ...` 是硬错误，**这是预期的终态** |
| **D2** | 搁置 S6 | phase 3 已激活（六面体贴体网格路线） |
| **D1 / D3 / D4** | 见原件 | ★ **未逐条核对** —— 如实记；phase 6 若要引用其中任何一条，**先回原件读，再在这里登记** |

★ **不假装承接了全部五条**：我核了 D5 与 D2（它们在后续文档里被引用过），D1/D3/D4 没有。
写"全部承接"会是一句无据的话。

## 3. 能力声明：**phase 6 必须产出自己的那一份**

`PHASE_TWO_CAPABILITY_BOUNDARY.md` 归档前是**唯一**陈述"求解器现在能做什么"的文件
（实测 phase 3/4/5 都没产出自己的）。它那份检查单说它是「最后一个可以移的」，
条件是「等后续产出**自己的**能力边界」。

⇒ **裁决把它变成参考，于是那个条件转成 phase 6 的一项交付**：

| | |
|---|---|
| 参考原件 | [`phases/p2/docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md`](../../phases/p2/docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md) —— 含五条指标的目标/实测/可达性、**14 条已排除路线**、未解释亏空、真实缺陷、证据完整性边界 |
| **八条产品指标的当前读数** | 在 [`progress.md`](progress.md) 的**产品指标追踪**表 —— 那是本阶段的活台账 |
| ★★ **phase 6 的交付** | **一份属于 phase 6 的能力边界**。在它产出之前，**引用能力声明就必须引用归档原件**，不许凭记忆转述 |

### 3.1 ★★★ 勘误 2026-09-06：**声明的马赫包线在 α = 0 上被实测推翻**

门审计（[合并报告](20260906-0100-gate-audit-synthesis.md) H11 / W1.5）实测：
CLAUDE.md 第 5 行写着 *"Target: wings at M∞ 0.3–0.87"*，而 `design.md`
§2 / §12 risk 2 早就登记了保守型全速势在 **M∞ ≈ 0.82–0.85、低升力**上的
Steinhoff–Jameson 非唯一性。**两句话之间的矛盾此前没有任何门去碰。**

**实测（NACA0012，α = 0，生产配方，冷启动）**：D05 自己的零升力带
（`|cl| ≤ 0.005`）在 **M 0.78 上守得住（0.012x）**，从 **M 0.80 起破**
（1.33x / 4.69x / 45x / **104x**）。而在 M 0.86 上，**只改初值**就得到
**四个各自收敛到 |R| ≈ 3e-13、零钳制的不同解**，cl 跨越 −0.0059 … +0.4360。

⇒ **能力声明必须分开写，不能只写一个 0.87**：

| 场景 | 可信上界 | 依据 |
|---|---|---|
| **有升力**（α ≠ 0，M6 / 翼身的生产工况） | 声明的 0.87 **未被本轮推翻** —— 本轮**没有在 α ≠ 0 上做同样的探针** | 未测，明写 |
| **α = 0 对称翼型** | **实测落在 M 0.78 与 0.80 之间** | `tests/D/test_D14`（能力边界门，红 = 好消息）+ `bench/studies/gate_audit_20260905/` |

★★ **这一条是 phase 6 那份能力边界文书的必备输入**：目前唯一在陈述能力的
是归档的 `PHASE_TWO_CAPABILITY_BOUNDARY.md`，而它**写不出上面这张表**
（那时既没有 D14，也没有这次的测量）。
★ **登记未做**：α ≠ 0 上的同款探针（固定一切、只改初值）。做了才谈得上
把有升力那一行从「未测」改成一个数。

## 4. ★★ 不许重开的封闭负结果（**引用归档，未重测**）

这一节存在的理由很实际：**重开一条已被测死的路要花真机时**。完整列表在归档的能力边界
（14 条）与 `docs/agent-rules.md` 的纪律 8。此处只钉最常被重提的几条：

- **G1.6 的修法**：h-refinement / recovery 调参 / Nitsche / 边界数据修正 —— **全部已排除**；
  P11 把 G1.6 **重新归因**为 h=0.08 上 P1 场的**固有能力**，不是壁面变分犯规。
  唯一通向字面判据的路是**等参 P2 壁面层**（route (b)，未采纳）。
- **M3a**（M6 七站 Cp RMS < 0.05）：★★ **测量为"无粘不可达"** —— 无粘模型地板 **0.0516–0.0707**，
  而参照是 Re 11.7e6 的**有粘**实验。⇒ 它是**粘性能力**的靶子，不是网格或离散的靶子。
- **超参数化（mapped-P1）曲面壁单元**：P11 实测负结果，**不要再提**。
- **B8 的约束侧翼尖治理**、**展向 Γ 平滑**：已测死。
- 「0.019 的差是分辨率」：口径仍是 ***strongly indicated, NOT earned***（2026-07-14 措辞裁决）。

## 5. 六件不许静默改掉的事

与立项文书 §5 同一份，此处不复制以免分叉 ⇒ 见
[立项文书 §5](20260824-0000-initiation.md#5-★★-六件不许静默改掉的事)。
★ 其中第 6 条（`usability.py` + `x_shock` 落带**不得**擅自接进 M1）是**裁决**，不是实现细节。

## 6. 东西都搬到哪了

| 原位置 | 现在 |
|---|---|
| `docs/dev_phase_two/`（7 个） | [`phases/p2/docs/dev_phase_two/`](../../phases/p2/docs/dev_phase_two/)（并入原有 85 个轮次文件） |
| `docs/dev_phase_three/`（68） | [`phases/p3/docs/dev_phase_three/`](../../phases/p3/docs/dev_phase_three/) |
| `docs/dev_phase_four/`（82） | [`phases/p4/docs/dev_phase_four/`](../../phases/p4/docs/dev_phase_four/) |
| `docs/dev_phase_five/`（61，含 `progress.md`） | [`phases/p5/docs/dev_phase_five/`](../../phases/p5/docs/dev_phase_five/) |
| `docs/inspection/` | ★ **已删除**。8 份 phase-1 报告在 [`phases/p1/docs/inspection/`](../../phases/p1/docs/inspection/)；2026-07-28 奠基审计 + exp1–exp6 在 [`phases/p2/docs/inspection/`](../../phases/p2/docs/inspection/)；2026-08-16 独立审计在 [`phases/p3/docs/inspection/`](../../phases/p3/docs/inspection/) |

`docs/` 现在只有 **7 个文件**：`overview.md`、`design.md`、`design_track_v.md`、
`agent-rules.md`，加 `dev_phase_six/` 的 4 个（立项、本文件、GS6.1-a 判定、`progress.md`）。

★ **轮次文件格式**：模板已归档 ⇒ 引用
[`phases/p2/docs/dev_phase_two/_TEMPLATE.md`](../../phases/p2/docs/dev_phase_two/_TEMPLATE.md)。
判定码不变：**PASS / FAIL / RECORDED / UNDEFINED**（收敛腿 < 2 报 UNDEFINED，不报方向）。

## 7. ★ 本文件自己的风险，写下来

**这份承接是"指向 + 少量实测"，不是"重新验证"。** §2 只核了 D5/D2，§4 完全是引用归档。
⇒ 任何一条被引用到会影响决定的地方，**回原件读**。
★★ 而 §1 是唯一**实测**的一节，所以它是最有资格被门保护的一节 —— **而门怎么设由使用者的设计指示决定**，
本文件只负责把事实摆在那里，**不预先决定门的形状**。
