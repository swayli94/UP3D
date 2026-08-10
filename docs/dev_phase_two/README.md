# `docs/dev_phase_two/` —— 目录名是历史的,**里面 6 个文件是 phase 3 的输入**

★ **先说清这个矛盾**:目录叫 `dev_phase_two`,但 2026-08-10 的整理已把 **85 个轮次文件**移进
[`phases/p2/docs/dev_phase_two/`](../../phases/p2/docs/dev_phase_two/)。**留在这里的 6 个,
没有一个是"phase two 的开发过程"** —— 它们是**当前计划 + 交接材料**。

**不改目录名是使用者 2026-08-10 的决定**:所有已提交文档与提交信息里的引用都不动,零链接改动、零风险;
代价就是这个名字要靠本文件解释。

## 这 6 个文件分别是什么

| 文件 | 是什么 | 现在是否**活的** |
|---|---|---|
| [roadmap.md](roadmap.md) | **phase two 唯一权威计划**:五条产品指标、八条工作原则、S0–S6、裁决 **D1–D5** | ★ **活的** —— D5(放弃 level-set)、D2(搁置 S6,phase 3 正要激活)现在还在管事 |
| [PHASE_TWO_CAPABILITY_BOUNDARY.md](PHASE_TWO_CAPABILITY_BOUNDARY.md) | ★★ **phase 3 接手先读这一份**:五条指标的目标/实测/为什么/可达性、14 条已排除路线、未解释亏空、真实缺陷、证据完整性边界 | ★★ **活的,且唯一** —— 目前只有它陈述"求解器现在能做什么" |
| [LEVELSET_DELETION_INVENTORY.md](LEVELSET_DELETION_INVENTORY.md) | **phase 3 第一件工作的清单**(AST 扫真实 import:库 9 文件/4624 行、18 整删 + 8 摘腿、七步顺序) | **活的** —— 删除还没做 |
| [progress.md](progress.md) | 71 轮台账 + 阶段进度概览 + 产品指标追踪 | 历史为主,两张表仍在用 |
| [DECISION-2026-08-02-precond.md](DECISION-2026-08-02-precond.md) | 预条件子/EW forcing 的决策记录,**含被排除选项与推翻条件** | ★ **活的** —— EW forcing 1e-6 是**当前默认** |
| [_TEMPLATE.md](_TEMPLATE.md) | 轮次文件格式。★ **预注册在它所管的测量之前提交** —— 这个先后就是要点本身 | 活的(若 phase 3 沿用) |

## ★★ 什么时候可以把它们移进 `phases/p2/docs/`

使用者 2026-08-10 问的就是这件事。答案是**能,但每个文件有自己的条件,不是一次性、也不是按日期**。
下表就是那份检查单 —— 将来那次移动**对着条件勾**,不要临时判断:

| 文件 | 可以移的条件 |
|---|---|
| `progress.md` | phase 3 有了自己的台账之后(纯历史) |
| `roadmap.md` | phase 3 的计划文件**已把仍然有效的裁决承接过去**。⚠ 在那之前不能移 —— 移走等于让 D5/D2 无处可查 |
| `LEVELSET_DELETION_INVENTORY.md` | **删除做完那一刻**(它自己就是"做完了"的标记) |
| `DECISION-2026-08-02-precond.md` | 那个决策被重新审视、或 EW forcing 不再是默认之后 |
| `_TEMPLATE.md` | phase 3 **不**沿用这个轮次格式时;沿用就别移 |
| **`PHASE_TWO_CAPABILITY_BOUNDARY.md`** | ★ **最后一个** —— 等 phase 3 产出**自己的**能力边界。早移就等于把当前能力声明归档了,而它是唯一那一份 |

★ 移动时**必须同时**:① 重指链接(本仓库有可复用的检查:全库 markdown 断链应为 0,唯一允许的
例外是 `[φ](TE) = Γ` 这个数学记号);② 更新 `CLAUDE.md` 的文档地图与 `PROJECT_STRUCTURE.md`;
③ 按收口仪式第 6 步 **grep 判据在 `tests/` 里的出现**。

## 相邻的两个面(不属于本目录,但常被一起问)

- **[`docs/inspection/`](../inspection/)** —— **持续存在的面,不属于任何阶段**:审查报告一律写在那里
  (使用者的长期要求)。里面现在是 **2026-07-28 的奠基审计 + 它的 exp1–exp6**,phase-two 的 roadmap
  把那份审计当作**全部判断的基线**;审计与它自己的实验是**一个单元**,所以没有拆去 p2。
  8 份 phase-1 的审查报告已在 [`phases/p1/docs/inspection/`](../../phases/p1/docs/inspection/)。
- **`bench/`** —— 反复重跑的**产品指标台**(见 `bench/README.md`);一次性的登记研究在
  `bench/studies/`。这条界 2026-08-10 才划清并写下来,在那之前 `bench/` 与 `cases/analysis/`
  是同一件事的两个家。
