# 20260824-0200 — **第六阶段第一轮（GS6.1-a）判定**：phase 2–5 文档归档 —— **移了 210 个文件，测试一个都没移**，而搬迁本身造成 **239 条断链**（已修回基线）

**性质**：判定。使用者指示 2026-08-24：**梳理 phase 2–5 的文档与测试代码、以及
`docs/inspection/`，移入 `phases/p*/`**，参照 [`phases/README.md`](../../phases/README.md)。
**`pyfp3d/` 代码零改动 —— 用 AST 证明（§4）。**

---

## ★★★ 勘误 2026-08-24（同日，第二波）—— **§2 的结论被使用者裁决覆盖**

**本文件记录的是第一波。**同日使用者裁决:**phase 2–5 的文档全部移入 `phases/p*/`,
`docs/inspection/` 删除**,理由是「**我们现在本来就是要进行全面梳理,以前的分析和文档都只是参考**」。

⇒ §2 里被检查单挡下的那 **4 个文件**(`roadmap.md`、`PHASE_TWO_CAPABILITY_BOUNDARY.md`、
`DECISION-2026-08-02-precond.md`、`_TEMPLATE.md`)以及 `phases/p5/docs/dev_phase_five/progress.md`
**现在也已归档**。**本文件下面的正文一字不改** —— 那次勾选的读数是真的,
而**裁决可以覆盖一份检查单,这两件事不冲突**。

★★ 但**挡它们的两条理由是可测的事实**(EW forcing 1e-6 仍是库默认;它是唯一的能力声明),
所以归档的同时必须有一个**新的落点**,否则那不是归档而是丢失 ⇒
**[承接台账](20260824-0300-carried-forward.md)** 就是那个落点。

★ 本文件正文里出现的 `phases/p2/docs/dev_phase_two/README.md`、`phases/p5/docs/dev_phase_five/progress.md`
等路径,按 CLAUDE.md 文档地图的既有口径读作**历史名**:
它们现在分别在 `phases/p2/docs/dev_phase_two/` 与 `phases/p5/docs/dev_phase_five/`。

---

## 0. 一句话

按 `phases/p2/docs/dev_phase_two/README.md` **自己那份 2026-08-10 写下的逐文件检查单**勾选，
**phase-2 六件套里只有 2 个能移，4 个被条件挡住** —— 而**两条阻塞是关于当前库的可测事实**。
**测试代码一个文件都没有移动**，理由可测（§3）。
★★★ 而本轮最贵的发现是**搬迁自己造成的**：**239 条相对 markdown 断链**，
那一类**只有它自己的仪器看得见**，而 2026-08-10 那次整理的"已修三类"里**没有它**。

## 1. 移动了什么

| 目标 | 文件数 | |
|---|---|---|
| `phases/p2/docs/dev_phase_two/` | **+2** | `progress.md`、`LEVELSET_DELETION_INVENTORY.md`（条件已满足，见 §2） |
| `phases/p3/docs/dev_phase_three/` | **68** | phase 3 全部轮次文件 |
| `phases/p4/docs/dev_phase_four/` | **82** | phase 4（GS4.1）全部轮次文件 |
| `phases/p5/docs/dev_phase_five/` | **60** | phase 5 轮次文件 —— ★ `progress.md` **留下** |
| `phases/p2/docs/inspection/` | **13** | 2026-07-28 奠基审计 + `20260728-audit/`（exp1–exp6 + results）**整体一个单元** |
| `phases/p3/docs/inspection/` | **1** | 2026-08-16 独立审计（触发 GS4.0） |
| **合计** | **226** | |

`docs/` 现在只剩 **12 个文件**：overview / design / design_track_v / agent-rules、
`dev_phase_six/`（2）、`dev_phase_five/progress.md`、`dev_phase_two/`（4）、`inspection/README.md`。

## 2. ★★★ 六件套：**对着条件勾**，而检查单起作用了

`phases/p2/docs/dev_phase_two/README.md` 2026-08-10 就写了「**将来那次移动对着条件勾，不要临时判断**」。
勾选结果：

| 文件 | 勾选 | 依据（**两条是实测，不是判断**） |
|---|---|---|
| `progress.md` | ✓ **移** | phase 3/4/5/6 各有自己的台账 ⇒ 纯历史 |
| `LEVELSET_DELETION_INVENTORY.md` | ✓ **移** | `pyfp3d/wake/` 不存在（PR #26）⇒ 它自己就是"做完了"的标记 |
| `roadmap.md` | ★ **留** | 实测**没有任何后续计划文件承接 D1–D5**：`D5` 只出现在 `CLAUDE.md`（文档地图不是计划），`D2` 只在 phase-3 的**轮次**文件里被引用 |
| `DECISION-2026-08-02-precond.md` | ★ **留** | 实测 `newton.py`：`ew_eta0: float = 1e-6`、`ew_eta_max: float = 1e-6` ⇒ **仍是库默认** |
| `_TEMPLATE.md` | ★ **留** | phase 6 的 `progress.md` 正在沿用**并链接**它 |
| `PHASE_TWO_CAPABILITY_BOUNDARY.md` | ★★★ **留** | `find docs -iname "*capability*"` 只有它 ⇒ phase 3/4/5 **都没有产出自己的能力边界** |

★★ **值得单记**：若按"phase 2 已结束"这个**日期**理由一次性搬走，就会把**当前库默认值的决策记录**
和**唯一的能力声明**归档掉。**那份检查单正是为了拦住这件事才存在的，它拦住了。**

★ `phases/p5/docs/dev_phase_five/progress.md` 按**同一条理由**留下 —— phase 6 立项把它的 §A/§B 列为前置阅读。
可移条件已登记：**phase 6 的结论不再依赖那两张表时**。

## 3. ★★ 测试代码：**一个文件都没有移动**，而这是测量结论

使用者的指示包含"测试代码"。答案是**不移**，依据是 `phases/README.md` 第 59–61 行自己的规则：

> **「把一条活的收口门指进归档比没有这道门更糟」** —— 它会跑、会绿、并为一条没人维护的路线断言能力。

- 实测：**72 个测试文件在门控全集里全部有通过项**（604 passed / 1 skipped / 4 xfailed，
  1:44:32 @8 线程）⇒ **每一个都在检验活的 `pyfp3d/` 代码**；
- phase 1 的 18 个测试能归档，是因为**它们那条路（level-set）被删除了**；
  **phase 2–5 没有删除任何路线** ⇒ 按同一条规则，**没有一个测试文件可归档**。

★ 唯一例外是一条**该删不该归档**的：`test_v6_wake_sheet::test_ab_bit_identity_gate_free_library`，
2026-07-29 因**前提为假**而 `mark.skip`（松循环 run-to-run 不可复现）。
登记：**删除，不是移进归档** —— 归档一条已知前提为假的测试等于把它藏起来。

## 4. ★★ 库文件被改了 8 个，而它是**纯 docstring 改动 —— 用 AST 证明**

绑定文本也写在 `pyfp3d/` 的模块 docstring 里，所以重指改到了
`kernels/entropy.py`、`meshgen/structured.py`、`meshgen/wing3d.py`、`solve/newton.py`、
`viscous/closures.py`、`viscous/closures_2d.py`、`viscous/ibl3.py`、`viscous/strip2d.py`
（23 insertions / 23 deletions）。

判据**不是**"看起来都在 docstring 里"（我第一版的 grep 过滤器就漏了 docstring 的**续行**），
而是：**把每个 docstring 清空后比较 AST**，8 个文件对移动前基线**逐位相同** ⇒ 代码零改动。
基线取自 `git worktree add <tmp> HEAD` —— 未提交状态下 HEAD 就是移动前的树。

## 5. ★★★ 本轮的仪器，与它抓到的两个新断裂类

### 5.1 顺序就是方法：**先证明仪器有牙，再信它**

| 步 | 动作 | 读数 |
|---|---|---|
| ① | 移动**前**取基线 | A=12（我范围内的预存在死引用）· C=9 · D=139 |
| ② | `git mv` | — |
| ③ | ★★★ **移动后、改写前**跑仪器 | **A 12 → 152** ⇒ **G-TEETH 通过** |
| ④ | 改写 | — |
| ⑤ | 复验 | 见 §5.4 |

★★ 第 ③ 步是本轮的要点：**一个在这一刻报 0 的仪器是坏的，不是好的。**

### 5.2 ★★★ 断裂类 4：**根相对路径字面量**（"已修三类"里没有它）

活脚本里的**绑定文本行** —— `Binding gate text: phases/p1/docs/roadmap/track_v.md GV1.1(a)-(e)` ——
是 `bench/studies/` 的研究**指向管它的门**的审计链。2026-08-10 把 `phases/p1/docs/roadmap/` 移进了 `p1/`，
**这些绑定行一条都没重指**。实测**预存在** **227** 条死的 `docs/**` 文件字面量
（`phases/p1/docs/roadmap/track_v.md` 57、`phases/p1/docs/roadmap.md` 50 领头）+ **139** 处指向已不存在的 `docs/` 目录。
★ 本轮的前缀限定改写**顺带修掉 65 条**（都落在我搬迁的那些目录里）；
**其余是 phase-1 旧账，本轮不修** —— 把一次目录搬迁和一次两百多处的链接修复混在一个提交里，
**两件事都变得无法验证**。⇒ 登记为 **GS6.1 正式对象**，而它现在**有仪器了**。

### 5.3 ★★★ 断裂类 5：**相对 markdown 链接** —— 本轮真正的代价

`[x](../inspection/y.md)` 相对**被引用文件自己的目录**解析，
所以根相对的正则**在原理上看不见它**。实测：搬迁**造成 239 条断链**（基线 9 → 248）。
最大一块是 **`phases/p5/docs/dev_phase_five/progress.md` 留在原位而它的 60 个轮次文件搬走了**，
11 条兄弟链接当场断掉 —— **"留一个、搬其余"这个决定本身有链接代价，而它不在任何检查单上。**

修法：按**唯一 basename** 重指并重算相对路径（唯一命中才改，0 或 >1 报告不动）⇒
自动修 **231**；**12 条**因 `progress.md` / `roadmap.md` / `README.md` 这类通名重复而按意图手改。

### 5.4 复验：**按发现的集合差比，不按计数比**

复验时 B/C/D 的**计数上升**了 —— 因为**我为记录这些亏空而写的新文档，本身引用了那些死路径**
（`../inspection/y.md` 就是我自己举的例子）。⇒ 改成**集合差**，并排除本轮我新写的文件：

| 仪器 | 基线 | 现在 | **我造成的新断裂** | 顺带修掉 |
|---|---|---|---|---|
| B 死字面量 | 227 | 165 | ★ **0** | **65** |
| C 相对链接 | 9 | 8 | ★ **0** | **5** |

两条初看"新增"的，核后都是**同一条发现换了路径**（键是 `(文件, 路径)`，改名读起来像新增）：
`20260816-2200-independent-audit-zh.md → phases/p1/docs/roadmap.md` 在旧路径里就有 2 次；
`progress.md → ](TE)` 是 `[φ](TE) = Γ` 这个**唯一被明确允许的数学记号例外**。
⇒ **净新增断裂 = 0，另修掉 70 条旧账。**

## 6. ★★ 而第一版仪器是错的，且是被**抽样检查**抓住的

它把 `[phases/p1/docs/roadmap.md](phases/p1/docs/roadmap.md)` 这种 **markdown 显示文本**算成路径引用
（CLAUDE.md 的文档地图**明写**那个名字是历史的），于是报出 **277 条"断裂"**。
我准备把这个数写进结论**之前**去抽样看了几条，才发现它量的不是我以为的东西。

★★★ **这与 `phases/README.md` 自己记的"第一版闭包太松（把注释里的提及也算了）"是同一个错误的第二次出现。**
⇒ 口径记下来：**报一个数之前，先抽样看那个数里的东西是不是你以为的东西。**

## 7. 事前预测的兑现（**不改预测**）

| # | 预测 | 实测 | |
|---|---|---|---|
| 1 | 六件套大部分可移 | ★ **4 个被挡**，只 2 个可移 | **落空** |
| 2 | 测试里有一批可归档 | ★★ **一个都没有**（72/72 都在检验活代码） | **落空** |
| 3 | 断裂主要是根相对字面量 | ★★★ **相对链接才是大头（239 条）**，且我最初的仪器看不见它 | **落空** |
| 4 | `pyfp3d/` 不会被碰到 | 被改 8 个文件（docstring 里的绑定文本） | **落空** |

⇒ 四条**全部落空**，一次都没有事后改预测。

## 8. 实测口径

- 移动 **226 文件**（`git mv`，历史跟随）；改写 **142 + 33 + 13** 个文件；`pyfp3d/` **AST 证明零代码改动**；
- 主回归 `tests/A/test_A01_freestream_preservation.py` **6 passed**；全套见台账行；
- 仪器脚本在 scratchpad（一次性），**判据与读数在本文件与 `phases/README.md` §2026-08-24 节**；
- 未修：**165 条死字面量 + 139 处死目录引用**（phase-1 旧账）⇒ **GS6.1 的对象**。
