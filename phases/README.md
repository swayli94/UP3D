# `phases/` —— phase 1–5 的归档

**这里放的是已结束阶段的材料。** phase 6 的工作面在仓库的原位置
(`pyfp3d/`、`tests/`、`bench/`、`cases/`、`docs/`)。

★★★ **2026-08-24 使用者裁决:phase 2–5 的文档全部移入 `phases/p*/`,`docs/inspection/` 删除。**
「以前的分析和文档都只是参考」—— phase 6 在自己的文档里指向或跟踪它们。
⇒ `docs/` 现在只剩 **4 份活文档 + 当前阶段目录**(`overview` / `design` / `design_track_v` /
`agent-rules` + `dev_phase_six/`)。

★★ 这次裁决**覆盖了** `p2/docs/dev_phase_two/README.md` 那份逐文件检查单挡下的 4 个文件。
挡它们的两条理由是**可测的事实**(EW forcing 1e-6 仍是库默认;它是唯一的能力声明),
所以**归档的同时必须有一个新的落点**,否则那不是归档而是丢失 ⇒
**[`docs/dev_phase_six/20260824-0300-carried-forward.md`](../docs/dev_phase_six/20260824-0300-carried-forward.md)**
就是那个落点(实测的库默认值、仍在管事的裁决、不许重开的封闭负结果、东西都搬到哪了)。
★ 而**能力边界仍然只有归档那一份** —— 产出 phase 6 自己的那一份是**phase 6 的交付**;
在那之前引用能力声明**必须引用归档原件**,不许凭记忆转述。

**测试代码一个文件都没有移动 —— 见下面「为什么测试不归档」。**

整理于 **2026-08-10**,基准提交 **`d224223`**(phase two 的 PR #25 合并点)。

## 分类口径(量出来的,不是凭印象)

相位归属 = **文件首次提交的日期**,边界 **2026-07-28**(phase two 纲领建立日)。
按此,整理前的 1606 个被跟踪文件里 phase 1 约 1029、phase 2 约 524。

## ★ 什么留在原位:被保留代码引用的东西

规则是**可测的**:**如果一个仍然存在的测试或脚本读它,它就是 phase 3 仍然需要的,于是留在原位。**
按这条规则,`cases/demo` 的 33 个子目录里 **9 个留下**、`cases/analysis` 的 33 个里 **10 个留下** ——
例如 `cases/demo/p11_curved_walls/`,`tests/test_laplace_sphere.py` 要读它来跑 **G1.6 Option C 的活判据**。

留在原位的还有:`pyfp3d/` 全部、`cases/reference_data/`(外部真值,永不编辑)、`cases/meshes/`、
73 个测试文件(**全部 conforming 锚点**)、`bench/` 的 7 个指标脚本闭包 + `bench/gate_results/`
(能力边界按路径引用那些 CSV)、以及 `docs/` 的 overview / design / design_track_v / agent-rules /
inspection 与 6 个 phase-3 要用的 phase-two 文件。

## 第二轮再减(2026-08-10,同日)

第一轮之后 `docs`/`bench` 里仍剩很多,**量了一下发现两处是我保守多留的**:

| | 第一轮后 | 再减后 | 依据 |
|---|---|---|---|
| `bench/` | 255 | **19** | `gate_results/` 的 246 个里,**被保留代码点名的只有 7 个 CSV**,`capability/` 的 189 个**一个都没有** ⇒ 只留 9 个被代码或 phase-3 文档点名的 + `capability_summary.json`,其余 **244 个**入 `p2/bench/gate_results/` |
| `docs/` | 31 | **23** | `inspection/` 的 21 个里 **8 个是 phase-1 审查报告** ⇒ 入 p1;留下 2026-07-28 的**奠基审计**及其 `20260728-audit/` 的 exp1–exp6(roadmap 把它当作 phase two 全部判断的基线) |
| `cases/` | 403 | **403** | ★ **测出来必须留**:`onera_m6_wakefree` 与 `naca0012_wakefree_2.5d` 两个 LS 专用族**仍被 4 个保留测试引用**(`test_b1_cut_elements` / `test_b3_lifting` / `test_b45_farfield` / `test_b7_onera_m6` —— 删除清单里"摘腿"的那批),要等 phase 3 摘掉 LS 腿才能走;19 个保留的 demo/analysis 是自洽的证据单元 |

★ 移走 244 个 CSV 之后,能力边界文件里"每条都指到提交的产物"那句话**就地更正**了产物的位置
(两处都写明,并强调**全部仍被跟踪** —— 归档是目录重组,不是取消证据)。

## 目录

| | 内容 |
|---|---|
| `p1/` | phase 1(Track P/M/B/V/A):24 个 demo + 23 个 analysis 证据链、phase-1 的 roadmap 与 demo_report、`design_track_b.md`、`analysis/`、`archive/`、**18 个纯 level-set 测试**、8 份审查报告 |
| `p2/` | phase 2:45 个 bench 脚本 + `s1_duct/`、**92 个文件**的 `dev_phase_two/`(2026-08-24 补入全部 7 个:`progress.md`、`LEVELSET_DELETION_INVENTORY.md`、`roadmap.md`(裁决 D1–D5)、`PHASE_TWO_CAPABILITY_BOUNDARY.md`(★ **仍是唯一的能力声明**)、`DECISION-2026-08-02-precond.md`、`_TEMPLATE.md`(★ **轮次格式仍在用**)、`README.md`(那份检查单))、**2026-07-28 奠基审计 + exp1–exp6** |
| `p3/` | phase 3:**68 个轮次文件** + **2026-08-16 独立审计**(触发 GS4.0) |
| `p4/` | phase 4(GS4.1 二维 strip IBL / 剪应力滞后):**82 个轮次文件** |
| `p5/` | phase 5(M1 二维无粘标定):**61 个文件**,含 `progress.md`(★ 它的 **§A 仍然成立 / §B 已撤回**两张表是 phase 6 的前置阅读 —— 引用 phase 5 任何数字之前先读那两张表) |

## ★ 第三轮:15 个目录搬回原位(保留规则原来只算了一跳)

全套复跑抓到 **2 failed**,机制是**传递依赖**:
`tests/test_meshgen_rae2822.py` → `bench/studies/v5_2_rae2822/run.py`(保留)
→ `_load("gv3_run", "bench/studies/v3_loose_coupling/run.py")`(**第二跳,被移走了**)。

⇒ 保留规则改为**传递闭包**,并且**判据先校验**:第一版闭包按"目录名出现在文本里"算,报 22 个,
而全套实测只坏 1 个 ⇒ 太松(把注释里的提及也算了);改成只认**路径样字面量**
`cases/(demo|analysis)/<name>/`,并断言它**必须包含那个已被实测证明会坏的** `v3_loose_coupling`,
才拿来用。据此 **15 个目录搬回原位**(`cases/` 403 → 601,`phases/` 1023 → 825)。

搬回不只是挪文件:那些目录里被改过的路径表达式要**反向还原**(41 处深度 + 归档前缀),
oracle 复验 **26 个文件里 0 处失败**;`b6_transonic` 回到原位后,**归档文档里指向它的 3 处链接又断了**
—— 移动是双向的,链接也得双向修。

★ 同轮还修掉一处更隐蔽的破坏:**快层 `run_capability_locks.py` 有 2 条锁指向已归档的 LS 测试**。
处置是**删掉而不是指进归档** —— 按 D5 那条路线已放弃、phase 3 第一件事就是删它们,而
**把一条活的收口门指进归档比没有这道门更糟**:它会跑、会绿、并为一条没人维护的路线断言能力。
⇒ 快层现在 **5 组、482 s、5/5 green**;conforming 锚点一条不少。

## ★ 归档保持旧布局是故意的

活的一侧在 2026-08-10 把 `cases/analysis/` **合并进了 `bench/studies/`**(使用者裁决:
`bench/` 与 `cases/analysis/` 本是同一件事的两个家)。**归档不跟这次改名** ——
`phases/p1/cases/analysis/` 保持原样,因为归档是**历史快照**,重命名它会让"当年长什么样"失真。
所以看到两种路径不必奇怪:**活的在 `bench/studies/`,历史的在 `phases/p1/cases/analysis/`**。

## ★ 已知限制(2026-08-11 补记):归档测试的兄弟模块导入没修

**3 个归档测试文件在归档里收集不动**:`test_b2_multivalued`、`test_b31_tip_fringe`、
`test_b8_span_blend` —— 它们用**兄弟模块导入**(`mesh_utils` / `_tol` / `conftest`),
而那些 helper 留在 `tests/`。

⚠ **这是本次整理引入的、而当时没验证到的一类** —— 下一节声称"三类路径断裂已修并验证",
**归档测试的 helper 导入不在那三类里**。发现它是在 phase 3 核门控账的时候
(`git worktree` 回到删除前的提交去数归档测试的项数,3 个文件报 ImportError)。

**处置:如实记,不修。** 归档不用于运行;而下一节给出的可运行配方
`git worktree add ../up3d-prereorg d224223` 指向**整理之前**的提交,本来就绕过这个问题。

## ★★ 归档脚本的可运行性:如实说

移动会打断脚本找路的方式,已**机械修好并验证**的有三类:
1. **58 处 `parents[N]`** 深度(两种写法:`Path(__file__).resolve().parents[k]` 与
   `HERE.parents[k]`,二者相差 1);
2. **22 处指向已移动目录的路径字面量**(只改真的移动了的,`cases/meshes/` 与留下的 demo 不动);
3. **41 个归档 bench 脚本的 `gate_results` 指向** —— 那个目录**留在** `bench/`,所以归档脚本
   是**横向**去取,不是往自己下面找;每一个都验证过解析结果落在含 CSV 的目录上。

**没有逐一验证到"能真的跑起来"。** 残留的少数交叉引用依赖被 import 的兄弟脚本来设置 `sys.path`,
静态检查看不透。所以口径是:

- **证据是已提交的 CSV / PNG**(项目纪律:"没有提交产物的说法不算证据"),它们**都还在,且仍被跟踪**;
- **脚本是出处(provenance)**;若要一个**保证能跑**的 phase 1/2 树,用一条命令:

```bash
git worktree add ../up3d-prereorg d224223     # 整理之前的完整状态
```

## 为什么是跟踪而不是 gitignore

最初的要求是把归档 gitignore。**改为跟踪**(使用者 2026-08-10 裁决),理由是:
gitignore 意味着**新 clone 里没有这些文件** —— 而那会让 phase 1/2 的每一个结论在 HEAD 上
**失去提交产物**,与本项目的头号纪律直接冲突。目录整洁的效果不变,证据仍在。

---

## ★★★ 2026-08-24 这一轮:两个**本文件上一节没列**的断裂类,以及为什么测试不归档

### 为什么**测试代码一个文件都没有移动**

使用者要求梳理「phase 2–5 的文档**和测试代码**」。文档移了,**测试一个都没移**,理由是**可测的**:

- 本文件第 59–61 行自己写着:**「把一条活的收口门指进归档比没有这道门更糟」** ——
  它会跑、会绿、并为一条没人维护的路线断言能力。
- 而 **72 个测试文件在门控全集里全部有通过项**(604 passed / 1 skipped / 4 xfailed,
  1:44:32 @8 线程),即**每一个都在检验活的 `pyfp3d/` 代码**。
- phase 1 的 18 个测试之所以能归档,是因为**它们那条路(level-set)被删除了**。
  **phase 2–5 没有删除任何路线** ⇒ 按同一条规则,**没有一个测试文件可归档**。

★ 唯一的例外是一条**该删不该归档**的:`tests/test_v6_wake_sheet.py::test_ab_bit_identity_gate_free_library`
是 2026-07-29 因**前提为假**而 `mark.skip` 的(松循环 run-to-run 不可复现)。
登记给 phase 6:**删除,而不是移进归档** —— 归档一条已知前提为假的测试等于把它藏起来。

### ★★ 断裂类 4:**根相对路径字面量**(本文件"已修的三类"里没有它)

上一节列了三类已机械修好的断裂(`parents[N]` 深度、移动目录的路径字面量、`gate_results` 横向指向)。
**这一轮实测出第四类**:活脚本里的**绑定文本行** ——

```
Binding gate text: docs/roadmap/track_v.md GV1.1(a)-(e)
```

这是 `bench/studies/` 的研究**指向管它的门**的审计链。2026-08-10 整理把
`docs/roadmap/` 移进了 `p1/`,**而这些绑定行一条都没有重指**:实测**仍有 139 处**
指向已不存在的 `docs/` 目录(`docs/roadmap/` 94、`docs/discussion_notes/` 17、
`docs/analysis/` 14、`docs/demo_report/` 11、`docs/archive/` 3),
另有 **184 处**死的 `docs/**/*.md` 文件字面量(`docs/roadmap/track_v.md` 57、`docs/roadmap.md` 50 领头)。

★★★ **本轮不修它们**,理由写下来:把一次目录搬迁和一次 184+139 处的链接修复混在一个提交里,
**两件事都变得无法验证**。⇒ 登记为 **phase 6 GS6.1 的对象**,并且它现在**有仪器了**(下节)。

### ★★ 断裂类 5:**相对 markdown 链接** —— 而这一类只有它自己的仪器看得见

`[x](../inspection/y.md)` 是**相对被引用文件自己的目录**解析的,
所以一个根相对的正则**在原理上看不见它**。实测代价:本轮搬迁**造成 239 条断链**
(移动前基线 9 → 移动后 248),其中最大一块是 `phases/p5/docs/dev_phase_five/progress.md`
**留在原位而它的 60 个轮次文件搬走了** —— 11 条兄弟链接当场断掉。

修法与验证:按**唯一 basename** 重指并重算相对路径(唯一命中才改,0 或 >1 报告不动)——
自动修 **231** 条,**12 条**因为 `progress.md` / `roadmap.md` / `README.md` 这类通名重复而按意图手改,
复验剩 **5 条,全部是 `[φ](TE) = Γ` 这个数学记号**,即本仓库唯一被明确允许的例外。

### ★ 这一轮的方法论,一句话

**先给搬迁做仪器,并先证明仪器有牙。** 顺序是:
① 移动前取基线(A=12 / C=9 / D=139);② `git mv`;
③ **移动后、改写前**跑仪器,**A 必须从 12 跳到 152** —— 这一步就是 G-TEETH,
一个报 0 的仪器在这里是坏的而不是好的;④ 改写;⑤ 复验回到基线。

★★ 而**第一版仪器是错的,且是被抽样检查抓住的**:它把
`[docs/roadmap.md](phases/p1/docs/roadmap.md)` 这种**markdown 显示文本**算成了路径引用
(CLAUDE.md 的文档地图明写那个名字是历史的),于是报出 277 条"断裂"。
**这与本文件上一节记的"第一版闭包太松(把注释里的提及也算了)"是同一个错误的第二次出现。**
⇒ 口径:**报一个数之前,先抽样看那个数里的东西是不是你以为的东西。**

### ★ 库文件被改了 8 个,而它是**纯文档改动 —— 用 AST 证明的,不是靠读 diff**

绑定文本也写在 `pyfp3d/` 的模块 docstring 里,所以重指改到了 8 个库文件
(`kernels/entropy.py`、`meshgen/structured.py`、`meshgen/wing3d.py`、`solve/newton.py`、
`viscous/closures.py`、`viscous/closures_2d.py`、`viscous/ibl3.py`、`viscous/strip2d.py`;23 insertions / 23 deletions)。
判据不是"看起来都在 docstring 里",而是:**把每个 docstring 清空之后比较 AST**,
8 个文件对移动前的基线**逐位相同** ⇒ 代码零改动。
(基线取自 `git worktree add <tmp> HEAD` —— 未提交状态下 HEAD 就是移动前的树。)
