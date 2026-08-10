# `phases/` —— phase 1 与 phase 2 的归档

**这里放的是已结束阶段的材料。** phase 3 的工作面在仓库的原位置
(`pyfp3d/`、`tests/`、`bench/`、`cases/`、`docs/`)。

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
| `p1/` | phase 1(Track P/M/B/V/A):24 个 demo + 23 个 analysis 证据链、phase-1 的 roadmap 与 demo_report、`design_track_b.md`、`analysis/`、`archive/`、**18 个纯 level-set 测试** |
| `p2/` | phase 2:45 个 bench 脚本 + `s1_duct/`、**85 个轮次文件** |

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
