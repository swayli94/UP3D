# 20260811-0100 — 判定:level-set 已删除,**只减不改**在两个层面都算平了

计划 [20260810-1900](20260810-1900-ls-deletion-plan.md)(动代码之前提交,`98a0cc9`)。
裁决 **D5**;清单 [LEVELSET_DELETION_INVENTORY.md](../../../p2/docs/dev_phase_two/LEVELSET_DELETION_INVENTORY.md)。

## 0. 三句话

**① 删掉了 9 个库文件 / 4624 deletions**,与 AST 逐文件相加的实测**精确一致** ——
也就此确认 roadmap 原写的"约 2500 行"**低了 1.85 倍**(已就地更正)。

**② "只减不改"不是"看起来没坏",是两层账都算平**:

| | 删除前 | 删除后 | 差额的逐项解释 |
|---|---|---|---|
| **常开全套** | 538 passed | **457 passed** | 归档 6 文件 71+11 项、摘 5 条腿、搬 1 条 → 逐项吻合 |
| ★ **删库那一步单独看** | 457 passed | **457 passed** | **删 4624 行库代码,passed 一个没动**;+1 skipped = 新增的门控锁 |
| **门控全套** | 720 / 2 / 8(3:04:44) | **466 / 1 / 4(2:08:50)** | 208 + 47 + 5 − 1 = **259** ✓ |

★ 门控那 47 项来自 3 个**在归档里收集不动**的文件(见 §3),用 `git worktree` 回到删除前的
提交数出来的 —— 不是估的。

**③ 途中补了一条能力锁**,所以本轮**不是纯只减**,理由见 §2。

## 1. 归档 vs 摘腿:按"目的"判,不按名字判

清单写的是"18 整删 / 8 摘腿"。**跟随 helper 之后重算,实际是 23 整删 / 3 摘腿** ——
8 个"混合"文件里 **5 个 100% 依赖 LS**(b1 32/32、b45 6/6、b4 5/5、b3 4/4、b18 4/4)。
★ 而我**第一版分类器低报了**(b1 说 30/32、b4 说 2/5),因为它只看测试体内的关键字与作为参数的
fixture,**不跟随体内按名字调用的 helper**(b45 的 `_solve_b`);改成传递闭包后,又在
**最可能出错的那个文件上验证**:`test_b4_te_control_volume` 的清单描述是"两路径对照",分类器说
5/5 —— 去读代码,每个测试都建 `MultivaluedOperator`,**conforming 只是参考值**,5/5 是对的。

demo 与研究按同一条规则(**目的即跨模型 ⇒ 归档;LS 只是侧腿 ⇒ 摘腿**),逐个查 docstring 与门的判据:

| 归档(目的即跨模型) | 依据 |
|---|---|
| b18 / b9 demo | 头条就是两路径一致(0.4%/0.6%、8/8 含 cross-model 行) |
| b6 demo | docstring 第一句 "on the level-set path" |
| **p13 demo** | **G13.1 门本身**就是同网格 conforming-vs-LS A/B |
| **p14 demo** | V14.6 跨模型 + **G14.7 判据原文** "< 1% vs the level-set oracle" |
| a1 研究 | 自述"conforming-vs-level-set 成本对比"的共享机具 |
| a2 研究 | 两路径 TE 保真归因 |
| b28 研究 | 纯 LS flat-sheet harness |

| 摘腿(LS 是侧腿) | 摘掉什么 |
|---|---|
| `bench/run_capability_matrix.py` | 10 个 level-set cell + 4 个 helper(567 → 498 行);dispatch 的 LS 分支改成**显式 raise** |
| `bench/studies/v2_transpiration_channel/run.py` | leg 5(五条腿之一,不是门的判据) |
| `tests/test_v2_newton_rhs_channel.py` | 2 条 LS 腿 |
| `tests/test_a1_instrumentation.py` | 3 条 LS 腿 |

★ p13 我先按"摘腿"改了一半,**读到 G13.1 的判据才发现分类错**,用 `git checkout HEAD --`
撤回(不是裸 `git checkout --` —— 那会从可能被污染的 index 恢复,本仓库为此付过两次代价)。

## 2. ★★ 为什么本轮不是纯"只减":先补了一条能力锁

使用者问"**既然放弃 LS 了,为什么还要跨模型对照**" —— 这一问戳出的不是措辞问题,是**真缺口**:
b18 的 demo 除跨模型外还记着 **conforming 翼身 M0.84 / cl_p 0.2738 / 0 clamps**,而那条能力
**没有任何活的锁**(b9 的测试只到 M0.5;b31 锁的是 blend 代数;b18 自己的测试 4/4 全是 LS 已归档)。
⇒ 能力边界在陈述 M0.84,却**没有东西会因为它坏掉而变红**,而 phase 3 下一步换网格范式正是最可能
弄坏它的工作。

所以顺序改成**先补锁再归档**:新增 `tests/test_b32_wingbody_conforming_transonic.py`(门控),
锚点全部取 b18 已提交的值(`checks.csv` 的 GB18.1 行 + `cl_vs_mach.csv`),配方是那个 demo 的
`conf_ramp` 逐字。**实测 1 passed / 652 s**。
docstring 里明写它是**漂移锁而非正确性声明**(翼身升力无外部真值),clamp 按 GS1.4 单独断言,
并声明生产 taper 的 **−1.3% cl_p 偏差**。

## 3. ★ 一处如实补记的限制:归档测试的兄弟模块导入没修

核门控账时发现 **3 个归档测试文件在归档里收集不动**
(`test_b2_multivalued`、`test_b31_tip_fringe`、`test_b8_span_blend`):它们用兄弟模块导入
(`mesh_utils` / `_tol` / `conftest`),而那些 helper 留在 `tests/`。
⚠ **这是目录整理引入的、而我当时没验证到的一类** —— `phases/README` 声称"三类路径断裂已修并验证",
归档测试的 helper 导入**不在那三类里**。

**处置:如实记,不修。** 归档不用于运行(证据是已提交的 CSV/PNG,脚本是出处),而
`phases/README` 给出的可运行配方 **`git worktree add ../up3d-prereorg d224223`**(整理**之前**的提交)
本来就绕过这个问题。已在 `phases/README` 就地补上这条限制。

## 4. 顺带必须标注的一个数字(否则会被悄悄改掉含义)

能力边界 M5 那条引用 `run_capability_matrix` 的"**13 配置 × 阶梯 = 78 点、72/78、clamp 0/78**" ——
那 78 点**含 10 个 level-set cell**。committed CSV 作为**历史读数**仍然有效,但该脚本此后只测
conforming ⇒ 三个数都已加上 **"pre-deletion"** 限定。
这正是纪律 #11 要防的:**删除会让一个仍在被引用的数字改变含义,而它自己不会报错。**

## 5. 删完之后只剩的一条尾迹路线

**conforming**:`mesh/wake_cut.py`(尾迹切割)+ `constraints/te_pressure.py`(P14 压强相等 Kutta)
+ `solve/newton.py` / `solve/picard.py`,配 `tip_taper`(B31/B32,带 **−1.3% cl** 的模型偏差)。
`pyfp3d/wake/` 整个包已不存在 —— 它只为承载 LS。

## 6. 交给下一轮

1. **六面体贴体网格**(phase 3 第 2 件工作)—— ★ 开工那一轮**先写 kill criterion**(原则 8,
   不许事后补);依据在能力边界 §7,逐条标了测量 vs 推断。
2. **D5 的代价现在生效了**:P14 跨模型 0.17%/0.36%、A1 四驱动器成本、B5 远场域稳健性、
   B18 两路径 —— 全部只剩已提交 CSV/PNG,**不再可重跑**。此后唯一的外部裁判是
   `cases/reference_data/`。
