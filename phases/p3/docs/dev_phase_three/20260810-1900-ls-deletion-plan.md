# 20260810-1900 — phase 3 第 1 件工作:删除 level-set(计划,执行前提交)

依据裁决 **D5**(2026-08-09,使用者):**放弃 level-set 路线,未来只开发 conforming**;
删除是 phase 3 的开篇。清单 =
[LEVELSET_DELETION_INVENTORY.md](../../../p2/docs/dev_phase_two/LEVELSET_DELETION_INVENTORY.md)。

**本文件在动任何代码之前提交。** 它不是预注册(删除不是一次测量,没有判据要事前钉),
它要固定的是**范围**、**顺序**和**每一步的验证**,以及**两条被我核实过的前提**。

## 1. ★ 范围澄清:删"活的",不删归档

清单写于 2026-08-10 的目录整理**之前**,所以它的第 1、2 步(删纯 LS 的 demo/bench、
删 18 个纯 LS 测试)在今天读起来是**已经做过了** —— 那些文件在整理中进了
`phases/p1/`,而使用者当时明确裁决**归档仍然跟踪**,理由是 gitignore 会让 phase 1/2 的
每个结论在 HEAD 上失去提交产物。

⇒ **本轮删除的对象是"活的" level-set 代码**:9 个库文件、8 个混合测试里的 LS 腿、
2 个 LS 专用网格族的生成器,以及 `post/unified.py` 的塌缩。
**`phases/p1/` 里的 LS 材料保留** —— 归档存在的意义正是这个。
若使用者要连归档一起删,那是**另一个**决定,不在本轮。

## 2. 实际要做的六步(清单第 3–7 步 + 一步核实)

| # | 内容 | 验证 |
|---|---|---|
| **0** | ★ **删除前的绿锚点** | 快层 `bench/run_capability_locks.py`(5 组);**已在跑** |
| **1** | 摘 **8 个混合测试**的 LS 腿,**一个一个摘,每摘一个跑一次该文件** | 单文件 pytest 绿 |
| **2** | 塌缩 `pyfp3d/post/unified.py`(全库唯一 import 进 LS 的地方,它是**设计好的两路径统一入口** ⇒ 塌缩成 conforming 那一半,不是删掉) | 常开全套 |
| **3** | 删 **9 个库文件 / 4624 行** | 常开全套 + 快层 |
| **4** | 删 **2 个 LS 专用网格族**的生成器与已提交 stats(`naca0012_wakefree_2.5d`、`onera_m6_wakefree`) | 常开全套 |
| **5** | **门控全套**(3:04:44 @8 线程) | 0 failed |
| **6** | 收口:六面 + 快层;★ 并**更正 roadmap 的"约 2500 行"为 4624**(实测差 1.85 倍) | — |

## 3. ★★ 两条前提,我核实过了(不照抄清单)

**① 清单第 4 步的前提已经满足,不需要新写测试。** 清单写着"**先给 V2 的 transpiration 通道
补 conforming 等价测试**,再动 `b_base`" —— 那是按**文件名**扫出来的担心。实际读
`tests/test_v2_newton_rhs_channel.py`:**8 个测试里 6 个已经是 conforming Newton 腿**,
且 LS 腿断言的那两条性质**在 conforming 上已经各有一条**:

| LS 腿断言 | conforming 上已有的对应 |
|---|---|
| `test_ls_wall_rhs_zero_bit_identical`(ṁ = 0 逐位相同) | `test_newton_zero_external_rhs_bit_identical` |
| `test_ls_wall_rhs_nonzero_live`(非零 ṁ 真正进入求解) | `test_newton_external_rhs_solves_consistently` |

⇒ 摘掉那 2 条 LS 腿,**粘性侧的 transpiration 通道不会失去任何一条锁**。
★ 这条前提**必须核实而不是继承** —— 按清单的写法我会去写一个已经存在的测试。

**② 4624 行是实测,不是估计**:逐文件相加确认(`newton_ls` 1196 / `picard_ls` 1044 /
`multivalued` 907 / `cut_assembly` 556 / `cut_elements` 333 / `levelset` 225 /
`surface_ls` 186 / `schur_ls` 165 / `wake/__init__` 12)= **4624**。roadmap 的"约 2500"
要在第 6 步更正。

## 4. 会永久失去的东西(D5 已记的代价,这里再点一次名)

摘掉那 8 个混合测试的 LS 腿之后,下面这些结论**只剩已提交的 CSV/PNG,不再可重跑**:

| 结论 | 数字 | 承载它的腿 |
|---|---|---|
| P14 两条尾迹模型跨模型一致 | **0.17% / 0.36%** | `test_p14_te_pressure` 侧 + demo |
| B9 翼身两路径一致 | **0.4% / 0.6%** | `test_b9_wingbody_ls`(**已归档**) |
| A1 四驱动器成本对比 | 三维两条 Newton 都是预条件子瓶颈(~40% 墙钟) | `test_a1_instrumentation` |
| **B5 远场域稳健性** | Γ 在 R ∈ {15,…,120}c 内 **0.45% / 1.09%** | `test_b45_farfield` |
| B18 翼身跨声速两路径对照 | — | `test_b18_wingbody_transonic` |
| 切割单元结构锁 | B25 `inboard_clip` 等 | `test_b1_cut_elements` |

⇒ **phase 3 之后的能力边界必须在没有交叉验证的前提下陈述**,唯一的外部裁判是
`cases/reference_data/`。这不是本轮引入的代价,是 D5 当时就记下的。

## 5. 本轮不做

- **不删归档**(§1)。
- **不动 conforming 的任何数值路径** —— 删除必须是"只减不改":conforming 的答案
  一个数字都不许动。第 2/3 步之后的常开全套若有**任何** conforming 断言变化,
  就是删错了,回退重来。
- **不趁手改别的** —— 4624 行的机械删除已经足够大,任何顺路改动会让"是删除造成的
  还是顺路改动造成的"无法分辨(这正是本项目反复付过代价的那类混淆)。
