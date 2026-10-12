# 20261011-1600 — **预注册:代码面收口(死代码处置)—— 删除全库零调用方的定义,翼身生成器族立测试策略**

**性质**:预注册(判据与 kill criterion 事前写死),随代码提交。
**对象**:`pyfp3d/` 的 10 个全库零调用方定义 + `meshgen/wingbody.py` 族的测试策略 +
`pyfp3d/__init__.py` 的无消费者属性。
**由来**:2026-10-10 独立审计的可达性分析(联合种子 tests/+bench/+cases/):
**10 个定义 / 301 行全库零调用方**;另有 `meshgen/wingbody.py` 整族 9 个
(含 443 行的 `onera_m6_wingbody_mesh`)只被两个网格生成脚本调用、**零测试执行**。
**纪律依据**:phase 3 的先例(D5)——「被杀路线的旋钮随路删除」;
「一个已知缺陷属于报告,不属于邀请调用方去扳的选项」。

---

## 1. 目的

把「哪些代码没有任何东西读它」从一份审计清单变成**仓库状态本身**:死代码要么删,
要么有活的调用方 + 测试。归档(phases/)里的引用**不算**调用方 —— 归档是
历史快照,明确不可运行(phases/README.md)。

**预先写下的失败标准**:
- 删除任何**活**路径可达的定义 ⇒ FAIL(可达性以联合种子重算为准);
- 删除后全套出现任何红 ⇒ FAIL(删的必须是死的);
- 翼身族"补测试"若只是 import 冒烟(不生成、不比对)⇒ FAIL(那是自欺的绿)。

## 2. 改动内容(事前登记)

### 2.1 删除(10 个定义 / 301 行,全库零调用方)

| 定义 | 处置理由 |
|---|---|
| `mesh/metrics.py::precompute_face_normals` | 零调用方(2026-07 旧审查已标记) |
| `meshgen/structured.py::airfoil_o_grid_2d` + `_outward_normals_closed` | 唯一调用方在归档 `phases/p3/bench/`(不可运行);O-grid 路线的生产路径是 `airfoil_hybrid_2d` |
| `physics/isentropic.py::speed_of_sound_squared`、`upwind_factor` | 零调用方(phase-1 审查写明 "no production callers") |
| `post/surface.py::_node_to_element_patches` | 零调用方 |
| `post/vtk_out.py::write_nodal_field_heatmap`、`export_error_heatmap` | 零调用方 |
| `solve/timing.py::records_to_columns` | 零调用方 |
| `viscous/coupling.py::make_picard_nonlifting_driver` | 唯一调用方在归档 `phases/p1/`;**且它在 `__all__` 里** ⇒ 从 `__all__` 一并移除(导出表不该卖不存在的东西) |

★ 每个删除点**在提交信息里列全**;文档若引用这些名字,顺手清理(预期只有
归档文档,不动归档 —— 归档引用死名字是历史真相)。

### 2.2 `pyfp3d/__init__.py` 的 `NOJIT` 属性

**无消费者**(各 kernel 自己重读环境变量;审计实测)。属性保留但注明
「便捷读法,不是开关源」⇒ 实际上更干净的是**删掉**这个误导性属性。
选删除;`PYFP3D_NOJIT` 环境变量语义不变。

### 2.3 翼身生成器族(`meshgen/wingbody.py` 9 个,443 行主函数)的测试策略

**不补全量测试**(生成器是全量级的重活,投入产出不配);改为**一条生成-比对
冒烟门**:跑 `onera_m6_wingbody_conforming` 的 **coarse** 生成(最快的一档),
断言 (a) 节点/单元数落在已提交 coarse 网格的同档带内,(b) 生成的网格与已提交
网格的拓扑指纹一致(manifest/sha 口径 —— 与 F01 的 mesh manifest 机制同族)。
G-TEETH:改坏生成器一个几何参数 ⇒ 门红。
**登记**:这**不是**覆盖,是「生成器改动会响」的最小保险;全覆盖属 W3.8 之后的事。

### 2.4 顺手(明确列出)

- `post/__init__.py` docstring 提到不存在的 `artifacts.py` ⇒ 删该提及。
- pyflakes 报的 6 个未用 import / 7 个未用局部变量:**本轮只删 import**
  (零风险);局部变量留给各自模块的下一次触诊(登记)。

## 3. 判据(事前写死)

| # | 判据 | 期望 |
|---|---|---|
| C1 | 删除清单与预注册 §2.1 逐字一致;`git diff` 无删除清单外的行 | 核对 |
| C2 | 删除后:全仓(不含归档)grep 这些名字 = 零调用方;`python -c "import pyfp3d"` 与各包 import 正常 | 实测 |
| C3 | 翼身冒烟门 G-TEETH:改坏几何 ⇒ 红;恢复 ⇒ 绿 | 实测 |
| C4 | 全套收口:**810 + Δ**(删除不增不减测试;A55 的 9 条不动;翼身门 +新增) | 0 failed,账闭合 |
| C5 | 覆盖率跑(A+B+F,NOJIT)的数字**如实**进判定文件(不承诺阈值) | 记录 |

**G-DOMAIN**:C4 的相反结果 = 红了 ⇒ 删到了活的 ⇒ 恢复该定义并在报告里写
「可达性分析漏了它」(仪器缺陷,不是删除决策错)—— 区分这两者。

## 4. kill criterion

若删除清单里任何一项在**活代码**(非归档、非测试)里其实有调用方
(审计的联合种子漏了)⇒ 该项从清单摘掉,照常登记;不许为了凑整删。

## 5. 明确不做

- 不动归档(phases/)里的任何文件。
- 不改任何判据/证据。
- 不做全量翼身生成器测试(§2.3 写明了为什么不配)。
- 不重构被删函数所在的模块(只删定义本身)。
