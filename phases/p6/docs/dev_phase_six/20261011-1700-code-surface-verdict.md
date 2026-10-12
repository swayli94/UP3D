# 20261011-1700 — **判定:代码面收口 —— 删 10 个零调用方定义(−338 行),翼身生成器立冒烟门;kill criterion 正确摘除一项**

轮次文件名:`20261011-1700-code-surface-verdict.md`
阶段 / 门:无门(代码面梳理)—— 预注册:
[20261011-1600-code-surface-prereg.md](20261011-1600-code-surface-prereg.md)
提交:分支 `claude/phase-six-audit`(本文件随判定提交)

---

## 1. 目的

把「哪些代码没有任何东西读它」从审计清单变成仓库状态本身。判据 C1–C5 与
kill criterion 见预注册,事前写死。

## 2. 改了什么

- **删除 10 个定义 / 338 行**(逐项二次 grep 核对零活调用方后删):
  `mesh/metrics.py::precompute_face_normals`、`meshgen/structured.py::airfoil_o_grid_2d`
  + `_outward_normals_closed`、`physics/isentropic.py::speed_of_sound_squared` +
  `upwind_factor`、`post/surface.py::_node_to_element_patches`、
  `post/vtk_out.py::write_nodal_field_heatmap` + `export_error_heatmap`、
  `solve/timing.py::records_to_columns`、`viscous/coupling.py::make_picard_nonlifting_driver`
  (`__all__` 除名)。
- **顺手**:pyflakes 6 个未用 import 清零(删后 pyflakes pyfp3d/ 零未用 import);
  `post/__init__.py` 的不存在模块提及清理;`PROJECT_STRUCTURE.md` 两处
  被删名字的引用清理;`kernels/upwind.py` 一处指向被删函数的 docstring 提及清理。
- **新增 `tests/A/test_A56_wingbody_generator_smoke.py`**:翼身生成器
  (9 个定义、含 443 行 `onera_m6_wingbody_mesh`,此前零测试执行)的
  生成-比对冒烟门 —— coarse 档生成 **9.6 s**、**位确定**(sha256 逐字复现
  committed manifest);指纹含写盘 sha256(计数看不见排序/符号翻转,注释写明);
  第二条测试用扰动副本喂比对函数(F06 式行为验证)。
- 预注册 §2.2 的 **NOJIT 属性:摘除未删**(§4.2)。

## 3. 怎么测的

```bash
pytest tests/A/test_A56_wingbody_generator_smoke.py   # 2 passed;G-TEETH 见 §4.1
pytest tests/A/ tests/B/ tests/F/ -q                   # 646 passed / 3 skipped / 0 failed
python -c "import pyfp3d" + 各被删模块逐个 import      # 冒烟全过
# 收口全套(16 线程纪律,BLAS=1)
```

## 4. 结果

### 4.1 判据逐条

| 判据 | 结果 |
|---|---|
| C1 删除清单与预注册逐字一致 | ✅(唯一偏离 = §4.2 的摘除,按 kill criterion) |
| C2 删后零残留 + import 正常 | ✅ grep + 冒烟 |
| C3 翼身门 G-TEETH | ✅ 改坏 committed sha256 一位 ⇒ 红在该断言;恢复复绿。★ 门还**抓到了它自己**:初版扰动是空操作(sha 以 0 开头,`"0"+sha[1:]` 不变)—— 行为验证在干净态假红,这正是那条测试存在的理由 |
| C4 全套收口 | 见基线行 |
| C5 覆盖率数字 | 第一轮测量因工作树中途被本轮删改而作废(已杀),提交后重测 —— 归 W3.8 自己的轮次 |

### 4.2 ★★ kill criterion 的正确使用:NOJIT 摘除

审计说「`pyfp3d.NOJIT` 属性无消费者」—— **错了**:`tests/A/test_A01:19` 读它
(`assert pyfp3d.NOJIT is False`)。审计的「无消费者」是对**库内**成立、
对**测试种子**不成立。预注册的失败标准(删了会让 A01 红)与 kill criterion
在删除前拦下了它 ⇒ **摘除,属性保留**。
★ 这正是这族纪律的用途:**静态可达性名单是线索,不是处决令**;
删除前的逐条 grep + kill criterion 是最后一道。

### 4.3 登记(不修,列清单)

- `M_CRIT_DEFAULT` / `UPWIND_CONST_DEFAULT`(isentropic.py:25-26)随
  `upwind_factor` 的删除而无读者 —— 模块级公开常量,不在预注册清单,
  删不留待裁决(登记);
- pyflakes 未用局部变量余 6 处(`wing3d.py:448`、`wake_cut.py:577`、
  `coloring.py:162`、`coupling.py:575`、`closures.py:615`、`dirichlet.py:150`),
  留待各模块下一次触诊;
- 翼身生成器族的「全覆盖」仍不在 —— A56 是**最小保险**(改动会响),
  不是覆盖;覆盖归 W3.8 之后。

## 5. 判定

**PASS** —— C1–C3 全过,C4 见基线行;kill criterion 触发一次(NOJIT)
并照章摘除 —— 这正是它存在的目的。

一句话:**pyfp3d/ 少了 338 行死代码,翼身生成器第一次有了「改动会响」的
最小保险;审计的一个结论(NOJIT 无消费者)被删除前的核查推翻 —— 流程赢了。**

## 6. 下一步

1. **W3.8 覆盖率实跑**(提交后干净树上,NOJIT;A+B+F 起步,按耗时决定扩不扩)。
2. 算例面收口(`cases/demo` 分类:6 个现存 demo 的真实读取面只有 p11 ——
   README 的 "9 of 33" 早已不实)。
3. 小尾巴:W2.5 / W2.6 / D08 footer / D03 xfail 退役 / INDEX:297 / C05 precond。
4. 裁决项:`M_CRIT_DEFAULT`/`UPWIND_CONST_DEFAULT` 删留(§4.3)。

## 基线行

**首轮收口全套红了 7 条证据锁 —— 已归因,不是回归**:触发器是删除引发的 numba
缓存重建(编译/缓存来源敏感性),与删除内容无关(逐行 diff 无运算改动;
换回旧 pyfp3d 绿 / 换回新的也绿 / 双 A/B 齐备)。落定后全套 = **818 passed /
2 failed / 32 skipped / 2 xfailed @1679.39 s**,两条红 = D08/D13(IBL 松耦合腿对
**进程历史**敏感,单跑绿且复现、全套红且两次全套互差 1e-5 —— 值锁在全套里
原理上守不住)。
**完整归因:[20261011-1800-process-history-determinism.md](20261011-1800-process-history-determinism.md);
D08/D13 锁形状的处置是该文件 §5 的裁决项(裁决①的前提已过期)。**
