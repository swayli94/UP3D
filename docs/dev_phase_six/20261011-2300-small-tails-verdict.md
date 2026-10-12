# 20261011-2300 — **判定:工程小尾巴批(T1–T4)—— D03 xfail 退役、D13 cf 口径更正、D08 PNG 前提已自愈、W3.8 覆盖率第一次实测**

轮次文件名:`20261011-2300-small-tails-verdict.md`
阶段 / 门:无门(工程/卫生)—— 预注册:
[20261011-2100-small-tails-prereg.md](20261011-2100-small-tails-prereg.md)
提交:分支 `claude/phase-six-audit`(本文件随判定提交)

---

## 1. 目的

收掉四条登记的小尾巴。判据逐项见预注册 §1;**不改任何判据的数值**。
★ 执行中**发现并修复了一个仪器缺陷**(T4 的第一次测量被它污染):
`PYFP3D_NOJIT` 的交换只覆盖 kernels/viscous 的 7 个模块,**metrics / isentropic /
coloring 三个有 njit 的模块不响应** ⇒ 它们的 njit 函数体在覆盖率下失明
(metrics 读成 10 % —— 仪器失明,不是覆盖事实)。本轮把交换扩展到这三模块
(`@numba.njit(` → `@_njit(`,与 kernels/entropy.py 同款),再重测。

## 2. 改了什么

- **T2**:`tests/D/test_D03_naca0012_m080_shock.py` 的非 strict xfail **退役** ——
  标记覆盖的两件事都已闭合(① seed 失败 2026-08-05 修;② 线程翻面 = BLAS,
  H30 归因 + 已钉死);门控实测(2026-10-10 全集 XPASS + 本轮复跑 **4 passed**)。
  该腿回到常规门控断言:再红就是真回归。
- **T3**:`cases/gates/INDEX.md` D13 条目的 cf 叙述 —— 数字更正到现行 committed
  CSV 的站位口径(转捩区 x=0.10:XFOIL 0.0069 / 我们 0.0042;尾缘 x=1.00:
  我们 0.0025 / XFOIL 0.0010),**两个定性结论都成立**;注明 recon 时代的旧数是
  松散引用。
- **T1**:D08 PNG footer —— **前提已自愈**:第 44 轮的暖态刷新已把第 36 轮的
  footer 字符串渲进 PNG;本轮重跑导出腿验证:PNG 与 committed **逐位一致**、
  summary.csv 逐位一致 ⇒ 无事可做,核销登记。
- **T4**:W3.8 覆盖率实测(A+B+F,`PYFP3D_NOJIT=1`,`coverage --branch`)——
  数字见 §4.2。

## 3. 怎么测的

- T1:暖态导出腿复跑 + 与 HEAD 的逐位比;
- T2:门控跑该文件(`PYFP3D_TRANSONIC_GATES=1`);
- T3:对照 committed CSV 逐数;
- T4:`PYFP3D_NOJIT=1 coverage run --source=pyfp3d --branch -m pytest tests/A tests/B tests/F`。
- ★ 自伤一次:T1 的第一次重跑把 D08 CSV 的 `n_threads` 写成 `unset` ——
  **我在第 36 轮纪律化的那条规则(「证据腿的环境变量一律 export 或同行前缀,
  不许 `&&` 后裸赋值」)被我自己又违反了一次**。物理列逐位不动(numba 线程
  无关早已是实测),出处列错;正确 export 重跑后逐位复原。**第二次了:
  这条规则写进了判定文件还是会被违反 ⇒ 它值得进 CLAUDE.md 的操作危害清单
  (和 `pgrep` 那族同一层)。**

## 4. 结果

| 项 | 结果 |
|---|---|
| T1 D08 PNG/CSV | ✅ 与 HEAD 逐位一致(前提是旧的,第 44 轮已渲过) |
| T2 D03 xfail 退役 | ✅ 门控 4 passed(退役后常规断言绿);再红 = 真回归 |
| T3 INDEX cf 更正 | ✅ 每数可指(站位、列) |
| T4 覆盖率实测 | 见 §4.2(只测,不设阈值) |

### 4.2 W3.8 覆盖率(A+B+F,NOJIT,branch;命令:`PYFP3D_NOJIT=1 coverage run --source=pyfp3d --branch -m pytest tests/A tests/B tests/F`)

**TOTAL 80 %**(8737 语句 / 2400 分支;632 passed / 17 skipped / 0 failed @1:15:47)。
★ **仪器修复前后**:metrics 10 %→**86 %**、isentropic 36 %→**72 %**、
coloring 62 %→**83 %** —— 三个差值全部是仪器失明,不是覆盖缺口。

强:kernels(entropy 96 / jacobian 98 / residual 96 / upwind 79)、
viscous/closures 98、ibl3 94、solve/newton 90、solve/picard 90、wake_cut 92。
真低(读时注意口径 —— 本轮只跑 A+B+F):post/shock 15 %、solve/wall_correction 19 %、
meshgen/planar 26 %、meshgen/wing3d 26 %(这两族由 cases/ 生成器消费,不在测试收集)、
post/section_cut 48 %、viscous/coupling 56 %、viscous/tight_driver 60 %
(这几族的消费者在 C/D 层 —— 低读数是**范围**属性,不必然是缺口;
全套 NOJIT 覆盖率是更全的图,成本数小时,登记未做)。

## 5. 判定

**PASS** —— 四项各自的判据全过;全套收口见基线行。

一句话:**小尾巴清零;D03 的非 strict xfail 退役后该腿回到常规断言;
W3.8 的覆盖率第一次是实测数字而不是印象。**

## 6. 下一步

1. **W2.5**(C06 的 CL_CONTRACTION 落在网格实现噪声里 —— 测噪声带 + 阈值重定,
   阈值必须测出来)与 **W2.6**(hex 族求解交叉门 —— hex 网格读进来是
   15480 个 tet,同拓扑两组,可直接求解)—— 下一轮,带预注册。
2. 裁决项仍开着:C05 precond(direct vs amg)、`M_CRIT_DEFAULT`/
   `UPWIND_CONST_DEFAULT` 删留。
3. 研究项:xtr030 多吸引子;M0.50 cl 阶 1.24 归因。

## 基线行

**820 passed / 32 skipped / 2 xfailed / 0 failed @1559.57 s @16 线程(BLAS=1,暖机后)**
—— 与第 44 轮逐项相同;本轮改动(D03 退役、INDEX、A01 skip、NOJIT 扩展三模块)
未移动任何答案。
