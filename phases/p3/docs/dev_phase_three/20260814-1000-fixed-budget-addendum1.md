# addendum #1 — ★★★ **那四条全尺度腿是 `tip_cap="flat"`**;上一轮的头条要按项目自己的标准规则降级

**在写任何代码之前提交**(GV5.2/GV5.4 的 addendum 惯例)。
修订对象:[20260814-0900 预注册](20260814-0900-fixed-budget-allocation-prereg.md);
被更正对象:[20260814-0700 判定](20260814-0700-p2-reassessment-verdict.md)。

## 1. 事实(读原文,两行代码)

- `phases/p2/bench/run_tip_allscales.py:49` → `from run_le_response import SCRATCH, build, le_geometry`;
- `bench/run_le_response.py:53` → `BASE = dict(r_far=15.0 * MAC, tip_cap="flat", embed_wake=True)`。

⇒ **那四条腿(`Tm1_coarse` / `Tp5_h020` / `T0` / `T1`)全部建在平截帽(flat)网格上。**

而生产族是**圆角帽**:`cases/meshes/onera_m6/generate_onera_m6.py:164` →
`tip_cap=("flat" if level.endswith("_flat") else "round")` ⇒ **`medium.msh` 是 round**,
我前三轮的 LE 读数(池化 0.2422 / η0.90 0.3306)都在 round 上。

## 2. ★★★ 项目自己的标准规则直接管住它(CLAUDE.md,逐字)

> **Round tip is the default; flat must be asked for by name.** … because P13/G13.3 measured the flat
> cap's sharp convex edge **DIVERGING under refinement** (peak-Mach exponent p = +0.321), so
> **any refinement-based claim on a flat-cap mesh has a false premise**. Measured consequence:
> the LE-band convergence order is **0.37 on flat against 0.87 on round**. Only P13 and M5, which
> exist to MEASURE flat, use the `_flat` levels.

⇒ 上一轮我把那四条腿的 **−54.9 %** 当作"**唯一被测到能降低对实验误差的杠杆**"，
而它是一条**基于加密的论断**、建在**平截帽**上 ⇒ **按这条规则,它的前提为假。**

★ 而我上一轮写的限定是"**不许接到 HEAD 上(跨时间)**" —— 那条对,但**不是最要紧的那条**:
真正的问题是**跨网格族(flat vs round)**,而我没查。**第 5 问我问了四项,漏了"同一种几何"。**

## 3. 什么留下、什么不留下(分开写,不含糊)

**留下(仍是事实)**

- 四条腿**都收敛**(0 lim / 0 flr,|R| 1e-13…1e-15),LE 上 RMS **单调下降** 0.3899→0.2704→0.2358→0.1758;
- GS2.1 §7.2 那条事前预测(**单调** + coarse 明显更差)**在它自己的平截帽框架内被满足**;
- 那条义务的**成本**确实被我记错了(4 011 s vs 实测 7.2 s),**推迟四轮**这件事不因帽型而改变(勘误 1 成立)。

**不留下(必须降级)**

- ★★ **不得**把它当作**生产(round)**的加密/分配行为 —— 那是 **UNMEASURED**;
- ★ 我派生的逐段阶 **0.902 / 0.476 / 0.725** 必须标注 **flat-cap**;★ 注意它们**跨在**已记录的
  **flat 0.37 / round 0.87** 两侧,所以既不能说"和 flat 一致"也不能说"和 round 一致" —— **就是未测**;
- ★★ 上一轮 §6 的**理由 2**("唯一被测到能降误差的杠杆不是 P2")**削弱**:它现在只在 flat 上成立。
  ⇒ **R2 的结论不变**,但**支撑改了**:仍成立的是**理由 1**(A1:前提未成立,证据来自两个都不是 M6
  且彼此不一致的算例)与**理由 3**(A3 穿透 4/4 + A0 明文禁止项 + design.md 自列为另一个阶段)——
  这两条与帽型无关。★ 而**诚实的后果是**:证明"分配/加密有效"的证据**变弱了**,
  于是"提高阶数"这条路**更开而不是更关** —— 这一点必须说,即使它对我上一轮的结论不利。

## 4. 本轮设计的改动(这就是 addendum 的用途)

1. ★★ **本轮全部网格用 `tip_cap="round"`**(HEAD 默认、生产族)——
   **不复用 `run_le_response.build`**(它写死 flat);本轮自己建 builder,并**断言** `tip_cap == "round"`;
2. ★ **参照曲线只用本轮 round 上重测的两点** —— 预注册 §3 已写"不接 2026-08-01 的
   0.3899/0.2358",现在它有了**第二个、更强的**理由(不只是跨时间,是跨族);
3. ★★ **新增守卫 G-C(帽型)**:脚本必须打印并落盘每条腿的 `tip_cap`,且断言全为 `"round"`
   —— 这一条守卫**本轮之前不存在**,而它正是这次漏掉的东西;
4. **判据、band、kill criterion 一字不改**(B1/B2/B3、r ≥ 1/3、±10 %、质量与预算守卫全部保留)。
   ★ 唯一变的是**参照曲线的来源**,而它本来就要求 HEAD 重测。

## 5. ★ 附带:这也解释了一处我以为是巧合的事

上一轮我注意到"四级阶梯 T0 = 350 718 tets **恰好等于** committed medium 的 350.7 k",
并据此说"T0 就是 committed medium(把 LE/TE 场拆开的控制)"。
⇒ **那句话是错的**:两者**帽型不同**(flat vs round),单元数相近是巧合而非同一张网格。
★ 教训与本阶段记过的同族:**数字对上不等于东西一样**(第 5 问)。
