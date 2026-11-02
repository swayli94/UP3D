# 20261102-2030 — **预注册：第 5 轮 = tet medium r_far 20/30 不收敛归因（S4）**

**性质**：研究轮预注册；必须先于正式测量提交。**不改库、不改判据、不改默认值、
不刷新任何 committed gate 证据。**
**由来**：phase 7 立项 §2 S4 / §3 轮 5；原始观察见
`phases/p6/docs/dev_phase_six/20261012-1000-w25-w26-verdict.md` §4.2 与 §6c，
立项出处见 `20261012-1300-c10-hex-cylinder-prereg.md` §2.3。
**分支**：`kimi/phase7-charter`。

---

## 1. 问题与已有证据边界

问：**NACA0012 tet medium 在 M0.80/α1.25、D05 生产配方下，r_far=20/30
重生成变体撞 400 步帽（|R|~2e-7）的机制是什么？**候选不是互斥标签：

1. 远场半径改变带来的物理/边界敏感性；
2. 重生成同时改变离散网格，在 medium 跨声速敏感腿上触发有界停滞或极限环；
3. 冷启动/初值路径选择；
4. 熵修正及其冻结路径与上述网格实现的交互。

原始脚本与逐腿 CSV 未提交，当前只有判定文书的一条观察，**所以本轮先复现，
不把 |R|~2e-7 当面值，更不先把撞帽命名为极限环。**

迁移机现状单独登记：新环境完整/定向测试有 11 项既有证据锁或散布断言红；
D05 最新差异是 M0.80 medium `cp_rms_upper` 相对差 2.43e-6，并非上述撞帽。
本轮不借 S4 重基线、不修改这 11 项，也不把新机器与旧机器的单次差直接称为
换吸引子。

## 2. 固定出处与测量台

- 环境：`PYTHONNOUSERSITE=1`，`NUMBA_NUM_THREADS=16`，`OMP_NUM_THREADS=16`，
  `OPENBLAS_NUM_THREADS=16`；BLAS 由 `pyfp3d` 运行时钉 1，脚本断言实际值。
  全部求解在同一暖缓存进程模式下执行；记录 Python/NumPy/SciPy/Numba/Gmsh 版本、
  OS/CPU、线程数、缓存模式与每条腿的墙钟。
- 几何：只动 `r_far ∈ {15,20,30}`，level ∈ {coarse,medium}，其余参数逐字取
  `cases/meshes/naca0012_2.5d/generate_naca0012.py::LEVELS`。用
  `naca0012_wake_2d` + `extrude_single_layer` 在研究目录的 gitignored 工作区生成，
  **不覆盖 `cases/meshes/`**。参数、节点/单元/边界组计数与 SHA-256 入 CSV。
- r=15 两级新生成网格必须与已上传 committed mesh 的 manifest SHA-256 一致；
  每个变体同参数生成两次，SHA-256 必须一致。网格本身不是提交证据，生成脚本、
  参数与指纹 CSV 才是可复算出处。
- 流场：M0.80/α1.25，`solve_newton_lifting`，D05 `RECIPE` 原样
  （direct、400 步上限、5 次 Picard seed、原冻结设置）。每条腿保存残差/F/Gamma/
  clamp/步长/线性求解与熵修正诊断，以及 cl/cd/x_shock/Cp。

## 3. 实验矩阵与判据（事前写死）

### C1 — 基线复现与指纹（6 腿）

跑 coarse/medium × r=15/20/30 共 6 条冷启动生产腿。

| 判据 | G-DOMAIN（相反结果落在哪） |
|---|---|
| 6 腿都有完整摘要、网格指纹和逐步历史；r=15 指纹与 manifest 相合，所有变体重复生成 SHA 相合。原观察若复现，r=20/30 medium 应撞帽而 r=15 medium 不撞帽。 | 原观察不复现也是合格结果：记 **RECORDED**，说明它在当前 Gmsh/机器/出处态不可复现；不得改网格或配方去追撞帽。指纹不合则停止流场测量，判 **UNDEFINED（仪器出处不可判）**。 |

### C2 — 失败分类器及其牙齿

分类器只看事先定义的逐步量，不用 cl 倒推机制：

- `converged`：驱动器 `converged=True` 且零 clamp；
- `clamped`：末态 `n_limited+n_floored>0`；
- `bounded_stall`：至少 201 个残差点、末值 >1e-10，最后 200 点首尾各 20 点
  中位数之比在 [0.5,2]，且 tail max/min <10；
- `periodic_limit_cycle`：满足 `bounded_stall`，且周期 p∈[2,8] 的按相位折叠
  RMSE ≤0.25×tail 极差；p=1 只报 fixed-point stall，不冒充极限环；
- 其余为 `nonconverged_other`。GMRES stall、freeze、churn、accept_reason 全部分列，
  不合并成一个“机制分数”。

| 判据 | G-DOMAIN |
|---|---|
| `--self-test` 用合成收敛、clamp、有界周期、非周期停滞历史验证分类；真实腿全部有分类与原始量。 | 有界但不周期 ⇒ 只报 `bounded_stall`；周期但发散/ clamp ⇒ 不报本类极限环。分类器无法区分时记 other，不靠肉眼改名。 |

### C3 — 同网格诊断臂（只攻原登记的 r=20/30 medium）

对每条未收敛的 r=20/30 medium 基线腿，各至多加两臂；r=15 若意外失败，只登记
出处异常，不追加诊断臂：

1. **continuation**：同一网格、同一目标工况和最终严格容差，走既有
   `solve_newton_transonic(m_start=0.70,dm=0.05,dm_min=0.01)` 向上 Mach 延续；
2. **entropy-off**：同一网格、同一冷启动生产路径，仅 `entropy_correction=False`。

两臂都是诊断，不是生产修复，不比较谁“调绿”。

| 判据 | G-DOMAIN |
|---|---|
| continuation 收敛而冷启动不收敛 ⇒ 支持“初值/延续路径敏感”，不是远场半径本身使离散方程无解；entropy-off 收敛而 entropy-on 不收敛 ⇒ 支持“熵修正路径参与失败”。两臂均须报最终物理解、收敛状态、成本和警告。 | 臂不收敛是有效负结果；收敛到不同 cl/Cp 结构不得称为同一吸引子，也不许据此换生产默认值。两臂都收敛只说明两条路径都参与可行域，按各自结果分别归因，不虚构唯一原因。 |

### C4 — 可归因语句与不可判边界

| 结论形状 | 所需最小证据 |
|---|---|
| **medium 网格实现触发的求解路径敏感** | r=20/30 medium 未收敛、对应 coarse 基线收敛，且至少一个同网格诊断臂改变收敛状态；明确是哪条臂，不把半径与网格实现混写成单一“r_far 效应”。 |
| **冷启动/continuation 敏感** | 同一失败网格 continuation 收敛、冷启动不收敛。 |
| **熵修正参与** | 同一失败网格 entropy-off 收敛、entropy-on 不收敛；只声称参与，不声称已证明具体冻结单元。 |
| **UNDEFINED** | 只有 r 扫描差异，或同网格两臂均不改变失败状态，或出处指纹不合；半径和网格实现仍混杂时不作机制claim。 |

**G-DOMAIN**：若 coarse 也在 r=20/30 全部失败，不能自动改判“纯远场物理效应”；
只有上述同网格对照才有资格移动机制结论。所有结论只适用于本离散族、M0.80/α1.25、
本配方和已记录出处，不外推到 M6、fine 或其他跨声速工况。

### C5 — 代价与证据

每条腿报 wall_s、n_newton、n_gmres_total/stalled、n_refactor、最终残差、
cl/cd/x_shock；continuation 另报各 Mach level 成本。证据落
`bench/studies/tet_rfar/`：run.py、README、mesh_manifest.csv、summary.csv、
逐腿 history CSV、probe_levels.csv（若有）、分类自测 CSV、残差/Cp PNG。
目录登记进 `bench/studies/README.md` 一次性研究名单，phase 7 收口归档。

| 判据 | G-DOMAIN |
|---|---|
| 每个判定数字可在提交 CSV 查到，命令可复算；无数据格写 RECORDED，不编数。 | 超时/中断腿保留已有历史并记 RECORDED，不删除失败证据。 |

## 4. 时间盒与 kill criterion

- **最多 10 条 fresh 顶层腿**：6 基线 + r=20/30 medium 各至多 2 诊断臂；
  continuation 内部 level 如实另计成本，不借此增加新工况。
- 求解总墙钟上限 **60 分钟**；到时不再启动新腿。只允许 coarse/medium，不跑 fine，
  不新增 r 值，不扫 seed、freeze_tol、EW、direct_refactor_every 或 pyamg 参数。
- 任一 r=15 指纹不合、网格重复生成不合或分类自测失败 ⇒ 停止正式测量，先修仪器；
  无法在不改生产库/网格的前提下修复 ⇒ **UNDEFINED**。
- 超时间盒仍只有混杂证据 ⇒ **UNDEFINED/RECORDED**，不得延长时间盒换答案。

## 5. 明确不做与收口

- 不改 `pyfp3d/`、`tests/`、`cases/gates/`、`cases/reference_data/`、生产网格和
  `environment.yml`；不处理 S5，不顺手修迁移机 11 项红。
- 不建永久门，故无新 G-CADENCE；本轮只属于研究轮，默认套件测试数不变。
- 收尾跑研究脚本 `--self-test`、语法检查和 `tests/F/test_F07_repo_invariants.py`；
  无库改动，不重跑贵重全套，不宣称迁移机全套转绿。判定写 PASS/FAIL/RECORDED/
  UNDEFINED，台账一行并如实列出当前机器测试边界。
