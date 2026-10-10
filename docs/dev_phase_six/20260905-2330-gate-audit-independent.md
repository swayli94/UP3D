# 20260905-2330 — **独立审计**:gate 体系的合理性、数据齐全性与覆盖递进性

**性质**:分析 note,不建立任何门、不动任何判据。
**动机**:为后续 debug 提供指导与验证基线 —— 现有 gate 体系是否设计合理、数据是否齐全、
能否**全面、递进、有针对性**地检验求解器各功能在各场景下的正确性。
**与既有文档的关系**:`20260824-0400-gate-taxonomy-analysis.md` 是重编号**之前**的归类分析,
当时 C 类只有 3 文件、D 类只有 3 文件、CFL3D 完全没有。此后它建议的算例(C05–C09、D05–D13)
已全部落地。本审计是对**落地后现状**的独立复核,不复读旧结论;凡引用旧裁决处显式标注。

---

## 0. 审计范围与方法

四个盘点(全部逐文件做,非抽样):

1. **tests/** 全量:A 组 54 文件(~451 测试)、B 6/35、C 9/56、D 13/100、E 3/5、F 7/59,
   含支撑件 `_tol.py` / `_gate_evidence.py` / `conftest.py` / 各共享 case 模块。
2. **bench/** 全量:12 顶层脚本 + `gate_results/` 81 个跟踪文件 + `studies/` 7+1 目录。
3. **cases/**:`gates/` 20 门目录、`reference_data/` 8 数据集、`meshes/` 网格族、`demo/` 6 个。
4. **pyfp3d/** 功能面:逐子包列能力清单,与测试覆盖交叉核对。

外部基准(网络调研,用于对照找缺口):
- AIAA G-077 / Oberkampf & Trucano 的 V&V 分层(code verification → solution verification →
  validation)([NASA GRC V&V 教程](https://www.grc.nasa.gov/www/wind/valid/tutorial/overview.html)、
  [Oberkampf 2002](https://www.osti.gov/servlets/purl/793406));
- Roache/Roy 的 MMS 与阶次验证纪律([Veluri & Roy, FV 码 MMS 全覆盖](https://www.aoe.vt.edu/content/dam/aoe_vt_edu/people/faculty/cjroy/Publications-Articles/VerificationJournal-Revised.pdf));
- [NASA Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/) 的湍流验证算例族
  (ZPG 平板 / bump / hump / NACA0012 / RAE2822);
- SU2 的回归基线实践(TestCases + 期望值比对,见
  [Running Regression Tests](https://su2code.github.io/docs/Running-Regression-Tests/))。

**亲验声明**:§5 的 P0 项与 §4 的若干关键判断由我本人直接读码/grep 复核
(D10、D07-m050、C09、A01、`_tol.py` 的 identity marker、`bench/run_g82_anchor_check.py`);
其余来自逐文件盘点记录,均带路径与行号。**未重跑全套测试**(ungated 基线按 README 为
479+12+2 @8t);本 note 的全部结构性结论不依赖重跑。

---

## 1. 总评(三句话)

1. **体系设计是合理的,且在三个方面超前于常见 CFD 项目实践**:三层参照按模型层级对齐
   (无粘→Euler、有粘→RANS、实验只记偏置)、"参照物分辨底噪先行"(|SST−SA| 带、参考自身
   网格收敛差、实验采样间距)、证据锁的四层防静默守卫。A→B→C→D→E→F 的**递进结构成立**。
2. **数据基本齐全**:8 个 reference 数据集全部在仓、可再生;新克隆除 ONERA M6 网格
   (~30 s 生成脚本)外不需要任何外部物。证据 CSV/PNG 齐全度高(20 门目录,14+ 门有锁)。
3. **但存在 1 个实测的 P0 一致性缺陷(干净 clone 会硬红)、若干"平时跑覆盖静默消失"
   的结构性风险、以及 6 个明确的覆盖缺口**(§4)。其中两个缺口(任意方向自由流、全局
   质量守恒)是**便宜且高价值**的补强。

---

## 2. 现状盘点(一张表)

| 组 | 文件/测试数 | 角色 | 真值来源 | 证据 | 常跑? |
|---|---|---|---|---|---|
| A | 54 / ~451 | 机制层(装配/通量/迎风/熵修正/网格/IBL 闭包/耦合机构/驱动性质) | 解析+FD+位相同+性质断言 | 按设计无 artifact(判据即断言) | 全常跑 |
| B | 6 / 35 | 后处理(表面提取/恢复/Cp 度量/截面 cl) | 解析二次场/合成/AGARD(B05) | 无 | 全常跑 |
| C | 9 / 56 | **验证组:解析解 + 阶次** | MMS/球/圆柱/KT/Ringleb/喷管精确解/PG 渐近 | 8/9 有 committed 锁;**C04 无** | 常跑 |
| D | 13 / 100 | **外部参照组** | Hess–Smith/PG-KT/CFL3D Euler+RANS/实验/XFOIL | 9 门有锁;**C04 状:D04/D09/D10 无** | 常跑(D03/D04 重腿 gated) |
| E | 3 / 5 | 提交值锚点(自述"中转站") | 代码自身历史值 | — | **全 gated** |
| F | 7 / 59 | 仓库卫生(provenance/环境/仪器/spec/invariants/CFL3D 生成器) | 行为式 | — | 常跑 |

支撑机制:`_tol.py` 的 REL_TOL=1e-12 纪律;`_gate_evidence.py` 的新鲜计算 vs committed CSV
锁(1e-6,四层守卫:缺列必响/重复键必响/n>0/逐列从未被比较必响);`PYFP3D_GATE_FIGURES=1`
才写盘;`PYFP3D_TRANSONIC_GATES=1` 开重腿;`PYFP3D_NOJIT=1` 第二车道。

---

## 3. 设计合理性 —— 做得对的(对照业界)

1. **三层参照模型对齐(2026-08-24 裁决三)**:FP→Euler 设门、FP+IBL→RANS 设门、
   FP+IBL→实验只记偏置。这正是 AIAA G-077 "在正确模型层级上比较"的思想,而且本项目
   比常见实践多走一步:**每道门的容差不小于参照自身的分辨底噪**(RANS 带 = |SST−SA|
   分歧,D08;激波可分辨性 = 参考自身 L2→L3 差,D05/D06 的镜像断言;实验采样间距
   0.040–0.050 c 使 5/7 站位激波不可分辨,D07)。这与 NASA TMR "先定底噪再定容差"
   一致,而多数项目根本没做。
2. **C 组解析解的选型是教科书级的**:MMS(Laplace 阶 1.9+)、球与圆柱(壁面精度+阶)、
   Kármán–Trefftz(**唯一能把 Kutta 自选 Γ 对精确值检验**)、Ringleb(**唯一 2-D 跨声速
   精确解,且是全速势方程本身的解**,含剂量-响应与 C=0 必不收敛的反例腿)、准一维喷管
   (**唯一带激波的精确解**,激波位置亚单元精度 + "收敛到机器零却错 40 格"的反例腿)。
   五个机制轴(不可压/可压/升力/激波/Kutta)各有一个精确锚。Ringleb+喷管两条是
   20260824-0400 点名的最高价值项,已建成且判据立在"位置"上(符合当时预注册的告警)。
3. **证据锁机制优于 SU2 式基线比对**:SU2 是"代码对期望值";本项目额外防"只写不读"
   (2026-08-28 实测 19 个 CSV 中 16 个只写不读后补齐锁)与"值==来源恒真"等四类
   静默通过模式(`_gate_evidence.py` docstring 记了四个实测踩过的坑)。
4. **能力边界门(D09/D10)锁"失败的形状"**:松耦合第 0/1 次外迭代即死被断言
   (M²<0 的站位数、mdot_max 量级),而不是含糊的 skip 或 xfail。这是 Oberkampf 的
   "domain of applicability" 显式化,业界少见。
5. **诚实性机制是体系级资产**:`failure_modes.classify_failure` 五模式分类器(永不报
   裸 conv=False);收敛旗标线程依赖的显式处理(D05–D08 不断言 converged,CSV 带
   `n_threads` 列——20260905-0100 记录过"差点把线程数焊进门里");三条 strict xfail
   全部是**记录在案的能力限制**(C02 的 G1.3、C03 的球 2%、E01 的熵修正 medium 不可锚)。
6. **F 组卫生超出一般 CFD 项目**:provenance 防篡改、网格 manifest sha256、
   repo invariants(门↔证据目录双向一致、活代码 import 可解析)、CFL3D 生成器行为式
   测试(F08,不跑二进制)。

---

## 4. 缺口分析 —— 对照业界标准的差距

### 4.1 Code verification(Roache/Roy 口径)

- **G1 · 非线性全速势算子没有 MMS**。C01 只覆盖 Laplace(线性)算子;密度律、人工密度
  迎风、熵修正、透射通道全部没有源项级制造解验证。Ringleb(C08)只测 R(exact) 残差阶
  +生产配置误差,不等价于 MMS(MMS 能逐层剥离每个项的实现错,见 Veluri & Roy)。
  当前替代物是 A03/A06/A41 的 FD/JVP 单元级锁 —— 覆盖"导数对",不覆盖"算子在非平凡
  场上的整体阶"。**工程量中等,价值高,建议 P1 后期或 P2。**
- **G2 · 没有全局质量守恒门**。求解的是守恒形式 ∇·(ρ∇φ)=0,但没有一道门断言
  "远场进出质量流之差 = 内部源(透射/尾迹片)"。A11 是单元**几何体积**守恒,不是质量。
  喷管(C05)与任意翼型 case 都可以便宜地加一条质量平衡残差断言。**便宜,建议 P1。**
- **G3 · 自由流保持只测 φ=x 单方向**(亲验 `test_A01:116` 只构造 `phi = nodes[:,0]`)。
  业界惯例(SU2/CFL3D 回归套件)是**任意方向均匀流**(旋转来流/旋转网格)——
  它抓的是"度量/装配只在坐标轴对齐时正确"那族 bug,正是非结构网格项目的常见病。
  **极便宜,建议 P1 第一项。**
- **G4 · 解验证没有 GCI 纪律**。参照侧(CFL3D)每个数据集带 `grid_convergence.csv`
  (implied order/asymptotic/error_bar),且 F08 已经有一个经合成数据验证的
  `implied_order` 估计器;但**我们自己的腿**只有"收缩比 <0.7 / ≥1.5×"这类经验界,
  不报告 observed order、不查渐近区间。Roache 的 GCI 是解验证的标准口径。
  **建议:把 F08 的估计器复用到 C04/D 类阶梯上(工具已存在,只差接线),P1。**

### 4.2 粘性/湍流验证(NASA TMR 口径)

- **G5 · 没有 ZPG 湍流平板的端到端验证**。TMR 算例族的第一条(平板 Cf vs Re_x)是
  湍流闭包最基础的外部锚。现状:层流侧有 Blasius(A35/A36 锁与 Blasius 的**偏差本身**,
  A37 对 Falkner–Skan ODE oracle);湍流侧只有闭包常数锚(A39:XFOIL 常数、ZPG 固定点
  2.590433、H 发展锚点)——全是**点锚**,没有"strip2d marcher 在湍流平板上推进的
  Cf 曲线对关联式(White / TMR 数据)"的端到端门。GS4.1 第 15 轮已记录"E-CF 关联式
  无出处 ⇒ 降为 RECORDED"的教训,所以这次补门要**先用有出处的关联式**(White 的
  Cf≈0.455/ln²(0.06·Re_x) 或 TMR 公开数据)。**建议 P1。**
- 转捩只有强制转捩,e^N 未实现 —— 这是**能力边界**(已记录),不是测试缺口。
- D13 对 XFOIL 的 (c)(d)(e) 锁 RECORDED 负方向(加密变差、上下表面分不开,根因指向
  湍流闭包 H 族)——处理正确,符合"参照不确定度无边界就不能归因"的纪律。

### 4.3 功能 × 场景覆盖矩阵

**场景轴**(实测覆盖):Mach 0 / 0.3 / 0.5 / 0.72 / 0.725 / 0.730 / 0.75 / 0.778 / 0.80 /
0.803 / 0.8395 / 0.84 —— 好;α 含负攻角与近失速(实验 12.86° 记录在 D11)—— 好;
几何:翼型×2、球、圆柱、喷管、Ringleb 流道、M6、翼身 —— 好。

**功能轴上的薄/无测试路径**(逐条核对过 tests/ 引用):

| 功能 | 现状 | 处置建议 |
|---|---|---|
| `precond="ilu"`(newton.py:1336 实装) | tests/ 零引用 | 补一条 A 类装配-求解冒烟,或声明弃用 |
| `n_kutta_polish` | 仅 D04 一处 gated 调用 | 薄;A32 只测 intermediate_tol。可接受,低优先 |
| 闭体 BoR 松耦合端到端 | A42 只有合成 icosphere 接线;GV3.3 的真实案例已归档 | 若要保留 BoR 能力声明,补活门;否则删能力声明 |
| 3-D IBL 横流(B/Ψ)耦合态 | 单元级有锁(A34/A35);耦合态只有 bench 记录(max\|B/A\|≤0.072) | 登记,等紧耦合复活时一并立判据 |
| 紧耦合翼尺度失败模式(IBL 地板/无二次盆) | 机构正确性有锁(A44–A47);收敛性结论全在归档 bench studies | **建议把 GV5.x 的三条实测失败形状各立一条"断言即证据"门(D09/D10 同款形状)** |
| `wake_transpiration` 端到端效应 | A51 锁源构造;GV6.2 效应测量已归档 | 效应不显著已裁决;保留 plumbing 锁即可 |
| wake-free 翼身族 / `make_inboard_clip` | LS 路线删除后消费者消失,孤儿代码 | 删除或标注,避免审计噪声 |
| NL7301 双元素(多尾迹) | legacy 资产,无门;wake_cut 只支持单尾迹链 | 明确"不支持多段翼型"写进能力边界 |

**design.md 声称但未实现**(文档-实现偏差,不是测试缺口):侧滑 β
(`dirichlet.py` 自述 quasi-2D beta=0)、C_m(`wall_force_coefficients` 无此项)、
Trefftz 诱导阻力(grep 零命中)、非升力 Newton 入口(G10.1 开项)。
⇒ **建议修正 design.md 的能力表,或把它们列入路线图**;现状会让审计/新人误判能力面。

### 4.4 结构性风险("平时跑"覆盖静默消失族)

1. **gated 层默认不跑**(A21/A28/A30 各一腿、D03/D04 重腿、E01/E02 全部):
   `bench/run_capability_locks.py` 有 NOT_COVERED 清单打印(好),但**没有任何机制
   保证 gated 层的新鲜度** —— 它们可以在两次阶段边界之间静默腐烂。建议:阶段收口时
   强制跑(已有惯例),并考虑给 gated 腿加"上次通过日期"记录。
2. **NOJIT 车道下紧耦合 FD 覆盖静默消失**(A42/A44/A45/A46/A47/A51 重腿 skip)。
   记录在案,但跑 NOJIT 全套的人容易误以为覆盖等价。
3. ★ **C09 的物理断言整个挂在 `skipif(not gate_figures_enabled())` 下**(亲验
   `test_C09:28-77`):平时跑只比对 committed CSV(4 个数),**PG 2% 判据本身平时不执行**。
   与其他 C 门形状相反(它们判据常跑、图 opt-in)。而计算其实已在共享 fixture 里 —
   **把那两条 assert 移出 skipif 几乎零成本**,建议 P1。
4. **C03 活跃门依赖 `cases/demo/p11_curved_walls/results/` 的 sweep CSV,缺则 skip**
   (`test_C03:216`):证据文件丢失的 clone 上活跃 G1.6 门静默消失(只剩 xfail+sanity)。
5. **E01/E03 的升格未执行**:E01 自述"参照到位后升格 D 类,E 是中转站";
   D05–D07 的 CFL3D 参照 2026-09-05 已到位,升格尚未发生。E 组继续存在会让
   "自锚定"长期化。

---

## 5. 数据齐全性与一致性缺陷

### 5.1 数据齐全性(结论:基本齐全)

- `cases/reference_data/` 8 个数据集全部在仓、各带 README 与生成脚本;再生成需要
  gitignored 的 `tools/cfl3d_seq` / `tools/xfoil` 二进制(README 给了配方),**使用**不需要。
- 新克隆可用性:全部 reference_data 与 cases/gates 证据在仓;唯一需要生成的是
  ONERA M6 网格(coarse+medium ~30 s)。设计上缺失即 skip —— **但有例外,见 P0**。
- `cases/gates/`:20 门目录;D03 是唯一**没有 PNG** 的设门目录(可恢复出图,未恢复);
  INDEX.md 未覆盖 C04/D04(D04 的证据在 `cases/demo/p5_onera_m6/results/`)。
- C04/D04/D09/D10 无证据 CSV。C04/D04 是 conftest 明文承认的;**D09/D10 是
  "断言即证据"的能力边界探针** —— 形状可辩护,但无 CSV 意味着失败形状没有
  图/历史可读,且 F07 的门↔证据双向一致检查因它们不含 `gate_evidence_dir` 而豁免,
  属孤儿风险。建议补最小 summary.csv(把断言里的量写盘)。
- 悬空 gitignore 条目(onera_m6_wakefree、naca0012_wakefree_2.5d、cases/analysis/b23)
  对应目录已不存在,可清理。

### 5.2 P0(实测,干净 clone 会硬红/仪器已坏)

1. ★★ **`tests/D/test_D10_rans_onera_m6.py:82` 直接 `read_mesh(cases/meshes/onera_m6/coarse.msh)`,
   全文件无 exists/skip 守卫**(亲验)。M6 网格被 gitignore ⇒ **干净 clone 上 D10 硬
   FileNotFoundError,不是 skip**。与 A16/A17/A43/B05/D04/E03 的 skip 约定直接冲突。
2. ★ **`tests/D/test_D07_euler_onera_m6.py` 的 m050 腿同样无守卫**(`_one_m050`,
   :409-410,亲验)。主腿经 `_m6_case`(:127)间接被 E01:159 的 skip 保护,但 m050 腿
   不经过它 ⇒ 干净 clone 上该腿硬红。
3. ★ **`bench/run_g82_anchor_check.py:36` 已坏**:`from test_p8_newton import _m6_case`,
   而该文件在 2026-08-24 重编号中已拆分为 `tests/A/test_A52…` + `tests/E/test_E01…`
   (亲验:裸 import,`ModuleNotFoundError`)。**仪器盲区**:F07 的 import 检查只认
   `tests.*` 前缀,抓不到裸模块名 import;同类问题 run_m3_budget 修过(并催生了 F07),
   此脚本漏修。修法:`from tests.E.test_E01_p8_newton_anchors import _m6_case`。

### 5.3 P1 级陈旧/断链(不红但误导)

4. `bench/gate_results/capability_locks.csv` 内容仍是重编号**前**的节点名(落后一轮);
   `capability_matrix.csv` 是 19 列旧 schema 且含已删 level-set 行 —— 现脚本真跑会被
   自己的 schema 守卫拒写(设计如此,但该矩阵在 HEAD 不可复现/不可续跑);
   `capability/` 是空目录而 docstring 声称其产物 TRACKED。
5. `bench/gate_results/ibl_vs_xfoil_recon.csv` **无生产者**(study 的 run.py 不写此 CSV),
   证据链断;`bench/studies/ibl_xfoil_recon/` 无 PRE_REGISTRATION/VERDICT,不符合
   studies 纪律;`bench/studies/README.md` 的"被谁加载"表是重编号前测试名,且
   v6_1_wake_sheet 在 HEAD 已无测试加载其 run.py。
6. `tests/D/test_D05/D06` 模块 docstring 陈旧(写"读 L3",代码已读 L4;D05 的错误消息
   仍写 "no L3 rows")—— 2026-09-05 L4 刷新没扫干净的注释。
7. `tests/_tol.py:11` 承诺的 `@pytest.mark.identity` **零使用且未在 pyproject 注册** —
   docstring 描述了一个不存在的机制。
8. `bench/README.md`(2026-08-10)与 `README_ROLES.md`(2026-08-24)并存,前者脚本表
   已过时;`PROJECT_STRUCTURE.md` 的 bench 段与目录树整体过时(仍列已删的 LS 遗物);
   `CLAUDE.md` 三处 bench 引用失效。**F07 的路径检查只扫 .py 不扫 .md**,这族文档漂移
   不会被抓 —— 建议把 .md 纳入或定期人工扫。
9. `tests/B/test_B04:112` cell-field roundtrip 只断言 `len>0`,与文件声称的
   "verify all data identical" 不符。
10. `m3_budget_head.csv`(现默认产物)无 committed 对应物;`m3_budget_head_medium*.csv`
    经 2026-08-26 平滑勘误后在 HEAD 不可复现(已登记,知情即可)。
11. F02 缺号(F01→F03),轻微。

---

## 6. 对 debug 指导性评估 —— gate 能不能当 debug 的地图用

**能,而且已经是半个地图**:failure_modes 五模式分类器、每门 docstring 里的根因叙事、
证据锁的 refresh_hint,都是 debug 友好设计。**缺的是一张总索引** —— 症状 → 先跑哪几门 →
各门隔离什么。建议补一份(可直接放进 `docs/` 或 README),草案:

| 改动/症状 | 第一梯队(分钟级) | 第二梯队 | 各门隔离什么 |
|---|---|---|---|
| 改 kernels(residual/jacobian/gradient) | A01 → A02 → A21 | C01 → C03 | 自由流 / 装配位相同 / Jacobian FD / 阶 / 壁面精度 |
| 改迎风/熵修正 | A05 → A06 → A31 | C08 → D03 → D05 | no-op 位相同 / 熵修正 31 锁 / 伪根 / 剂量响应 / 激波带 / 对 Euler |
| 改尾迹/Kutta | A24 → A25/A26 → A54 | C06 → C07 → D01 | 切割 / 估计器 / 对称 / 规定 Γ / 自选 Γ / 两路互检 |
| 改远场/边界 | A01 → A54 | C09 → D05(α=0 绝对带) | 自由流 / 对称 / PG / 零升 |
| 改网格生成 | A13–A19 | A08/A09 → B06 | 生成器几何 / 邻接着色 / 截面 cl 交叉 |
| 改 IBL 闭包/推进 | A34/A37/A38/A39 → A36 | D13 | 闭包导数 / FS-oracle / Blasius / XFOIL |
| 改耦合(松/紧) | A40/A41 → A42 → A44–A47 | D08 | 通道守恒 / 松环 / 紧系统 FD / 对 RANS 带 |
| 改后处理 | B01–B04 → B03 | B05/B06 | 恢复精度 / 度量 / 锯齿压制 / 截面 cl |
| 收敛失败(任何) | `classify_failure`(bench/failure_modes.py) | A30/A31 → run_bench --compare | 五模式归因 / 种子回退 / 漂移表 |
| 改求解器驱动 | A52 → A32/A53 | E01(gated) | 驱动性质 FD / 默认值门 / 锚点 |

另两条 debug 基建建议:
① gs40/task3 系列(29+10 个 CSV)里的判别器(种子非唯一性、近奇分区、捕获选择等)是
debug 的宝库,但生产脚本已归档 —— 建议只把**判别器函数**提为活工具或写索引,不必复活脚本;
② 每门 FAIL 时的"下一步去哪查"指针,目前散在各 docstring,质量参差,建议在
`tests/_gate_evidence.py` 的报错消息里统一带 `refresh_hint` 之外的 `debug_hint`。

---

## 7. 建议优先级汇总

**P0(一致性/正确性,先修,都便宜)**
1. D10 与 D07-m050 补网格缺失 skip(§5.2-1/2,亲验)。
2. 修 `run_g82_anchor_check.py` import;**同时把 F07 扩到裸模块名 import**(§5.2-3)。
3. 修 D05/D06 的 L3→L4 陈旧 docstring 与错误消息(§5.3-6)。

**P1(覆盖补强,按价值/成本排序)**
4. A01 扩任意方向自由流(3 个旋转角 ×  lifting/nonlifting,§4.1-G3)。
5. 全局质量守恒门(§4.1-G2,先立在 C05 喷管与一个翼型 case 上)。
6. C09 的 PG 2% 断言移出 skipif(§4.4-3,零成本)。
7. ZPG 湍流平板 strip2d 端到端锚(§4.2-G5,用有出处的关联式)。
8. F08 的 implied_order 估计器复用到自身阶梯(GCI 纪律,§4.1-G4)。
9. E01/E03 升格评估(D 类参照已到位,§4.4-5);D09/D10 补最小 summary.csv(§5.1)。
10. ILU 路径补冒烟或声明弃用;清理孤儿(make_inboard_clip、wake-free 族、悬空 gitignore、
    capability/ 空目录)(§4.3、§5.3-4)。

**P2(可选/工程量大/需求驱动)**
11. 非线性全速势算子 MMS(§4.1-G1)。
12. 椭圆/椭球算例(分离壁面几何误差与 TE 奇点;20260824-0400 已列,仍有效,低优先)。
13. design.md 能力表与实现对齐(侧滑/C_m/Trefftz/非升力 Newton)(§4.3)。
14. 紧耦合三条实测失败形状立"断言即证据"门(§4.3)。
15. `.md` 路径引用扫描入 F07 或定期人工扫;`@pytest.mark.identity` 要么实现要么删
    docstring(§5.3-7/8)。

---

## 8. 一句话结论

**gate 体系本身健康、分层正确、证据纪律强于业界常见实践;数据齐全;真正要紧的是
两件事 —— 修掉干净 clone 会红的 D07-m050/D10 网格守卫(P0),以及补上任意方向自由流、
全局质量守恒、湍流平板三条便宜的外部锚(P1)。** 其余缺口都有记录在案的替代物或明确的
能力边界裁决,不构成 debug 盲区。
