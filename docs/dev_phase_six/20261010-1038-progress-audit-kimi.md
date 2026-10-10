# 20261010-1038 — **独立进展审计:复测 Claude 审计全部成立,补上它留的洞(BLAS 档位 = 默认 24),另抓到两个同族新缺陷**

轮次文件名:`20261010-1038-progress-audit-kimi.md`
阶段 / 门:无门(工程 / 调查)—— 使用者指示「独立开展全面审计,不要轻信 Claude 的
[20261009-2200 审计](20261009-2200-progress-audit.md)」
提交:分支 `claude/phase-six-audit`,HEAD `f931cf5`(+ 工作区未提交的 A52 ilu 测试与 progress.md 一行)

**性质**:只读审计 + 实测。本轮没有改 `pyfp3d/`、没有改任何判据、没有改任何测试。
方法 = **独立复测 Claude 的每一条主结论**(不同探针、不同入口)+ **扫它没扫的面**
(bench/ 健康度、判据恒真性、tests 之外的可达性、依赖声明、git 拓扑)。

---

## 0. 四句话

1. ★★★ **Claude 审计的主结论我逐条独立复现,全部成立**:未门控全套 **795 passed /
   6 failed / 32 skipped / 2 xfailed**(我实测 2062.73 s,它 1705 s;6 条失败门与
   首不符列**逐位一致**);BLAS 固定在 numpy 先导入时失效(我的探针直接读库);M1 门
   `n_newton_max=80` 仍在 `bench/run_m1_gate.py:187`。
2. ★★★ **它留下的最大开放问题我测出了答案:已提交证据锁在 BLAS = 24(机器默认档)**。
   它测了 {1, 8, 16} 三档全对不上,推断"锁在了一个没记下来的档位上";我测了它没测的
   第四种 —— 外部**不设** `OPENBLAS_NUM_THREADS`(本机 24 核 ⇒ OpenBLAS 默认 24):
   **D05/D06/D11/D12 四条锁 4/4 全绿**(437.7 s)。这同时解释了 76b70c0 当时为什么
   "六条判据全过"、以及为什么按纪律设 16 线程跑现在是红的。
3. ★★ **它漏了两个同族新缺陷**(都是这个项目反复登记过的族):
   (a) `bench/run_m3_budget.py:261` 读 `r.get("sigma_converged", r.get("entropy_converged"))`
   —— **两个键都不存在**(真键 `sigma_transport_converged`,`newton.py:1658`)⇒ R2 守门员的
   transport 支路 `(sconv is None or bool(sconv))` **按构造恒真**,与 `m_last_converged`
   先例(GS4.0)同族;
   (b) `tests/F/test_F05_usability_anchor.py:108` 裸 `import usability`,**只靠
   `bench/failure_modes.py:77` 的 `sys.path.insert` 在 D11/D12 收集期的副作用才能解析**
   ⇒ **F05 单跑即红**(实测 `ModuleNotFoundError`),全套里靠收集顺序碰巧绿。
4. ★★ **它的不可达代码名单方向对、口径太宽**:13 项里 5 项其实是活的生产路径
   (cases/meshes 生成器与 capability_matrix 的调用方都在),不能按"疑似死代码"处理;
   同时**漏报了 12 个 tests 不可达定义**,包括全库最大的零测试代码块
   `meshgen/wingbody.py::onera_m6_wingbody_mesh`(443 行)—— tests/ 里它只出现在
   `test_A28` 的 **docstring** 里(A28:8),从不被调用;Claude 自己写明"提到即算可达"
   的口径,正是被 docstring 提及骗过的那一类。

## 1. 目的

回答:(a) Claude 的审计哪些经得起独立复测、哪些不全面;(b) 现在真正要先解决什么;
(c) 哪些代码没有任何测试够得到。

**事前写下的失败标准**:与 Claude 相同 —— 任何结论建立在没有实测的读数上即失败;
对 Claude 的结论,凡标"复现"的都由我自己跑出数字,凡我没跑的都明写"未复测"。

## 2. 改了什么

- 新增本文件。**未改** `pyfp3d/`、`tests/`、`bench/`、`cases/`、`docs/` 的任何既有内容。
- ★ 全套跑完后 D01/D02 的 3 个已跟踪证据 CSV 被测试改写(与 Claude §4.8 同一现象,
  独立复现),已按纪律用 `git checkout HEAD -- <paths>` 还原;工作区回到审计前状态
  (只剩 Claude 留下的两样:未提交的 A52 ilu 测试、未跟踪的 `layer_attribution_D.csv`)。
- 探针与日志全部在 `/tmp`(不入库):`/tmp/kimi_audit_suite_full.log`、
  `/tmp/kimi_audit_blas24.log`、`/tmp/up3d_audit/reach.py`(可达性脚本)。
  ★ 与 Claude 同一纪律问题:本轮日志未落成已提交 CSV,数字抄自实跑输出,命令在 §3 可复算。

## 3. 怎么测的

```bash
conda activate up3d; export PYTHONNOUSERSITE=1
# (1) 全套基线(开跑时 load 0.39,24 核)
NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python -m pytest tests/ -q -rfEsx -p no:cacheprovider
# (2) BLAS 探针:ctypes 直调 libopenblasp-r0.3.33.so 的 openblas_get_num_threads()
#     读库实际线程数,四种 import/环境组合
# (3) Claude 没测的档位:外部完全不设 BLAS/OMP(本机默认 = 24)
env -u OPENBLAS_NUM_THREADS -u OMP_NUM_THREADS NUMBA_NUM_THREADS=16 \
  python -m pytest -q tests/D/test_D0{5,6}_*.py::TestCommittedEvidenceIsLoadBearing \
                      tests/D/test_D1{1,2}_*.py::TestCommittedEvidenceIsLoadBearing
# (4) 修复方向验证:numpy 先 import 后,运行时 openblas_set_num_threads(1) 是否仍有效
# (5) F05 自包含性:pytest tests/F/test_F05_usability_anchor.py 单跑
# (6) 静态面:四个并行子审计(bench/ 健康度、tests/ 结构与判据质量、docs/cases 一致性、
#     pyfp3d/ 可达性),全部要求给 文件:行号 证据,关键条目我抽查复核
```

- 参考值全部是**防漂移的已提交证据**(`cases/gates/*/summary.csv`),不是外部真值。
- 全套跑的时候工作区已含未提交的 A52 ilu 测试 —— 与 Claude 全套的条件相同,可比。

## 4. 结果

### 4.1 Claude 审计的逐项复核

| Claude 结论 | 复核判定 | 我的手段 |
|---|---|---|
| §4.1 全套 795/6/32/2,6 条全是证据锁 | ✅ **逐位复现**:同样的 6 门、同样的首不符列、同样的相对差(D05 3.34e-05、D06 7.02e-05、D08 6.54e-03、D11 2.05e-06、D12 3.37e-06、D13 1.43e-04) | 自己跑全套 2062.73 s |
| §4.2 BLAS 固定在 pytest 里不生效 | ✅ 成立,且**多测出第四档**:外部不设时库实际跑 **24**(nproc),不只是 {1,8,16} | 自己的 ctypes 探针,4 种组合 |
| §4.3 已提交证据对不上 {1,8,16} | ✅ 成立但**不全面**:它没测默认档 ⇒ 见 §4.2 补充,**24 档 4/4 全对上** | 实跑 4 条锁,437.7 s |
| §4.4 M1 门仍跑 `n_newton_max=80` | ✅ `bench/run_m1_gate.py:187` 实测仍在;★ 另见 D05:73 / D11:62 注释自称「与 bench/run_m1_gate.py 逐字一致」—— **这两句注释现在已经是假的**(它们自己已改 400) | 读码 |
| §4.5 W 台账(W2.1/W2.2/W2.5/W2.6/W3.1/W3.2/W3.5/W3.8 等) | ✅ 抽查 10 项全部一致;**一处数字修正**:「约 10 个提交没有轮次行」实测是 **15 个**(2026-09-06 全部提交) | 独立 grep + git log |
| §4.6 不可达 13 定义 / 626 行 | ◐ 方向对、口径宽:13 项从 tests 单种子**全部成立**(行数 626 交叉吻合);但 5 项有活的 bench/cases 调用方;**漏报 12 个**(§4.4) | 独立 AST,不同种子与边口径 |
| §4.7 文档欠账 | ✅ 全部成立,另有 6 处补充(§4.6) | 逐面核对 |
| §4.8 D01/D02 证据被默认跑改写 | ✅ 独立复现(我的全套改写了同样 3 个文件);**它漏了 2 处同类**:D02:213 的 PNG `savefig` 也无守卫;D03:344 写 `G4.3/summary.csv`(未跟踪文件、仅门控跑时) | 自己跑 + 全面扫写盘点 |

### 4.2 ★★★ 证据锁定的档位 = 机器默认 24(本轮最重要的新测量)

| 外部 BLAS 设置 | D05 | D06 | D11 | D12 | 来源 |
|---|---|---|---|---|---|
| =1 | 红(3.11e-05) | 红(3.33e-03) | 红 | 红 | Claude §4.3(我未复测这三档) |
| =8 | 红 | 红 | 绿 | 红 | Claude §4.3(同上) |
| =16(项目纪律档) | 红(3.34e-05) | 红(7.02e-05) | 红 | 红 | **我复测,逐位吻合** |
| **不设(本机默认 24)** | **绿** | **绿** | **绿** | **绿** | **我实测,437.7 s,4/4** |

- **推断闭合**:2026-09-06 刷证据(`76b70c0`)时,跑法的外部环境**没有设**
  `OPENBLAS_NUM_THREADS`(pytest 里 numpy 先导入 ⇒ 固定无效 ⇒ 库跑默认 24)。
  这解释了"六条判据全过"(自洽环境自我复现)与"现在红"(纪律跑法设 16)。
- **numba 线程数仍不是自变量**:我这次用的是 numba=16,已提交 CSV 记的是 `n_threads=8`,
  四条全对上 ⇒ Claude 的 F1(numba 无关)在我的读数上再次成立。
- ★ **注意它与项目纪律的张力**:AGENTS.md/CLAUDE.md 的线程纪律是 16 含 BLAS,
  而证据实际锁在 24 ⇒ 按纪律跑永远红。**修法不是把证据留在 24**,而是先把固定
  做成真的(§6-1),再在真正的 BLAS=1 下重刷证据。
- ★ **对修复方案的一个新约束**(Claude 的建议 1 没有写到这层):运行时
  `openblas_set_num_threads(1)` 我已验证对已加载的库有效(16→1)—— 但它一旦进
  `pyfp3d/__init__.py`,**快层 `bench/run_capability_locks.py` 的子进程会被强制成
  BLAS=1**(它现在靠"进程启动前设 env=8"实现确定性,运行时覆盖优先于 env)。
  G8.2 的墙钟预算是在 8 线程标定的,位级答案也可能动 ⇒ 修固定时快层要么显式走
  逃生舱 `PYFP3D_ALLOW_BLAS_THREADS=1`,要么接受 BLAS=1 并重刷 `capability_locks.csv`。

### 4.3 ★★ 新缺陷:两个"读不存在的键/靠副作用存活"的同族实例

1. **`bench/run_m3_budget.py:261`**:`sconv = r.get("sigma_converged", r.get("entropy_converged"))`
   —— 全 `pyfp3d/` grep:`sigma_converged` 只是 workspace **属性**(newton.py:268/563/998),
   返回 dict 里从来没有这两个键;真键是 `sigma_transport_converged`(:1658)。
   ⇒ R2 validity gate 的 transport 支路**永远为真、永不触发**,CSV 的
   `transport_converged` 列恒为 None。**与 GS4.0 的 `m_last_converged` 事故同族**
   (删除路线时 bench/ 没人 grep)。修 = 一行;但它自 2026-08-24 归档以来就这样,
   期间 R2 的读数都在半瞎状态下产出。
2. **`tests/F/test_F05_usability_anchor.py:108` 裸 `import usability`**:
   实测单跑 **1 failed**(ModuleNotFoundError)。机制:仓库根在 sys.path ⇒
   `from bench.usability import ...`(F05:14)没问题;裸名需要 `bench/` 本身在
   sys.path —— 由 `bench/failure_modes.py:77` 的模块级 `sys.path.insert` 提供,
   而它只在 D11/D12 收集期被 import ⇒ **F05 的绿色依赖目录字母序收集顺序**
   (D < F)。F 层一道"仪器锁"自身不自包含,与 F07 想防的那族(函数内 import 骗过
   作用域检查)是邻居。

### 4.4 ★★ 不可达代码:修正后的名单

独立 AST(种子 = tests/ + bench/ + cases/ + 模块级;边 = 剔除 docstring 后的定义体提及):

- **真死代码(两个种子集都不可达,10 定义 / 301 行)**:
  `mesh/metrics.py::precompute_face_normals` · `meshgen/structured.py::airfoil_o_grid_2d`
  与其私有 helper `_outward_normals_closed` · `physics/isentropic.py::speed_of_sound_squared`
  与 **`upwind_factor`(Claude 漏报,23 行,kernels/upwind.py:471 只在 docstring 提它)**
  · `post/surface.py::_node_to_element_patches` · `post/vtk_out.py::write_nodal_field_heatmap`
  与 `export_error_heatmap` · `solve/timing.py::records_to_columns` ·
  **`viscous/coupling.py::make_picard_nonlifting_driver`(Claude 漏报;它还在 `__all__` 里,
  唯一调用方在归档的 phases/p1/)**。
  注:`airfoil_o_grid_2d` 在归档的 `phases/p3/bench/run_hex_q1_donor.py` 有调用方,
  删除前斟酌;其余 8 个全仓库(含归档外)零调用方。
- **tests 不可达但有活调用方(15 定义 / 1183 行,Claude 名单里 5 个属于这类)**:
  5 个网格生成器入口(`airfoil_hybrid_2d`、`cylinder_o_grid_2d`、`cylinder_annulus_2d`、
  `naca0012_wake_2d`、`make_inboard_clip` —— 它的调用方 `bench/capability_matrix.py:53`
  是活脚本)+ **`meshgen/wingbody.py` 整族 9 个(Claude 全部漏报)**:
  `onera_m6_wingbody_mesh`(443 行)、`te_polyline`、`_fillet_junction`、`_junction_curves`、
  `_sheet_piece_in_body`、`_surface_on_fuselage`、`_wing_edge_curves`、`_wingbody_asserts`、
  `_wingbody_wake_asserts` —— 只被两个网格生成脚本调用;**翼身网格生成器 = 全库最大的
  零测试执行代码块**(tests/ 只在 A28:8 的 docstring 里提到它,我复核属实)。
- 差集涉及的 9 个 bench/cases 脚本**全部实测可 import**,无死脚本(抽查验真)。

### 4.5 Claude 没扫的面:bench/、依赖声明、git 拓扑

**bench/ 健康度**(12 顶层 + 9 studies,README_ROLES 说 7 —— 过期):
- 死键 1 处(§4.3-1);配方漂移 1 处(run_m1_gate.py:187 的 80,外加 D05/D11 两句
  「逐字一致」假注释)。
- **无生产者的证据 2 件**:`bench/gate_results/ibl_vs_xfoil_recon.csv`(无生产者、
  无运行时读者,只有 D13:21 docstring 提及)与未跟踪的 `layer_attribution_D.csv`
  (gate_audit 目录 4 个脚本都不写它,列结构也对不上 ⇒ **仓库内无法重新生成它**,
  "补生成脚本"需要先重建方法,否则只能删)。
- **断链证据 1 件**:phases/p4 判定文档称 `bench/gate_results/gs40k_subcell.csv`
  「判定前已落盘」—— 全仓(含归档、含 git 历史)都找不到这个文件 ⇒ **大概率当年
  只落盘未提交**,那条 G-SUB 判定的证据链是断的(纪律 3 的又一实例)。
- 真孤儿 CSV 2 个(`task3_le_registration_stations.csv`、`task3_m6_triage_16t.csv`,
  全仓零提及);`bench/gate_results/capability/` 空目录仍由 `capability_matrix.py:69`
  import 副作用反复重建。
- `tools/`(整个被 gitignore):`cfl3d_mpi` **全仓零引用,死二进制**;`cfl3d_seq` 只被
  手工生成器引用(不在任何 cadence 上,属设计);`cfl3d_work/` 有 9 月的运行残留。
- `bench/README.md`(旧版)点名 4 个已不存在/已改名的脚本。

**依赖声明**:`pyproject.toml` 的 `scipy>=1.10`、`pyamg>=4.2` 是**不安全的下界** ——
代码要求 scipy≥1.12(`rtol=`,environment.yml 有 81 failed 的实测记载)、pyamg≥5.2
需要 scipy≥1.12 的 `eye_array`;pip 按 pyproject 解析会装出必坏的组合。
`environment.yml` 的 pandas 注释指向不存在的测试路径(现为 A19,且 tests/ 无人 import
pandas;真正用户是 bench 脚本)。★ 附带实测:A19 在无 `PYTHONNOUSERSITE=1` 时红 2 条
(pip-user pandas 2.3.3 遮蔽 conda pandas 3.0.5)—— CLAUDE.md 已把该变量定为
load-bearing,这是它又多盖住的一处。

**git 拓扑**:`main == origin/main @ 8404548`(phase-5 收口);当前分支
`claude/phase-six-audit` **领先 64 个提交未合入** —— phase 6 的全部工作(W0/W1/H30/
W2.4/W3.7 + 两份审计)都不在 main 上。按 AGENTS.md 纪律这需要一个 PR 决策;
分支越大,与 main 的语义距离越远(尽管 main 没动,暂无冲突)。

### 4.6 文档欠账(在 Claude §4.7 基础上的增量)

| 面 | 增量发现 |
|---|---|
| `docs/agent-rules.md` | Current phase 停在 GV5.6(2026-07-25)—— 同 Claude;**另:头部基线写 652、尾部 Baseline 段写 651,同一文件自相矛盾** |
| `tests/conftest.py` docstring | 说「14 道全部转换完毕」—— **实际 19 道**(C 类 8 + D 类 11),D05–D08/D14 是后来补的,docstring 没跟 |
| `bench/README_ROLES.md` / `bench/studies/README.md` | 「7 个 studies」—— 实际 **9 个**(ibl_xfoil_recon 与 gate_audit_20260905 未登记;后者还违反 README 自定的 study 契约:无 PRE_REGISTRATION、硬编码 `/home/lrz/codes/UP3D`、无 `__main__` 守卫) |
| `cases/demo/README.md` | 「9 of 33 仍被读取」—— 现存 **6** 个,测试真读的只有 **p11** 一个;`tests/A/test_A20_tip_taper.py:28` 指向已归档的 p13 目录(失效指引) |
| `cases/gates/INDEX.md` | **缺 D14 节**(目录、PNG、summary.csv 都在,INDEX 零提及);D09/D10 只有散文无证据落点(设计如此但 INDEX 没写清) |
| `docs/dev_phase_six/progress.md` | 无轮次行的提交是 **15 个**(Claude 写"约 10 个");概览表仍写「轮次 0」 |
| `pyfp3d/` 卫生 | pyflakes:6 个未用 import、7 个未用局部变量;`pyfp3d.NOJIT` 属性**无消费者**(各 kernel 自己重读环境变量);`post/__init__.py` docstring 提到不存在的 `artifacts.py` |
| F09 自身 | 同一文件内 :38 说「强制覆盖」、:108-109 注释说「setdefault 所以外部值被保留」—— **机制描述自相矛盾**(:114 的 `blas_seen=="1"` 证明实际走的是强制覆盖;blas_ext 能存活是因为那是另一个 pyfp3d 不碰的变量,不是 setdefault)。守固定的门带着过期的机制描述 |

### 4.7 门控全集的账

- 最后一次门控全套实测 = **2026-08-24(604/1/4/0,1:44:32 @8)** —— 距今 7 周,
  期间未门控套件从 595 长到 835。D04 的重型 gated 腿(45–75 min)同样 7 周未跑。
- 快层 `capability_locks.csv` 最近刷新 = 2026-09-06(`7a81a48`)—— **早于** H30 层 1/2
  落地的 `1f4a49b`。快层自己在子进程启动前设 env=8 ⇒ 不受固定失效影响,读数仍自洽;
  但 H30 三层之后它没再跑过,重跑一次(~10–15 min)才能把"快层仍 5/5"说成现状。

## 5. 判定

**RECORDED**(只读审计,不对任何门作判定)。

一句话:**Claude 的审计主结论全部经得起独立复测,但它把最大的那个问题留了一半 ——
证据锁在 BLAS 默认 24 档(我实测 4/4 复现);另外它漏了两个同族缺陷(run_m3_budget
死键、F05 不自包含)、漏报了最大的零测试代码块(翼身网格生成器族),而 phase 6 的
64 个提交至今未合入 main。**

## 6. 下一步

**建议顺序**(每项独立一轮、带预注册;与 Claude §6 的差异在 1b/3/4):

1a. **让 BLAS 固定真正生效**(同 Claude 建议 1):运行时 `openblas_set_num_threads(1)`
    (我已验证对已加载库有效)+ conftest 顶部固定;**F09 的前提断言改读库的实际线程数**
    (ctypes/threadpool 风格),并新增一条 **numpy-先导入顺序**的载体腿。
    G-TEETH:撤掉固定,F09 必须红。
1b. ★ **新增,Claude 没写**:修固定前先决定快层的处置 —— 运行时固定会让
    `run_capability_locks.py` 从 BLAS=8 静默滑到 BLAS=1(§4.2 末条),要么显式逃生舱、
    要么重刷 `capability_locks.csv`;这个决定写进预注册。
2. **在真正的 BLAS=1 下重刷 D05/D06/D11/D12 证据**(纪律 11 勘误清单),证据 CSV 加
   **库实际 BLAS 线程数**列(从库读,不从 env 读 —— 现在的 `n_threads` 列记的是
   已证明无关的 numba 旋钮)。先单查 D06 为什么 BLAS 敏感度大两个数量级。
   ★ 本轮已证明 24 档可复现 ⇒ 若选择暂不重刷,至少把"证据锁在 BLAS=24、纪律跑法会红"
   记进各门 docstring,避免每个新会话重新撞一次。
3. ★ **修 `run_m3_budget.py:261` 死键**(一行:真键 `sigma_transport_converged`),
   并 grep `bench/` 里所有 `.get(` 读返回键的地方对一遍 newton.py 的写入清单 ——
   同族缺陷不能只修看到的这一处。
4. ★ **修 F05:108 裸 import**(改 `from bench.usability import ...` 与同文件 :14 对齐),
   并给 F 层加一条"每个文件单跑能收集"的机制检查(F07 的自然延伸)。
5. M1 门 `80 → 400` 重跑(同 Claude 建议 3),顺带把 D05:73/D11:62 的假注释更正。
6. D01:101 / D02:105 / D02:174 / **D02:213(PNG)** 四处写盘补 `gate_figures_enabled()`
   守卫(Claude 漏了第四处);D03:344 的 G4.3 写盘一并纳入。G-TEETH:不设标志跑一遍,
   `git status` 必须干净。
7. 文档面:agent-rules Current phase 到第六阶段(顺手统一 652/651)、progress.md 补
   15 个轮次行、INDEX.md 补 D14、bench 两个 README 的 studies 数、pyproject.toml 的
   scipy/pyamg 下界改到代码真实要求。
8. 然后才轮到合并分析排第一的 **W2.1**(我方误差棒)。

**提交决策(留给使用者)**:
- `claude/phase-six-audit` 领先 main 64 个提交,含两份审计与 W0–W3.7 全部工作,
  需要一个 PR 决策;工作区里 A52 ilu 测试(已通过,可随提交)与
  `layer_attribution_D.csv`(仓库内无生产者 ⇒ **建议删**,除非有人能补出生产脚本)。
- 本审计文件与 Claude 的 `20261009-2200-progress-audit.md` 目前都未跟踪。

**登记、不在本轮修**:D08/D13 证据锁形状(任何扰动下都红,IBL 松耦合,需使用者裁决);
H11 处置、W3.5 翼身 CFL3D 数据、W2.2 是否提前(Claude 清单上的三件裁决仍然开着);
真死代码 10 定义的删/留;`onera_m6_wingbody_mesh` 族的测试策略(生成器无 cadence
是现状设计,但 443 行零执行值得一条登记)。
