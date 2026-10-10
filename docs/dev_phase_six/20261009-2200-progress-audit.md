# 20261009-2200 — **进展审计:套件是红的,H30 的 BLAS 固定在 pytest 里不生效**

轮次文件名:`20261009-2200-progress-audit.md`
阶段 / 门:无门(工程 / 调查)—— 使用者指示「审计本仓库的进展,看现在需要解决什么问题,
还有哪些东西没有测试和修改」
提交:分支 `claude/phase-six-audit`,HEAD `f931cf5`(+ 工作区未提交的 A52 ilu 测试)

**性质**:只读审计 + 实测。**本轮没有改 `pyfp3d/`、没有改任何判据、没有改任何测试**。
上一份计划书是 [合并分析 20260906-0100](20260906-0100-gate-audit-synthesis.md)。本文件
对着它的台账逐项核对完成情况,并补上 2026-09-06 之后无人测过的现状。

---

## 0. 四句话

1. ★★★ **套件现在是红的**:未门控全套 **795 passed / 6 failed / 32 skipped / 2 xfailed**
   (1705 s,numba 16 线程,负载约 8)。6 条全是证据锁(D05 D06 D08 D11 D12 D13)。
   上一轮收口(`76b70c0`)记的是「tests/D 100 passed / 0 failed」,**这个状态现在复现不出来**。
2. ★★★ **H30 的层 2(把 BLAS 固定成单线程)在 pytest 进程里不生效** —— 已实测(§4.2)。
   F09 看不出来,因为它只测了固定**会生效**的那种 import 顺序,而它的前提断言读的是
   固定代码自己写进去的 `os.environ`,**按构造恒真**。
3. ★★ **M1 门仍跑 `n_newton_max = 80`**,也就是 H30 认定不够的那个预算,而那条腿就是 M1 的靶子。
   已提交 CSV 里有 4 条 medium 腿正好停在 80 步、0 钳制(§4.4)。**这是线索,未重跑验证。**
4. ★★ **`docs/agent-rules.md` 的 "Current phase" 停在 GV5.6(2026-07-25)**,CLAUDE.md 通过 `@`
   引入它,所以每个会话开工时都会被告知错误的当前阶段。

---

## 1. 目的

回答三个问题:**(a)** 合并分析 §5 的 W0–W3 计划做到哪了;**(b)** 现在的真实状态
(实测,不抄文档)是什么、要先解决什么;**(c)** 哪些代码没有任何测试够得到。

**事前写下的失败标准**:本轮只要有一个结论建立在**没有实测**的读数上,就算失败。
所以每个判断都带上它的手段;没测的(§4.4 的 M1)明写「线索,未验证」。

## 2. 改了什么

- 新增本文件,以及 `progress.md` 的一行。
- **未改** `pyfp3d/`、`tests/`、`bench/`、`cases/`。
- ★ 全套跑完后工作树多了 3 个改动的已跟踪文件(D01/D02 的历史子目录证据,见 §4.8),
  **不是有意的修改**,已按 CLAUDE.md 用 `git checkout HEAD -- <paths>` 还原。
- 工作区原有的两样东西**没动**:`tests/A/test_A52_newton_driver_properties.py` 的未提交改动、
  未跟踪的 `bench/studies/gate_audit_20260905/results/layer_attribution_D.csv`(处置见 §6)。

## 3. 怎么测的

```bash
conda activate up3d; export PYTHONNOUSERSITE=1
# 全套(开跑时负载 0.37,跑的过程中约 8.6;24 核)
NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python -m pytest tests/ -q -rfEsx -p no:cacheprovider
# BLAS 判别式:同四条证据锁,只改外部线程设置(每次约 4 min)
NUMBA_NUM_THREADS={16,8} OMP_NUM_THREADS={1,8} OPENBLAS_NUM_THREADS={1,8} \
  python -m pytest -q "tests/D/test_D0{5,6}_*.py::TestCommittedEvidenceIsLoadBearing" \
                      "tests/D/test_D1{1,2}_*.py::TestCommittedEvidenceIsLoadBearing"
```

- import 顺序探针:用 ctypes 直接调 `libopenblasp-r0.3.33.so` 的 `openblas_get_num_threads()`,
  读的是**库实际在用的线程数**,不是环境变量。
- 静态可达性:AST 列出 `pyfp3d/` 的全部顶层定义,以 `tests/` 里出现过的标识符
  加上各模块的模块级代码为种子,按「定义体里提到的名字」做传递闭包。**这个口径偏宽**
  (提到即算可达),所以「不可达」是可靠的,「可达」不等于「被执行过」。
- ★ 参考值全部是**防漂移的已提交证据**(`cases/gates/*/summary.csv`),不是外部真值。
- ★ **本轮的日志与探针输出没有落成已提交的 CSV**(纪律 3)。数字全部抄自实跑输出,
  命令如上可复算;若要把它们当证据引用,先补一份 CSV。

## 4. 结果

### 4.1 全套现状

| 量 | 文档写的 | 实测 2026-10-09 |
|---|---|---|
| 收集到的测试数 | — | **835** |
| 未门控全套 | 674 / 33 / 2(CLAUDE.md,08-28)· 479 / 12 / 2(README,08-16) | **795 passed / 6 failed / 32 skipped / 2 xfailed** |
| 墙钟 | — | 1705 s,numba 16 线程,负载约 8 |

失败的 6 条(全是 `TestCommittedEvidenceIsLoadBearing`;锁在**第一个**不符处就抛,
所以下表每门只看得到一列,其余列是否也动了**不知道**):

| 门 | 第一个不符的列 | 相对差 |
|---|---|---|
| D05 | `M0.80_a1.25\|medium` cp_rms_upper | 3.34e-05 |
| D06 | `M0.725_a2.55\|medium` cp_rms_upper | 7.02e-05 |
| D11 | `M0.803_a-0.10\|medium` cl_p | 2.05e-06 |
| D12 | `Case7_M0.725_a2.55\|medium` mach_max | 3.37e-06 |
| D08 | `xtr005\|medium` cp_rms_upper_sst | 6.54e-03 |
| D13 | `medium/upper/x=0.05` dstar_ibl | 1.43e-04 |

★ 未提交的 A52 ilu 测试(W3.7 / H20)在这次全套里**通过**。

### 4.2 ★★★ BLAS 固定在 pytest 里不生效

`pyfp3d/__init__.py:53` 在 import 时执行 `os.environ["OPENBLAS_NUM_THREADS"] = "1"`。
但 OpenBLAS **只在 numpy 首次加载时读一次**这个变量。外部设 `OPENBLAS_NUM_THREADS=16` 实测:

| import 顺序 | `os.environ` | 库实际线程数 |
|---|---|---|
| `import pyfp3d, numpy` | `1` | **1** |
| `import numpy, pyfp3d` | `1` | **16** |
| `import numpy, scipy.sparse.linalg, pyfp3d` | `1` | **16**(openblas + scipy 的两个 BLAS 包装) |

- **测试走的正是第二种顺序**:D05/D06/D08 第 58/43/77 行都是先 `import numpy as np` 再
  `from pyfp3d...`,A02/A03 也一样,而 pytest 第一个加载的就是 A 类。⇒ **整个 pytest 进程的
  BLAS 线程数由外部 shell 决定**;外部不设时由 OpenBLAS 自己的默认值决定。
- **F09 为什么看不见**:它的子进程脚本第 55 行写着
  `import pyfp3d  # <- 必须最先,它负责固定 BLAS`,也就是只测了固定能生效的那种顺序。
  它的前提断言(第 73/115 行)读 `os.environ.get("OPENBLAS_NUM_THREADS")`,而这个值正是
  固定代码自己刚写进去的。⇒ 合并分析 §2.3 批评的「**值等于它自己的来源**」,又出现了一次。
- **同一族错误的前身**:CLAUDE.md 记过「`git status` 不是验证」和「一个前提与它所声称的结论
  不符的检查」。这次是「**环境变量不是库的状态**」。

### 4.3 ★★ 已提交的证据对不上任何一个 BLAS 档位

同样四条锁(D05 D06 D11 D12),只改外部线程设置:

| 门 · 列 | 已提交值 | 外部 BLAS 16 | 外部 BLAS 1 | 外部 BLAS 8 |
|---|---|---|---|---|
| D05 cp_rms_upper | 1.573486174e-01 | 1.573538717e-01(3.34e-05) | 1.573535037e-01(3.11e-05) | 1.573458060e-01(1.79e-05) |
| D06 cp_rms_upper | 1.437568286e-01 | 1.437669133e-01(7.02e-05) | 1.432788297e-01(**3.33e-03**) | 1.436960928e-01(4.22e-04) |
| D11 cl_p | −4.231428646e-02 | −4.231437324e-02(2.05e-06) | −4.231436051e-02(1.75e-06) | **绿** |
| D12 mach_max | 1.340916705e+00 | 1.340921229e+00(3.37e-06) | 1.340709302e+00(1.55e-04) | 1.340889437e+00(2.03e-05) |
| 结果 | — | 4 红 | 4 红 | 3 红 / 1 绿 |

> **勘误标注(2026-10-10,第 34 轮)**:上表「已提交值」列已于当日重刷到真 BLAS=1
> (运行时固定生效后),且证据 CSV 新增 `blas_threads` 出处列;本表是重刷前的历史
> 快照,数值勿再引用。见
> [20261010-1345-evidence-relock-blas1-verdict.md](20261010-1345-evidence-relock-blas1-verdict.md)。
> 另:本审计 §4.3 的「未决:当时的档位」已由同日独立审计实测闭合 = **默认 24 档**
> (外部不设时四条锁 4/4 复现)。

- **numba 线程数不是变量**:外部 BLAS 1 时,numba 8 和 16 给出**逐位相同**的四个数。
  H30 的 F1 在这四条腿上也成立。
- **推断**(只验证了一半):2026-09-06 刷证据(`76b70c0`)时固定同样没生效,于是证据锁在了一个
  **没有记下来**的 BLAS 档位上;而那次「tests/D 100 passed」是在同一个环境里自我复现的。
  支持这个推断的有两点:BLAS=1 对不上(若固定当时生效,应该对上);{1, 8, 16} 三档没有
  一档能把四条一起对上。没验证的是:当时的档位**究竟是多少**。
- **证据 CSV 的 `n_threads` 列记的是 numba 线程数**(D05/D06 已提交行均为 `8`),
  这个旋钮已证明不影响结果;真正影响结果的 BLAS 线程数**没有列**。⇒ 出处列量错了旋钮。
- **D06 的 BLAS 敏感度比 D05 大两个数量级**(BLAS 1 对 8 差 3.3e-03,H30 在 D05 上测的是
  1.8e-05)。这条腿可能有自己的问题(预算、或近非唯一),**未查**。
- ★ 未跟踪的 `layer_attribution_D.csv` 里 D06 cp_rms_upper 那一行是 `4.22e-04`,
  与本轮 BLAS 8 的读数**逐位一致** ⇒ 上一轮已经测到过这件事,只是没提交、没写进文档。
- **D08/D13 是另一回事**:IBL 松耦合本来就不能逐次复现(纪律 12),`76b70c0` 已登记为
  开放议题,BLAS 固定治不了。

### 4.4 ★★ M1 门的 80 步预算(线索,未验证)

> **勘误标注(2026-10-10,第 35 轮)**:本条线索已重跑闭合 —— 预算假设**被否**:
> 80 步停在 limit_cycle 的 medium seed-0 C=3.0 腿在 400 步**仍 limit_cycle**
> (|R| 4.7e-07,cl 不动)⇒ M1c seed-0 的覆盖失败**不是预算缺陷**。
> 见 [20261010-1500-m1-gate-budget-verdict.md](20261010-1500-m1-gate-budget-verdict.md)。

- `bench/run_m1_gate.py:187` 仍是 `n_newton_max=80`。H30 层 1 按预注册 §2
  **有意不改 `bench/`**,理由是「A 类机制测试依赖预算耗尽的行为」。但 M1 门是**产品指标门**,
  断言的是答案,按那份预注册自己的划分标准应属于「答案门」。
- 已提交的 `m1_gate_default.csv` / `m1_gate_default_3seed.csv` 里有 **4 条 medium 腿
  `n_newton = 80`、`n_limited = n_floored = 0`、`accept_reason` 为空**,正是预算耗尽的特征
  (另有若干条是真钳制,那些与预算无关)。
- ⇒ M1c 在 seed 0 上那个「覆盖失败」(medium 只有一条收敛腿),**可能有一部分是同一个预算缺陷**。
  要确认只能重跑 M1 门;**本轮没跑**。

### 4.5 合并分析计划的完成台账

| 项 | 状态 | 核对手段 |
|---|---|---|
| W0.1–W0.7 | ✅(附录 A.1) | — |
| W0.5 的两项延期 | ◐:`capability_locks.csv` 已重生成(5 组,新路径);**README 基线仍是 479(08-16)**;`bench/gate_results/capability/` 空目录仍在 | `head` / `ls -A` |
| W1.1–W1.5 | ✅(附录 A.2) | — |
| W1 验收的一项 | ✗ CLAUDE.md 第 5 行仍写 "M∞ 0.3–0.87",而 D14 已实测 M ≥ 0.82 不唯一 | `grep` |
| H30 层 1 / 层 3 | ✅(`1f4a49b`、`76b70c0`) | — |
| **H30 层 2** | **✗ 在 pytest 里不生效(§4.2)**,且由它刷出的证据要重刷(§4.3) | 实测 |
| **W2.1** 我方误差棒 | ✗ 未开始:`tests/D` 里没有一个文件用到 `xcoarse`、观察阶或 GCI | `grep` |
| W2.2 G10.1 | ✗ `newton.py` 只有 `solve_newton_lifting` / `solve_newton_transonic`;C05/C08 仍走 `tests/_nozzle_case.py::DuctSystem` | `grep` |
| W2.3 非唯一性探针 | ◐ 由 D14 以能力边界形状承担 | — |
| W2.4 | ✅(`c991e76`) | — |
| W2.5 网格实现噪声 | ✗ C06 的 `CL_CONTRACTION = 2.0` 未动(余量落在噪声里,登记未修) | 读 :69 |
| W2.6 hex 网格族求解 | ✗ `naca0012_hex_2.5d` / `cylinder_hex_2.5d` 只有 A18 引用,A18 里 0 次求解 | `grep -c solve_` |
| W3.1–W3.6 | ✗ 未开始;W3.2 现状不变(`design.md` §12 没有一条风险指向守它的门) | — |
| W3.7 | ◐ H20(未提交,测试已通过)、H23 ✅;**H21** C03 仍从 `cases/demo/p11_curved_walls/` 读证据;**H22** `cases/gates/` 下没有 D09/D10 目录、D03 没有 PNG、`INDEX.md` 里 C04/D04 出现 0 次;**H27** 门控层没有新鲜度机制 | `grep` / `ls` |
| W3.8 覆盖率 | ✗ 未跑。★ 但 **`coverage` 7.16.0 已经装在 `up3d` 里** ⇒ §7 第 2 条「允许装 coverage 吗」**不再需要裁决** | `import coverage` |

### 4.6 没有任何测试够得到的代码

339 个顶层定义里有 **13 个、共 626 行**(全库 18 007 行)从 `tests/` 出发完全不可达。

| 定义 | 行数 | 其他调用方 | 读法 |
|---|---|---|---|
| `meshgen/structured.py::airfoil_hybrid_2d` | 193 | `cases/meshes/naca0012_hex_2.5d/generate_naca0012_hex.py` | 只被生成脚本用,无测试 |
| `meshgen/structured.py::cylinder_o_grid_2d` | 56 | `cases/meshes/cylinder_hex_2.5d/generate_cylinder_hex.py` | 同上 |
| `meshgen/planar.py::cylinder_annulus_2d` | 54 | `cases/meshes/cylinder_2.5d/generate_cylinder.py` | 同上 |
| `meshgen/planar.py::naca0012_wake_2d` | 27 | `cases/meshes/naca0012_2.5d/generate_naca0012.py` | 同上 |
| `meshgen/structured.py::airfoil_o_grid_2d` | 82 | 无 | 疑似死代码 |
| `meshgen/structured.py::_outward_normals_closed` | 17 | 无 | 同上 |
| `post/vtk_out.py::write_nodal_field_heatmap` | 26 | 无 | 同上 |
| `post/vtk_out.py::export_error_heatmap` | 52 | 无 | 同上 |
| `mesh/metrics.py::precompute_face_normals` | 43 | 无 | 同上 |
| `physics/isentropic.py::speed_of_sound_squared` | 15 | 无 | 同上 |
| `post/surface.py::_node_to_element_patches` | 14 | 无 | 同上 |
| `solve/timing.py::records_to_columns` | 8 | 无 | 同上 |
| `meshgen/fuselage.py::make_inboard_clip` | 39 | `bench/capability_matrix.py:53` **只 import 不调用** | **确认死代码**:它服务的 `CutElementMap` 已随水平集路线在 phase 3 删除(D5),全仓只剩它 docstring 里的一处提及 |

★ 这是函数级的「完全没碰到」。函数里**哪些分支**没跑到,要靠 W3.8 的覆盖率,
并且必须配 `PYFP3D_NOJIT=1`,否则 `@njit` 函数体整体读成未覆盖。

### 4.7 文档欠账(五面收口)

| 面 | 现状 |
|---|---|
| `docs/agent-rules.md` "Current phase" | **停在 GV5.6(2026-07-25)**,即 phase 1 的 Track V。CLAUDE.md 通过 `@` 引入,**每个会话开工都读到它** |
| `progress.md` | 最后一行是第 32 轮(2026-09-05);之后约 10 个提交(两份审计与合并、W0/W1、H30 全链、W2.4、W3.7)**没有轮次行**;概览表还写着「轮次 0」 |
| 测试基线 | README 停在 479(08-16),CLAUDE.md 和 PROJECT_STRUCTURE 停在 674(08-28),实测 795 / 6 / 32 / 2 |
| 能力边界 | 第六阶段自己的能力边界文书**还没写**,而 CLAUDE.md 把它列为第六阶段的交付物 |

### 4.8 ★ 默认跑会改写已跟踪的证据文件

全套跑完后,`git status` 显示 3 个已跟踪文件被改写:

| 文件 | 写入位置 | 改动 |
|---|---|---|
| `cases/gates/D01_naca0012_incompressible_panel/G2.3/summary.csv` | `test_D01…py:101` | medium 行末位(如 4.785760988e-01 → …987e-01) |
| `cases/gates/D02_naca0012_m05_panel/G3.2/summary.csv` | `test_D02…py:105` | 末位 |
| `cases/gates/D02_naca0012_m05_panel/G3.2/v3_2_cp_curves.csv` | `test_D02…py:174` | 8 行末位 |

- 这两处写盘**没有** `gate_figures_enabled()` 守卫,而同一文件里的出图腿有,其 skip 理由写着
  「平时不写,避免每次 pytest 都脏工作树(2026-08-24)」。⇒ 这条约定在 D01/D02 的历史子目录上**没有兑现**。
- 后果不只是工作树变脏:谁顺手 `git add -A`,就会把一次随便哪个环境下跑出来的末位漂移
  当成证据提交。这正是 CLAUDE.md 记过的「`git add -A` 之前先核对」那类风险。
- 本轮已用 `git checkout HEAD -- <paths>` 还原,**没有提交任何改写**。

## 5. 判定

**RECORDED**(只读审计,不对任何门作判定)。

一句话:**套件是红的(795 / 6 failed),主因是 H30 的 BLAS 固定在 pytest 里从未生效、
用它刷出的 D 类证据锁在了一个没记下来的 BLAS 档位上;另外 M1 门仍跑那个已知不够的
80 步预算,`agent-rules.md` 还停在 phase 1 的当前阶段。**

## 6. 下一步

**建议顺序**(每一项都是独立的一轮,带预注册):

1. **让 BLAS 固定真正生效**:在 `tests/conftest.py` 最顶部(任何 numpy import 之前)固定,
   和/或在 `pyfp3d/__init__.py` 里于运行时调用 `openblas_set_num_threads(1)`(对已加载的库也有效)。
   ★ **F09 要同时改**:前提断言读**库的实际线程数**,子进程脚本**先** `import numpy`。
   G-TEETH:把固定撤掉,F09 必须红。
2. **在真正的 BLAS = 1 下重刷 D05/D06/D11/D12 的证据**,走纪律 11 的勘误清单;
   证据 CSV 加一列**实际 BLAS 线程数**(由库读出,不由环境变量读出)。
   ★ 先单独查 D06 那条腿为什么对 BLAS 敏感两个数量级。
3. **重跑 M1 门**(`run_m1_gate.py` 的 80 → 400),看 M1c 在 seed 0 上的覆盖失败还剩多少。
   会移动产品指标 ⇒ 记进 `progress.md` 的产品指标表。
4. **文档**:`agent-rules.md` 的 "Current phase" 改到第六阶段;三处基线改成**实测值**;
   给 2026-09-05 之后的轮次补 `progress.md` 行;CLAUDE.md 第 5 行的 M∞ 上界改成 D14 实测的可信值。
5. **D01/D02 的证据写盘补上 `gate_figures_enabled()` 守卫**(§4.8)。G-TEETH:不设标志跑一遍,
   `git status` 必须干净。顺手扫一遍 `tests/` 里其他不带守卫的 `open(..., "w")`。
6. 然后才是合并分析排第一的 **W2.1**(我方的误差棒)。

**工作区里悬着的两样东西**:
- A52 的 ilu 测试已通过,可以提交(W3.7 / H20)。
- `layer_attribution_D.csv`:内容与本轮实测吻合(§4.3),实际上已经回答了上一轮挂着的
  「只有层 1 的对照」那个问题。但它**没有生成脚本,也没有读者**。要么补生成脚本、写进文档再提交,
  要么删掉 —— 和 §2.3 处理 `ibl_vs_xfoil_recon.csv` 是同一个二选一。

**登记、不在本轮修**:
- D08/D13 的证据锁在任何扰动下都会红(IBL 松耦合),证据锁该改成什么形状,需要使用者裁决。
- `make_inboard_clip` 及 §4.6 里「无调用方」那 8 个定义:删,还是补测试。
- `bench/gate_results/capability/` 空目录(纯本地噪声)。

**需要使用者裁决**(合并分析 §7 还开着的三件,外加一件新的):
1. H11 的处置:(a) 缩包线、(b) 改激波处理,还是 (c) 两者都做。
2. W3.5:要不要为翼身构型补一批 CFL3D 数据。
3. W2.2:是否按「让精确解经过生产驱动器」这个新理由提前做 G10.1。
4. D08/D13 证据锁的形状(上一条「登记」里那件)。
