# 20261104-0524 — **审计（第 10 轮）：UP3D 开发流程审计与 openspec-devflow 采用评估**

轮次文件名：`20261104-0524-devflow-adoption-audit.md`
阶段 / 门：phase 7 / 无门（工程/调查）—— 不是预注册研究轮，是采用决策的审计文书
审计对象：`080f291`（分支 `kimi/phase7-ibl-ap-projection`，已于本轮经 PR #34 合入 `main` = `9b4f425`）
证据：`bench/studies/devflow_adoption_audit/`（devflow 指纹 / JUnit 判定最小实验 + 12 红定向重跑原始日志）

## 0. 轮次定位（按模板）

- **目的**：判断能否、以及怎样在不改变求解器行为与科学证据的前提下采用 openspec-devflow；列出正式切换前的阻塞。
- **改了什么**：库 / 测试 / 默认值 / 参考数据 / committed gate 证据**零改动**。新增本文件、`bench/studies/devflow_adoption_audit/`，以及 adoption Change 草案 `openspec/changes/adopt-devflow/`（只是文书，devflow 尚未安装）。
- **使用者裁决（2026-11-04）**：① 合并第 6–9 轮研究分支（PR #34）；② 缓存新鲜度采用 §5 P0-2 **方案 A**（Run 自证暖态）；③ skip/xfail 采用 §7.3 **方案 A**；④ 本审计提交入库、起草 adoption Change。
- **判定**：**RECORDED** —— 采用可行，但 §一 的四个前置阻塞中阻塞 1（12 红归因）未闭合、阻塞 3/4 的机制已裁决但前提待实测；阻塞 2 已解除。
- **下一步**：读完 baseline-provenance 四臂结果（§12，运行中）后另起一轮对 12 红逐项分类；adoption Change 在阻塞 1 处置前不 `start`。


- 日期：2026-11-04
- 版本：第二次修订。初版审计 → 独立复核后合并重写（第一次修订）→ 第二轮外部审计的五项结论经逐条复测后并入（第二次修订，见附录 B）。初版原文备份在 `/tmp/up3d_audit/original-20261104.md`，第一次修订稿在 `/tmp/up3d_audit/revision1-20261104.md`（均为机器本地）
- 对象提交：UP3D `080f291`（分支 `kimi/phase7-ibl-ap-projection`）；openspec-devflow `cc0a9f1`（tag `v0.2.0`）
- 本文件最初写在被 gitignore 的 `bench/results/`，本轮提交入库。★ 文中多数数字是**审计读数**而非 committed 证据：只有 `bench/studies/devflow_adoption_audit/` 里的产物是已提交的；44 分钟全套的 console 日志未保留（§3.3），作为证据引用时以该目录与后续 baseline-provenance 产物为准（纪律 3）。

## 目的

在不改变 UP3D/pyFP3D 的求解器行为、科学证据、参考数据和历史档案的前提下，审计仓库当前的代码、测试、文档、依赖、证据链与协作状态，并评估如何基于 `~/code/openspec-devflow` 重新配置开发流程。

目标不是把项目改造成另一套文档体系。UP3D 已有的预注册、负结果、不可判结果和已提交科学证据纪律要保留；要补上的是：

- 机器可校验的依赖图与任务状态；
- 证据新鲜度；
- 独立复核与用户批准；
- Change 级 worktree 隔离；
- 可恢复交接。

# 一、结论

建议采用 openspec-devflow，但只能做**治理与证据层的叠加**：不重写仓库，也不全量反向补写 OpenSpec 规格。

UP3D 现有的研究流程在许多方面已经比普通 OpenSpec 流程更严格：

- 预注册先于代码；
- G-DOMAIN、G-TEETH、G-CADENCE、G-ATTRACTOR；
- 允许 FAIL / RECORDED / UNDEFINED；
- 只有已提交的 CSV/PNG 才算科学证据。

devflow 的真正增量有五点：

- 把 phase、研究问题、Change 和 task 的依赖关系变成机器可读图；
- 把验证 Run 绑定到提交、平台、环境、命令和实际输入指纹；
- 把实现、任务完成、当前验证、独立复核、用户批准、归档、合并和发布分开；
- 用 Change 级 branch/worktree 和资源锁规范 Kimi/Claude 并行；
- 用 handoff 保存可恢复的执行状态，而不依赖会话记忆。

## 正式切换前的四个前置阻塞

**1. 当前机默认全套不是绿线。**
- 实测 **12 failed / 827 passed / 35 skipped / 2 xfailed**。
- 第二次、测试顺序不同的定向重跑（暖缓存态）得到**同一组 12 个失败节点**，报告表中 3 位有效数字的差值一致。这**不是**位级复现：全套那次的日志没有保留，断言只打印 `.9e`，两跑都是暖态（§3.3）。
- 输入网格全部由 git 跟踪或由代码生成，已排除。
- 测试顺序与运行间散布只是**弱排除**；首要候选是**跨机器的编译/平台变量**（§3.4）。
- 必须先分类和处置。不能通过改阈值、批量刷新证据或把必需检查改成 optional 来“采用框架”。

**2. 集成分支未对齐。** ★ **已解除（2026-11-04）**：使用者裁决合并，研究分支经 PR #34 合入 `main`（merge commit `9b4f425`），本地 `main` 已 fast-forward 对齐 `origin/main`。
- 撰写时：本地 `main` 落后 `origin/main` 12 个提交；当前研究分支另有 24 个未合入提交。

**3. 两类出处 devflow 默认看不见：被 gitignore 的 M6 族网格，以及 Numba 编译出处（有效 codegen + 缓存态）。**
- 网格侧：须配置 Case raw（或 `extraInputs`）输入指纹；范围比初版所写的窄，12 红所用的网格全部被 git 跟踪，不在此列。
- 编译侧：**devflow v0.2.0 目前没有闭合这条路**。非 milestone 默认复用旧 Run，复用时项目 wrapper 根本不执行，平台指纹只有 OS 名加 provider 返回的版本串，archive 只读现有 Run 状态、不重跑。机制已由使用者裁决（2026-11-04）：**方案 A，每个 Run 自证暖态**（§5 P0-2）。尚待完成的是：实测它的两个前提，并用“通过后清缓存，再执行 status/archive”的测试和 G-TEETH 测试证明它生效。

**4. Check 判定语义与项目的 skip/xfail 习惯直接冲突。**
- pytest 的 JUnit 把 xfail 也记成 `<skipped>`；devflow 对声明了 `data` 的 Check 只要见 skip 就判 `skip`，而 archive 拦截 skip。
- 若按“网格 Check 声明 data + 保留 strict xfail + verify 不开图导出”三条同时配置，这些 Check 将**永远无法通过**（§7.3）。
- 机制已由使用者裁决（2026-11-04）：**§7.3 方案 A**，在 Check 定义里就不收集预期会 skip 的测试，xfail 随整个门文件成一个不声明 `data` 的 Check。

# 二、审计范围与实测方法

本轮未修改生产代码、测试、参考数据或既有文档：

- 执行过 `git fetch origin` 更新远端引用；
- 测试只产生被忽略的缓存；
- devflow 的 `install --check` 与测试只在 `/tmp` 副本里运行；
- 两个仓库的工作树最终均干净。

审计覆盖：

- 现行治理文档、phase 1–6 归档索引、phase 7 立项/台账/研究判定；
- `pyfp3d/`、`tests/`、`bench/`、`cases/` 的结构和静态依赖；
- 全部已跟踪 Python 文件的 AST 统计；
- 活跃树 CSV/JSON 结构扫描、Markdown 本地链接扫描、Git diff 空白检查；
- UP3D 默认全套、12 红所在门文件的独立定向重跑、实际运行环境、BLAS 钉线程行为、wheel 内容和元数据；
- 网格出处：全部 `.msh` 的跟踪状态，以及 24 份 manifest 的 SHA-256 逐一比对；
- openspec-devflow 的协议、brownfield 指南（`devflow/REFACTORING.md` §2）、配置/指纹/判定实现（逐条定位到源码行）、最小模板 `install --check`、全量测试；
- pytest JUnit 对 xfail/skip 的实际写法（最小实验）。

**未运行：**

- 重型 `PYFP3D_TRANSONIC_GATES=1` 全集、fast capability tier、贵重 demo；
- NOJIT 覆盖率；
- Windows/macOS 测试；
- 冷缓存臂；
- 任何 `NUMBA_CPU_NAME=generic` 或旧机复跑臂；
- 任何全精度（`repr`/位级）结果保存：全套那次只留下了本报告 §3.3 表中的数字，console 日志未保留。

历史归档脚本按项目规定视为 provenance，未逐条重跑。wheel 内容与 Markdown/CSV 扫描只做了一次，未独立复核。

# 三、当前基线

## 3.1 Git 与交接状态（已复核；撰写时的状态，阻塞 2 已于本轮解除）

| 项 | 状态 |
|---|---|
| 当前提交 | `080f291`，分支 `kimi/phase7-ibl-ap-projection`，无 upstream |
| `origin/main` | `6c405b2`（PR #33 合并），已含 phase 7 第 1–5 轮 |
| 当前分支 vs `origin/main` | 0 behind / 24 ahead；第 6–9 轮未合入主线 |
| 本地 `main` | `dfac561`，落后 `origin/main` 12 个提交，不能据它直接 `devflow start` |
| 24 个提交的内容 | 研究脚本、CSV/PNG 和轮次文档；相对 `origin/main`，`pyfp3d/`、`tests/`、`cases/reference_data/`、`environment.yml`、`pyproject.toml` **零差异** |

因此 12 个失败不能归因于第 6–9 轮修改了生产代码或测试。

另外，`AGENTS.md` 现行规则是“Claude 在 `main` 上直接修改和提交，Kimi 走 `kimi/*` 分支 + PR”；push、开 PR、合并和归档都需用户单独确认。采用 devflow 必须改写这条规则（§8 阶段 1 第 5 步）。

## 3.2 运行环境（已复核）

当前机：`aliyun`，Intel Xeon Platinum 8369B @ 2.90GHz，32 逻辑核，Linux 5.10 x86_64。

`up3d` Conda 环境：

- Python 3.11.18；pytest 9.2.0；
- NumPy 2.4.6、SciPy 1.17.1、Numba 0.66.0、PyAMG 5.3.0；
- meshio 5.3.5、matplotlib 3.11.2、pandas 3.1.0；
- `python -m pip check`：无损坏依赖；
- 无论先导入 `pyfp3d` 还是先导入 NumPy，`blas_thread_count()` 都返回 1。

这些版本与 phase 7 第 6、9 轮提交的 `environment.json` 一致。

★ **这一致性只说明迁移以来当前机的环境没变，说明不了产生已提交证据的那台机器。** 那两份 json 正是在这台迁移机上记录的（platform 字段为同一 `al8` 内核、同一 Python 构建）。

关于旧机，能读到的只有：

- 门证据的提交时区为 +0100；
- p6 的 `phases/p6/bench/studies/c05_amg_investigation/README.md:75` 记录为 Python 3.11.15 / NumPy 2.4.6 / SciPy 1.17.1 / PyAMG 5.3.0；
- 旧机的 **Numba 版本、CPU 型号未被记录**；
- 门 `summary.csv` 只带 `blas_threads` / `n_threads` 出处列，不带包版本或 CPU。

## 3.3 默认全套与独立复现

**全套（初版执行）：**

```bash
conda run -n up3d env \
  PYTHONNOUSERSITE=1 \
  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python -m pytest tests/ -q
```

```text
collected 876 items
12 failed, 827 passed, 35 skipped, 2 xfailed, 11 warnings
2668.05 s (0:44:28)
```

历史开发机账目为 839 / 0 / 35 / 2（phase 7 第 3 轮 `831 + 8（D15）= 839`）。总项数一致，`839 − 12 = 827`，账目闭合。

**定向重跑（复核执行）：**

- 同环境、同线程，只跑 12 红所在的 11 个门文件，测试顺序与全套不同，numba 缓存为暖态；
- 结果 `12 failed, 96 passed, 12 skipped`，1051 s（日志已提交：`bench/studies/devflow_adoption_audit/targeted_12red_20261104_console.txt`）；
- 失败节点与下表 12 项一一对应，下表的差值（3 位有效数字）与定向重跑一致。

★ **比较精度的边界**（第二次修订收窄）：
- 全套那次的 console 日志**没有保留**，现存只有下表的 3 位有效数字；
- 定向重跑的断言消息只打印 `.9e`（`tests/_gate_evidence.py:145-148`），且每个测试只报**第一个**越界列，其余列未比较；唯一例外是 D13 不对称度（`1.2022292349421009`，打印了 `repr`）；
- 两跑都是暖态。
⇒ 能说的只是“**两个已测暖态运行中，失败节点一致，显示精度的差值一致**”。要声称位级复现，须对两跑保存全精度结果（逐列 `repr` 或 `.npz`）并做位级比较。

| 门 | 失败节点 | 输入网格 | 本次 vs 已提交 |
|---|---|---|---|
| C05 | `TestCommittedEvidenceIsLoadBearing::test_fresh_run_reproduces_the_committed_summary` | 代码生成 | `x_shock` 相对差 1.07e-04 > 1e-06 |
| C07 | 同上 | 代码生成 | `gamma` 相对差 5.10e-06 |
| C08 | 同上 | 代码生成 | `rms` 相对差 4.63e-04 |
| D05 | `…::test_matches_committed_summary` | naca0012_2.5d（跟踪） | `cp_rms_upper` 相对差 2.43e-06 |
| D06 | 同上 | rae2822_2.5d（跟踪） | `cp_rms_upper` 相对差 3.78e-03 |
| D08 | 同上（散布带层） | naca0012_2.5d（跟踪） | `cp_rms_lower_sst` 绝对差 5.43e-03 > 带 3.89e-03 |
| D11 | `…::test_fresh_run_reproduces_the_committed_summary` | naca0012_2.5d（跟踪） | `mach_max` 相对差 4.87e-06 |
| D12 | 同上 | rae2822_2.5d（跟踪） | `mach_max` 相对差 1.76e-04 |
| D13 | `TestRecorded::test_surface_asymmetry_is_missing` | naca0012_2.5d（跟踪） | coarse 上/下不对称度 1.202229 > 1.15 |
| D13 | `…::test_fresh_run_reproduces_the_committed_summary`（散布带层） | 同上 | `dstar_ibl` 绝对差 1.20e-04 > 带 8.87e-05 |
| D14 | `…::test_matches_committed_summary` | naca0012_2.5d（跟踪） | `cl` 相对差 2.62e-04 |
| D15 | 同上 | naca0012_2.5d / naca0012_hex_2.5d（跟踪） | `bubble_max_cp`（tet medium 锚）相对差 1.14e-06 |

`docs/dev_phase_seven/20261102-2030-tet-rfar-prereg.md:25` 已登记“迁移机有 11 项既有证据锁或散布断言红”。本表给出了完整节点清单：11 项证据锁/散布锁，加上 D13 的一条“缺陷仍在”记录断言。

★ **D13 的不对称度断言不应与散布锁归为一类。** 该断言锁的是“我们的上/下 δ\* 比不足（历史实测 0.90–1.12）”这一已登记缺陷，红可能是好消息，也可能是状态移动。本次读数 1.202 比历史上界高约 0.08，**超出该量自己登记的散布带**（`asym_ibl` 5.6e-2，`tests/D/test_D13_ibl_vs_xfoil.py:329`）。所以它不能归入已标定的缓存出处散布，更像是 D13 的耦合末态在新机器上整体移动了。D08 同理：5.43e-3 是其散布带的 1.4 倍、实测散布包络（1.94e-3）的 2.8 倍。

已确认通过：F07 仓库不变量、F09 BLAS 确定论，以及其余 827 项。根因尚未确认（见 §3.4）。

## 3.4 根因的已排除项与剩余候选

| 候选 | 状态 | 依据 |
|---|---|---|
| 第 6–9 轮改了库/测试 | **排除** | §3.1 零差异 |
| 测试顺序 / 进程内前序状态 | 弱排除 | 全套与定向重跑顺序不同，失败节点与 3 位有效数字差值一致；另有 p6 第 44 轮实测“单跑 vs 先跑 A+B 层再单跑 = 逐位相同”（`phases/p6/docs/dev_phase_six/20261011-1800-process-history-determinism.md:78-83`） |
| 当前机运行间散布 | 弱排除 | 只比到显示精度，且两跑都暖态（§3.3）；位级结论需全精度保存 |
| 输入网格不同 | **排除** | D 类网格全部被 git 跟踪且工作树干净；C05/C07/C08 网格由代码生成 |
| BLAS 线程 | **排除** | 运行时实测 1 |
| 包版本（NumPy/SciPy/PyAMG） | 弱排除 | 与 p6 记录相同 |
| Python 补丁版本 3.11.15 → 3.11.18 | 未排除 | 旧机仅此一条记录 |
| Numba/LLVM 版本 | **未排除** | 旧机版本未记录 |
| **CPU 指令集相关的代码生成** | **未排除，首要候选** | 内核为 `@njit(cache=True, fastmath=True, parallel=True)`（`pyfp3d/` 内 33 处 `fastmath`），按宿主 CPU 编译；fastmath 允许 FMA 收缩与重排，不同 ISA 下舍入不同；敏感腿（折叠区、多吸引子邻域、IBL 地板）会放大这类差异 |
| 冷/暖缓存出处 | 未测 | 纪律 11 已登记 cache-load 与 fresh-compile 不逐位一致；p6 第 44 轮实测旧机上“冷编译跑”有 8 条证据锁红、“暖缓存跑”复现到 ≤ 4.5e-10；定向重跑为暖态，且复现了同一组 12 红（全套那一跑的缓存态没有日志可查），所以 12 红**不能**归为冷编译扰动 |

★ **初版提出的“同 commit、同网格 SHA、同包版本，在当前机做冷/暖缓存 + 测试顺序单变量 A/B”判别不了首要候选**：其中测试顺序已被排除，剩下的一臂（冷/暖缓存）只能回答缓存问题。要判别跨机器变量，至少需要以下之一：

1. **ISA 臂**：`NUMBA_CPU_NAME=generic` 配独立的空 `NUMBA_CACHE_DIR`。实测机制：设 `generic` 时 `NUMBA_CPU_FEATURES` 默认取空串（`numba/core/config.py:306-308`），所以它确实是去特性的通用目标；Numba 的缓存键本身含 CPU 名与 feature 串（`codegen.py:1298-1303` 的 `magic_tuple`），所以独立缓存目录不是为了防止误载，而是为了让冷/暖态明确、不污染 native 缓存。
   ★ **独立空缓存意味着 generic 首跑同时也是冷编译**，而项目已实测冷/暖态本身就会移动这些数（纪律 11、p6 第 44 轮）。所以**不能**拿 generic 首跑与已有的 native 暖态比较并归因于 ISA。对照关系须事前写死：
   - **冷态 ISA 敏感性**：native 冷（清空 native 缓存目录后首跑）对 generic 冷（空缓存首跑）；
   - **暖态 ISA 敏感性**（与 12 红同态，是真正要回答的问题）：在上一对之后各自再跑一次，比 native 暖对 generic 暖；
   - 两对都要保存全精度结果（§3.3），否则比较精度不够。
   ★ **结论范围**：`NUMBA_CPU_NAME` 只控制 Numba/LLVM 的 JIT 代码生成，不控制 OpenBLAS、SciPy、PyAMG 等**预编译本地库**的运行时 CPU 分派（例如 OpenBLAS 按 CPU 选 kernel，由它自己的 `OPENBLAS_CORETYPE` 控制）。所以 ISA 臂只能检验“答案对 Numba JIT 代码生成敏感”，不能排除预编译库的 ISA 敏感。它也**不会**复现旧机的逐位结果。
2. **旧机复跑臂**：如果旧机仍可访问，用同一 commit 在旧机重跑这 11 个文件，这是唯一的决定性判别。
3. 若两者都做不到，处置只能是登记为“跨平台证据差、根因 UNDEFINED”，**由用户裁决**是否生成平台专用证据或再基线；不得借 adoption 夹带。

无论哪种结论，都要为今后的证据补上出处列，否则下一次迁移会重演同一问题。只写“宿主 CPU 型号 / `numba.config.CPU_NAME`”**不够**：当前环境实测 `CPU_NAME = None`、`CPU_FEATURES = None`，实际编译用的是 LLVM 检测到的宿主 CPU（`icelake-server`）和 114 项 feature 组成的 feature 串（字符串长 939 字符）（`numba/core/codegen.py:1322-1330`）。应记录**有效** codegen：
- numba / llvmlite / LLVM 版本（当前 0.66.0 / 0.48.0 / 22.1.8）、target triple；
- 有效 CPU 名（`CPU_NAME` 覆盖值，或 `llvmlite.binding.get_host_cpu_name()`）；
- 有效 feature 串（`CPU_FEATURES` 覆盖值，或 `get_host_cpu_features()` 经 `NUMBA_ENABLE_AVX` 处理后的结果，即 `numba/core/codegen.py:1427-1445`）；
- 影响 codegen 的配置覆盖：`NUMBA_CPU_NAME`、`NUMBA_CPU_FEATURES`、`NUMBA_ENABLE_AVX`、`NUMBA_OPT`，以及 SVML 是否可用；
- Python 构建、BLAS 库及其 coretype。
否则虚拟机的 feature mask 变化，或只改 `NUMBA_CPU_FEATURES`，都不会让旧证据或旧 Run 失效。

# 四、仓库与科学能力审计

## 4.1 架构和规模（已复核）

| 范围 | 规模 |
|---|---|
| `pyfp3d/` | 51 个 Python 文件、21,920 行、460 个函数 |
| `tests/` | 119 个 Python 文件，其中 100 个 `test_*.py`（A 57 / B 6 / C 10 / D 15 / E 3 / F 9），786 个 AST 测试函数；参数化后收集 876 项 |
| `bench/` | 26 个 Python 文件、11,368 行 |
| `cases/` | 30 个 Python 文件、11,730 行 |
| `phases/` | 284 个归档 Python 文件、70,527 行（历史快照） |
| `cases/gates/` | 57 个已跟踪证据文件，24 个 `summary*.csv` |
| `cases/meshes/` | 39 个 `.msh`：**19 个被 git 跟踪**（29 MB，全部 2.5-D 族、圆柱、球壳、表面网格）、**20 个被 gitignore**（ONERA M6 四族 17 个、fuselage_bor 2 个、naca0012 fine 1 个）；**24 份 `.msh.manifest.json`**，覆盖全部 20 个被忽略网格及 4 个被跟踪网格，SHA-256 **全部与当前文件一致** |
| `cases/reference_data/` | 152 个已跟踪参考文件，当前分支零修改 |

子包静态依赖图无环：

```text
constraints → mesh
kernels     → mesh, physics
meshgen     → mesh
post        → _compat, mesh, physics
solve       → constraints, kernels, physics
viscous     → physics, post, solve
```

没有必要为采用 devflow 重排求解器。

大模块：`pyfp3d/solve/newton.py`（2,236 行）、`pyfp3d/viscous/ibl3.py`（1,642 行）、`pyfp3d/viscous/closures.py`（1,169 行）。它们是后续垂直重构的候选，不应混入框架采用 Change。

## 4.2 科学能力不能过度声明（已复核）

现行能力边界 `phases/p6/docs/dev_phase_six/PHASE_SIX_CAPABILITY_BOUNDARY.md`：

- 无粘内核对同层级 Euler 参照已有定量闭合，但不是全包线可信；
- α=0 的实测可信上界只有 M0.78–0.80；
- 有粘耦合在举力机翼上仍跑不起来（D09/D10）；
- `converged ∧ 零钳制` 不是可用性的充分条件；
- M1a ✓、M1b ✓（两种子）、M1c ◐（seed 5 PASS / seed 0 覆盖 FAIL）、M3b ◐（medium 达标、无锁）、M2/M4 未动、M5 仅一条子句达成。

phase 7 的 AMG 方向 2 在 C05 尺度研究中 PASS，但生产库未改，生产尺度代价未测。S4/S5 及后续 IBL 研究（第 5–9 轮）仍有 UNDEFINED 项。

当前项目应描述为活跃的研究求解器，而不是已完成产品化验收的通用 CFD 软件。

## 4.3 现有优势必须保留

- `docs/agent-rules.md:68-79` 的轮次制度（预注册、G-TEETH、四值判定、kill criterion 照章）与 devflow 的 proposal/design/tasks 很匹配；
- `tests/_gate_evidence.py:60` 的 `assert_matches_committed` 与 `assert_within_spread` 已形成“新鲜复算对 committed 证据”的两层锁（出处层 1e-6 + 散布带层）；
- F01/F07/F09/F10 覆盖网格出处、路径、BLAS 和风险↔门索引；
- strict xfail 和“缺陷仍在，红=好消息”式门保留负结果，而不是强迫换绿；
- 参考数据不可修改、贵重证据不重算、归档脚本不冒充可运行现状，边界写得清楚。

这些应成为 devflow 的配置和复核规则，而不是被通用模板替代。§3.4 也说明：这套锁在迁移时**如实变红**，这正是它该做的事。

# 五、主要工程问题与优先级

## P0：采用前必须处理

### 1. 基线复现未闭合

当前 12 红不能写进“已知问题”后继续宣称绿线。需要先做独立的 `baseline-provenance` Change：

- 固定相同 commit、包版本、BLAS、线程；网格已由 git 与 manifest 固定，只需记录 SHA；
- 按 §3.4 补做**冷缓存臂**和 **ISA 臂**（`NUMBA_CPU_NAME=generic` + 独立缓存目录，按事前写死的冷/冷、暖/暖两对比较），可行时加**旧机复跑臂**；所有臂保存全精度结果；测试顺序臂不作为主臂（已弱排除）；
- 逐项区分：平台/编译敏感（期望大多数）、缓存出处敏感、真实回归；
- 输出已提交的比较 CSV 和逐项处置，同时为证据出处列补上 Numba 版本与 CPU 信息；
- D13 不对称度断言单独处置：它是缺陷记录断言，不是散布锁，变红需要读物理意义（§3.3）；
- 不放宽带、不删除测试、不批量重刷证据。

若处置需要新的证据基线（含“平台专用证据”），须另走用户裁决和全文档勘误（纪律 10）；若是回归，单独修复。

正式 adoption 的验收不能以未处置的红色必需检查为背景。

### 2. 被忽略输入与 Run 新鲜度

**devflow 的实际行为**（`devflow/fingerprint.py`）：

- 内容指纹读 HEAD 的 `git ls-tree`（:137-143），脏检查用不带 `--ignored` 的 `git status`（:51-73）；所以声明目录里的被忽略文件不进指纹；
- 若在 `sourcePaths`（或 `environmentFiles`）里直接声明一个 HEAD 里没有、checkout 里有的文件，它算 uncommitted（:150-159、:434-438），并在 `verification.py:448-456` 阻止 Run；
- 两处框架文档在这一点上**一致**（第一次修订称其“有张力”，有误）：README.md:121 同一段既说“declare that file itself in `sourcePaths`”，也说“a declared input that is untracked or ignored by Git counts as uncommitted”；`openspec/README.md:44` 说的是同一件事。也就是说，声明被忽略文件的效果就是逼你提交它。

**被忽略文件有三条路，不是一条**（第一次修订写“只能走 Case raw”，是事实错误）：

| 路径 | 指纹 | 脏检查 | 新 Change worktree |
|---|---|---|---|
| `sourcePaths` 直接声明 | — | 判 uncommitted，**阻止 Run** | — |
| Case 的 `extraInputs` | 直接读**工作树**文件算 SHA-256（`fingerprint.py:246-250`） | 进入 `git status` 检查（:396），但不带 `--ignored`，被忽略文件的编辑**不报 dirty** | 文件不在时哈希为 `None`（:249），**不报错**，Run 照常进行 |
| Case raw `<dataRoot>/<case_id>/raw/<localName>` | 读 raw 文件算 SHA-256（:233、:244） | raw 落在本 checkout 或主 checkout 内时进入 `git status`（:405-420），同样不带 `--ignored` ⇒ **被忽略的 raw 也不报 dirty** | 外部 `dataRoot` 由主 checkout 与各 worktree 共享 |

独立复测（`bench/studies/devflow_adoption_audit/exp_extra_inputs.py`，输出 `exp_extra_inputs_output.txt`；用 devflow 自身测试夹具的同款构造，含对照臂）：

```text
arm1 extra-input hash changed: True | dirty after ignored edit: []
arm1 ignored in-checkout raw: hash changed: True | dirty: []
arm2 (tracked control) dirty after edit: ['cases/meshes/m6.msh']
arm3 absent ignored extra-input: {'cases/meshes/m6.msh': None} | dirty: []
```

⇒ 被忽略文件不是“只能走 raw”，raw 放在主 checkout 内也不必然被拦截，只要它被 Git 忽略即可。另须知道：`extraInputs` 在 README.md 与 `openspec/README.md` 里**都没有文档**，只出现在 `registry.py:192-221`、`fingerprint.py` 和单元测试中，属于未承诺稳定的接口。

**受影响范围：** 只有 §4.1 列出的 20 个被忽略网格（M6 四族为主）。19 个被跟踪网格由 git 定址，不需要额外机制。

**推荐方案（设计选择，不是框架硬约束）：**

- 使用外部 `UP3D_DATA_DIR`（`dataEnvironmentVariable` 会覆盖 `dataRoot`，`devflow/project.py:22-38`）和 Case raw manifest；
- raw 文件放在 `<dataRoot>/<case_id>/raw/<localName>`（这是 raw 路径的唯一布局，`fingerprint.py:233-244`），manifest 每项要有 `localName`；
- 推荐把 `dataRoot` 放在主 checkout 与各 Change worktree **之外**。理由是跨 worktree 共享，以及不依赖 `.gitignore` 规则是否覆盖到 raw 目录，而不是“放在内部会被拦截”。被忽略的内部 raw 同样可以工作；
- 不推荐用 `extraInputs` 指向被忽略的网格：新 worktree 里文件缺失只会得到 `None` 哈希而不会报错（正是验收条件 2 要防的“依赖主 checkout 里的偶然文件”），且该接口未文档化；
- Check 前由项目 wrapper 校验 raw SHA 与已提交的 `.msh.manifest.json`（现有 24 份正好覆盖全部被忽略网格），再物化到 worktree 的 `cases/meshes/`；
- 网格有意变化时，在同一提交里更新 raw、manifest、相关证据和勘误；
- 成本须知：raw 文件每个实例都要全量重算 SHA-256（缓存只在单实例内，`fingerprint.py:103,207-218`），M6 大网格每次 Run 都要付这笔 I/O。

不要把大型网格强行提交，也不要只依赖当前目录里的偶然文件。

**Numba 编译出处同样是答案出处，而 devflow v0.2.0 下它没有闭环**（第二次修订重写本段）。

第一次修订写了“清缓存或换机不能复用旧通过”（§8 阶段 2、§9 条 4），但只给了目标，没有机制。实测/读码的框架行为：

- 非 milestone profile 默认复用：`reuse = not force and profile not in MILESTONE_PROFILES`（`verification.py:806`）；
- 复用判定只比内容指纹、OS 名和 `packages`（`verification.py:153-187`，环境比较在 :170-176）；平台指纹就是 `system` 加 provider 返回的版本串（`fingerprint.py:301-303`）；
- 清 numba 缓存既不改内容指纹，也不改任何包版本 ⇒ 旧 Run 照常复用；
- **复用时不执行 Check 命令**，所以项目 wrapper 根本不运行，wrapper 里的暖缓存前置或拦截在复用路径上一律不生效；
- milestone 只在 **verify** 时强制新鲜（:806）；archive 经 `archive_blocker` 读 `status["verified"]`，即现有 Run 的有效性（`execution.py:1098-1113`、`archive.py:73-77`），**不重跑**。所以“milestone verify → 清缓存 → archive”不会被拦截；
- 注册表里没有 per-check 强制新鲜的字段，复用只由 profile、`--force` 和 `manual` 决定。

这里其实有两类事件，需要两套机制：

**(1) 换机 / 换 codegen**（CPU、feature mask、numba/llvmlite/LLVM 版本、覆盖变量）。框架原生可做：在 `environmentPackages` 里加一个伪包（例如 `numba-codegen`），由 provider 以子进程进入真实的 `up3d` 环境，返回 §3.4 所列有效 codegen 字段的摘要。Run 的环境在执行**之后**记录（`verification.py:771` 的 `_environment()`），复用时整体比较。代价是任何变化都让本平台全部 Run 失效，对换机来说这正是想要的。

**(2) 清缓存 / 缓存态。** ★ **使用者裁决（2026-11-04）：采用方案 A。** 方案 B、C 保留在下面，作为已考虑但未采用的备选，并写明理由，供 adoption design.md 引用：

- **方案 A（已采用）：每个 Run 自证暖态。** 项目自己对规范态的定义是“测量进程内不发生编译”（p6 第 44 轮勘误 §1）。具体做法：
  - wrapper 先在**独立进程**里白跑一遍（纪律 11）；
  - 敏感证据 Check 在测量后读取 pyfp3d 各 dispatcher 的 `stats.cache_misses`（`numba/core/dispatcher.py:167,952`）；任一未命中 ⇒ 该 Check 判 **fail**（不是 skip），命中/未命中计数写进 step trace。
  
  这样一个 pass 本身就证明了“在规范暖态、在记录的 codegen 下测得”。之后清缓存并不会让它变错，因为它陈述的是“代码 + codegen ⇒ 答案”，不是缓存目录的状态。相应地，§9 条 4 改写为“冷态测量不能产生 pass”，而不是“清缓存使旧 Run 失效”。待实测的前提有两个：dispatcher 统计能覆盖 pyfp3d 的全部 JIT 入口（`cache=True` 共 61 处，部分经 `numba.njit(*args, **kwargs)` 包装生成）；“完全暖缓存进程复现出处态”在旧机上只测过两次。
- **方案 B（未采用）：缓存态伪包。** provider 返回 pyfp3d 缓存索引的摘要，清缓存后摘要改变 ⇒ Run 失效。缺点有两个：`packages` 是整体比较，任何 Check 编译出新特化都会让本平台**全部** Run 失效；索引内容跨重编译是否确定未测。最坏情况是过度失效（安全，但代价高）。
- **方案 C（未采用）：框架扩展。** 加 per-check always-fresh（复用路径对该 Check 恒判不可复用），或让 archive 对指定 Check 重跑。这是最直接的语义，但要改 devflow，须走 devflow 自己的 Change，采用要等它发版。

**不推荐**只把敏感门放进 milestone：archive 不重跑，填不上 verify 与 archive 之间的缺口。

按方案 A，§8 阶段 2 必须有一条实测：**通过后清缓存，再执行 status/archive**。期望结果写进 design.md：status 仍为 verified，这是设计如此，因为那次 pass 已自证是暖态下测得的；但紧接着的冷态新鲜运行（测量进程内发生编译）必须判 **fail**。

方案 A 的两个前提须在阶段 2 **先实测**，再依赖：
- dispatcher 统计能覆盖 pyfp3d 的全部 JIT 入口（G-TEETH：故意让一个内核在测量进程内编译，断言必须变红）；
- 暖态在当前机上可复现：同一 codegen 下连续两次暖跑，全精度逐位一致。
若任一前提不成立，回到本节重新裁决，不得静默退到 B 或 C。

只有已证明与编译出处无关的轻量检查（schema、F07/F10、`--analyze-only`）才默认复用 Run。

### 3. Git 切换点

- 先完成当前 24 个研究提交的复核/合入，或明确冻结其交接；
- 再在授权后 fast-forward 本地 `main` 到 `origin/main`；
- 不改写现有研究历史，也不把旧轮次伪造成已完成的 devflow Change。

## P1：采用过程中必须显式配置

### 4. 两套 Python 环境必须分离

UP3D 固定 Python 3.11；devflow 要求 Python ≥3.12（`pyproject.toml` 的 `requires-python`）。应新建独立的 `up3d-devflow` 工具环境，而不是改变求解器环境：

- Python 3.12、Node.js ≥22.12（README.md:161、模板 `package.json:5`）；
- devflow 固定 v0.2.0 / `cc0a9f1`；
- OpenSpec CLI 固定 1.16.0，**必须**提交 `package.json` 与 `package-lock.json`：1.16.0 只在模板里精确钉死，`archive.py` 运行时不校验版本，会依次取 `node_modules/.bin/openspec` 或 PATH 上的 `openspec`；
- `devflowCommand` 指向该工具环境；
- 每个求解 Check 显式通过 `conda run --no-capture-output -n up3d ...` 进入 Python 3.11 环境。

### 5. 科学判定与执行状态不能混同

研究的 FAIL/RECORDED/UNDEFINED 可以是按预注册成功交付，不等于 devflow task 失败。proposal/design 必须预先定义：kill criterion 触发或结果不可判时，交付完整证据和判定即算完成任务。

`tasks.md` 只管任务文本和 checkbox；`progress.md` 保留为科学结果索引，不再成为另一套任务执行状态库。

### 6. Check 的 skip、xfail 与 milestone 语义（见 §7.3，属阻塞级）

devflow 的 JUnit 判定（`verification.py:544-583`）：

| 情形 | 判定 |
|---|---|
| 零个测试 | `fail` |
| 有失败/错误 | `fail` |
| 全部 skip | `skip` |
| 声明 `data` 且有任一 skip | `skip` |
| 未声明 `data` 且部分 skip | `pass` |

archive 把 `skip` 与 `fail` 一样拦下（`execution.py:1113`）。

milestone 不复用 Run（`verification.py:805`）；manual Check 的 operator 证据在内容与环境未变时，在任何 profile 下都复用（:810）。

由此推出：

- 所有 pytest Checks 都必须输出 JUnit，不能只看 console 绿色；
- **pytest 把 xfail 写成 `<skipped type="pytest.xfail">`**（最小实验实测：strict xfail + 普通 skip ⇒ `skipped="2"`；`bench/studies/devflow_adoption_audit/exp_junit_xfail.py`），所以 xfail 在 devflow 眼里就是 skip；
- milestone 会新鲜执行所有非 manual Checks，所以不能把一次性贵重研究测量注册成全局自动 Check；它们应作为带已提交 CSV/PNG 的 manual Check；自测、schema 检查和 `--analyze-only` 才适合作为轻量自动 Check。

### 7. 依赖元数据不一致（已复核）

- `pyproject.toml:26` 声明 `scipy>=1.10`，但 README.md:21 与 `pyfp3d/solve/linear.py` 的 `rtol=` 调用要求 `>=1.12`；构建出的 wheel 元数据同样写成 `>=1.10`（初版实测，未独立复核）；
- `dev` extra 未声明 pandas，但 `bench/studies/v5_2_rae2822/run.py:46` 在 `import pandas`，且被 `tests/A/test_A19_meshgen_rae2822.py`、`tests/D/test_D12_rae2822_experiment_bias.py` 加载（`environment.yml:88` 有 pandas，pip 路径没有）；
- wheel 共 56 个文件，只含 `pyfp3d` 包和 dist-info，没有 tests/bench 泄漏（初版实测）；
- `up3d` 环境中没有 ruff、mypy、build，仓库没有 CI（无 `.github/`），不能声称这些检查有常规 cadence。

依赖修复应单独成为 Change，并在干净环境中实际安装和运行，不与 devflow adoption 混合。

### 8. 文档与证据治理漂移（已复核）

- `PROJECT_STRUCTURE.md:563` 起仍描述已不存在的 `cases/analysis/`；
- `PROJECT_STRUCTURE.md:634-665` 仍是旧测试规模（88 文件 / 585 条、D 10 文件、“D = gated 全集”）；
- `CLAUDE.md:22` 写 99 个测试文件，实际 100（多出的是 2026-10-12 加入的 D15）；
- `docs/overview.md:21-22` 仍说 `tests/D` 是门控全集，与 `CLAUDE.md:670` 起的勘误（门控由环境变量定义）不一致；
- `cases/gates/INDEX.md` 有 D15 节而无 D14 节，但 `cases/gates/D14_alpha0_transonic_nonuniqueness/` 存在（phase 7 第 3 轮已登记未修）；
- `bench/gate_results/capability_matrix.csv` 表头 19 列，70 行 19 列、64 行 24 列；写入端已在 `bench/capability_matrix.py:137-151` 拒绝继续追加，但既有的混合 schema 文件仍需数据兼容性处置。

Markdown 扫描未发现经人工确认的失效本地链接，问题主要是内容、计数和语义漂移。通用 CSV 扫描有 17 个结构告警，其中包含注释前言和自由文本格式，不能一概判为损坏；上述混合 schema 是已确认的具体问题。

# 六、推荐的目标模型

| UP3D 现有概念 | devflow/OpenSpec 中的落点 |
|---|---|
| phase 7 范围与 S1–S6 | `openspec/plan/roadmap.json` 的 phase/lane/item 和依赖 |
| 一个独立的科学问题或工程切片 | 一个 Change；不把整个 phase 或每次 commit 机械地变成一个 Change |
| 轮次预注册、判据、kill criterion | `proposal.md` + `design.md`；G-DOMAIN/G-TEETH/G-CADENCE/G-ATTRACTOR 为项目必填约束 |
| 仪器、自测、测量、判定、台账 | `tasks.md` + `execution.json` 的依赖、Checks、重试和 checkpoint |
| 单元/回归/能力门 | 项目拥有的 `openspec/checks/` |
| 被忽略的大型网格（M6 族等） | 外部 `dataRoot` + `openspec/cases/` raw manifest（布局见 §5.2） |
| 被跟踪的网格 | 照旧由 git 定址，不进 Case raw |
| 贵重研究测量、人工视觉复核 | `manual: true` 的 Check，列入 `requirements.manualChecks`，用 `devflow record-manual` 记录，指向当前提交中的 CSV/PNG |
| committed CSV/PNG | 长期科学证据，继续按现有规则提交 |
| devflow Run/JUnit/console | 机器本地执行证据，放 `runs/`，不替代 committed 科学证据 |
| 独立审计 | 真正独立的 reviewer；产品变化后重新复核，不能只换个 reviewer 名字 |
| 用户裁决 | checkpoint、archive approval；与 merge/push/release 分开 |
| OpenSpec archive | Change/acceptance 的权威归档，留在 `openspec/changes/archive/` |
| `phases/pN/` | 科学阶段档案和一次性研究证据；与 OpenSpec archive 通过 ID 链接，不互相搬运 |

采用 tooling-only Change 时，`.openspec.yaml` 设 `skip_specs: true`（`devflow/archive.py:17-22` 匹配顶层该键）。不为整个仓库反向编写规格；只有当前或目标可观察行为发生变化时，才为对应 capability 写 spec delta。不要把内部函数、历史轮次和未实现的愿望复制成 specs。

# 七、Check 与配置蓝图

## 7.1 初始 Check 分层

devflow 只有四个 profile：`PROFILES = ("docs", "task", "affected", "milestone")`（`verification.py:25`）。“manual”不是第五层，而是 Check 的属性，可以出现在任何层。

| profile | 内容 |
|---|---|
| `task` | A01 freestream、F07 仓库不变量、F10 风险索引、当前研究脚本自测 |
| `affected` | A/B 及受影响的 C/D/E/F、研究 `--analyze-only`、schema/证据一致性检查 |
| `docs` | F07/F10，加上文档计数、INDEX 覆盖、测试分层语义等专项检查 |
| `milestone` | fast capability locks、门控全集（含默认项）、环境/参考数据前置检查，全部新鲜执行 |

manual Check（贯穿各层）：贵重研究测量和人工视觉复核，记录必须指向当前提交中的 CSV/PNG。

求解器/内核改动仍先跑 A01。普通研究 Change 只运行受影响的检查，不把“每轮都重跑全套或贵重 demo”设成新规矩；只有 phase 边界的 milestone Change 才要求完整门控。

## 7.2 `openspec/project.json` 的关键选择

以下 12 个键均已在 devflow 源码中确认存在（另有 `dataRoot`，默认 `data`）：

- `integrationBranch: main`，`integrationRemote: origin`；
- `devflowCommand` 指向独立的 `up3d-devflow` 环境；
- `evidenceRoot: runs/evidence`，`resourceStateRoot: runs/devflow`；
- `dataEnvironmentVariable: UP3D_DATA_DIR`，要求机器本地设置外部绝对路径；
- `evidenceDirectoryVariable` 和 `stepTraceVariable`，供项目 wrapper 写 JUnit 和步骤证据；
- `environmentFiles`：求解与工具环境、构建文件、OpenSpec CLI lockfile、项目配置和环境 provider；
- `environmentPackages`：Python（含构建）、NumPy、SciPy、Numba、llvmlite、PyAMG、meshio、gmsh，**以及一个有效 codegen 伪包**（§3.4 所列字段的摘要：有效 CPU 名、有效 feature 串、LLVM 版本、相关 `NUMBA_*` 覆盖、SVML、BLAS coretype）。只写 `numba.config.CPU_NAME` 不够：默认值是 `None`，换 CPU 或改 feature mask 都不会改变它。provider 本身在 devflow 的 Python 3.12 环境里被导入（`project.py:41-66`），所以必须以子进程查询真实的 `up3d` 环境，而不是读 devflow 环境；
- CFL3D/XFOIL 只在需要再生成时强制要求，记录二进制路径、版本和 SHA；普通使用 committed 参照数据的检查不要假装需要外部二进制；
- `productPaths` 覆盖真实的活跃代码、测试、文档、证据和 workflow 配置，不用窄声明掩盖 uncovered changes；
- 初期 `resources: {"solver-cpu": {"slots": 1}}`，串行重型验证，防止并行 Change 超订 CPU、内存和共享数据。

`runs/` 与 `node_modules/` 加入 `.gitignore`（当前均未加）。active Runs 不靠提交传输，跨机器用 handoff。

## 7.3 skip/xfail 与 JUnit：机制（已裁决：方案 A）

**冲突的三个事实：**

1. 21 个 C/D 门文件引用了 `gate_figures_enabled()`。本次定向重跑的 11 个文件每个都至少有一条默认 skip 的图导出腿（共 12 skipped）；按纪律，正式 verify **不得**设 `PYFP3D_GATE_FIGURES=1`。
2. 默认全套有 2 条 xfail（含 G1.6-literal、G1.3 这类 strict 负结果锁），JUnit 里它们记为 skipped。
3. 网格专项 Check 若声明 `data`（为防缺网格时静默变绿），任一 skip 都会让它判 `skip`，进而被 archive 拦下。

所以初版“wrapper 校验预期 skip/xfail 清单”这一句不够：devflow 自己解析 JUnit，wrapper 只看 exit code 改变不了判定。★ **使用者裁决（2026-11-04）：采用方案 A。** 方案 B 保留为未采用的备选。

- **方案 A（已采用，改动在 Check 定义）：**
  - 每个 Check 用 `--deselect` / `-m` 显式排除图导出腿，使默认 verify 下的 Check 不含预期 skip；
  - 含 strict xfail 的门文件（默认套件：`tests/C/test_C02_cylinder_wall_correction.py:235`、`tests/C/test_C03_laplace_sphere.py:128`；门控集另有 `tests/E/test_E01_p8_newton_anchors.py:111`）各自以**整个文件**（排除图导出腿后）成一个**不声明 `data`** 的 Check；
  - ★★ **不能把 xfail 节点单独拎成 Check**（第二次修订后自查修正）：devflow 读 `<testsuite>` 的 `tests`/`skipped` 计数（`evidence.py:199-216`），而“全部 skip ⇒ `skip`”不论是否声明 `data`（`verification.py:577-579`）。只含 xfail 的 Check 是 `tests=skipped`，同样永远卡在 `skip`。随整个文件跑，同文件的通过项保证 `tests > skipped`，于是落进“未声明 `data`、部分 skip ⇒ `pass`”。启用前须逐文件确认这一点成立，否则给该 Check 配一个同门的通过项；
  - wrapper **只能往红里判**：它对照该 Check 预期的 xfail 节点名单与 skip 节点名单（应为空），任何不符、任何 XPASS 都以非零退出码返回，devflow 据此判 `fail`（`verification.py:562-564`）。wrapper 从不改写 JUnit，也不把 skip 变成 pass；
  - 缺网格的保护由 wrapper 的 raw/manifest 前置校验承担（同样以非零退出码判 `fail`），而不是靠 `data` 声明。
- **方案 B（未采用，改动在 wrapper）：** wrapper 在比对预期 skip/xfail 清单后，产出过滤后的 JUnit 供 devflow 读取，原始 JUnit 一并存档。代价是 devflow 看到的已不是 pytest 原始输出，复核者必须能读到两份。

无论选哪种，都要实测：缺网格不能变绿；预期外的 skip 不能变绿；strict xfail 意外通过（XPASS strict ⇒ fail）仍然是红。

## 7.4 项目 wrapper 的必要职责

- 进入真实的 `up3d` 环境，设置 `PYTHONNOUSERSITE=1`、线程上限和 headless 后端；
- 校验/物化 Case raw 网格，核对 committed mesh manifest；
- 检查实际 BLAS（运行时读库，而不只是看环境变量）；
- 记录有效 codegen 字段与缓存目录；对敏感门在独立进程里白跑暖缓存，再测量，并按 §5 P0-2 的方案 A（已裁决）自证暖态：测量进程内零 `cache_misses`，否则判 fail；
- ★ wrapper 只在 Check **真正执行**时运行，复用路径上不运行。所以任何“wrapper 会拦下”的保证，都只覆盖新鲜运行；复用的正确性只能靠指纹（§5 P0-2）。
- 按 §7.3 方案 A（已裁决）输出原始 JUnit，并校验预期 xfail 名单、预期 skip 为空、无 XPASS、零覆盖保护；任何不符以非零退出码返回；
- 不在正式 verify 中设置 `PYFP3D_GATE_FIGURES=1`。证据刷新是独立的测量/再基线动作：先生成并提交证据，再做只读验证。

# 八、分阶段实施计划

## 阶段 0：基线和交接处置

1. 建立独立的 `baseline-provenance` Change：按 §3.4 / §5.1 补做冷缓存臂、ISA 臂（按冷/冷、暖/暖两对事前写死，可行时加旧机臂），所有臂保存全精度结果，逐项分类 12 红，形成已提交的比较证据，并补齐有效 codegen 出处列。
2. 完成当前研究分支的复核/合入，或明确冻结；授权后对齐本地 `main` 与 `origin/main`。
3. 明确：任何再基线、平台专用证据或真实修复都需要独立的用户裁决，不夹带在 adoption 中。

## 阶段 1：tooling-only adoption

遵循 devflow 的 brownfield 指南（`devflow/REFACTORING.md` §2）：先只读评估，在未修改的 commit 上复现行为基线，再把框架接入作为一个独立的 tooling Change，不夹带产品重构；Check 定义提交后才能作为干净的指纹输入。

1. 新建 `up3d-devflow` 工具环境，固定 devflow/OpenSpec 版本。
2. 执行 `devflow install --host claude --host kimi`（两个 host 都受支持，`install.py:23`），提交 schema、适配器和安装记录。注意：最小模板 `install --check` 报 `outdated`，需补 6 个 schema 文件；加两个 host 后另需 18 个 SKILL.md。
3. 新增 `openspec/config.yaml`、`project.json`、初始轻量 Checks、`package-lock.json` 和 ignore 规则。
4. 把 `AGENTS.md` 设为协作与状态转换的单一权威，`CLAUDE.md` 保持薄入口，不复制第三份协议。
5. 改写 `AGENTS.md` 的现行规则：Claude 不再“在 `main` 上直接修改和提交”，同样走 Change 分支；Kimi 分支使用 `kimi/<change-id>`；push/PR 仍需确认。
6. adoption 本身设 `skip_specs: true`，不改 `pyfp3d/`、测试、参考数据和 phases。
7. 验证 `list-checks`、`list-cases`、plan-only verification、OpenSpec strict validation，以及一次真实的 worktree smoke。

## 阶段 2：出处与验证适配

1. 实现外部 Case raw 数据（按 `<dataRoot>/<case_id>/raw/<localName>` 布局）、网格 manifest 校验和 worktree 物化。
2. 按 §7.3 方案 A（已裁决）实现 Check 定义中的 `--deselect`/`-m` 与 xfail 门文件的 Check 划分，以及 wrapper 的 xfail/skip 名单校验；按 §5 P0-2 的方案 A（已裁决）实现暖缓存前置与暖态自证，并先实测方案 A 的两个前提。
3. 实现真实环境（含有效 codegen 伪包）/外部工具版本 provider 和资源锁。
4. 逐条证明：
   - 改 raw 网格会使 Run 失效；
   - 只改 `NUMBA_CPU_FEATURES`（或 `NUMBA_CPU_NAME=generic`）会让 codegen 伪包变化、旧 Run 判环境不符；
   - **通过后清缓存，再执行 status 和 archive**：status 仍为 verified（方案 A 的设计语义），紧接着的冷态新鲜运行判 fail；
   - 故意让一个 pyfp3d 内核在测量进程内编译，暖态自证断言变红（G-TEETH）；
   - 冷态下（测量进程内发生编译）的新鲜运行不能产生 pass；
   - 缺网格不能 skip 成绿；
   - strict xfail 意外通过会红；
   - 含 xfail 的门文件 Check 在正常状态下判 `pass`（不是 `skip`）；
   - 有人新加一条默认 skip 的测试而没有更新 Check 定义时，wrapper 判红；
   - 新 worktree 不依赖主 checkout 里的偶然文件。

## 阶段 3：一次完整的真实试点

建议以当前 IBL 内层收敛问题作为一个边界清楚的 Change，因为它能实际检验“研究结论为 UNDEFINED 也可按规程交付”：

1. proposal/design 预先写死判据、时间盒和 kill criterion；
2. start 创建独立的 branch/worktree；
3. 仪器和合成自测一个 task，测量与已提交证据一个 task，判定/台账一个 task；
4. 每个 task 先提交声明的输入，再 formal verify 和 complete；
5. 独立 reviewer 只读复核最终内容；
6. 用户分别决定 archive 和 push/PR/merge，不把任何一步当作另一步的授权。

## 阶段 4：phase 收口与后续 Change

- 一次性研究脚本归档到 `phases/pN/` 前，同步退役或重指其 Check 定义，避免 registry 指向不存在的路径；
- OpenSpec 的 archived Change 不搬进 phases；phase README 链接其 Change ID 和科学证据；
- 文档漂移、SciPy/pandas 元数据、混合 schema CSV、INDEX 缺 D14 分别作为独立 Change；
- 本地完整生命周期稳定后再加 CI。通用 CI 只跑轻量、无外部数据的检查；重型数值门保留给有资源锁和数据保障的自托管或 operator milestone。若 CI 机 CPU 与证据机不同，按 §3.4，数值证据锁在 CI 上不能期望逐位通过，需要先有阶段 0 的结论；
- `newton.py`、`ibl3.py` 等大模块只在有稳定的行为 oracle 后做垂直切片重构。

# 九、采用验收条件

以下条件全部满足，才能称“devflow 已在 UP3D 可用”：

1. 当前 12 红逐项有已提交的分类和处置，不存在被隐藏的红色必需检查；
2. 干净 clone 和新 Change worktree 能按文档取得相同的网格、工具和输入；
3. raw 网格变化、环境变化（含有效 codegen：CPU 名、feature 串、numba/llvmlite/LLVM 版本、`NUMBA_*` 覆盖）、声明输入变化都会使旧 Run 失效；
4. 敏感证据门的 pass 只能在规范暖态下产生（方案 A：测量进程内零编译，否则 fail）；清缓存不使已自证暖态的旧 pass 失效，但之后的冷态运行不能产生 pass。这些由“通过后清缓存，再执行 status/archive”的测试和 G-TEETH 测试实证。换机由条件 3 覆盖；
5. 缺网格、零测试、全跳过、非预期 skip 和 strict XPASS 都不能形成绿色必需 Check，预期的图导出 skip 和 strict xfail 也不会让必需 Check 永远卡在 `skip`；
6. 一次真实研究 Change 完整走通 propose → start → task commits → verify → independent review → user approval → archive；
7. FAIL/RECORDED/UNDEFINED 的科学结论能按预注册完成交付，不需要换绿；
8. 产品变化后独立复核失效；archive 不授权 merge/push/release；
9. 一次性脚本归档后，所有 Check/sourcePaths 仍然有效；
10. adoption 自身没有改变求解器 API、默认值、参考数据、已提交科学结果或历史档案语义。

# 十、明确不做的事

- 不做全仓库规格反向转写，不重排目录，不改写历史；
- 不把 `runs/` 当科学证据提交，也不删除既有 committed CSV/PNG；
- 不通过放宽阈值、删除测试、改 optional、增加 skip 来制造绿色；
- 不在正式验证中自动刷新 committed gate 证据；
- 不借 adoption 或 baseline-provenance 夹带平台专用再基线；
- 不让 Claude 继续以“直接 main 开发”绕过 Change 的验证和复核；
- 不把独立 reviewer 仅当作另一个用户名；
- 不把一次性重型研究全部塞进全局 milestone；
- 不把 archive、merge、push、发布合并成一个授权动作。

# 十一、openspec-devflow 本体验证与限制

实测框架环境：Conda `openspec-devflow-dev`，Python 3.12.16、Node.js 26.12.1、OpenSpec 1.16.0、devflow 0.2.0；QuickJS 可用，`pip check` 正常。

| 项 | 结果 |
|---|---|
| 版本 | HEAD `cc0a9f1`（`Merge pull request #4 … release/v0.2.0`），tag `v0.2.0`，`devflow/__init__.py:4` 写 `__version__ = "0.2.0"` |
| 全量测试 | 409 collected / 409 passed / 0 failed / 0 skipped；两次独立运行 24.65 s 与 23.3 s |
| 框架自身 `install --check` | current |
| 最小模板 `templates/minimal-project/` | `install --check` 报 `"status": "outdated"`，需创建 6 个文件（`schema.yaml` + templates 下 5 个）；README 原句：“It does not contain the schema and adapters that `devflow install` adds” |
| `init` / `doctor` | 不存在（`cli.py:25-123`）；README.md:13 标 “Status: experimental”，并说明初始化与诊断尚未完成 |

因此采用时必须做真实的生命周期验证，不能只复制模板或照文档假定。干净消费项目的跨平台完整验收和 Windows CI 尚未完成；本报告只实测了 Linux。

# 十二、建议的下一步

1. **先做 `baseline-provenance`**：对 11 个门文件跑四次，全部保存全精度结果：native 冷 → native 暖 → generic 冷（独立空缓存）→ generic 暖。比较 native 冷对 generic 冷、native 暖对 generic 暖，并将 native 暖与已提交证据比较。暖跑单次约 17.5 min / 16 线程，冷跑另加编译时间，合计约 70 min 以上；跑之前先查机器负载。能联系旧机就加旧机臂。结论的范围限于 Numba JIT（§3.4）。结论连同有效 codegen 出处列补丁一起提交。
   ★ **已启动（2026-11-04 05:16，detached）**：runner `bench/results/baseline_provenance_20261104/run_arms.sh`，在固定于 `080f291` 的 detached worktree `/home/lrz/code/up3d-baseline-prov` 里运行（不受主 checkout 并发修改影响）。全精度值由只记录、不改判定的 pytest 插件 `prov_plugin.py` 落盘，四臂结束后 `compare.py` 自动生成 `compare_summary.md` / `compare_values.csv`。插件同时记录每臂 pyfp3d dispatcher 的 `cache_hits/misses`，作为裁决 A 前提 1 的试点读数（只记录，不作判定）。 ★ **结果见第 11 轮**：[20261104-0640-baseline-provenance-verdict.md](20261104-0640-baseline-provenance-verdict.md)（RECORDED：不是单一根因；证据已提交到 `bench/studies/baseline_provenance/`）。
2. 同时处理当前研究分支的交接（复核/合入或冻结），授权后对齐本地 `main`。 ★ **已完成（2026-11-04）**：PR #34 合入，本地 `main` = `origin/main` = `9b4f425`；合并前 F07 + F10 当前机 11 passed。
3. 然后创建 tooling-only 的 adoption Change。**第一步是提交 proposal/design/check 清单供复核**，其中写明 §7.3 方案 A（已裁决）的 Check 划分与 deselect 清单；不要直接安装后批量改仓库。 ★ **草案已起草（2026-11-04）**：`openspec/changes/adopt-devflow/`（proposal / design / tasks / change.json / `.openspec.yaml`）；devflow 尚未安装，草案待复核，在阻塞 1 处置前不 `start`。
4. 任何实际的 git 写操作、安装、提交、push、PR、archive，都按现有授权边界另行执行。

# 附录 A：第一次修订（相对初版）

标 †B 的行已被附录 B 部分推翻或收窄，原文保留以便追溯。

| # | 初版表述 | 修订 | 依据 |
|---|---|---|---|
| 1 †B | 12 红待做“冷/暖缓存和测试顺序单变量复核” | 测试顺序、运行间散布、网格已排除；首要候选改为跨机器的 CPU/ISA 代码生成（fastmath + 宿主 CPU 编译）；增加 ISA 臂与旧机臂 | 独立定向重跑逐位复现；网格跟踪状态；`pyfp3d/` 的 33 处 `fastmath` |
| 2 | 包版本与第 6、9 轮 `environment.json` 一致 ⇒ 不能简单归因于包版本 | 那两份 json 记录的是迁移机本身，不代表旧机；旧机 Numba/CPU 未被记录 | json 的 platform 字段；p6 README:75 |
| 3 | D13 非对称度断言“是否同一出处状态的另一表现，尚未闭合” | 保留未闭合，但补充：它是缺陷记录断言，读数超出该量自己的散布带，不属于已标定散布 | `test_D13_ibl_vs_xfoil.py:106,262-283,329` |
| 4 | 网格专项 Checks 声明 data；wrapper 校验预期 skip/xfail 清单 | 升为前置阻塞 4：xfail 在 JUnit 中即 skipped，与 data 规则及“不开图导出”冲突；给出方案 A/B | pytest JUnit 最小实验；`verification.py:577`；`execution.py:1113` |
| 5 | Check 分层含 `manual/operator` 一层 | devflow 只有 4 个 profile；manual 是 Check 属性 | `verification.py:25`；`record-manual` |
| 6 | “UP3D 的大型 `.msh` 正处于这个边界” | 只涉及 20 个被忽略网格（M6 族为主）；19 个被跟踪网格不受影响；12 红的输入全在被跟踪或代码生成一侧 | `git ls-files` / `git check-ignore` |
| 7 | `cases/gates/`：……24 份网格 manifest | manifest 位于 `cases/meshes/`，覆盖全部被忽略网格，SHA 全部一致 | 逐份 SHA-256 比对 |
| 8 | 声明 data 的 Check 有 skip “不会算作通过” | 准确状态是 `skip`（不是 `fail`），archive 对两者同样拦截 | `verification.py:577-579` |
| 9 | OpenSpec CLI 固定 1.16.0 | 补充：只有模板精确钉死，运行时不校验版本，故 lockfile 是硬要求 | `archive.py`；模板 `package.json:8` |
| 10 †B | 外部 `UP3D_DATA_DIR` + raw manifest | 补充强制布局 `<dataRoot>/<case_id>/raw/<localName>`、全量重哈希成本，以及框架两处文档的张力 | `fingerprint.py:103,207-244`；README.md:121 vs `openspec/README.md:44` |
| 11 †B | `environmentPackages` 列 Python/NumPy/…/gmsh | 增加 CPU 型号与 `numba.config.CPU_NAME` | §3.4 |
| 12 | 阶段 1 未引用 brownfield 指南 | 补充 `devflow/REFACTORING.md` §2 的顺序要求 | 同左 |

# 附录 B：第二次修订（第二轮外部审计，逐条复测后并入）

第二轮审计提出五项意见。按 agent-rules 轮次制度第 6 条，五项都由本次修订者独立复测，没有按面值收下。五项全部成立，另外复测中自查出两处问题（B6、B7）。

| # | 第一次修订的表述 | 第二次修订 | 复测依据 |
|---|---|---|---|
| B1（高） | “清缓存或换机不能复用旧通过”（§5 P0-2、§8 阶段 2、§9 条 4） | 只有目标没有机制：非 milestone 默认复用；复用不执行 wrapper；清缓存不改指纹也不改包版本；archive 不重跑。改为区分换 codegen（伪包，框架原生）与清缓存（方案 A 自证暖态 / B 缓存态伪包 / C 框架扩展；**使用者 2026-11-04 裁决采用 A**），并新增“通过后清缓存，再执行 status/archive”的测试 | `verification.py:153-187,771,806`；`fingerprint.py:301-303`；`execution.py:1098-1113`；`archive.py:73-77`；`numba/core/dispatcher.py:167,952` |
| B2（中） | “被忽略的大文件只能走 Case raw”；raw“必须”放在 checkout 之外 | 事实错误。`extraInputs` 直接读工作树算 SHA，被忽略文件不报 dirty；被忽略的内部 raw 同样不报 dirty。外部 raw 降为推荐的设计选择。另加：`extraInputs` 在缺文件时得 `None` 不报错，且未文档化，所以不推荐用于网格 | 最小实验三臂 + 对照臂（`bench/studies/devflow_adoption_audit/exp_extra_inputs.py`）；`fingerprint.py:233,244,246-250,396,405-420` |
| B3（中） | ISA 臂 = generic + 独立缓存，与现有结果比 | generic 首跑也是冷编译，不能与 native 暖态比较后归因于 ISA。改为事前写死冷/冷、暖/暖两对；结论范围限于 Numba JIT，不含 OpenBLAS 等预编译库 | 纪律 11；p6 第 44 轮勘误；`numba/core/config.py:306-308`；`codegen.py:1298-1303` |
| B4（中） | “差值逐位相同”“运行间散布已排除” | 证据不足：全套日志未保留；断言只打印 `.9e`，且每个测试只报第一个越界列；两跑都是暖态。改为“两个已测暖态运行中，失败节点一致、显示精度差值一致”，测试顺序与运行间散布降为弱排除 | `tests/_gate_evidence.py:145-148`；全套日志未保留；定向日志已提交为 `bench/studies/devflow_adoption_audit/targeted_12red_20261104_console.txt` |
| B5（中） | 出处列 / `environmentPackages` 写“宿主 CPU 型号 / `numba.config.CPU_NAME`” | 不够：实测 `CPU_NAME = None`、`CPU_FEATURES = None`，实际用的是 LLVM 检测到的宿主 CPU 名与 feature 串。改为记录有效 codegen 字段（含 `NUMBA_ENABLE_AVX`、`NUMBA_OPT`、SVML、BLAS coretype） | `numba/core/codegen.py:1322-1330,1427-1445`；当前机实测 `icelake-server`、114 项 feature |
| B6（自查） | README.md:121 与 `openspec/README.md:44` “有张力” | 两处文档一致：声明被忽略文件即判 uncommitted | 重读 README.md:121 全段 |
| B7（自查） | §3.4 冷/暖缓存“未测” | 补充：暖态的定向重跑复现了 12 红，所以它们不能归为 p6 第 44 轮记录的“冷编译扰动”一类 | `phases/p6/docs/dev_phase_six/20261011-1800-process-history-determinism.md:76-87` |
| B8（自查，裁决后） | §7.3 方案 A：strict xfail 门“单独成一个 Check” | 只含 xfail 的 Check 会被判 `skip`（全部 skip 规则不看 `data`）。改为随整个门文件成 Check，wrapper 只以非零退出码判红 | `evidence.py:199-216`；`verification.py:562-579`；pytest JUnit 实测 `tests="3" skipped="2"`（1 pass / 1 skip / 1 xfail） |

**使用者裁决（2026-11-04）**：§5 P0-2 缓存新鲜度采用方案 A（Run 自证暖态）；§7.3 skip/xfail 采用方案 A（在 Check 定义里排除预期 skip，xfail 随门文件成 Check）。

第二次修订**没有**新增求解运行，没有改动任何库、测试或已提交文件；新做的实测只有 devflow 最小实验（`/tmp`）与 Numba 配置读取。

未独立复核、沿用初版的实测：44 分钟全套中 12 红以外的 827 passed、wheel 56 文件与元数据、Markdown 链接扫描、17 个 CSV 结构告警、`git diff` 空白检查；以及第一次修订称“两跑均为暖态”中**全套那一跑**的缓存态（没有日志可查）。
