# Design

## Context

动机见 proposal.md；全部事实与源码定位见审计
[docs/dev_phase_seven/20261104-0524-devflow-adoption-audit.md](../../../docs/dev_phase_seven/20261104-0524-devflow-adoption-audit.md)（下称“审计”）。
这里只列决定设计形状的约束：

- devflow v0.2.0 判断 Run 能否复用，只比内容指纹、OS 名和 `packageVersionProvider` 返回的版本串；
  非 milestone 默认复用，**复用时不执行 Check 命令**（所以 wrapper 不运行）；archive 只读现有 Run 状态，不重跑（审计 §5 P0-2）。
- devflow 读 JUnit `<testsuite>` 的 `tests`/`skipped` 计数；pytest 把 xfail 计入 `skipped`；
  “全部 skip ⇒ `skip`”不看 `data`，“声明 `data` 且有 skip ⇒ `skip`”，archive 把 `skip` 与 `fail` 同样拦下（审计 §7.3）。
- 求解答案对 numba 编译出处敏感：p6 第 44 轮实测冷编译跑与暖缓存跑是两个确定但不同的数值状态；
  迁移机暖态下 12 条证据锁红，首要候选是跨机器的 codegen（审计 §3.4）。
- 20 个网格被 Git 忽略（M6 族为主），24 份已提交的 `.msh.manifest.json` 覆盖全部；19 个被跟踪网格由 git 定址。
- 求解环境固定 Python 3.11；devflow 要求 Python ≥ 3.12。

## Goals / Non-Goals

**Goals:**
- 一次旧的“通过”只有在它的出处（内容、环境、有效 codegen、输入网格、暖态）仍然成立时才被复用。
- 必需 Check 不能因缺网格、零测试、全跳过、意外 skip 或 strict XPASS 而变绿；预期的图导出 skip 与 strict xfail 也不能让它永远卡在 `skip`。
- 科学判定 FAIL/RECORDED/UNDEFINED 能按预注册完成交付，不需要换绿。

**Non-Goals:**
- 不为整个仓库反向编写 specs；不重排目录；不改写历史；不把旧轮次伪造成已完成的 Change。
- 不处置 12 红（那是前置的 baseline-provenance 工作）；不刷新任何 committed 证据。
- 不加 CI（本地完整生命周期稳定后另起 Change）。
- 不重构 `newton.py` / `ibl3.py` 等大模块。

## Decisions

**D1 叠加层，tooling-only。** 只新增 `openspec/`、`tools/devflow/`、锁文件与工具环境；`skip_specs: true`。
备选“全量反向补写 specs”被否：会把内部函数、历史轮次和未实现的愿望复制成规格（审计 §6）。

**D2 两套环境分离。** 工具环境 `up3d-devflow`（Python 3.12、Node ≥ 22.12、devflow 0.2.0、OpenSpec 1.16.0 经 lockfile）；
`devflowCommand` 指向它。每个求解 Check 显式经 `conda run --no-capture-output -n up3d` 进入求解环境。
备选“把 up3d 升到 3.12”被否：会改变求解环境，违反 adoption 不动求解器的边界。

**D3 被忽略网格走外部 Case raw。** `UP3D_DATA_DIR`（经 `dataEnvironmentVariable`）下 `<case_id>/raw/<localName>`；
wrapper 先以已提交的 `.msh.manifest.json` 校验 SHA-256，再物化到 worktree 的 `cases/meshes/`；校验失败以非零退出判 `fail`。
备选 `extraInputs` 指向被忽略文件被否：它读工作树、新 worktree 缺文件时得 `None` 而不报错，且未文档化（审计 §5 P0-2 实测 arm3）。
放在 checkout 之外是**设计选择**（跨 worktree 共享、不依赖 `.gitignore` 是否覆盖 raw 目录），不是框架硬约束。

**D4 换机 / 换 codegen ⇒ 环境指纹失效。** `environmentPackages` 列 Python（含构建）、NumPy、SciPy、Numba、llvmlite、PyAMG、
meshio、gmsh，外加伪包 `numba-codegen`：有效 CPU 名、有效 feature 串摘要、LLVM 版本、target triple、`NUMBA_CPU_NAME`/
`NUMBA_CPU_FEATURES`/`NUMBA_ENABLE_AVX`/`NUMBA_OPT` 覆盖、SVML、BLAS 库与 coretype。provider 在 devflow 的 3.12 环境里被导入，
所以必须以子进程查询 `up3d`。只写 `numba.config.CPU_NAME` 不够：默认是 `None`（审计 §3.4、附录 B5）。

**D5 缓存态 = 每个 Run 自证暖态（使用者裁决 2026-11-04，审计 §5 P0-2 方案 A）。**
敏感证据 Check：wrapper 先在独立进程里白跑一遍（纪律 11），测量进程结束时读取 pyfp3d 各 dispatcher 的
`stats.cache_misses`，任一未命中 ⇒ `fail`（不是 skip），计数写入 step trace。语义：一次 pass 陈述的是
“代码 + 有效 codegen ⇒ 答案，在规范暖态下测得”；清缓存**不**使它失效，但之后的冷态测量不能产生 pass。
**两个前提必须先实测再依赖**：(a) dispatcher 统计覆盖 pyfp3d 全部 JIT 入口（G-TEETH：故意让一个内核在测量进程内编译，必须红）；
(b) 同一 codegen 下连续两次暖跑全精度逐位一致。任一不成立 ⇒ 停下回使用者重新裁决，不得静默退到方案 B（缓存态伪包）或 C（框架扩展）。

**D6 skip/xfail = 在 Check 定义里不收集预期 skip（使用者裁决 2026-11-04，审计 §7.3 方案 A）。**
- 每个 Check 以 `--deselect` / `-m` 排除默认 skip 的图导出腿（正式 verify 不设 `PYFP3D_GATE_FIGURES=1`）；
- 含 strict xfail 的门文件（默认：C02、C03；门控：E01）各自以**整个文件**成一个不声明 `data` 的 Check ——
  **不得**把 xfail 节点单独成 Check（`tests == skipped` ⇒ 永远 `skip`）；同文件的通过项保证落进“部分 skip ⇒ pass”；
- wrapper **只往红里判**：对照该 Check 预期的 xfail 名单与 skip 名单（应为空），任何不符、任何 XPASS、零测试 ⇒ 非零退出 ⇒ `fail`；
  不改写 JUnit，不把 skip 变成 pass。缺网格由 D3 的前置校验判 `fail`。
备选方案 B（wrapper 产出过滤后的 JUnit）被否：devflow 看到的不再是 pytest 原始输出，复核者要读两份。

**D7 Check 分层。** 四个 profile（`task` / `affected` / `docs` / `milestone`）按审计 §7.1；manual 是 Check 属性不是层。
一次性贵重研究测量注册为 `manual: true` 的 Check，记录指向当前提交中的 CSV/PNG；不进全局 milestone。

**D8 科学判定与执行状态分开。** 预注册的 kill criterion 触发或结果 UNDEFINED 时，交付完整证据与判定即算 task 完成；
`tasks.md` 只管任务与 checkbox，`progress.md` 保留为科学结果索引，不成为第二套执行状态库。

**D9 协作规则单一权威。** `AGENTS.md` 写协作与状态转换；Claude 与 Kimi 都走 Change 分支（`claude/<change-id>`、`kimi/<change-id>`）；
push / PR / merge / archive 仍各需使用者确认。`CLAUDE.md` 只留指针，不复制第三份协议。

**D10 资源锁。** 初期 `resources: {"solver-cpu": {"slots": 1}}`，重型验证串行，防并行 Change 超订这台共享机（纪律 1）。

## Risks / Trade-offs

- [devflow 标注 experimental，无 `init`/`doctor`] → 固定 `cc0a9f1`；每一项保证都以本 Change 的实测任务证明，不照文档假定。
- [D5 前提 (a) 不成立：部分 JIT 入口经 `numba.njit(*args, **kwargs)` 包装生成，模块级扫描可能漏掉] → tasks 4.4 先做 G-TEETH；不成立即回使用者。
- [D5 前提 (b) 不成立：暖态本身不唯一] → 同上；p6 第 44 轮只测过两次“暖缓存复现出处态”。
- [D4 伪包让本平台全部 Run 一起失效] → 对换机恰当；日常不会变化。
- [D6 的 deselect 清单随测试演化而过期，新加的默认 skip 测试会让 Check 判红] → 这是想要的行为（红而不是静默），wrapper 报错信息指向要更新的清单。
- [Case raw 每次 Run 全量重算 SHA-256，M6 大网格有 I/O 成本] → 记录；只在声明了该 Case 的 Check 上发生。
- [迁移机证据仍红] → 本 Change 不 `start`，直到阻塞 1 处置。

## Migration Plan

1. 阻塞 1 处置完成后 `devflow start adopt-devflow`（独立分支 + worktree）。
2. 按 tasks.md 分组推进，每组自带验证。
3. 回滚：本 Change 只新增文件并改写 `AGENTS.md`/`CLAUDE.md`/`.gitignore`，回滚即 revert 其提交；求解器与证据不受影响。

## Open Questions

- CI：通用 CI 只跑轻量、无外部数据检查，留到本地生命周期稳定后另起 Change（不影响本 Change 的任务划分）。
