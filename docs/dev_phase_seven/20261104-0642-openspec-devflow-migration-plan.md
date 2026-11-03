# UP3D 迁移到 openspec-devflow：清场重启方案

日期：2026-11-04
状态：待执行
唯一执行顺序：`bootstrap-devflow` → `restart-solver-core` → 后续单能力 Change

## 1. 决策与边界

立即停止旧的 phase 7 开发方式。UP3D 从现在起使用 openspec-devflow 重新开发，不以修复旧求解器、保持旧 API、让旧测试转绿或复现旧结果作为前提。

本次是用户明确授权的产品重启，不是行为保持型重构：

- 允许删除旧代码、旧测试和旧流程文件，不建立第二套 `legacy/` 活跃目录。
- 不处置旧红，不修复旧 bug 后再切换，不跑旧基线或新旧一致性实验。
- 不迁移旧 phase、gate、任务台账、测试账目和历史证据，不反向补写旧系统规格。
- 新实现只接受新 specs、活跃 Change 和新 Checks 的约束；不得用 skip、xfail 或 optional 掩盖新代码的必需检查失败。

迁移后，当前开发权威只有：

1. 新的 `AGENTS.md`；`CLAUDE.md` 只作为指向它的入口。
2. `openspec/specs/` 中的当前产品行为。
3. 活跃 Change 的提案、设计、spec delta、`change.json` 和 `tasks.md`。
4. `devflow protocol execution` 发布的执行协议。

旧文档、旧测试和旧证据不产生当前任务、Check 或验收条件。本方案只指导下列三步迁移，也不成为后续开发的并行流程或进度台账。

## 2. 统一的 devflow 生命周期

### 2.1 每次创建 Change worktree 之前

除第一步唯一豁免 `devflow start` 外，每个 Change 都先完成 propose，再在主 checkout 执行：

```bash
git switch main
git fetch origin
git merge --ff-only origin/main
git status --short
devflow start --change <change-id> --create
```

执行 `start --create` 前必须确认：

- 本地 `HEAD` 与 `origin/main` 一致，不只是本地 `main` 可自行 fast-forward。
- 除新 Change 目录 `openspec/changes/<change-id>/` 的未跟踪文件外，主 checkout 没有其他修改或未跟踪文件；不得为通过检查而清理他人工作。
- `openspec/project.json` 使用 `"integrationBranch": "main"` 和 `"integrationRemote": "origin"`。不得设为 `null` 绕过远程同步检查。
- Change 声明独立 branch/worktree；Kimi 的分支使用 `kimi/<topic>`，不在 `main` 直接开发或提交。

`devflow start` 不代替上述 fetch 和 fast-forward 对齐。进入其报告的 worktree 后，先运行 `devflow resume` 或 `devflow next`，再选择任务。

### 2.2 每个任务和整个 Change 的顺序

所有 Change，包括 bootstrap，都遵守同一顺序：

```text
propose
  → start --create（仅 bootstrap 改为手工建独立 branch/worktree）
  → devflow next / devflow run --task
  → 实现当前任务并做项目快速检查
  → 提交当前任务及全部声明输入
  → devflow verify --change <id> --task <task-id>
  → devflow complete --change <id> --task <task-id>
  → 提交 tasks.md、执行记录和本次 Run 元数据所形成的状态变更
  → 下一任务，直至全部完成
  → devflow verify --change <id> --profile affected
  → 独立复核并关闭全部发现
  → 用户明确批准 archive
  → devflow archive --change <id> --approval "<用户批准>"
  → 提交 spec 同步（如有）、acceptance 和归档结果
  → 分别取得用户对 push、PR、merge 的明确授权
```

具体约束：

- 快速检查只用于及时反馈，不能代替正式 `verify --task`。正式 Run 所覆盖的声明输入必须先提交。
- `devflow complete` 会更新任务复选框和执行记录，因此完成后还要提交任务状态；不能手工勾选复选框。
- 根目录 `/runs/` 中的本地证据、worktree 和资源状态不提交；Change 内需要持久化的任务状态、执行记录和 Run 元数据按 devflow 生命周期提交。
- 每个任务至少声明一个必需 Check。affected verification 的必需 Run 必须对当前内容有效；失败或过期时先修正当前任务，不借助旧基线豁免。
- 独立复核必须由未参与实现的其他上下文以只读方式完成。复核后的产品改动会使既有复核失效，必须重新验证和复核。
- archive 不自动执行 `git commit`，也不等于 push、PR 或 merge。归档产生的 specs 变更、`acceptance.json`、归档目录和执行记录必须检查后再提交。
- push、开 PR、merge 是三件事，分别授权；任何一项都不由任务完成、archive 或本方案自动授权。

## 3. 三步迁移

### 第一步：`bootstrap-devflow` —— 最小工具接入

这一步只让框架可用，不删除旧产品代码，不接入旧测试，也不实现求解器功能。

**唯一流程例外：**不运行 `devflow start --create`。先按第 2.1 节同步并确认 `main` 干净，再从该提交手工创建独立 branch/worktree：

```bash
git worktree add -b kimi/bootstrap-devflow runs/worktrees/bootstrap-devflow HEAD
```

随后只在该 worktree 工作。按 devflow schema 建立 Change，并如实填写 `change.json`：`id`、上述分支、实际 worktree 相对路径和所属机器必须一致；`impact.capabilities` 为空，`impact.checks` 至少包含 `framework-smoke`，并声明真实平台和 `independentReview: true`。不得保留模板占位值或虚构不存在的分配。

它是 tooling-only Change，因此必须同时满足：

- `.openspec.yaml` 设置 `skip_specs: true`。
- `proposal.md` 的 New Capabilities 和 Modified Capabilities 都为空，`change.json` 的 `impact.capabilities` 也为空。
- 不创建任何 spec delta，不为通过校验而虚构产品能力。

除不用 `devflow start` 外，没有第二个流程豁免：任务提交、`verify --task`、`complete`、任务状态提交、affected verification、独立复核、用户批准、archive 和归档后提交全部照常执行。由于当前工具没有 `devflow init`，按最小项目模板手工编写初始配置是安装方式，不是跳过上述门槛的理由。

工作范围固定为：

1. 固定 devflow `v0.2.0` / `cc0a9f1`、OpenSpec CLI `1.16.0`；devflow 使用独立 Python 3.12 环境，OpenSpec 由项目 `package.json` 和 lockfile 固定。
2. 运行 `devflow install --host claude --host kimi`，提交 schema、agent 适配器和安装记录；新建最小 `openspec/config.yaml` 与 `openspec/project.json`。
3. 只注册一个自动 Check：`framework-smoke`。它检查安装完整性、配置可读取以及一个不依赖求解器的新 pytest 冒烟测试，并如实声明全部输入路径。
4. 将 `AGENTS.md` 重写为简短的新 devflow 协作规则，`CLAUDE.md` 只保留入口指针；删除 `openspec/changes/adopt-devflow/` 旧草案，以本 Change 为唯一接入 Change。
5. `.gitignore` 忽略根目录 `/runs/` 和 `node_modules/`。不提交机器本地证据，不增加旧流程台账。
6. 每个任务严格执行第 2.2 节；全部任务结束后运行 affected verification、独立复核，取得用户 archive 批准，归档并提交归档结果。

明确不做：旧测试接入、旧 bug 修复、Numba/codegen provider、缓存包装、外部网格 Case、skip/xfail 适配、CI、roadmap HTML 和历史规格迁移。确有新需求时由后续 Change 单独提出。

第一步完成条件：

- `devflow install --check` 报告 `current`，`devflow list-checks` 和 `devflow list-cases` 正常。
- `framework-smoke` 真实通过，OpenSpec strict validation 通过，affected verification 和独立复核均有效。
- 用户已批准 archive，归档结果已提交；随后经分别授权的 push、PR、merge，使 `main` 包含最小可用框架。

### 第二步：`restart-solver-core` —— 清场并重建最小闭环

第一步归档并合入后，先通过 OpenSpec propose 新行为、写新的 spec delta 和任务，再按第 2.1 节同步 `main` 并执行：

```bash
devflow start --change restart-solver-core --create
```

该 Change 在一个独立 worktree 中完成以下工作：

1. 删除旧 `pyfp3d/`、`tests/`、`bench/`、`phases/`，再建立全新的最小 `pyfp3d/` 包和 `tests/` 布局；重写当前 README 和依赖清单。
2. 删除旧 `cases/gates/`、`cases/demo/`、旧网格生成与 case 运行脚本，以及其余旧 case 配置。只保留明确列为只读资料的 `cases/reference_data/` 和原始论文，不复制旧 case 到新框架。
3. 删除旧 phase/gate 工作文档、旧审计报告和旧流程规则；新开发只由新 `AGENTS.md`、specs、活跃 Change 和 devflow protocol 约束。
4. 只实现一个最小纵向闭环：从一个新的极小四面体网格入口读取输入，施加均匀无黏来流，从真实入口输出速度场和自由流保持误差。新 spec 明确输出契约和误差判据。
5. 添加新的 API 测试和端到端测试，注册对应的新必需 Checks。测试只证明新 spec，不引用旧期望值、旧测试账目或历史性能基线。
6. 每个任务均按“实现/快速检查 → 提交 → `verify --task` → `complete` → 提交任务状态”执行；最后完成 affected verification、独立复核、用户批准、archive 和归档后提交。

清场后的边界必须保持清楚：

- 旧文档、旧测试、旧 case 和旧证据即使仍能从 Git 历史读取，也不生成当前任务、Check 或验收条件。
- `cases/reference_data/` 和原始论文只读，不自动进入 Check 输入，不构成数值真值或验收门槛；不得为迁移重算或补交历史证据。
- 本 Change 不追求跨声速、升力、Kutta、IBL、AMG、复杂边界条件或完整 CFD 功能。

第二步完成以新最小闭环的真实 API/端到端验证和完整 devflow 生命周期为准，不以任何旧测试是否通过为准。归档结果提交后，再分别授权 push、PR、merge。

### 第三步：后续开发 —— 一次增加一个可观察能力

`restart-solver-core` 归档、提交并合入后，只规划下一个 Change。每个 Change 必须能用一条新或修改后的 spec、一个真实入口测试和独立 worktree 说明其可观察能力及回滚边界。

后续可以从最小势流离散、亚声速边界条件与求解、一个新定义的外部校验算例中选择一个开始，但不预设大型路线图，不批量启动 Change。每个都从 **propose → 同步并检查 main → start --create → next** 开始，再重复第 2.2 节的任务与验收顺序。

## 4. 回滚与迁移完成判据

回滚不涉及恢复旧流程：第一步之前保留旧 `main`；第一步之后、第二步之前只回滚框架接入提交；第二步之后如需读取旧实现，只能从迁移前 commit 建只读 worktree。把旧代码恢复为当前产品必须另开带 spec delta 的行为 Change，不直接复制旧目录，也不要求先修复旧系统。

只有同时满足以下条件，才称迁移完成：

- 最小框架已合入 `main`，干净 clone 能安装并运行 `framework-smoke`。
- `restart-solver-core` 已在独立 worktree 中走通 propose、任务提交、逐任务 verify/complete、任务状态提交、affected verification、独立复核、用户批准、archive 和归档后提交，并经授权合入。
- 当前产品只由新 specs、新代码和新测试描述；旧红、旧文档、旧 case 和历史证据不再阻塞新开发。
- 下一个 Change 能从 **propose → start --create → next** 开始，无需阅读历史会话、旧 phase 台账或旧审计。

迁移完成不表示求解器功能完整或没有 bug，只表示 UP3D 已能在 openspec-devflow 下持续进行可验证的新开发。
