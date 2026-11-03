# UP3D 迁移到 openspec-devflow：清场重启方案

日期：2026-11-04  
状态：待执行  
唯一执行顺序：`bootstrap-devflow` → `restart-solver-core` → 后续单功能 Change

## 1. 决策

立即停止旧的 phase 7 开发方式，UP3D 从现在起改用 openspec-devflow 管理全部新开发。

这次迁移的目标不是修复旧求解器、保持旧测试全绿或保留旧 API，而是尽快在新框架下重新做出可验证的产品。旧代码可以删除，旧科学结果只作为历史资料；需要重新开发的功能，通过新的 OpenSpec spec 定义。

两项旧审计仍然正确定位了旧体系的问题：数值答案受平台、JIT 和缓存出处影响，旧测试、证据、文档与流程彼此缠绕。但它们回答的是“如何无损保留旧系统”，不是本次“如何重新开始”。因此：

- 不再处置或再基线旧 12 红，不跑四臂出处实验，不修复旧 bug 后才允许切换。
- 不迁移旧 phase、轮次、progress 台账、门编号、测试账目或证据体系。
- 不反向为旧仓库补规格，不追求新旧结果一致。
- 不把旧测试先登记为 devflow 必需 Check 再豁免；它们随旧实现退出活跃产品树。

这是对重构指南中“保持旧行为”的明确例外：本次是用户授权的**产品重启**，不是行为保持型重构。新框架仍须如实验证新实现，不能用删除测试、skip 或 optional 掩盖新代码的失败。

## 2. 目标状态

迁移完成后，仓库只保留一套当前事实来源：

| 内容 | 唯一事实来源 |
|---|---|
| 当前产品行为 | `openspec/specs/` |
| 正在做什么、为什么做 | 活跃 Change 的 `proposal.md`、`design.md`、spec delta |
| 任务与完成状态 | 活跃 Change 的 `tasks.md`，由 devflow 更新 |
| 分支、worktree、验证、复核、批准、归档 | devflow |
| 产品代码与新测试 | `pyfp3d/`、`tests/` |
| 历史实现与旧研究 | Git 历史中的只读归档，不参与当前开发 |

日常开发只有一条规则：**一个 Change 做完当前任务就验证，验证通过后提交并完成该任务；不再写另一份进度台账。**

push、PR、merge、archive 仍需用户分别明确批准；本地任务提交按 Change 授权执行。

## 3. 三步迁移

### 第一步：`bootstrap-devflow` —— 最小接入，不碰旧产品

这是唯一允许手工引导的 Change，因为框架尚不存在，不能假装它已经能管理自己的安装。

只做以下事项：

1. 固定工具版本：devflow `v0.2.0` / `cc0a9f1`，OpenSpec CLI `1.16.0`；devflow 使用独立的 Python 3.12 环境，不改变旧求解环境。
2. 运行 `devflow install --host claude --host kimi`，提交 schema、agent 适配器和安装记录。
3. 新建最小 `openspec/config.yaml`、`openspec/project.json`、`package.json` 和 lockfile；`integrationBranch` 为 `main`。
4. 只注册一个自动 Check：`framework-smoke`，检查安装完整性、配置可读取以及一个无求解依赖的 `python -m pytest` 冒烟测试。
5. 将 `AGENTS.md` 重写为简短的 devflow 协作规则，`CLAUDE.md` 只保留入口指针；删除当前 `openspec/changes/adopt-devflow/` 旧草案，以本方案为准重新建立最小 Change。
6. `.gitignore` 增加 `runs/`、`node_modules/`；运行证据不提交。

明确不做：旧测试接入、Numba/codegen provider、暖缓存 wrapper、外部网格 Case、skip/xfail 适配、资源锁、CI、roadmap  HTML、历史规格迁移。需要时由后续 Change 单独增加。

完成条件：

- `devflow install --check` 为 `current`；`devflow list-checks` 和 `devflow list-cases` 正常。
- `framework-smoke` 真实通过，OpenSpec strict validation 通过。
- 独立复核后，由用户批准归档和 PR/merge；主分支包含可用的最小框架。

### 第二步：`restart-solver-core` —— 首个真实 devflow Change

这是迁移的关键验收：从主分支用 `devflow start --change restart-solver-core --create` 创建独立 branch/worktree，完整走新生命周期。

在同一个行为变更 Change 中：

1. 删除旧 `pyfp3d/`、`tests/`、`bench/`、`phases/` 和旧 phase 工作文档；旧内容已由 Git 历史保存，不复制成第二套 `legacy/` 活跃目录。`cases/reference_data/` 和原始参考文献可保留为只读参考资料，但不再自动构成验收基线。
2. 建立全新的最小包与测试布局，重新编写当前 README 和依赖清单。
3. 只实现一个最小纵向闭环：读取一个极小四面体网格，施加均匀无黏来流，输出速度场与自由流保持误差。暂不求解跨声速、升力、Kutta、IBL、AMG 或复杂边界条件。
4. 为该行为写新的 spec delta 和新测试；至少包括一个 API 测试和一个从真实入口读取输入、产生输出的端到端测试。
5. 每个任务依次执行：实现 → 项目快速检查 → 提交 → `devflow verify --task` → `devflow complete`。
6. 最后运行 affected verification，由未参与实现的独立 agent 复核，用户批准后再 archive；push、PR、merge 另行批准。

第一步只负责让框架可用，第二步才证明新开发方式可用。不要在第一步顺手删除旧代码，也不要把第二步扩成完整 CFD 重建计划。

### 第三步：后续开发 —— 一次只增加一个可观察能力

`restart-solver-core` 归档并合入后，再确定下一个 Change。每个 Change 必须能回答：

1. 用户或调用方能观察到什么新能力？
2. 哪条 spec 和哪个真实入口测试证明它？
3. 能否在一个独立 worktree 内小步完成并回滚？

建议顺序是：最小势流离散 → 亚声速边界条件与求解 → 一个外部校验算例 → 单个更复杂能力。这里只是依赖顺序，不是预先承诺的大型路线图；每次只规划下一个 Change。

## 4. 新框架下的固定工作方式

每个新 Change 都遵循：

```text
propose
  → devflow start --create
  → devflow next / run --task
  → 实现一个任务
  → 提交
  → devflow verify --task
  → devflow complete
  → affected verify
  → independent review
  → 用户批准 archive
  → 用户另行批准 push / PR / merge
```

约束只有六条：

1. 不在 `main` 直接开发；一个 Change 一个 branch/worktree。
2. 一个任务只做一件事，完成当前任务后再选下一个。
3. 新代码的必需 Check 必须真实通过；失败可以记录，但不能伪造成通过。
4. 科学结论可以是 FAIL、RECORDED 或 UNDEFINED；这不等于 devflow 任务失败，交付约定证据即可完成任务。
5. 独立复核之后再有产品改动，必须重新验证和复核。
6. 不恢复旧 phase 文书体系，不新增与 `tasks.md` 并行的进度状态。

## 5. 回滚与迁移完成判据

回滚很简单：第一步之前继续使用旧 `main`；第一步之后、第二步之前只回滚框架接入提交；第二步之后旧实现仍可从迁移前 commit 建立只读 worktree。任何恢复旧代码为当前产品的行为都必须另开行为变更 Change，不能直接把旧目录复制回来。

只有同时满足以下条件，才称迁移完成：

- 最小框架已合入主分支，干净 clone 能安装并运行 `framework-smoke`。
- `restart-solver-core` 已在独立 worktree 中完整走通 propose、任务提交、verify、独立复核、用户批准和 archive。
- 当前产品只由新 specs、新代码和新测试描述；旧红、旧证据与旧流程不再阻塞新开发。
- 下一个开发任务能够从 `devflow next` 开始，而不需要阅读历史会话或旧 phase 台账。

迁移完成不表示求解器功能完整或没有 bug，只表示**UP3D 已经能够在 openspec-devflow 下持续开发**。
