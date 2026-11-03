# UP3D 迁移到 openspec-devflow：清场重启方案

日期：2026-11-04
状态：待执行
执行顺序：`bootstrap-devflow` → `restart-solver-core` → 后续单能力 Change

## 1. 决策与边界

立即停止旧的 phase 7 开发方式，使用 openspec-devflow 重新开发。本次是用户明确授权的产品重启，不是行为保持型重构：

- 允许删除旧代码、旧测试和旧流程文件，不保留第二套 `legacy/` 活跃目录。
- 不先修旧 bug，不追求旧测试转绿，不跑旧基线、新旧一致性实验或历史性能比较。
- 不迁移旧 phase、gate、任务台账、测试账目和历史证据，不反向补写旧系统规格。
- 新实现只接受新 specs、活跃 Change 和新 Checks 的约束；不得用 skip、xfail 或 optional 掩盖必需检查失败。

迁移后的开发权威只有新 `AGENTS.md`、`openspec/specs/`、活跃 Change artifacts 和 `devflow protocol execution`；`CLAUDE.md` 只作为指向 `AGENTS.md` 的入口。旧文档、旧测试和旧证据不产生当前任务、Check 或验收条件。本方案只指导下列三步，不成为另一套流程或进度台账。

## 2. 统一的 devflow 生命周期

### 2.1 普通 Change 的进入方式

先在主 checkout 同步干净 `main`，再 propose。`start --create` 前必须满足：本地 `HEAD` 与 `origin/main` 一致；除新 Change 目录 `openspec/changes/<change-id>/` 的未跟踪文件外没有其他改动；`openspec/project.json` 使用 `"integrationBranch": "main"` 和 `"integrationRemote": "origin"`。不得设为 `null` 绕过远程同步检查，也不得清理他人工作。

```bash
devflow start --change <change-id> --create
cd <start 报告的 worktree>
devflow next --change <change-id>
```

Change 必须声明独立 branch/worktree；Kimi 使用 `kimi/<topic>`，不在 `main` 直接开发或提交。第一步只用第 3 节明确的一次性初始化替代上述进入方式，之后的任务和验收没有豁免。

### 2.2 每个任务和整个 Change 的顺序

字面顺序固定为：`next` 先于实现，提交先于正式验证：

```text
devflow next --change <id> / devflow run --change <id> --task <task-id>
→ 实现当前任务并做快速检查
→ 提交当前任务及全部声明输入
→ devflow verify --change <id> --task <task-id>
→ devflow complete --change <id> --task <task-id>
→ 提交 tasks.md、执行记录和本次 Run 元数据形成的状态变更
→ 下一任务，直至全部完成
→ devflow verify --change <id> --profile affected
→ 独立复核，修正并重新验证，关闭全部发现
→ 用 devflow review 记录当前内容的批准结论
→ 用户明确批准 archive
→ devflow archive --change <id> --approval "<用户批准>"
→ 提交 spec 同步（如有）、acceptance 和归档结果
→ 分别取得用户对 push、PR、merge 的明确授权
```

约束如下：

- 快速检查不能代替正式 `verify --task`。Check 的 `sourcePaths`、环境文件和 Case 输入必须真实存在并被 Git 跟踪；未提交、未跟踪或被忽略的声明输入不能作为有效 Run 的依据。
- 每个任务至少声明一个必需 Check。失败或过期时修正当前内容，不借助旧基线豁免。`complete` 更新任务复选框和执行记录后必须再次提交，不能手工勾选。
- 根目录 `/runs/` 中的本地证据、worktree 和资源状态不提交；Change 内需要持久化的任务状态、执行记录、评审和 Run 元数据按生命周期提交。
- 独立复核必须由未参与实现的其他上下文以只读方式完成。复核后的产品改动会使复核失效，必须重新验证和复核。
- archive 不自动提交，也不等于 push、PR 或 merge；后三者分别授权。

独立复核发现全部关闭后、请求用户批准前，必须显式记录批准，不能只写评审文字：

```bash
devflow review --change <id> --reviewer <独立复核者> --approved
```

有 findings 时，在同一条最终复核命令中增加 `--findings <findings.json>`，逐项保留原 finding `id`，并如实记录 `"status": "closed"`；确由产品改动取代的项才使用 `"superseded"`。不得省略仍开放的发现。最终复核必须针对当前产品内容且没有开放发现，否则不得请求 archive 批准。

## 3. 三步迁移

### 第一步：`bootstrap-devflow` —— 最小工具接入

这一步只让框架可用，不删除旧产品代码，不接入旧测试，也不实现求解器功能。由于没有 `devflow init`，且 `next` 需要先加载项目配置和 Check，严格采用以下**一次性初始化顺序**：

1. **同步干净 `main`：**在主 checkout 执行下列命令，确认 `HEAD` 与 `origin/main` 一致且 `git status --short` 为空。

   ```bash
   git switch main
   git fetch origin
   git merge --ff-only origin/main
   git status --short
   ```

2. **手工建立独立 branch/worktree：**本 Change 不运行 `devflow start`；其他 Change 均使用 `start --create`。

   ```bash
   git worktree add -b kimi/bootstrap-devflow runs/worktrees/bootstrap-devflow HEAD
   cd runs/worktrees/bootstrap-devflow
   ```

3. **安装固定工具并安装项目适配器：**从已固定的 devflow `v0.2.0` / `cc0a9f12704f` 安装到独立 Python 3.12 环境；用项目 `package.json` 和 lockfile 固定 OpenSpec CLI `1.16.0`，然后运行：

   ```bash
   devflow install --host claude --host kimi
   ```

4. **写最小项目配置和 `framework-smoke`：**新建 `openspec/config.yaml`（选择 `devflow` schema）和最小 `openspec/project.json`。核心字段如下，不照抄模板中不存在或当前被忽略的 `src`、`tools`：

   ```json
   {
     "formatVersion": 1,
     "integrationBranch": "main",
     "integrationRemote": "origin",
     "evidenceRoot": "runs/evidence",
     "resourceStateRoot": "runs/devflow",
     "productPaths": ["pyfp3d"],
     "taskLayerPaths": ["pyfp3d"],
     "environmentPackages": [],
     "resources": {}
   }
   ```

   `pyfp3d` 当前真实存在，第二步清场后也必须在同一路径重建。省略 `packageVersionProvider`，不引入任何环境 provider。

   新建并 Git 跟踪 `openspec/framework_smoke.py` 和 `openspec/checks/framework-smoke.json`。脚本本身只导入 Python 标准库，通过固定 devflow CLI 检查安装完整性、配置可读取和 registry 可加载；不导入 `pyfp3d`，不读取 `tests/`。Check 使用 `python openspec/framework_smoke.py`，`required: true`，只声明 `console` evidence，不使用 pytest 或 JUnit。`sourcePaths` 逐项覆盖脚本、`openspec/checks`、项目配置、安装记录、schema 和两类 agent adapters；所有声明路径必须真实存在并被 Git 跟踪。

5. **写 tooling-only Change artifacts：**按已安装 schema 建立 `openspec/changes/bootstrap-devflow/`。`.openspec.yaml` 设置 `skip_specs: true`；proposal 的 New/Modified Capabilities 均为空；`change.json` 的 `impact.capabilities`、`impact.cases` 为空，`impact.checks` 只含 `framework-smoke`，并声明真实平台、匹配的 branch/worktree/machine 和 `independentReview: true`；`tasks.md` 至少有一个任务。不创建 spec delta，不留模板占位值。

6. **首次选择任务：**上述最小配置、Check 和 Change artifacts 写好后再运行：

   ```bash
   devflow next --change bootstrap-devflow
   ```

进入首个任务后，补齐本 Change 的简短 `AGENTS.md`、只含入口指针的 `CLAUDE.md`，删除 `openspec/changes/adopt-devflow/` 旧草案，并让 `.gitignore` 忽略根目录 `/runs/` 和 `node_modules/`。首个任务在快速检查后、正式 `verify --task` 前，仍必须提交当前任务及**全部声明输入**，包括安装生成的 schema/adapters、安装记录、最小配置、冒烟脚本、Check 和 Change artifacts；随后严格按第 2.2 节执行。初始化写在 `next` 之前不构成提交、验证或完成门槛的豁免。

第一步不做旧测试接入、旧 bug 修复、Numba/codegen provider、缓存包装、旧网格或外部网格 Case 适配、skip/xfail 适配、CI、roadmap HTML 和历史规格迁移。完成以 `devflow install --check` 为 `current`、`list-checks`/`list-cases` 正常、`framework-smoke` 真实通过、OpenSpec strict validation 和 affected verification 通过为前提；再记录有效独立复核，取得用户 archive 批准，提交归档结果，并经分别授权的 push、PR、merge 合入 `main`。

### 第二步：`restart-solver-core` —— 清场并重建最小闭环

第一步归档并合入后，按第 2 节从同步干净 `main` 开始，propose 新行为、写新的 spec delta 和任务，再执行：

```bash
devflow start --change restart-solver-core --create
cd <start 报告的 worktree>
devflow next --change restart-solver-core
```

该 Change 在一个独立 worktree 中完成：

1. 删除旧 `pyfp3d/`、`tests/`、`bench/`、`phases/`，在新的最小 `pyfp3d/` 包中重建产品代码和新的测试布局，重写 README 与依赖清单。`productPaths` 保持 `["pyfp3d"]`。
2. 删除旧 `cases/gates/`、`cases/demo/`、旧网格生成与 case 运行脚本及其余旧 case 配置。只保留只读资料 `cases/reference_data/` 和原始论文，不复制旧 case。
3. 删除旧 phase/gate 工作文档、旧审计报告和旧流程规则。保留第一步的 `openspec/framework_smoke.py` 与 `framework-smoke` Check；它们不依赖旧 `tests/`，因此清场后仍必须真实通过。
4. 只实现一个最小纵向闭环：从一个新的极小四面体网格入口读取输入，施加均匀无黏来流，从真实入口输出速度场和自由流保持误差。新 spec 明确输出契约和误差判据。
5. **此时才为新 API/E2E 产品测试选择并固定测试框架**，添加新测试及对应必需 Checks。测试只证明新 spec，不引用旧期望值、旧测试账目或历史性能基线；bootstrap 的 `framework-smoke` 不因此改用该框架。
6. 全部任务、验证、复核记录、用户批准、archive 和归档后提交均按第 2 节执行。

`cases/reference_data/` 和原始论文不自动进入 Check 输入，不构成数值真值或验收门槛。本 Change 不做 Numba、旧网格或外部网格、缓存适配，也不追求跨声速、升力、Kutta、IBL、AMG、复杂边界条件或完整 CFD 功能。完成只以新最小闭环的真实 API/E2E 验证和完整 devflow 生命周期为准；归档提交后，再分别授权 push、PR、merge。

### 第三步：后续开发 —— 一次增加一个可观察能力

`restart-solver-core` 归档、提交并合入后，只规划下一个 Change。每个 Change 用一条新或修改后的 spec、一个真实入口测试和独立 worktree 说明可观察能力及回滚边界，不批量启动 Change，也不预设大型路线图。

后续可以从最小势流离散、亚声速边界条件与求解、一个新定义的外部校验算例中选择一个开始。每次都按 **同步干净 main → propose → start --create → next** 进入，再重复第 2.2 节；不恢复旧基线，不要求阅读历史会话、旧 phase 台账或旧审计。

## 4. 回滚与迁移完成判据

回滚不恢复旧流程：第一步之前保留旧 `main`；第一步之后、第二步之前只回滚框架接入提交；第二步之后如需读取旧实现，只能从迁移前 commit 建只读 worktree。把旧代码恢复为当前产品必须另开带 spec delta 的行为 Change，不直接复制旧目录。

只有同时满足以下条件，才称迁移完成：

- 最小框架已合入 `main`，干净 clone 能安装并运行不依赖 pytest/JUnit 的 `framework-smoke`。
- `restart-solver-core` 已在独立 worktree 中走通完整生命周期并经授权合入，新最小闭环的 API/E2E 验证真实通过。
- 当前产品只由新 specs、新代码和新测试描述；旧红、旧文档、旧 case 和历史证据不再阻塞新开发。
- 下一个 Change 能从 **同步干净 main → propose → start --create → next** 开始，无需旧基线或复杂适配。

迁移完成不表示求解器功能完整或没有 bug，只表示 UP3D 已能在 openspec-devflow 下持续进行可验证的新开发。
