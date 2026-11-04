# Design

## Context

See `proposal.md` — Why. 当前 `devflow` 环境是 v0.2.0 的非 editable VCS 安装，仓库清单只记录 `engineVersion: 0.2.0`。v0.3.1 tag 已解析为 `1c7ab6f223d4301aea8ea1deeb4dde2b246de37a`；只读兼容性检查已用新版严格 `Registry.load` 成功读取现有项目、4 个 Checks 和空 Case 集合。

对上游 tag 的实际比较显示，`devflow/assets/` 没有变化，因此项目 schema、模板和两套 agent adapters 不需要内容迁移。v0.3.1 的关键行为是所有非 install/init 命令校验清单版本，真实 verify/archive 还会校验已安装 distribution 的完整提交。升级环境后，主 checkout 的旧清单会有意与新版引擎不匹配，所以切换必须局限在本 Change worktree。

## Goals / Non-Goals

**Goals:**

- 用真实 tag 安装和真实 `devflow install` 建立版本、提交、清单三者一致的开发环境。
- 让持久 Checks 直接验证安装身份、严格配置读取、新协议和证据门前的拒绝行为。
- 在不实现任何 CFD 能力的前提下调整静态计划优先级，并在升级收尾恢复 p1-gradient 为下一详细能力。
- 保留从 v0.2.0 启动的 Change 执行历史，不在切换前生成旧引擎 Runs。

**Non-Goals:**

- 不使用 `devflow init` 重建已有配置，不引入 `hold`、optional Check 或新 Case。
- 不升级 OpenSpec CLI、Node、Python 或产品依赖，不修改上游引擎源码。
- 不迁移或修复归档 Runs，不让旧归档证据重新成为当前证据。
- 不在主 checkout 运行新版 `devflow install` 来暂时消除版本不匹配。

## Decisions

### 1. 用旧引擎启动，用新引擎继续同一个 Change

规划、`devflow start --change upgrade-devflow-0-3-1 --create` 和首次 `next/run` 仍在 v0.2.0 下完成，避免升级前出现版本守卫拒绝。随后通过 pip 的 direct URL 语法从 `v0.3.1` 非 editable 重装 `openspec-devflow`，先核对 distribution 版本、`direct_url.json` 的 requested revision 和完整 commit，再在 Change worktree 执行新版 `devflow install`。

此后所有 devflow 命令都在该 worktree 中执行，直到合入后主 checkout 获得新清单。备选方案是另建一套临时 conda 环境同时保留两套引擎，但这会偏离用户要求更新的 `devflow` 环境并增加切换状态，因此不采用。

### 2. 安装器只负责其拥有的文件

运行 `devflow install --host claude --host kimi`，不使用 `--force`，也不手工写入清单。预期 24 个 schema/adapter 文件全部 unchanged，安装器仅更新 `openspec/devflow-install.json` 的 `engineVersion` 和 `engineCommit`。若出现 conflict 或意外资产更新，停止并审查本地改动，不强制覆盖。

`workflow-contract` 增加对清单版本/完整提交及 AGENTS 版本锚点的直接读取，并把清单纳入自身 `sourcePaths`。`framework-smoke` 继续覆盖实际 distribution 和安装资产，避免只检查文档而漏掉机器环境。

### 3. 将新版身份与拒绝路径变成持久 smoke

`framework_smoke.py` 固定 `v0.3.1`、0.3.1 和完整 commit，验证：

- distribution 非 editable、direct URL 来自要求的 tag 且 commit 完全匹配；
- 清单 `engineVersion`/`engineCommit` 与实际安装一致，`devflow install --check` 返回 current；
- 新版严格 Registry 能读取项目，`list-checks` 和 `list-cases` 成功；
- `protocol execution` 包含安装身份与执行语义契约，`protocol configuration` 可读取完整配置参考；
- 在 runs 下隔离的临时项目中，真实 CLI 对错误 `engineVersion` 在 Registry 加载前拒绝，对错误 `engineCommit` 的真实 verify 在 Run/execution 写入前拒绝；负例不修改真实清单或证据。

正向路径由本 Change 的真实 `verify --task` 和 affected verification 经过提交预检。拒绝测试只构造最小临时输入并核对目录快照，不重装不同版本、不伪造通过 Run。

### 4. 静态计划先表示当前升级，再恢复后续能力

第一任务把详细 `nextChange` 改为 `upgrade-devflow-0-3-1`，其 observable/entry/tests/rollback 对应真实 devflow 安装与 smoke；原 p1-gradient 变成一条未展开 backlog 主题，不创建其 proposal 或 spec。

升级安装和持久契约完成后，收尾任务恢复原有 p1-gradient 详细字段，并移除升级的临时详细项。这样最终合入内容不会把已完成的框架升级继续展示为下一能力，也不丢失用户原先的产品开发顺序。

### 5. 不承接旧引擎证据作为本 Change 验证

在首次 `next/run` 后才升级环境，之前只产生 start 等执行历史，不产生本 Change 的旧引擎验证 Runs。v0.3.1 的自动化执行语义版本和扩展执行定义将用于所有新 Runs；归档数据保持只读兼容，不重算或改写。

每个任务仍按项目约定先提交实现和全部声明输入，再运行 `workflow-contract`、`framework-smoke`、`solver-core-api`、`solver-core-e2e`，随后 `complete` 并提交任务状态与 Run 元数据。全部任务后运行 affected verification、OpenSpec strict validation 和独立复核；archive、push、PR、merge 分别等待明确授权。

## Risks / Trade-offs

- [全局 conda 环境升级后，主 checkout 的旧清单暂时拒绝新版 devflow] → 只在 Change worktree 执行后续 devflow；合入前回滚时重装 v0.2.0，不能通过改主清单绕过。
- [tag 移动或安装来源不符导致提交不一致] → 安装前后核对 tag 解析、requested revision 和 40 位 commit，失败则停止，不手工修补清单。
- [永久负例误写真实执行状态] → 只在 runs 下临时目录构造错误清单，比较写入前后文件集合；真实项目只运行正向验证。
- [v0.3.1 严格校验暴露旧配置问题] → 已预先确认当前项目可加载；若归档读取或后续 Change 元数据失败，按实际字段做最小兼容修复，不放宽上游校验。
- [计划收尾遗漏导致伪当前 nextChange] → 将恢复 p1-gradient 设为依赖升级任务的独立收尾任务，并由 `workflow-contract` 和最终复核检查。

## Migration Plan

1. 用 v0.2.0 完成 Change start 和首次任务选择，确认主 checkout 除已移入 worktree 的规划 artifacts 外保持干净。
2. 在 `devflow` 环境安装 v0.3.1，核实非 editable direct URL 与完整 commit；在 Change worktree 运行真实 install 并提交清单、流程文档、Checks 和临时计划优先级。
3. 通过任务级四项 Checks 后完成任务；再恢复 p1-gradient 详细计划，执行新版拒绝路径、正向验证和完整产品回归，提交收尾状态。
4. 运行 affected verification 与独立复核，向用户分别呈现实现、证据、复核和归档影响。只有明确批准后才 archive；push、PR、merge 仍分别授权。
5. 合入前放弃升级时，删除未合入 Change 的流程改动并按原 tag `v0.2.0`（commit `cc0a9f12704fdce38ae98d5750c32b338207790f`）重装环境。合入后不直接执行此回退，另开新 spec Change。
