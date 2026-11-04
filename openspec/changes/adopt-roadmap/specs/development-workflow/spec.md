# Spec Delta

## ADDED Requirements

### Requirement: 唯一协调 graph

仓库 SHALL 以 `openspec/plan/roadmap.json` 作为唯一协调 graph，并在迁移完成后移除 `development-plan.json`。graph 驱动开发顺序和 Change 关联，但 MUST NOT 取代主 specs 或活跃 Change 的产品行为与验收契约。

#### Scenario: 查看当前协调计划

- **WHEN** 开发者需要确定 roadmap 顺序、依赖或 Change 关联
- **THEN** 只需读取 canonical graph，无需比较或同步第二份静态计划

#### Scenario: 行为契约与协调计划分离

- **WHEN** graph 与主 specs 或活跃 Change 的行为要求不一致
- **THEN** 以主 specs 和活跃 Change 为产品行为依据，并通过规划流程修正 graph

### Requirement: 最小初始 graph 与逐步展开

初始 graph SHALL 只包含四个已归档历史节点、`roadmap-adoption`、`tetrahedron-p1-gradient` 和一个 `cell-field-reconstruction` deferred 项。P1 的详细计划 SHALL 保存在 plan brief；后续主题 SHALL 保持未展开，并在相应能力完成后逐个添加或提升，不得批量生成 Change artifacts。

#### Scenario: 初次采用路线图

- **WHEN** 初始 graph 被读取
- **THEN** 四个历史节点关联各自归档 Change，adoption 关联本 Change，P1 依赖 adoption，后续只有 `cell-field-reconstruction` 为 `deferred` 且属于 `unscheduled`

#### Scenario: 完成一个后续能力

- **WHEN** 已详细展开的 roadmap 项完成并归档
- **THEN** 维护者才添加或展开下一个主题，不同时详细规划多个后续能力 Change

### Requirement: graph 存储状态与 Change 关联

roadmap item SHALL 只使用 `planned`、`draft`、`deferred` 或 `done` 作为存储状态。每个 Change SHALL 关联恰好一个 item；active、ready、archived 和 merged SHALL 由 Change、依赖、归档记录和集成分支推导，不得手工写入存储状态。

#### Scenario: Change 开始执行

- **WHEN** 一个 item 关联的 Change 处于活跃目录
- **THEN** 该 item 显示为 active，启动前其 `deps` 必须全部完成

#### Scenario: Change 完成并进入主线

- **WHEN** 维护者已将 item 标记为 done 且其 Changes 全部归档
- **THEN** item 显示为 archived，并仅在归档记录存在于集成分支时于 live 视图显示 merged

### Requirement: 官方路线图视图

仓库 SHALL 通过 devflow 生成 `internal/roadmap.html` 确定性快照和 `runs/views/roadmap.html` 本机 live 视图。live 视图 SHALL 保持机器本地且不提交；HTML MUST NOT 作为修改 graph、任务、Run、review 或 approval 的入口。可选 `globalMap` 与 `parallelPlan` MUST NOT 成为首版生成或验收的前提。

#### Scenario: 生成普通视图

- **WHEN** 开发者从相同已提交输入重复生成普通视图
- **THEN** 输出为官方模板的确定性 HTML，不包含本机 worktree、资源锁或任务计数

#### Scenario: 生成本机 live 视图

- **WHEN** 开发者生成 live 视图
- **THEN** 页面包含当前 worktree、任务进度、资源锁和 merged 状态，输出位于被忽略的 runs 路径

### Requirement: graph 维护与视图刷新

graph、brief 和 item 存储状态 SHALL 通过正常规划流程维护；任务状态 SHALL 仅通过 devflow 执行流程更新。普通视图 SHALL 在规划、Change、PR、merge 或 archive 边界刷新，live 视图 SHALL 在本机任务、worktree、提交、锁或 merged 状态变化后刷新，不得假定安装过程提供自动 Git hook 或 CI。

#### Scenario: 修改计划或完成状态

- **WHEN** 维护者改变 graph、关联、存储状态或 brief
- **THEN** 通过相应规划流程保存变更并重新生成普通视图，需要本机状态时同时重新生成 live 视图

#### Scenario: 执行状态发生变化

- **WHEN** 任务、worktree、资源锁或集成分支状态发生变化
- **THEN** 不手改 HTML 或 graph 来伪造进度，而是重新生成本机 live 视图

## MODIFIED Requirements

### Requirement: 流程改动的真实验证覆盖

主 specs、canonical graph、plan briefs、官方普通视图和永久工作流检查脚本 SHALL 纳入产品指纹及相应 Check 输入。流程、graph 或文档改动 SHALL 由工作流与 roadmap Checks 直接验证，不能只依靠产品 API/E2E 退出码推定合规。

#### Scenario: 修改流程文档或计划

- **WHEN** 一个任务修改 AGENTS、主 specs、canonical graph、plan briefs、官方普通视图或工作流检查脚本
- **THEN** 工作流与 roadmap Checks 直接读取这些内容，验证结构、引用、状态推导、必需规则与视图确定性，且所有声明输入真实存在

#### Scenario: 验证依赖和本机状态

- **WHEN** roadmap 依赖或 live 状态需要验收
- **THEN** Checks 在真实临时输入上验证未完成依赖的启动拒绝、完成后的 ready 推导，以及 live 视图的任务与 merged 信息，不依赖本机备份或历史 Run

## REMOVED Requirements

### Requirement: 简洁静态能力队列

仓库 SHALL 维护一个静态计划，至多详细描述一个 `nextChange`，包含其真实入口、可观察结果和回滚边界。backlog SHALL 只含未展开主题，不得批量生成 Change artifacts；完成详细项后 SHALL 清空或替换它，不得遗留伪当前计划。

#### Scenario: 查看下一个能力

- **WHEN** 当前没有正在执行的后续能力 Change
- **THEN** 计划只给出一个可独立 propose 的下一能力，其他主题没有 proposal、design、spec delta 或任务文件

#### Scenario: 完成计划中的能力

- **WHEN** 计划中详细描述的 Change 已完成
- **THEN** 计划清空或替换该详细项，不把它继续展示为待执行的当前下一项

**Reason**: 自定义静态队列不能被 devflow 读取，也无法提供 Change 关联、依赖推导和官方视图。

**Migration**: 将现有 P1 nextChange 转为 roadmap item 与 plan brief，将唯一 backlog 转为 deferred/unscheduled item；验证 graph 和视图后移除 `development-plan.json`。
