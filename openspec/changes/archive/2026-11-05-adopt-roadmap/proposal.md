# Proposal

## Why

当前仓库只有项目自定义的 `development-plan.json` 静态队列，devflow 无法据此形成规范 graph、推导依赖状态或生成官方路线图；一次性审计 HTML 也不能承担长期追踪。需要显式采用 devflow 的协调计划和视图，同时保持 UP3D 的单 Change、真实验证与逐步展开纪律。

## What Changes

- **BREAKING（仅规划接口）**：以 `openspec/plan/roadmap.json` 作为唯一协调 graph，完成迁移后退役 `development-plan.json`，不长期双写两份计划；产品 API、CLI 与数值行为不变。
- 建立最小初始 graph：四个已归档历史 Change、`roadmap-adoption`、`tetrahedron-p1-gradient`，以及唯一未展开的 `cell-field-reconstruction` deferred 项；后续节点在相应能力完成后再逐步添加。
- 将现有 P1 梯度的入口、可观察结果、输入输出、测试与回滚边界迁入 plan brief，不把详细范围塞进 graph item，也不批量创建后续 Change artifacts。
- 首版不维护可选的 `globalMap` 或 `parallelPlan`；以 `phases / lanes / items / deps` 和 Change 关联表达真实协调关系。
- 生成可提交、确定性的 `internal/roadmap.html`，以及不提交的本机 `runs/views/roadmap.html` live 视图；明确只有 live 视图包含 worktree、任务进度、资源锁和 merged 状态。
- 更新 AGENTS、development-workflow 契约、项目视图配置和永久 Checks，直接验证 graph、引用、状态推导、依赖启动门及两类官方视图。
- 在独立 worktree 和非 main 分支实施，经任务级验证、affected verification、独立复核和用户批准后归档；push、PR、merge 仍分别授权。

## Capabilities

### New Capabilities

无。

### Modified Capabilities

- `development-workflow`：以 canonical roadmap graph 取代自定义静态能力队列，规定唯一计划源、最小图维护、Change 关联、状态推导、官方视图与验证要求。

## Impact

- 规划与治理：`openspec/plan/`、`openspec/specs/development-workflow/`、AGENTS、永久工作流检查与 `openspec/project.json`。
- 本 Change：proposal、design、spec delta、tasks、Change 元数据和执行元数据。
- 生成视图：`internal/roadmap.html` 与机器本地、被忽略的 `runs/views/roadmap.html`。
- 不影响：`pyfp3d/` 产品代码、现有 API/E2E 语义、`cases/reference_data/`、只读备份和已归档科学证据。
- 回滚边界：合入前可移除本 Change 的 graph、brief、视图和治理改动并恢复静态队列；合入后如需回退，另开带新 spec 的 Change，不直接恢复旧流程。
