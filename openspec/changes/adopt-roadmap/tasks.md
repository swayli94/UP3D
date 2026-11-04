# Tasks

## 1. 唯一计划源、官方视图与验收

- [x] 1.1 在首次 `next/run` 后，按 design 建立仅含四个历史节点、`roadmap-adoption`、`tetrahedron-p1-gradient` 和 `cell-field-reconstruction` 的 `openspec/plan/roadmap.json`，将现有 P1 契约完整迁入 plan brief；原子移除 `development-plan.json`，更新 AGENTS、`openspec/project.json`、`openspec/workflow_check.py` 与 `workflow-contract` 输入，保持 adoption 自关联和 P1 依赖可解析，不创建后续 Change artifacts、`globalMap` 或 `parallelPlan`。运行 Change OpenSpec strict validation，提交全部声明输入后执行 `devflow verify --change adopt-roadmap --task 1.1`，以 framework-smoke、workflow-contract、solver-core-api、solver-core-e2e 验证 graph 迁移和现有产品回归。
- [ ] 1.2 在 1.1 完成后，新增 `roadmap_contract.py`、`roadmap_view_check.py` 及对应 Check 定义，在真实临时 Git 项目中验证依赖拒绝/放行和状态推导；生成并提交官方 `internal/roadmap.html`，生成机器本地 `runs/views/roadmap.html`，同提交更新产品路径、`framework_smoke.py` 六项清单和最终 execution 映射。确认普通视图双生成字节一致且无机器状态、live 视图含当前 worktree/tasks/locks 与历史 merged，运行 Change OpenSpec strict validation 和 `devflow verify --change adopt-roadmap --task 1.2`，以全部六项 Checks 验证 graph、两类视图和产品回归。
- [ ] 1.3 在 1.2 证据齐全后将 `roadmap-adoption` 存储状态改为 `done`，保持 P1 为 planned、后续仅一个 deferred/unscheduled 项，重新生成普通与 live 视图并核对 archive 前 adoption 仍为 active、P1 尚未 ready；运行 Change OpenSpec strict validation 和 `devflow verify --change adopt-roadmap --task 1.3`，再次以全部六项 Checks 验证最终 graph、brief、视图和 39 项产品测试，为 affected、milestone 与独立复核固定当前内容。

## Workflow follow-up

- 每项任务均遵循 `next/run → 实现与快速检查 → 提交任务和全部声明输入 → verify --task → complete → 提交任务状态与 Run 元数据`；不得手工勾任务、提前引用未注册 Check 或在主 checkout 落地 graph/视图。
- 全部任务后运行 affected verification 和 milestone 新鲜验证，再由独立上下文复核；复核后任何 graph、proposal、design、spec delta 或任务文字变化均重新验证并重新复核。
- archive 须取得用户明确批准，归档时同步主 spec 并刷新普通视图；提交归档结果后，push、PR、merge 分别另行授权，merge 后重新生成普通和 live 视图。
- 合入前回滚仅移除 graph、brief、内部视图、新 Checks 和治理改动并恢复 `development-plan.json`；合入后回退另开带新 spec 的 Change，不直接恢复旧流程。
