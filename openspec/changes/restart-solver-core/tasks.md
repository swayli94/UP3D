# Tasks

## 1. 清场并建立最小自由流闭环

- [ ] 1.1 在首次 `next` 后，以同一提交删除旧 `pyfp3d/`、`tests/`、`bench/`、`phases/` 并建立设计规定的最小包、真实 JSON fixture、全新 API/E2E pytest、README、`environment.yml` 与 `pyproject.toml`；创建并配置 `pyfp3d-core`，注册 `solver-core-api`/`solver-core-e2e`，加入 JUnit 插件和标准库版本 provider，更新 `openspec/project.json`，且把三个任务的 `execution.json` 都扩展为 `framework-smoke` 加两个产品 Checks；先运行真实 CLI、API/E2E 和 smoke 快速检查，再提交全部声明输入并以三个必需 Checks 正式 verify，确认非零测试数和有效 JUnit。
- [ ] 1.2 在同一任务提交中删除 `cases/gates/`、`cases/demo/`、`cases/meshes/`、`docs/dev_phase_seven/`、`docs/agent-rules.md`、`docs/overview.md`、`docs/design.md`、`docs/design_track_v.md`、`PROJECT_STRUCTURE.md` 与 `setup.py`，保留 `cases/reference_data/`、`docs/migration-plan.md`、LICENSE 和框架文件；用路径 diff 核对无过时声明，并以 `framework-smoke`、`solver-core-api`、`solver-core-e2e` 正式 verify。
- [ ] 1.3 清理 `.gitignore` 中已删除 case/gate/bench 的专属规则与过时注释，保留根级 `/reference/`、`/docs/references/`、`/tools/` 及仍适用通用规则，不提交论文、工具或本地证据；用 `git check-ignore`/`git status` 和保留路径核对结果，并以 `framework-smoke`、`solver-core-api`、`solver-core-e2e` 正式 verify。

## Workflow follow-up

- 每项任务均遵循 `next/run → 实现与快速检查 → 提交任务和全部声明输入 → verify --task → complete → 提交任务状态与 Run 元数据`；不得手工勾任务或先提交纯删除。
- 全部任务后执行 affected verification，检查路径覆盖、OpenSpec strict validation 与当前产品内容，再由未参与实现的上下文独立复核并关闭全部发现。
- archive 须取得使用者明确批准；归档提交后，push、PR、merge 分别另行授权。实际 worktree 只由后续 `devflow start --create` 创建。
