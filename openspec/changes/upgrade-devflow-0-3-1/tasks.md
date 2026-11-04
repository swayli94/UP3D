# Tasks

## 1. 框架切换与持久契约

- [x] 1.1 在首次 `next/run` 后将 `devflow` 环境按 tag `v0.3.1` 非 editable 升级并核对完整 commit，在本 worktree 运行真实 `devflow install`；更新 AGENTS、`framework_smoke.py`、`workflow_check.py` 与 `workflow-contract` 的清单输入，覆盖 distribution/清单身份、严格 Registry、新协议及隔离负例，并把静态 `nextChange` 暂改为本升级、p1-gradient 转为未展开 backlog；提交全部声明输入后运行 `devflow verify --change upgrade-devflow-0-3-1 --task 1.1`，以 workflow-contract、framework-smoke、solver-core-api、solver-core-e2e 四项 Checks 验证安装和回归行为。
- [ ] 1.2 在 1.1 完成后确认主 checkout 的 v0.2.0 清单未被改写、worktree 的 v0.3.1 清单与实际安装一致，执行隔离的版本守卫和 commit 预检负例并确认没有真实 Run/execution 写入；随后恢复原 `tetrahedron-p1-gradient` 详细计划并移除临时升级项，运行 Change OpenSpec strict validation 和 `devflow verify --change upgrade-devflow-0-3-1 --task 1.2`，再次以四项 Checks 验证最终计划、工具链和现有 39 项产品测试，不创建 p1-gradient 产品 artifacts。

## Workflow follow-up

- 每项任务均遵循 `next/run → 实现与快速检查 → 提交任务和全部声明输入 → verify --task → complete → 提交任务状态与 Run 元数据`；不得手工勾任务或在主 checkout 运行新版 install。
- 全部任务后运行 affected verification 和独立复核；复核后任何 proposal、design、spec delta 或任务文字变化均重新复核。
- archive 须取得用户明确批准，归档提交后 push、PR、merge 分别另行授权。合入后的框架回退只能另开带新 spec 的 Change。
