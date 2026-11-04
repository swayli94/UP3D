# Tasks

## 1. 长期工作流与迁移收尾

- [x] 1.1 在首次 `next` 后重写 AGENTS 的长期权威、`pyfp3d-core` 环境和单 Change 生命周期，新增 `openspec/plan/development-plan.json`、标准库 `openspec/workflow_check.py` 与 `workflow-contract` Check，按设计把 specs/plan/脚本加入 `openspec/project.json` 真实指纹；注册后把两个任务的 `execution.json` 都扩展为现有三个 Checks 加 `workflow-contract`，以四项 Checks、39 项 API/E2E 回归和 OpenSpec strict validation 验证流程改动，不实现 p1-gradient。
- [x] 1.2 按设计先执行备份禁推、项目根外精确导入、origin tag、bootstrap `b1db684` 与 restart `5fbed29` 的 fetch/merge-base 和两份 archive acceptance、旧路径清场、现有主 spec `tetrahedron-freestream --type spec --strict` 与 Change strict、39 项测试和 smoke；在 runs 内分别完成当前提交的 clean clone smoke 及无其他活跃 Change 的隔离 clone 中合成 tooling-only Change 的真实 propose/start/next，记录各命令、退出码、关键输出和证据边界；全部其余第 8 节判据通过后才在同一提交中写入 `migration-acceptance.md` 并删除 `docs/migration-plan.md`，确认 AGENTS/plan 无其现行引用、临时 clone/cache 已清理，再以四项 Checks 正式 verify。

## Workflow follow-up

- 每项任务均遵循 `next/run → 实现与快速检查 → 提交任务和全部声明输入 → verify --task → complete → 提交任务状态与 Run 元数据`；不得手工勾任务或在验收前删除迁移计划。
- 全部任务后运行 affected verification 和独立复核；复核后任何 proposal、design、spec delta 或任务文字变化均重新复核。
- archive 须取得用户明确批准，归档提交后 push、PR、merge 分别另行授权。后续只从静态计划独立 propose 一个 Change，不批量生成 backlog artifacts。
