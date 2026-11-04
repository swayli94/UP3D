# Proposal

## Why

最小产品、主 specs 和 devflow 生命周期已经建立，但现行 AGENTS 仍引用迁移入口，仓库还没有面向未来的单 Change 开发契约与简洁能力队列。现在需要把开发流程转为长期规则，并在核实全部迁移判据后移除临时迁移计划。

## What Changes

- 重写 AGENTS 的开发流程部分：权威仅限 AGENTS、主 specs、活跃 Change 和 devflow protocol；产品环境明确为 `pyfp3d-core`，标准入口固定为同步干净 main → propose → `start --create` → `next/run`。
- 新增 `development-workflow` 行为契约，规定一次最多一个活跃 Change，每个 Change 必须有 spec delta、真实入口测试、独立 worktree、可观察能力、回滚边界、任务级 Checks、独立复核及分级审批。
- 新增简洁静态 `openspec/plan/`：只详细描述下一个 `tetrahedron-p1-gradient` 能力 Change，backlog 只保留未展开主题，不批量生成 Change artifacts，也不实现该能力。
- 新增标准库 `workflow-contract` Check，并将主 specs、静态计划和检查脚本纳入真实产品指纹；现有 `framework-smoke`、API、E2E 继续作为每项任务的回归门。
- 通过一次性命令和 Change 内验收记录核实备份禁推与精确导入、origin tag、bootstrap `b1db684` 和 restart `5fbed29` 的 merge-base/归档 acceptance、现有主 spec strict、39 个产品测试及 runs 内 clean clone smoke；另在隔离临时 clone 中对合成 tooling-only Change 真实执行 propose/start/next，全部其余判据通过后才删除 `docs/migration-plan.md`。
- 不恢复旧 phase/round 文书、旧代码或旧测试，不新增 CFD 功能，不重设计 devflow 工具本身。

## Capabilities

### New Capabilities

- `development-workflow`: 定义仓库如何以单 Change、规格驱动、真实入口验证、独立复核和静态能力队列持续开展小步开发。

### Modified Capabilities

无；本 Change 不改变 `tetrahedron-freestream` 的产品行为。

## Impact

- 流程与计划：更新 `AGENTS.md`，新增 `openspec/plan/development-plan.json` 和新的主 spec capability；归档后不保留对临时迁移计划的现行引用。
- 验证配置：新增 `openspec/workflow_check.py`、`openspec/checks/workflow-contract.json`，更新 `openspec/project.json` 的 `productPaths` 和本 Change 的 `execution.json`。
- 回归范围：每个任务最终运行 `workflow-contract`、`framework-smoke`、`solver-core-api`、`solver-core-e2e`；OpenSpec strict validation 继续必需。
- 迁移收尾：在 Change 内记录一次性验收命令及原始结果，最后删除 `docs/migration-plan.md`；永久 Checks 不依赖备份路径、旧 conda 环境或本机 tag 状态。
- 下一步但不在本 Change 实现：规划 `tetrahedron-p1-gradient` 的公共 API/CLI 输入输出、真实测试与独立回滚边界。
