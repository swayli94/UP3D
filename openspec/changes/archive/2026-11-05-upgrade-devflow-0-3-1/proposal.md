# Proposal

## Why

当前仓库和 `devflow` conda 环境固定在 `openspec-devflow` v0.2.0，无法使用 v0.3.1 的安装身份保护、严格配置校验和新执行协议。用户要求优先升级框架，因此需要把当前 `tetrahedron-p1-gradient` 详细计划暂缓为后续能力，并以一个可验证、可回滚的工具链 Change 完成切换。

## What Changes

- 将 `devflow` 环境中的 `openspec-devflow` 从 tag `v0.2.0`（commit `cc0a9f12704fdce38ae98d5750c32b338207790f`）升级为 tag `v0.3.1`（commit `1c7ab6f223d4301aea8ea1deeb4dde2b246de37a`），保持 Python 3.12、非 editable VCS 安装和 `PYTHONNOUSERSITE=1`。
- 在本 Change 的独立 worktree 中运行真实 `devflow install`，把 `openspec/devflow-install.json` 切换到 `engineVersion: 0.3.1` 并记录完整 `engineCommit`；不使用 `--force`，不手工伪造安装清单。
- 更新 AGENTS 与 `framework-smoke`，使版本、tag 对应提交、安装清单身份和新版执行协议成为持久契约；增加对 v0.3.1 严格配置读取、版本守卫和验证前安装提交预检的真实验证。
- 暂时把静态计划的详细 `nextChange` 改为本框架升级，并将 `tetrahedron-p1-gradient` 保留为未展开 backlog 主题；升级工作完成后再恢复它为唯一的详细下一能力，不创建其产品 artifacts。
- 保持 OpenSpec CLI 1.16.0、现有 Cases/Checks、`pyfp3d` API/CLI 和全部求解器行为不变；不采用与本升级无关的新功能，不修改引擎源码，不重写历史 Runs 或归档记录。

## Capabilities

### New Capabilities

无。

### Modified Capabilities

- `development-workflow`: 增加开发框架 tag/commit 固定、安装清单身份、版本守卫与升级回路边界，并在不改变单 Change 规则的前提下调整当前能力优先级。

## Impact

- 工具环境：`devflow` conda 环境中的 `openspec-devflow` 非 editable 安装升级到 v0.3.1；主 checkout 的清单仍保持 v0.2.0，只有 Change worktree 在升级期间使用新清单。
- 流程文件：更新 `AGENTS.md`、`openspec/devflow-install.json`、`openspec/framework_smoke.py`、`openspec/workflow_check.py`、相应 Check 声明及 `openspec/plan/development-plan.json`；`devflow install` 预计保持现有 schema 和 agent 适配器内容不变。
- 验证：真实执行新版 `devflow protocol execution|configuration`、严格 Registry 读取、`install --check`、版本/提交预检、现有 `workflow-contract`、`framework-smoke`、API 和 E2E Checks，以及 OpenSpec strict validation。
- 产品：不改变 `pyfp3d` 数值行为、公共入口、测试输入或参考数据。
- 回滚边界：合入前可丢弃本 Change 的流程文件改动并重装 v0.2.0 tag；升级不迁移活跃产品状态或历史证据。合入后的任何回退必须另开带新 spec 的 Change，不能直接恢复旧流程。
