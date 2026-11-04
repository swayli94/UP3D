# Proposal

## Why

UP3D 需要统一的开发流程、真实的工具链检查和可复核的 Change 元数据。当前仓库尚无已安装的 devflow 项目配置；本 Change 只建立最小框架接入，不延续旧 phase 7 的开发方式。

## What Changes

- 以 `package.json` 和 `package-lock.json` 精确固定 `@fission-ai/openspec` 1.16.0，并安装主 checkout 使用的本地 CLI。
- 安装 devflow schema、Claude/Kimi adapters 和安装记录，建立最小项目配置及 `framework-smoke` 必需 Check。
- 建立本 Change 的 tooling-only artifacts，以标准库 Smoke 检查 Python、devflow 提交、Node、OpenSpec、安装状态及 registry 列表命令。
- 重写协作入口规则，删除作废的 `adopt-devflow` 草案，更新本地运行路径的忽略规则，并把当前迁移计划移到 `docs/migration-plan.md`。
- 不删除旧产品代码，不接入旧测试，不修复旧 bug，不实现求解器功能。

## Capabilities

### New Capabilities

无。本 Change 不改变 spec 级产品行为，`.openspec.yaml` 设置 `skip_specs: true`。

### Modified Capabilities

无。

## Impact

影响开发工具、agent 协作规则和迁移文档；`pyfp3d` 仅作为真实存在的最小交付路径声明，不由本 Change 修改。框架检查不导入 `pyfp3d`，不读取旧测试，也不声明任何 Case 或产品任务层覆盖。
