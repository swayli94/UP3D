# Spec Delta

## ADDED Requirements

### Requirement: 固定开发框架安装身份

仓库 SHALL 将 `devflow` 环境的 `openspec-devflow` 固定到 tag `v0.3.1` 对应提交 `1c7ab6f223d4301aea8ea1deeb4dde2b246de37a`，并要求非 editable VCS 安装。安装清单 SHALL 记录相同版本与完整提交，永久检查 SHALL 验证实际安装、清单和项目配置的一致性。

#### Scenario: 当前工具链通过检查

- **WHEN** 开发者在声明的 `devflow` 环境中运行框架 smoke
- **THEN** 检查确认版本 0.3.1、完整安装提交、`install --check` 为 current，且新版严格配置读取成功

#### Scenario: 引擎版本与清单不一致

- **WHEN** 运行引擎版本与安装清单的 `engineVersion` 不同
- **THEN** 除安装修复入口外的 devflow 命令在工作流写入前拒绝执行，并提示重新安装匹配资产

#### Scenario: 安装提交与证据门不一致

- **WHEN** 实际安装提交与清单的 `engineCommit` 不同并进入真实 verify 或 archive
- **THEN** devflow 在产生新 Run、验收或其他工作流写入前拒绝执行

### Requirement: 框架升级的隔离与回滚

框架升级 SHALL 通过单一 Change 和独立 worktree 完成，安装清单 SHALL 由真实安装命令生成且不得用手工编辑或强制覆盖代替。合入前回滚 SHALL 同时还原流程文件与旧 tag 环境；合入后回退 SHALL 另开带新 spec 的 Change。

#### Scenario: 切换期间隔离主 checkout

- **WHEN** 升级 Change 已在独立 worktree 安装新版引擎资产但尚未合入
- **THEN** 新版清单和流程改动只存在于该 worktree，主 checkout 不被直接改写

#### Scenario: 合入前放弃升级

- **WHEN** 使用者决定放弃尚未合入的升级
- **THEN** 可丢弃该 Change 的流程文件改动并重装旧 tag，不迁移或重写产品状态和历史证据

#### Scenario: 合入后需要回退

- **WHEN** 升级版已经合入后需要恢复旧框架行为
- **THEN** 另开带新 spec 的 Change 重新推导验证与回滚边界，不直接恢复旧流程文件
