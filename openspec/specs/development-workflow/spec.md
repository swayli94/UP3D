# development-workflow Specification

## Purpose
为 UP3D 建立可持续的单 Change 开发方式，使每个小步能力都有明确规格、真实入口验证、独立复核和回滚边界，同时让静态计划与永久检查不依赖迁移历史或特定机器。

## Requirements

### Requirement: 统一开发权威与产品环境

仓库 SHALL 仅以 AGENTS、主 specs、活跃 Change artifacts 和 devflow execution protocol 为开发权威。AGENTS SHALL 将 `pyfp3d-core` 记录为产品环境，并 MUST NOT 把临时迁移文档或已完成迁移 Change 作为未来工作入口。

#### Scenario: 开始常规开发

- **WHEN** 开发者读取 AGENTS 准备常规产品工作
- **THEN** 能找到当前环境、主 specs 和标准 Change 入口，且无需查阅临时迁移文档或历史 phase/round 记录

### Requirement: 单 Change 标准入口

仓库 SHALL 一次最多规划或执行一个活跃 Change。每个 Change SHALL 从同步且干净的 main 依次经过 propose、`start --create` 和 `next/run`，在声明的独立 worktree 与非 main 分支中实施。

#### Scenario: 启动一个能力 Change

- **WHEN** 主 checkout 与 origin/main 同步且除该 Change 规划目录外干净
- **THEN** 开发者通过标准入口建立独立 worktree，并在其中选择任务后开始实现

### Requirement: 完整 Change 交付契约

每个 Change SHALL 包含一个新的或修改的 capability spec delta、对应真实入口测试、可观察能力与回滚边界。每项任务 SHALL 声明必需 Checks；交付 SHALL 经过当前内容验证、独立复核、用户 archive 批准，并对 push、PR、merge 分别取得授权。

#### Scenario: 审阅能力交付

- **WHEN** 一个 Change 请求进入 archive
- **THEN** 其规格、真实入口证据、任务 Checks、回滚说明和独立复核均可检查，且后续发布操作仍需分别授权

### Requirement: 简洁静态能力队列

仓库 SHALL 维护一个静态计划，至多详细描述一个 `nextChange`，包含其真实入口、可观察结果和回滚边界。backlog SHALL 只含未展开主题，不得批量生成 Change artifacts；完成详细项后 SHALL 清空或替换它，不得遗留伪当前计划。

#### Scenario: 查看下一个能力

- **WHEN** 当前没有正在执行的后续能力 Change
- **THEN** 计划只给出一个可独立 propose 的下一能力，其他主题没有 proposal、design、spec delta 或任务文件

#### Scenario: 完成计划中的能力

- **WHEN** 计划中详细描述的 Change 已完成
- **THEN** 计划清空或替换该详细项，不把它继续展示为待执行的当前下一项

### Requirement: 流程改动的真实验证覆盖

主 specs、静态计划和永久工作流检查脚本 SHALL 纳入产品指纹及相应 Check 输入。流程或文档改动 SHALL 由工作流契约 Check 直接验证，不能只依靠产品 API/E2E 退出码推定流程合规。

#### Scenario: 修改流程文档或计划

- **WHEN** 一个任务修改 AGENTS、主 specs、静态计划或工作流检查脚本
- **THEN** 工作流契约 Check 读取这些内容并验证其结构与必需规则，且所有声明输入真实存在

### Requirement: 永久检查可在干净环境运行

永久 Checks SHALL 只依赖仓库内容、声明的工具环境和真实临时输入，不得要求本机备份目录、备份 conda 环境、origin tag 或历史证据才能通过。clean clone 验证 SHALL 把 clone、Node 依赖和缓存放在 runs 下的临时位置。

#### Scenario: 在新机器运行永久检查

- **WHEN** 开发者在没有旧备份和历史运行记录的干净 clone 中安装声明依赖
- **THEN** 永久 Checks 能基于当前仓库内容运行，不访问旧备份路径或旧产品环境

### Requirement: 一次性迁移验收与收尾

备份禁推、精确备份导入、origin tag、bootstrap/restart 合入祖先关系与两份 acceptance、现有主 spec strict、旧路径清场和 39 项测试 SHALL 由一次性命令及原始记录核实，不得固化为永久 Check 的本机依赖。`docs/migration-plan.md` MUST 只在其余判据均核实后的收尾提交中删除。

#### Scenario: 完成迁移收尾

- **WHEN** 一次性验收、主 spec 与 Change strict validation、39 项测试、smoke、clean clone smoke 和标准入口实证均已通过
- **THEN** Change 记录各项命令、退出码和关键输出并删除临时迁移计划，AGENTS 与静态计划不再引用它，回滚改用只读备份或 origin tag

### Requirement: 标准入口一次性实证

迁移验收 SHALL 在没有其他活跃 Change 的隔离临时 clone 中，为一个非产品 synthetic smoke Change 真实执行 propose、`start --create` 和 `next`，记录同步状态、生成的 worktree 与 ready 任务；MUST NOT 在主仓库创建 p1-gradient artifacts。

#### Scenario: 隔离验证标准入口

- **WHEN** 临时 clone 的 HEAD 与其 origin/main 同步、工作树干净且没有其他活跃 Change
- **THEN** 合成 smoke Change 通过真实 propose/start/next 建立独立 worktree 并返回 ready 任务，不实现产品能力也不改变主仓库分支

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
