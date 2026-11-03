# Proposal

> 草案（2026-11-04，phase 7 第 10 轮）。依据：[docs/dev_phase_seven/20261104-0524-devflow-adoption-audit.md](../../../docs/dev_phase_seven/20261104-0524-devflow-adoption-audit.md)
> （下称“审计”）。devflow 尚未安装；本 Change 在审计 §一 阻塞 1（迁移机 12 红）处置前不 `start`。
> ★ install（tasks 2.1）之前，`openspec validate adopt-devflow --strict` 会报 `unknown schema 'devflow'`，这是预期的；
> 2026-11-04 在带 devflow schema 的临时副本里校验：`Change 'adopt-devflow' is valid`。

## Why

UP3D 的研究纪律（预注册、G-DOMAIN/G-TEETH/G-CADENCE/G-ATTRACTOR、四值判定、只认已提交证据）已经很严，
但依赖关系、验证新鲜度、独立复核/用户批准、并行隔离和交接仍靠散文与会话记忆。迁移机上 12 条证据锁变红
（审计 §3.3）说明了代价：**答案的出处（编译目标、缓存态、被忽略的网格）在现有流程里看不见**，
一次旧的“通过”无法说明它今天是否仍然成立。

## What Changes

把 openspec-devflow v0.2.0（`cc0a9f1`）作为**治理与证据层**叠加到仓库上，不重写仓库、不反向补写规格：

- 新增 `openspec/`：`config.yaml`、`project.json`、`checks/`、`cases/`、`plan/roadmap.json`；
  以及 `devflow install --host claude --host kimi` 写入的 schema 与两套 agent 适配器。
- 新增 `package.json` + `package-lock.json`，把 OpenSpec CLI 精确钉在 1.16.0（运行时不校验版本，lockfile 是唯一的钉）。
- 新增独立工具环境文件（Python 3.12 + Node ≥ 22.12），与求解环境 `up3d`（Python 3.11）分离。
- 新增 `tools/devflow/`：pytest Check wrapper、环境/有效 codegen provider、Case raw 网格校验与物化。
- 新增被忽略网格的 Case raw 定义（外部 `UP3D_DATA_DIR`），以已提交的 `.msh.manifest.json` 为准。
- `.gitignore` 增加 `runs/`、`node_modules/`。
- **BREAKING（协作流程）**：改写 `AGENTS.md` —— Claude 不再“在 `main` 上直接修改和提交”，与 Kimi 一样走
  Change 分支 / worktree；`AGENTS.md` 成为协作与状态转换的单一权威，`CLAUDE.md` 保持薄入口。
- **不改**：`pyfp3d/`、`tests/`、`cases/reference_data/`、committed gate 证据、`phases/`、求解器默认值与 API。

## Capabilities

### New Capabilities

（无。tooling-only Change：不改变任何可观察的求解器行为，`.openspec.yaml` 设 `skip_specs: true`。）

### Modified Capabilities

（无。）

## Impact

- **流程**：所有后续 Change 走 propose → start → task commits → verify → independent review → user approval → archive；
  archive、merge、push、release 是四个分开的授权动作。
- **代码**：只新增 `openspec/`、`tools/devflow/`、工具环境与 Node 锁文件；求解器、测试、证据零改动。
- **依赖**：新增工具环境（Python 3.12、Node ≥ 22.12、devflow 0.2.0、OpenSpec 1.16.0）；求解环境不变。
- **前置条件（不在本 Change 内）**：审计阻塞 1 —— 12 红的逐项分类与处置（独立的 `baseline-provenance` 工作，
  四臂运行于 2026-11-04 启动）。任何再基线或平台专用证据都需另行的使用者裁决，**不得**夹带进本 Change。
- **已裁决（2026-11-04）**：缓存新鲜度 = 审计 §5 P0-2 方案 A；skip/xfail = 审计 §7.3 方案 A（见 design.md D5、D6）。
