# Proposal

## Why

旧求解器、测试、阶段文书和历史门已不再构成当前产品的开发权威；继续在其上修补会阻碍按新规格建立可验证能力。现在需要按迁移计划第 2 步清场，并以一个真实、极小的四面体自由流闭环重新开始。

## What Changes

- **BREAKING**：在同一任务、同一提交中删除旧 `pyfp3d/`、`tests/`、`bench/`、`phases/`，重建可导入的最小 `pyfp3d` 包、全新 API/E2E 测试、README 和实际依赖清单；不保留纯删除中间提交。
- 新增唯一产品能力：从新的极小四面体网格 JSON 输入读取顶点、单元和均匀无黏来流，通过公共 API 与真实命令行入口输出逐单元速度场及自由流保持误差，并明确验证失败行为。
- 建立新产品环境 `pyfp3d-core` 及两个必需 pytest Checks；保留并继续运行独立的 `framework-smoke`，不把旧测试值、旧 gate、历史性能或 `cases/reference_data/` 作为验收依据。
- 按迁移计划删除其余旧 case、阶段文档和过时打包文件，保留 `cases/reference_data/` 只读资料、`docs/migration-plan.md`、LICENSE 与第 1 步工具链接入产物，并同步清理 `.gitignore`。
- 更新 devflow 的产品路径、任务覆盖路径、Check 输入指纹、JUnit 证据目录和包版本绑定，使验证对应当前新实现而非旧目录树。

## Capabilities

### New Capabilities

- `tetrahedron-freestream`: 从极小四面体网格输入施加均匀无黏速度场，并通过 API/CLI 提供结构化速度场、输入错误行为和自由流保持误差判据。

### Modified Capabilities

无；当前没有主 specs，本 Change 不修改既有 capability。

## Impact

- 产品代码与入口：重建 `pyfp3d/`，新增公共 API 和 `python -m pyfp3d` 命令行入口。
- 测试与 fixture：重建 `tests/`，新增 API/E2E pytest 与新的单四面体输入，不恢复旧测试树。
- 环境与依赖：新建 `pyfp3d-core`；核心只使用 Python 标准库，pytest 仅用于测试，不引入 NumPy、Numba、SciPy、网格生成、可视化或线性求解栈。
- devflow：影响 `framework-smoke`、新 `solver-core-api` 和 `solver-core-e2e` Checks，以及 `openspec/project.json`、Check 定义和仅标准库的测试/版本辅助模块；要求 Linux 验证和独立复核。
- 清理：影响迁移计划第 6.1 节列出的旧产品、case、文档、打包与忽略规则路径；不修改只读备份，不归入历史参考数据为数值真值。
