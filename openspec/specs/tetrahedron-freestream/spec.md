# tetrahedron-freestream Specification

## Purpose
为重启后的最小产品建立一个可独立验证的四面体自由流能力：调用者能从真实网格文件获得均匀无黏速度场和明确的保持误差，而不依赖旧求解器、旧测试或历史基线。

## Requirements

### Requirement: 四面体自由流输入文档

系统 SHALL 接受 UTF-8 JSON 对象，包含等于 `1` 的 `format_version`、至少四个三维有限 `vertices`、非空 `tetrahedra` 和三维有限 `free_stream_velocity`。每个单元 SHALL 按输入顺序给出四个互异、零基且有效的顶点索引。

#### Scenario: 合法单四面体输入

- **WHEN** 调用者提供满足文档契约的一个非退化四面体及非零自由流速度
- **THEN** 系统接受输入，并保持顶点和单元的输入顺序

### Requirement: 输入验证与失败行为

系统 SHALL 在输出结果前拒绝 malformed JSON、字段缺失或类型错误、非有限数值、维度错误、空单元、非整数或越界或重复索引，以及共面或重合导致的退化单元；布尔值 MUST NOT 作为整数索引接受。

#### Scenario: 文档结构无效

- **WHEN** 输入不可解析、缺少必需字段、数值非有限或数组维度不符合契约
- **THEN** 系统报告包含无效字段或输入原因的诊断，且不产生成功结果

#### Scenario: 单元拓扑无效

- **WHEN** 任一单元含布尔、非整数、越界或重复索引，或其四个顶点退化
- **THEN** 系统拒绝整个输入，且不输出部分速度场

### Requirement: 均匀无黏速度场

系统 SHALL 将输入 `free_stream_velocity` 的三个分量原样赋予每个四面体单元，不做时间推进、通量计算或近似求解。

#### Scenario: 多单元保持同一来流

- **WHEN** 合法输入含一个或多个单元及任意有限来流向量
- **THEN** 每个输出单元的三个速度分量分别与输入来流分量完全相等

### Requirement: 结构化输出契约

成功结果 SHALL 仅包含 `velocities` 与 `free_stream_preservation_error`。`velocities` SHALL 是按单元输入顺序排列的 `N × 3` 有限数值数组；误差 SHALL 是非负有限标量，API 与命令行成功结果 SHALL 具有相同语义。

#### Scenario: 输出形状和顺序可消费

- **WHEN** 系统成功处理含 `N` 个单元的输入
- **THEN** 输出仅含两个约定字段，速度数组有 `N` 行且每行三个分量，顺序与输入单元一致

### Requirement: 自由流保持误差判据

自由流保持误差 SHALL 定义为所有单元和分量上 `|velocity[i][j] - free_stream_velocity[j]|` 的最大值。满足本 capability 的成功结果 MUST 使该误差精确等于 `0.0`；验收 MUST NOT 使用旧测试容差或历史性能阈值替代此判据。

#### Scenario: 代数上无近似误差

- **WHEN** 系统按均匀无黏速度场要求生成结果
- **THEN** 最大绝对分量误差为 `0.0`，且报告值与由输出速度重新计算的误差一致

### Requirement: 真实文件入口

系统 SHALL 提供读取同一输入文件路径的公共 API 和命令行入口。命令行成功时 SHALL 在 stdout 写出一个符合输出契约的 JSON 对象；失败时 SHALL 返回非零状态、在 stderr 给出诊断且不写出成功 JSON。

#### Scenario: 命令行端到端读取

- **WHEN** 调用者以真实存在的合法输入文件路径调用命令行入口
- **THEN** 进程返回零状态，stdout 可直接解析为一个成功结果 JSON 对象

#### Scenario: 命令行拒绝无效文件

- **WHEN** 调用者以不存在或不符合输入契约的文件路径调用命令行入口
- **THEN** 进程返回非零状态，stderr 包含诊断，stdout 不包含成功结果 JSON
