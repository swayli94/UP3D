# P1 单元常梯度

> 本 brief 记录 `tetrahedron-p1-gradient` 的当前详细计划，只有正式提出对应 Change 后才成为活跃工作。

## Scope

- 从节点势函数计算每个四面体单元的常梯度。
- 保持现有自由流 API、CLI 和回归，不在本能力中实现流动控制方程、通量或时间推进。

## Contract

- 公共 API 为 `run_p1_gradient(path)`，命令行为 `python -m pyfp3d p1-gradient <input.json>`。
- 输入包含 `vertices`、`tetrahedra` 和与顶点一一对应的有限 `nodal_potential`。
- 输出为按单元输入顺序排列的 `N×3 gradients`。
- API 与命令行共享输入验证和结果语义，失败时不产生成功 JSON。

## Validation

- 用多个方向的仿射场验证精确梯度，并覆盖非轴对齐几何和多个单元。
- 用手工推导的非仿射节点值验证梯度，不以同一实现自证。
- 覆盖错误势函数长度、非有限值、布尔值、非法索引、退化单元和真实 API/E2E 失败行为。
- 保留现有自由流 API/E2E 回归。

## Rollback

- 只移除新势函数输入、梯度入口、对应 spec、测试和 Checks。
- 保留现有自由流 API、CLI、fixture 和全部既有回归。

## 待决事项

- 数值输出语义、极端尺度和非退化病态单元的验收方式。
- 新命令与自由流共享网格验证的边界及 `format_version` 语义。
- 新增子命令时旧 `python -m pyfp3d <input.json>` 裸文件入口的兼容方式。
