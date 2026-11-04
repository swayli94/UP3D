# UP3D / pyFP3D

当前产品只提供一个最小能力：从四面体 JSON 文件读取网格和均匀无黏来流，输出逐单元速度场及自由流保持误差。它不求解流动控制方程，也不包含网格生成、外部网格格式、时间推进或可视化。

## 创建环境

```bash
conda env create -f environment.yml
conda env config vars set -n pyfp3d-core PYTHONNOUSERSITE=1
```

所有命令都从仓库根目录运行。测试直接导入当前 checkout 的源码，不需要 editable 安装。核心能力只使用 Python 标准库，pytest 仅用于测试。

## 输入

输入是 UTF-8 JSON 对象，字段如下：

- `format_version`：整数 `1`。
- `vertices`：至少四个 `[x, y, z]` 有限数值点。
- `tetrahedra`：非空单元数组；每项是四个互异、零基、有效的顶点索引，且单元非退化。
- `free_stream_velocity`：三个有限数值分量。

合法整数和浮点数保留 JSON 解析后的原始数值语义，不强制转换为双精度。退化判断以来源数值构造的精确有理数标量三重积为准，因此大整数共面单元不会被误判为非退化，小尺度非退化单元也不会因下溢被拒绝。

仓库内的真实示例是 `tests/fixtures/single_tetrahedron.json`。缺失字段、错误形状、非有限数值、非法或重复索引、退化单元都会被拒绝。

## 命令行

```bash
conda run -n pyfp3d-core --no-capture-output python -m pyfp3d tests/fixtures/single_tetrahedron.json
```

成功时 stdout 只包含一个 JSON 对象：

```json
{"velocities":[[0.75,-1.25,2.5]],"free_stream_preservation_error":0.0}
```

`velocities` 按输入单元顺序排列。`free_stream_preservation_error` 是所有单元和速度分量的最大绝对来流偏差；本能力原样复制来流，因此成功结果的误差精确为 `0.0`。文件或验证失败时，进程返回非零状态，诊断写入 stderr，不输出成功 JSON。

## Python API

```python
from pyfp3d import run_freestream

result = run_freestream("tests/fixtures/single_tetrahedron.json")
```

API 与 CLI 返回相同的两个字段。文件错误以 `OSError` 报告，输入验证错误以 `ValueError` 报告； malformed JSON 使用其 `JSONDecodeError` 子类。

## 测试

```bash
conda run -n pyfp3d-core --no-capture-output python -m pytest tests/api tests/e2e
```

API 测试覆盖输出契约和非法输入分支；E2E 测试通过真实子进程读取 fixture 并检查 stdout、stderr 和退出状态。
