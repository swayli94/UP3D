# Design

## Context

见 `proposal.md` 的 Why。主 specs 为空，现有 devflow 注册表只有 `framework-smoke`；旧 `pyfp3d/`、`tests/`、`bench/`、`phases/` 仍在工作树中，但不提供新能力契约。`specs/tetrahedron-freestream/spec.md` 已固定唯一新能力的输入、输出、失败行为与精确零误差判据。

已核对 devflow v0.2.0 实现：`ExecutionEngine._metadata` 在每次 `next` 时校验所有任务的 Check id；`start --create` 末尾还会用启动前加载的注册表调用 `next`。Check 命令以参数数组直接启动，不经过 shell；JUnit 只从当前 Run 证据目录读取 `*.junit.xml`。这些行为决定下面的注册与证据设计。

## Goals / Non-Goals

**Goals:**

- 用一个任务提交同时完成旧核心/测试移除和最小核心、测试、文档、环境、Checks 的建立，任何提交都可导入并验证。
- 让 API 与 CLI 共用同一路径读取、验证和结果生成逻辑，E2E 通过真实子进程与真实 JSON 文件观察行为。
- 让正式 Run 绑定全部声明输入、产品解释器与包版本，并由 JUnit 防止零用例或全 skip 被误判为通过。
- 完成迁移计划第 6.1 节剩余清理，同时保留指定只读资料和框架文件。

**Non-Goals:**

- 不实现流动控制方程、通量、边界条件、时间推进、网格生成或外部网格格式。
- 不引入 Numba、SciPy、AMG、可视化、缓存、旧 gate、性能基线或旧测试期望值。
- 不把 `cases/reference_data/` 注册为 Case、自动输入或数值真值，不修改 `/home/lrz/code/up3d_backup`。
- 不在本 Change 中重设计后续开发流程，不删除 `docs/migration-plan.md`，不创建 CI 或完整 CFD 功能。

## Decisions

### 1. 最小包与真实输入入口

新包只保留 `pyfp3d/__init__.py`、`pyfp3d/freestream.py` 和 `pyfp3d/__main__.py`。公共 `run_freestream(path)` 返回含两个约定字段的普通 Python 字典；`python -m pyfp3d <input.json>` 调用同一函数并将结果序列化到 stdout。文件、JSON、验证错误在 API 中保留明确异常，在 CLI 中转换为非零退出和 stderr 诊断，不在 stdout 混入日志。

输入 fixture 新建为 `tests/fixtures/single_tetrahedron.json`，使用简单非退化单元和有正有负的来流分量。验证包括字段类型、布尔值排除、有限性、数组形状、索引范围与互异性。合法 `int/float` 不强制转换为双精度；以原始数值构造 `fractions.Fraction` 并精确计算标量三重积，避免大整数舍入及浮点上下溢。逐单元复制原始来流对象后，以 `Fraction` 计算输出相对原始输入的最大绝对分量差，因此精确零是由新规格和代数定义推导的，不继承旧数值。

备选方案是恢复 meshio 和旧网格入口，或以 NumPy 双精度执行验证；前者扩大格式与依赖范围并混入旧行为，后者不能保持大整数语义或精确判断所有有限输入，故均不采用。

### 2. 新环境和实际依赖

产品环境命名为 `pyfp3d-core`，不复用 `up3d`、`up3d-dev` 或 `devflow`。`environment.yml` 固定 `python=3.12`、conda-forge/`nodefaults` 和 `pytest>=8`；核心只使用 Python 标准库，pytest 是本 Change 新选定的测试框架。`pyproject.toml` 只保留最小项目元数据、空运行依赖和 pytest 配置；`README.md` 写明环境创建、`PYTHONNOUSERSITE=1`、CLI 与测试命令。无需 NumPy 或 SciPy，因为数值保存、精确退化判断和误差计算均由标准库完成。

apply 时按新清单创建环境并执行 `conda env config vars set -n pyfp3d-core PYTHONNOUSERSITE=1`，随后用 `conda run -n pyfp3d-core --no-capture-output ...` 验证；规划阶段不创建环境。测试直接从 Change worktree 导入源码，不做 editable 安装，避免旧环境的安装指向遮蔽当前树。

### 3. 原子重启与后续清场

任务 1.1 的一个提交同时包含：删除旧 `pyfp3d/`、`tests/`、`bench/`、`phases/`；建立新包、fixture、API/E2E 测试、README、`environment.yml`、`pyproject.toml`；新增两个 Check 及辅助模块；更新 `openspec/project.json` 和本 Change 的 `execution.json`。不得先提交纯删除，也不得让新 Check 指向未跟踪文件。

任务 1.2 才删除 `cases/gates/`、`cases/demo/`、`cases/meshes/`、列出的旧 docs、`PROJECT_STRUCTURE.md` 和 `setup.py`。任务 1.3 清理 `.gitignore` 中旧 case/gate/bench 专属规则与注释，保留根级 `/runs/`、`/node_modules/`、`/reference/`、`/docs/references/`、`/tools/` 和仍适用的通用规则。`cases/reference_data/` 保持只读且不成为任何 Check 的 `sourcePaths` 或 `data`。

### 4. 两阶段 Check 元数据，避免 `next` 死锁

新 Checks 定为 `solver-core-api` 和 `solver-core-e2e`，但它们在 apply 前不可能注册到 main 的 `openspec/checks/`。因此随规划提交的初始 `execution.json` 只能给三个任务声明已注册的 `framework-smoke`；这是启动元数据，不是减少最终验收。

`devflow start --create` 和首次 `devflow next` 成功后，任务 1.1 在 Change worktree 中创建 `openspec/checks/solver-core-api.json` 与 `solver-core-e2e.json`，并在同一原子提交前把 `execution.json` 扩展为每个任务均声明 `framework-smoke`、`solver-core-api`、`solver-core-e2e`。后续 CLI 进程会重新加载 worktree 注册表，首次正式 `verify --task 1.1` 和后续 `next` 因而都能解析这些 id。若未创建定义就写入新 id，`next` 会以 `Task <id> names unknown checks` 失败，必须避免这种中间状态。

每个任务的最终正式验证都包含真实 API 与 E2E Checks；不能用 `framework-smoke` 替代它们。`change.json` 可提前声明影响的新 Checks，因为 `next` 不解析 impact 选择，只有验证计划才要求这些 id 已注册。

### 5. pytest 命令与 JUnit 证据

两个产品 Check 的 `command` 都直接以 `conda run -n pyfp3d-core --no-capture-output python -m pytest -p openspec.pytest_evidence ...` 调用，分别只运行 `tests/api` 或 `tests/e2e`。二者均为 Linux、`required: true`、`profiles: ["task", "affected"]`、`evidence: ["console", "junit"]`，使用不同 `executionKey`；`framework-smoke` 保持标准库脚本和 console-only，不接入 pytest。

`openspec/project.json` 设置 `evidenceDirectoryVariable: "UP3D_EVIDENCE_DIR"`。由于 devflow 不做 shell 展开，新增仅使用 pytest 钩子与标准库的 `openspec/pytest_evidence.py`：通过 `-p` 显式加载，在 `pytest_load_initial_conftests` 阶段读取该变量，并把 `xmlpath` 设为当前证据目录下的 `pytest.junit.xml`。不得在 Check JSON 中写字面量 `$UP3D_EVIDENCE_DIR`，否则它不会被展开且 devflow 找不到报告。

新测试不得使用 skip/xfail。devflow 的 JUnit 判定继续作为硬门：无报告、零 testcase、计数不一致均失败，全 skip 为必需 Check 不可接受的 skip；控制台退出码不能单独构成通过。

### 6. 包版本 provider

`environmentPackages` 规划为 `["python", "pytest", "openspec-devflow"]`。`python` 表示产品环境解释器版本；pytest 是产品环境实际测试包。`openspec/package_version_provider.py` 只导入标准库：

- 对 `python`，经 `conda run -n pyfp3d-core --no-capture-output python -c ...` 查询 `platform.python_version()`；
- 对 `pytest`，经同类 conda 子进程查询 `importlib.metadata.version()`，并显式传递 `PYTHONNOUSERSITE=1`；
- 对 `openspec-devflow`，在 devflow 进程内读取 distribution 版本及 `direct_url.json` 中 `vcs_info.commit_id` 的完整提交，组合成版本字符串，缺失时失败而不静默返回空值。

这样平台指纹不会把 devflow 的解释器或包版本误记为产品环境版本，环境变化会使旧 Run 失效。

### 7. `productPaths`、`taskLayerPaths` 与 Check `sourcePaths`

`productPaths` 保留第 1 步仍存在的条目，并加入 `tests`、`README.md`、`environment.yml`、`pyproject.toml`、`openspec/package_version_provider.py`、`openspec/pytest_evidence.py`。`pyfp3d` 路径保留但内容重建；已删除的 `setup.py`、旧 bench/case/docs 等不得加入或残留。

`taskLayerPaths` 从空列表改为仅 `["pyfp3d", "tests"]`，因为两个新 task-profile Checks 真实覆盖这些路径；不把仅保留资料或未测试路径虚报为覆盖。默认环境文件指纹继续覆盖 `environment.yml`、`pyproject.toml`、`package.json`、`package-lock.json`。

两个新 Check 的 `sourcePaths` 只列各自真实输入：`pyfp3d`、对应 `tests/api` 或 `tests/e2e`、共享 `tests/fixtures`、`environment.yml`、`pyproject.toml`、`openspec/project.json`、两个辅助模块及各自 Check 定义。E2E 另列 `README.md`，因为其公开 CLI 示例是该测试的调用契约。所有路径必须与原子提交中的 tracked 文件一致；`cases/reference_data/`、`bench/`、`phases/`、旧 case 和旧 gate 不得残留。

### 8. 测试分工

API 测试覆盖合法单/多单元结果、输出键与顺序、精确零误差，以及 malformed/缺字段/非有限数/错误形状/非法索引/退化单元等失败分支。E2E 测试通过 `subprocess` 调用真实 `python -m pyfp3d`，读取 tracked fixture，解析 stdout JSON，并验证无效路径的非零状态、stderr 与空成功 stdout。所有预期值由新 spec 和 fixture 中的简单向量直接计算，不读取旧文件。

## Risks / Trade-offs

- [原子任务范围较大] → 只允许任务 1.1 形成一个包含全部声明输入的提交；提交前用快速 API/E2E 与 smoke 检查消除可导入性和配置错误，再进入正式 verify。
- [未来 Check id 使启动或 `next` 失败] → 严格按两阶段元数据顺序，初始只引用已注册 id，定义文件与最终映射在任务 1.1 同提交落地。
- [pytest 报告路径未展开导致假失败] → 用早期 pytest 钩子读取 devflow 注入的变量，并在任务 1.1 正式 verify 中确认 Run 含有效 JUnit 文件。
- [环境 provider 查询到错误环境] → 产品包全部经命名 conda 环境子进程查询，devflow 自身才在进程内查询；缺失版本或直接 URL commit 时失败关闭。
- [精确零误差被误解为完整求解器精度声明] → 它只证明规范规定的均匀赋值，不宣称任何通量、收敛或 CFD 精度；后续能力另立 Change。
- [清理误删保留资料] → 按第 6.1 节显式路径操作并检查 diff；不批量恢复备份，不删除 `cases/reference_data/`、迁移计划或第 1 步框架文件。

## Migration Plan

1. 使用者审阅规划后，在同步干净 main 运行 `devflow start --change restart-solver-core --create`；实际 worktree 由该命令创建，分支为 `kimi/restart-solver-core`。随后进入其报告的 worktree 并执行首次 `next`。
2. 按任务 1.1 原子重建、创建 `pyfp3d-core`、注册 Checks 并扩展执行元数据；提交全部声明输入后正式 verify、complete，再提交任务状态与 Run 元数据。
3. 依次完成任务 1.2 剩余路径删除和任务 1.3 `.gitignore` 清理，每项均遵循提交前快速检查、提交后正式 verify/complete 和状态提交。
4. 全部任务后运行 affected verification，确认无过时路径和未覆盖声明；再由独立上下文只读复核。规划或任务文字在复核后的任何修改同样要求重新复核。
5. archive 须另获使用者明确批准；archive 后提交规格同步与归档结果，push、PR、merge 分别另行授权。
6. 回滚只读参考为 `/home/lrz/code/up3d_backup` 或 `pre-devflow-restart`。不得把旧目录直接恢复为当前产品；任何旧行为恢复都须另开带新 spec 的 Change。
