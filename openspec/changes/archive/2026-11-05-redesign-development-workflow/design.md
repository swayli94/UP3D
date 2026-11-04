# Design

## Context

见 `proposal.md` 的 Why。当前主 specs 只有 `tetrahedron-freestream`，产品回归为 35 个 API 测试和 4 个 E2E 测试；已注册 Checks 为 `framework-smoke`、`solver-core-api`、`solver-core-e2e`。`restart-solver-core` 已归档并有 `verified: true` 的 acceptance，主 checkout 已合入其结果。

AGENTS 仍引用 `docs/migration-plan.md` 和已完成的 restart Change，`openspec/specs` 尚未进入 `productPaths`，也没有静态能力队列。已核对的 devflow 行为仍是：`next` 会校验所有任务 Check id，`start --create` 使用启动前注册表，因此新 workflow Check 必须沿用两阶段注册时序。

## Goals / Non-Goals

**Goals:**

- 把 AGENTS 转为长期开发规则，并以新的主 spec capability 固定单 Change 交付契约。
- 用一个小型 JSON 静态计划描述下一能力，以标准库 Check 验证流程内容和真实输入覆盖。
- 将一次性迁移判据与永久 Checks 分离，在 runs 内验证 clean clone 后再删除临时迁移计划。

**Non-Goals:**

- 不实现 `tetrahedron-p1-gradient` 或其他 CFD 功能，不修改现有自由流 spec/API/CLI。
- 不恢复旧 phase/round、任务台账、gate、旧代码或旧测试。
- 不修改 devflow/OpenSpec 引擎，不引入 roadmap HTML、并行 Change 调度或新的依赖包。
- 不把备份、旧环境、origin tag 或历史 Runs 变成永久 Check 输入。

## Decisions

### 1. AGENTS 的长期工作流

重写 AGENTS 的开发权威、环境和 Change 生命周期部分，保留已验证的线程、Git、进程、贵重计算和独立复核机器约束：

- 权威仅为 AGENTS、主 specs、活跃 Change artifacts 与 `devflow protocol execution`；静态计划只是队列，不是第五套行为权威。
- `devflow` 继续只运行工具链；`pyfp3d-core` 按 `environment.yml` 创建或更新，只运行当前 `pyfp3d` 和测试，核心只用标准库、pytest 只用于测试，两个环境均要求 `PYTHONNOUSERSITE=1`。
- 一次最多一个活跃 Change，入口固定为同步干净 main → propose → `start --create` → `next/run`。每 Change 恰好一个新或修改的 capability delta、真实入口测试、独立 worktree、可观察结果与回滚边界。
- 每任务遵循 `next/run → 实现与快速检查 → 提交全部声明输入 → verify --task → complete → 状态提交`；affected verification 后由独立上下文复核。archive、push、PR、merge 分别授权。
- 删除对 `docs/migration-plan.md` 和 `restart-solver-core` 的现行入口引用。回滚只读参考改为 `/home/lrz/code/up3d_backup` 或 origin 的 `pre-devflow-restart`，恢复任何旧行为仍须新 Change 和新 spec。

### 2. 静态计划采用独立 JSON，而非 devflow roadmap

新增 `openspec/plan/development-plan.json`，不创建 `roadmap.json`，避免引入阶段、lane、HTML 或并行调度语义。格式为：

```json
{
  "formatVersion": 1,
  "nextChange": {
    "id": "tetrahedron-p1-gradient",
    "capability": "tetrahedron-p1-gradient",
    "observable": "从节点势函数计算每个四面体单元的常梯度",
    "entry": "公共 run_p1_gradient(path) 与 python -m pyfp3d p1-gradient <input.json>",
    "input": "vertices、tetrahedra 和与顶点一一对应的有限 nodal_potential",
    "output": "按单元输入顺序排列的 N×3 gradients",
    "tests": "仿射场精确梯度、手工推导的非仿射节点值、输入错误及真实 API/E2E",
    "rollback": "只移除新势函数输入、梯度入口、spec、测试和 Checks，保留现有自由流 API/CLI 与回归"
  },
  "backlog": [
    {"title": "评估后续单元场重建能力"}
  ]
}
```

`nextChange` 允许为 `null` 或一个完整对象，不允许数组或多项并行规划。详细项完成后，完成它的 Change 在最终验证前将其置为 `null` 或由后续明确规划替换；backlog 项只允许 `title`，不得附带 proposal、design、tasks、spec delta 或预先注册的 Checks。

下一能力的具体边界：从线性四面体节点势函数计算单元梯度，优先复用当前精确有理数验证思路。仿射测试可从独立选择的线性函数推导期望梯度，非仿射节点值按 P1 公式独立手算；不引用旧测试或旧求解器。真实入口是新的公共 API 与 CLI 子命令，现有 `python -m pyfp3d <input>` 自由流入口保持不变。未来新 Checks 同样先完成注册再扩展 execution，不在本 Change 创建其 artifacts。

### 3. 标准库 `workflow-contract` Check

新增 `openspec/workflow_check.py`，命令为 `python openspec/workflow_check.py`，console-only、Linux、required、`task`/`affected` profiles，capability 为 `development-workflow`。脚本只读取静态、可指纹化内容并验证：

- AGENTS 含当前权威、`pyfp3d-core`、单 Change、标准入口、真实入口测试、回滚、独立复核及分级审批的稳定语义锚点；不含 `docs/migration-plan.md`、`restart-solver-core` 或旧 phase 路径的现行引用。
- 计划 `formatVersion` 正确，`nextChange` 为 `null` 或具有全部必需字段的单个对象，id/capability 为 kebab-case，backlog 每项只有非空 `title`。
- `openspec/project.json` 的 `productPaths` 包含 `openspec/specs`、`openspec/plan`、`openspec/workflow_check.py`；所有 project/Check 声明路径存在。
- 永久 Check 的 command、data、sourcePaths 不含 `/home/lrz/code/up3d_backup`、`up3d-backup`、`pre-devflow-restart` 或对历史 Runs 的依赖。
- `docs/migration-plan.md` 是否存在不属于该永久 Check 的判定条件，避免收尾前任务被人为制造失败。

Check 的 `sourcePaths` 为 AGENTS、`openspec/specs`、`openspec/plan`、`openspec/workflow_check.py`、`openspec/project.json` 和 `openspec/checks`。脚本不读取活跃 Change 目录，避免任务状态、Runs 和执行日志成为未声明或持续变化的隐藏输入。

备选方案是只用 pytest 产品回归或把一次性迁移判据塞入 framework smoke；前者不能验证流程，后者会让新机器永久依赖备份，均不采用。

### 4. 指纹与 Check 注册时序

`productPaths` 增加 `openspec/specs`、`openspec/plan` 和 `openspec/workflow_check.py`。`taskLayerPaths` 保持真实的 `pyfp3d`、`tests`，不虚报文档覆盖；workflow Check 自身的 `sourcePaths` 负责流程文件 affected 选择。

初始 `execution.json` 只能给任务 1.1、1.2 声明现有三个 Checks，确保 `start --create` 和首次 `next` 可解析。首次 next 后，任务 1.1 创建 `openspec/checks/workflow-contract.json`，再把两个任务都扩展为：

- `framework-smoke`
- `solver-core-api`
- `solver-core-e2e`
- `workflow-contract`

定义和扩展元数据在任务 1.1 同一提交中落地。每个任务的最终正式 verify 都运行四项 Checks，不用现有三项代替新流程验证。

### 5. 一次性迁移验收与完整证据链

任务 1.2 在其开始时已具有任务 1.1 的提交和状态提交。`migration-acceptance.md` 必须按迁移计划第 8 节逐项列出判据、命令、退出码、关键 stdout 和对应提交，不能只写结论；所有下列项通过前保留 `docs/migration-plan.md`。这些命令和记录只存在于当前 Change，不注册为永久 Check。

**A. 合入、远程同步与历史 acceptance。** 先取得远程引用，再证明两个合入提交均为当前 `origin/main` 祖先；记录 `origin/main` 的完整 oid 和两个 merge-base 的零退出码：

```bash
git fetch origin
git rev-parse origin/main
git merge-base --is-ancestor b1db684d10fee69c4c2c3c6e93252eef65baa898 origin/main
git merge-base --is-ancestor 5fbed29a5cb6d0db31d03709c944cdc1530715c4 origin/main
```

对 `bootstrap-devflow` 和 `restart-solver-core` 分别读取归档 acceptance，要求 change 名匹配、`verified is True`、至少一个 Check 且全部 pass，并记录 acceptance commit：

```bash
PYTHONNOUSERSITE=1 conda run -n devflow --no-capture-output python -c 'import json; from pathlib import Path
for name in ("bootstrap-devflow", "restart-solver-core"):
    path=Path(f"openspec/changes/archive/2026-11-04-{name}/acceptance.json")
    a=json.loads(path.read_text(encoding="utf-8"))
    assert a["change"] == name and a["verified"] is True
    assert a["checks"] and all(c["status"] == "pass" for c in a["checks"])
    print(name, a["commit"], [(c["check"], c["status"]) for c in a["checks"]])'
```

**B. 备份和 tag。** 备份导入必须在 `/home/lrz/code` 运行，使当前 worktree 不作为 `sys.path` 的项目根；同时禁止 user site 和字节码写入。断言不是模糊的父路径包含，而是解析路径精确等于备份入口：

```bash
test -d /home/lrz/code/up3d_backup
test "$(git -C /home/lrz/code/up3d_backup remote get-url --push origin)" = DISABLED
(
  cd /home/lrz/code
  PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 conda run -n up3d-backup --no-capture-output python -c 'from pathlib import Path; import pyfp3d; actual=Path(pyfp3d.__file__).resolve(); expected=Path("/home/lrz/code/up3d_backup/pyfp3d/__init__.py").resolve(); assert actual == expected, (actual, expected); print(actual)'
)
git ls-remote --exit-code --tags origin refs/tags/pre-devflow-restart
```

不得修改备份；该命令只读取导入路径，两个 Python 环境变量必须同时出现在实际执行记录中。

**C. 当前产品、主 spec 与回归。** 记录现有主 spec 和当前 Change 的独立 strict 结果，主 spec 命令不得被 Change delta validation 替代：

```bash
PYTHONNOUSERSITE=1 conda run -n devflow --no-capture-output node_modules/.bin/openspec validate tetrahedron-freestream --type spec --strict
PYTHONNOUSERSITE=1 conda run -n devflow --no-capture-output node_modules/.bin/openspec validate redesign-development-workflow --strict
```

用 `test ! -e` 分别确认 `bench`、`phases`、`cases/gates`、`cases/demo`、`cases/meshes`、`docs/dev_phase_seven` 和其他已清理旧路径不存在，用 `test -f` 确认主 spec、新 `pyfp3d` 和新 API/E2E 文件存在；`cases/reference_data` 必须存在且 `git diff --exit-code pre-devflow-restart -- cases/reference_data` 为零。另用 project/Check 路径检查证明永久声明不引用旧 case/gate/bench 或历史 Runs。

以现有 Check 命令运行 35 个 API 和 4 个 E2E，保存 JUnit 并记录 `tests=35`、`tests=4`、零 failed、零 skipped；运行 framework smoke 和新增 workflow Check。旧测试不运行，也不作为阻塞或验收依据。

**D. 当前提交的 clean clone 与全新工具链 smoke。** D–E 在同一带 `set -euo pipefail` 的验收 shell 中执行，先记录 `worktree_root="$PWD"` 并设置返回该目录后清理临时目录的 `trap`。从任务 1.1 的已提交 HEAD 创建独立 clone，在 clone 内新建 prefix devflow 环境，而不是借用现有 devflow 环境：

```bash
test ! -e runs/migration-clean-clone
git clone --no-hardlinks . runs/migration-clean-clone
clone_root="$worktree_root/runs/migration-clean-clone"
clean_devflow="$clone_root/.runtime-devflow"
mkdir -p "$clone_root/.cache/conda-pkgs" "$clone_root/.cache/pip" "$clone_root/.cache/tmp" "$clone_root/.cache/npm"
CONDA_PKGS_DIRS="$clone_root/.cache/conda-pkgs" conda create --prefix "$clean_devflow" -c conda-forge "python=3.12" "nodejs>=22.12"
conda env config vars set --prefix "$clean_devflow" PYTHONNOUSERSITE=1
TMPDIR="$clone_root/.cache/tmp" PIP_CACHE_DIR="$clone_root/.cache/pip" PYTHONNOUSERSITE=1 conda run --prefix "$clean_devflow" --no-capture-output python -m pip install "openspec-devflow @ git+https://github.com/swayli94/openspec-devflow@v0.2.0"
conda run --prefix "$clean_devflow" --no-capture-output npm ci --prefix "$clone_root" --cache "$clone_root/.cache/npm"
cd "$clone_root"
PYTHONNOUSERSITE=1 conda run --prefix "$clean_devflow" --no-capture-output python openspec/framework_smoke.py
PYTHONNOUSERSITE=1 conda run --prefix "$clean_devflow" --no-capture-output python openspec/workflow_check.py
```

这样同时验证第 2.2 节的 Python/Node/devflow 安装和 clone 内 Node 依赖。conda packages、pip 下载、编译临时目录、npm cache、环境和 clone 全部位于 runs；不得创建或修改命名 conda 环境。`framework_smoke.py` 必须从 clone 自身定位 main checkout，记录中的 OpenSpec 可执行文件和 Python 必须分别位于该 clone 的 `node_modules` 与 `.runtime-devflow`，不得指向主 checkout 或既有 devflow 环境。workflow Check 也必须通过。

**E. 标准入口的真实 startup smoke。** 返回 `worktree_root`，另建一个从实际 origin URL 的 main 克隆，不能直接使用含当前活跃 redesign Change 的任务 1.1 clone。startup 使用 D 创建的隔离工具链：

```bash
cd "$worktree_root"
remote_url=$(git remote get-url origin)
git clone --no-hardlinks --branch main --single-branch "$remote_url" runs/migration-startup-smoke
startup_root="$worktree_root/runs/migration-startup-smoke"
conda run --prefix "$clean_devflow" --no-capture-output npm ci --prefix "$startup_root" --cache "$clone_root/.cache/npm"
cd "$startup_root"
test "$(git rev-parse HEAD)" = "$(git rev-parse origin/main)"
git status --short
PYTHONNOUSERSITE=1 conda run --prefix "$clean_devflow" --no-capture-output "$startup_root/node_modules/.bin/openspec" list --json
```

要求 `git status --short` 在 propose 前为空，OpenSpec `changes` 为空。随后用同一隔离工具链执行 `openspec new change migration-startup-smoke` 完成真实 propose，为合成 tooling-only Change 写最小工件：`.openspec.yaml` 设 `skip_specs: true`，proposal 的 capabilities 为空，`change.json` 声明 `kimi/migration-startup-smoke` 和 `runs/worktrees/migration-startup-smoke`，一个任务 1.1，execution 只引用当时已注册 Checks。strict validation 后，在 clone 主 checkout 运行：

```bash
PYTHONNOUSERSITE=1 conda run --prefix "$clean_devflow" --no-capture-output devflow start --change migration-startup-smoke --create
cd "$startup_root/runs/worktrees/migration-startup-smoke"
PYTHONNOUSERSITE=1 conda run --prefix "$clean_devflow" --no-capture-output devflow next --change migration-startup-smoke
```

记录 start 返回的独立 worktree/branch、next 的 `state=ready`、任务 id 和 Checks。start 自身会强制核对 main 同步和“只有该新 Change 未跟踪目录”，因此以上输出同时证明标准入口的真实约束生效。

此 smoke 仓库没有任何其他活跃 Change，合成项是唯一 Change；它不是 p1-gradient，不实现产品、不运行 verify/complete、不 push，也不改变主仓库分支。它证明的是 propose/start/next 的真实工具行为，不证明 p1 的 spec 内容；p1 仍由后续独立 propose 定义。记录后返回 `worktree_root`，删除两个临时 clone、隔离工具链和全部临时 cache。

只有 A–E、39 项回归、两份 strict 和 smoke 全部通过，任务 1.2 才能在同一提交中写入完整验收记录并删除 `docs/migration-plan.md`。任一项失败都保留迁移计划、保留失败输出并修复当前内容，不放宽判据。

### 6. 任务划分

- 任务 1.1：更新 AGENTS，建立静态计划、workflow 脚本和 Check，更新 productPaths 并扩展 execution；以四个 Checks 和 strict validation 验证。
- 任务 1.2：先完成 A–E 一次性证据、39 项回归和两类临时 clone 验证，再写验收记录并删除迁移计划；提交后以四个 Checks 正式 verify。

任务 1.1 的提交是当前 clean clone 的确定输入；startup smoke 使用独立的 origin/main clone 避免与当前 redesign 冲突；任务 1.2 是唯一允许删除迁移文档的收尾提交。

## Risks / Trade-offs

- [静态计划变成另一套繁重文书] → 只保留一个 JSON、至多一个详细 nextChange 和标题 backlog；workflow Check 拒绝预生成工件。
- [流程 Check 依赖动态 Change 状态而频繁过期] → 只读取静态指纹文件，活跃 Change 生命周期继续由 devflow 负责。
- [备份导入误读当前 worktree] → 在 `/home/lrz/code` 执行并设置 `PYTHONDONTWRITEBYTECODE=1`、`PYTHONNOUSERSITE=1`，断言精确 `__init__.py` 路径。
- [clean clone 或 startup smoke 误用主 checkout 状态] → 两类 clone 分离，使用 runs 内全新 prefix 工具链和本地 Node 依赖，记录 HEAD/origin、Python/OpenSpec 路径、活跃 Change 列表和 start/next JSON。
- [一次性判据污染永久验证] → 永久 Check 扫描并拒绝备份/tag/历史 Run 依赖，一次性输出只存于本 Change 的验收记录。
- [提前删除迁移计划导致判据丢失] → 删除只位于任务 1.2，且 A–E 和全部回归必须先通过并逐项记录。
- [下一能力计划过细而被误认为本 Change 范围] → synthetic smoke 与 p1-gradient 明确分离，本 Change 不创建其产品代码、测试或 Change artifacts。

## Migration Plan

1. 使用者审阅后在同步干净 main 运行 `devflow start --change redesign-development-workflow --create`，进入其生成的 worktree 并执行首次 next。
2. 按任务 1.1 建立长期流程、计划和 workflow Check，提交后正式 verify/complete 并提交状态。
3. 按任务 1.2 完成 A–E、39 项回归、主 spec/Change strict 和两类临时 clone 验证；全部通过后才提交验收记录和迁移计划删除，再正式 verify/complete。
4. 运行 affected verification 和独立复核，规划或任务文字变化同样使复核失效。archive、push、PR、merge 分别取得授权。
5. 后续只从静态计划选择一项进行独立 propose。回滚读取只读备份或 `pre-devflow-restart`；不得直接恢复旧流程，任何旧行为恢复均另开带 spec 的 Change。
