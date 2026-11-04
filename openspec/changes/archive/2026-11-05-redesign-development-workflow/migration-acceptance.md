# 迁移完成一次性验收记录

- 执行时间：2026-11-04T23:46:34+08:00
- 验收执行提交：`6f8d1101a60c948e4b11b4338e8f8f91af207b3d`
- 执行分支：`kimi/redesign-development-workflow`
- 结论：A–E 全部真实命令通过，随后才写入本文件并删除 `docs/migration-plan.md`。
- 可读性变量：`WT=/home/lrz/code/UP3D/runs/worktrees/redesign-development-workflow`、`CLONE=$WT/runs/migration-clean-clone`、`STARTUP=$WT/runs/migration-startup-smoke`、`CLEAN=$CLONE/.runtime-devflow`。实际日志使用展开后的绝对路径。

首次脚本封装在 A 的 acceptance 调用前因 `PYTHONNOUSERSITE=1: command not found` 停止，未形成通过结论；修正为显式 `env` 后从 A 重新执行。本文件只接受第二次完整执行的证据，没有跳过任何失败判据。

## A. 远程同步、合入祖先与 archive acceptance

以下每条命令退出码均为 `0`：

```bash
gh repo view swayli94/UP3D --json url --jq .url
# stdout: https://github.com/swayli94/UP3D

gh api repos/swayli94/UP3D/git/ref/heads/main --jq .object.sha
# stdout: 5fbed29a5cb6d0db31d03709c944cdc1530715c4

git clone --no-hardlinks --branch main --single-branch https://github.com/swayli94/UP3D.git "$STARTUP"
git -C "$STARTUP" rev-parse origin/main
# stdout: 5fbed29a5cb6d0db31d03709c944cdc1530715c4

git -C "$STARTUP" merge-base --is-ancestor b1db684d10fee69c4c2c3c6e93252eef65baa898 origin/main
git -C "$STARTUP" merge-base --is-ancestor 5fbed29a5cb6d0db31d03709c944cdc1530715c4 origin/main
```

远程证据来自 GitHub HTTPS API 与公共 HTTPS clone，不使用 SSH、`git fetch` 或全局 Git 配置。两个 merge-base 都没有 stdout，零退出码分别证明 bootstrap 合入提交 `b1db684` 和 restart 合入提交 `5fbed29` 是实际远程 `origin/main` 的祖先。

分别读取两份 acceptance 的 Python 命令退出码为 `0`，关键 stdout 为：

```text
bootstrap-devflow fcb2cfee057476ba1a9c35c1ba21d5adfa93d98a [('framework-smoke', 'pass')]
restart-solver-core 2165f8a65a19ec63e9df40e7241d0c36812d7ded [('framework-smoke', 'pass'), ('solver-core-api', 'pass'), ('solver-core-e2e', 'pass')]
```

断言同时要求 change 名匹配、`verified is True`、至少一个 Check 且全部 pass。它们证明两份归档验收本身有效；merge-base 则独立证明其合入祖先关系。

## B. 只读备份、精确导入与 origin tag

以下命令退出码均为 `0`：

```bash
test -d /home/lrz/code/up3d_backup
test "$(git -C /home/lrz/code/up3d_backup remote get-url --push origin)" = DISABLED
```

备份导入在 `/home/lrz/code`（不是任何 pyfp3d 项目根）执行，实际命令和输出为：

```bash
cd /home/lrz/code
PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 conda run -n up3d-backup --no-capture-output python -c 'from pathlib import Path; import pyfp3d; actual=Path(pyfp3d.__file__).resolve(); expected=Path("/home/lrz/code/up3d_backup/pyfp3d/__init__.py").resolve(); assert actual == expected, (actual, expected); print(actual)'
# exit=0
# stdout: /home/lrz/code/up3d_backup/pyfp3d/__init__.py
```

tag 的远程核验使用 gh HTTPS API，退出码为 `0`：

```bash
gh api repos/swayli94/UP3D/git/ref/tags/pre-devflow-restart --jq .ref
# stdout: refs/tags/pre-devflow-restart
```

备份只被读取；执行前后其 Git 状态均为 `## main...origin/main`，无工作树改动。

## C. 清场、保留路径、主 spec、39 项回归与永久 Checks

逐个 `test ! -e` 的退出码均为 `0`：`bench`、`phases`、`cases/gates`、`cases/demo`、`cases/meshes`、`docs/dev_phase_seven`、`docs/agent-rules.md`、`docs/overview.md`、`docs/design.md`、`docs/design_track_v.md`、`PROJECT_STRUCTURE.md`、`setup.py`。

逐个 `test -e` 的退出码均为 `0`：`openspec/specs/tetrahedron-freestream/spec.md`、`pyfp3d/__init__.py`、`pyfp3d/freestream.py`、`tests/api/test_freestream.py`、`tests/e2e/test_cli.py`、`cases/reference_data`。

```bash
git diff --exit-code pre-devflow-restart -- cases/reference_data
# exit=0，无 diff
```

永久声明检查命令退出码为 `0`，其断言包括所有声明路径存在，且永久 Check 的 command/data/sourcePaths 不含旧 case/gate/bench、活跃或归档 Change、历史 evidence、备份环境或 tag。关键 stdout：

```text
declared paths exist: 27
portable checks: ['framework-smoke', 'solver-core-api', 'solver-core-e2e', 'workflow-contract']
```

worktree 起初没有 node_modules，使用声明 lockfile 安装，cache 位于 `$WT/runs/npm-cache`：

```bash
conda run -n devflow --no-capture-output npm ci --cache "$WT/runs/npm-cache"
# exit=0；added 57 packages, 0 vulnerabilities
```

主 spec 与当前 Change 分别验证，不能互相替代，退出码均为 `0`：

```bash
"$WT/node_modules/.bin/openspec" validate tetrahedron-freestream --type spec --strict
# stdout: Specification 'tetrahedron-freestream' is valid

"$WT/node_modules/.bin/openspec" validate redesign-development-workflow --strict
# stdout: Change 'redesign-development-workflow' is valid
```

实际 pytest 命令经 `conda run -n pyfp3d-core --no-capture-output python -m pytest -p openspec.pytest_evidence` 分别运行 `tests/api` 和 `tests/e2e`，并设置 `UP3D_EVIDENCE_DIR` 到 runs 临时目录。退出码均为 `0`：

```text
API: 35 passed
E2E: 4 passed
api {'tests': 35, 'skipped': 0, 'failed': 0, 'reports': ['pytest.junit.xml']}
e2e {'tests': 4, 'skipped': 0, 'failed': 0, 'reports': ['pytest.junit.xml']}
```

worktree framework smoke 与 workflow Check 退出码均为 `0`。framework smoke 的 Python 为 3.12.16、Node 为 v26.13.0、devflow commit 为 `cc0a9f12704fdce38ae98d5750c32b338207790f`，识别全部四个 Checks；workflow Check 通过 14 个 AGENTS 锚点、单个 nextChange 和四个可移植 Checks。worktree smoke 按其实现读取主 checkout 的 OpenSpec，D 另行证明 clone 自含路径。

## D. 当前提交的 clean clone 与全新工具链

以下命令退出码均为 `0`：

```bash
git clone --no-hardlinks . "$CLONE"
git -C "$CLONE" rev-parse HEAD
# stdout: 6f8d1101a60c948e4b11b4338e8f8f91af207b3d

mkdir -p "$CLONE/.cache/conda-pkgs" "$CLONE/.cache/pip" "$CLONE/.cache/tmp" "$CLONE/.cache/npm"
CONDA_PKGS_DIRS="$CLONE/.cache/conda-pkgs" conda create --prefix "$CLEAN" -c conda-forge "python=3.12" "nodejs>=22.12"
conda env config vars set --prefix "$CLEAN" PYTHONNOUSERSITE=1
TMPDIR="$CLONE/.cache/tmp" PIP_CACHE_DIR="$CLONE/.cache/pip" GIT_TERMINAL_PROMPT=0 PYTHONNOUSERSITE=1 conda run --prefix "$CLEAN" --no-capture-output python -m pip install "openspec-devflow @ git+https://github.com/swayli94/openspec-devflow@v0.2.0"
conda run --prefix "$CLEAN" --no-capture-output npm ci --prefix "$CLONE" --cache "$CLONE/.cache/npm"
```

新 prefix 环境安装 Python 3.12.16、Node v26.13.0、devflow 0.2.0；npm 安装 57 个包，全部缓存和工具链均位于 runs 临时 clone 内。pip 安装、npm 安装、clone 创建和 framework/workflow 均没有借用既有 devflow 环境。

路径与版本证据：

```text
Python: $CLONE/.runtime-devflow/bin/python
OpenSpec: 1.16.0 at $CLONE/node_modules/.bin/openspec
Node: v26.13.0
devflow commit: cc0a9f12704fdce38ae98d5750c32b338207790f
```

在 `$CLONE` 内以 `$CLEAN` 执行 `python openspec/framework_smoke.py` 和 `python openspec/workflow_check.py`，退出码均为 `0`。framework smoke 明确报告上述 clone 内 Python/OpenSpec 路径并识别四个 Checks；workflow Check 也在 clone 内通过。

## E. 无其他活跃 Change 的真实 startup smoke

startup clone 来自 A 的独立 origin/main HTTPS clone，HEAD 与 origin/main 均为 `5fbed29a5cb6d0db31d03709c944cdc1530715c4`，不是包含当前活跃 redesign Change 的 D clone。以下命令退出码均为 `0`：

```bash
conda run --prefix "$CLEAN" --no-capture-output npm ci --prefix "$STARTUP" --cache "$CLONE/.cache/npm"
cd "$STARTUP"
test "$(git rev-parse HEAD)" = "$(git rev-parse origin/main)"
git status --short
# stdout 为空
test -z "$(git status --porcelain)"
"$STARTUP/node_modules/.bin/openspec" list --json
# changes: []
```

随后真实执行 `"$STARTUP/node_modules/.bin/openspec" new change migration-startup-smoke`，创建唯一合成 tooling-only Change。其 `.openspec.yaml` 设置 `skip_specs: true`，proposal 无 capabilities，声明分支 `kimi/migration-startup-smoke`、worktree `runs/worktrees/migration-startup-smoke`，一个任务 1.1，execution 使用当时已注册的三个 Checks。没有 p1-gradient 或任何产品实现。

```text
openspec validate migration-startup-smoke --strict: exit=0
Change 'migration-startup-smoke' is valid
```

真实 `devflow start --change migration-startup-smoke --create` 退出码为 `0`，返回 `state=ready` 和独立 worktree：

```text
$STARTUP/runs/worktrees/migration-startup-smoke
```

对该 JSON 的解析断言也通过。实际分支读取和断言退出码均为 `0`：

```text
kimi/migration-startup-smoke
```

在生成的 worktree 中真实执行 `devflow next --change migration-startup-smoke`，退出码为 `0`，关键结果为：

```text
state=ready
task=1.1
attempts=0
checks=['framework-smoke', 'solver-core-api', 'solver-core-e2e']
waiting=[]
```

证据边界：该 synthetic Change 是隔离 clone 中唯一活跃 Change，只证明 propose/strict/start/next 的真实工具行为；未执行其产品实现、verify/complete 或 push，也未改变主仓库分支。

## 清理与最终允许删除迁移计划的条件

D–E 的 `set -euo pipefail` 脚本通过 EXIT trap 返回 `$WT`，并删除 `$CLONE`、`$STARTUP`、prefix 工具链、全部 cache、JUnit 和本次创建的 worktree node_modules。实际核对 `test ! -e` 这些路径均为 `0`；主 checkout 和备份的 `git status --short --branch` 均保持 `## main...origin/main`。

A–E 在上述真实结果全部存在后，才在同一批未提交改动中新增本文件并删除 `docs/migration-plan.md`。一次性外部依赖只记录于本 Change，没有写入永久 Check、AGENTS 或静态计划。
