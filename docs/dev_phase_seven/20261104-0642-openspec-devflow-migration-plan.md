# UP3D 迁移到 openspec-devflow：清场重启方案

日期：2026-11-04（第 4 版：按 Claude 审计 B1–B5 / H1–H5 / M1–M5 / L1–L2 与同日使用者裁决重写）
状态：待执行
执行顺序：第 0 步 `收拢、备份与切换` → 第 1 步 `bootstrap-devflow` → 第 2 步 `restart-solver-core` → 第 3 步 `重新设计开发流程与计划` → 后续单能力 Change

## 1. 决策与边界

### 1.1 使用者裁决（2026-11-04）

1. **立即停止旧的 phase 7 开发方式，清场重启。**这是使用者明确授权的产品重启，不是行为保持型重构。它与 devflow 自身对存量项目的建议相反（README 第 71、91 行：先复现基线、逐个纵向切片重构），因此本方案的依据是使用者裁决，不是 devflow 的推荐路径。
2. **先完整备份，再迁移。**当前项目整体（含全部未跟踪与被忽略文件）保存为与本项目同级的 `/home/lrz/code/up3d_backup`。备份前先清理所有分支、全部合入 `main`。实际操作是**重命名**当前目录，再从**同一 GitHub 仓库重新 clone** 到原路径 `/home/lrz/code/UP3D` 开始迁移（不新建 git 仓库）。
3. **框架迁移完成后，再重新设计开发流程和计划。**旧 phase 7 的计划、轮次、门、任务台账和测试账目不继承。
4. **旧测试按需从备份取回。**重新设计时若认为有必要复现以前的某个测试或结果，从 `up3d_backup`（或 tag `pre-devflow-restart`）取回，作为一个新 Change、在新 spec 下重新推导期望值并配新 Check 引入；不照抄旧字面值，不整批恢复旧基线。

### 1.2 由此确定的边界

- 允许删除旧代码、旧测试和旧流程文件，不保留第二套 `legacy/` 活跃目录。
- 不先修旧 bug，不追求旧测试转绿，不跑旧基线、新旧一致性实验或历史性能比较。
- 不迁移旧 phase、gate、任务台账、测试账目和历史证据，不反向补写旧系统规格。
- 新实现只接受新 specs、活跃 Change 和新 Checks 的约束；不得用 skip、xfail 或 optional 掩盖必需检查失败。

### 1.3 旧裁决与旧阻塞的处置

- 第 10 轮的两项裁决——缓存新鲜度**方案 A**（Run 自证暖态）与 skip/xfail **方案 A**——在本重启下**休眠**：新产品重新引入 numba，或重新引入旧门文件时，自动恢复适用，届时再实测其两个前提。
- 第 10 轮阻塞 1（迁移机 12 红，第 11 轮 RECORDED）属于旧产品，不再阻塞迁移；需要时在备份中复现。
- `openspec/changes/adopt-devflow/` 草案作废，在第 1 步删除。
- 本裁决须在第 0 步合入 `main` 的同一批提交中，在 `docs/dev_phase_seven/progress.md` 登记一行（第 12 轮），使其在仓库内有出处，而不只存在于本文件。

### 1.4 迁移后的开发权威

迁移后的开发权威只有新 `AGENTS.md`、`openspec/specs/`、活跃 Change artifacts 和 `devflow protocol execution`；`CLAUDE.md` 只作为指向 `AGENTS.md` 的入口。旧文档、旧测试和旧证据不产生当前任务、Check 或验收条件。

但**与产品无关、经本机实测得出的机器与 git 约束不是“旧流程”**，必须迁入新 `AGENTS.md`（清单见第 5 节第 1 个任务）。

本方案只指导第 0–3 步，不成为另一套流程或进度台账；迁移完成时删除（第 8 节）。

## 2. 仓库拓扑与工具链（迁移后固定）

### 2.1 唯一主 checkout

- `/home/lrz/code/UP3D` 是**唯一**的主 checkout，始终停在干净、与 `origin/main` 一致的 `main` 上，不在其中开发或提交。
- 原因（实测）：devflow 用 `main_worktree_root()` 定位主 checkout，从任何 linked worktree 调用都解析到同一目录。旧的 `UP3D_kimi` 是 `UP3D` 的 linked worktree，在它里面运行 `start --create` 时，与 origin 的同步检查被**静默跳过**（`execution.py:649` 只在“当前目录就是主 checkout 且在集成分支上”时比较），新分支会从 kimi 分支长出。因此：
  - **`UP3D_kimi` 退役**，不再重建；
  - Claude 与 Kimi **都不在 `main` 上直接提交**；Change 一律经 `devflow start --create` 在 `runs/worktrees/<change-id>` 中开发，分支命名 Kimi 用 `kimi/<topic>`、Claude 用 `claude/<topic>`；
  - `devflow start` 只在主 checkout 中运行；
  - devflow 之外不创建 linked worktree（独立复核用的只读 worktree 见 3.3）。

### 2.2 两套 Python 环境

| 环境 | 用途 | 要求 |
|---|---|---|
| `devflow` | 运行 devflow 与 OpenSpec | Python 3.12；从 tag 非 editable 安装 `openspec-devflow`（`v0.2.0` = `cc0a9f12704f`）；`nodejs >= 22.12`；`PYTHONNOUSERSITE=1` |
| 新产品环境 | 运行新 `pyfp3d` 与其测试 | 第 2 步按新依赖清单新建，名称在该 Change 中确定；`PYTHONNOUSERSITE=1` |

- **不复用 `openspec-devflow-dev`**：它是 devflow HEAD（比 `v0.2.0` 多 17 个提交）的 editable 安装，却同样报告 `__version__ = "0.2.0"`，`devflow-install.json` 无法区分二者。安装后以 `importlib.metadata` 读取 `direct_url.json` 的 `commit_id`，必须等于 `cc0a9f12704f…`。
- **Node**：OpenSpec 1.16.0 声明 `engines.node >= 22.12.0`，系统 node 为 v20.18.0（实测），不满足。Node 装进 `devflow` 环境；`node_modules/.bin/openspec` 通过 `#!/usr/bin/env node` 取 PATH 上的 node，所以 devflow 命令**必须在激活 `devflow` 环境后运行**。
- **`PYTHONNOUSERSITE=1`** 用 `conda env config vars set` 写进两套环境：devflow 传给 Check 子进程的环境变量不含它，而它在本机是必需的（`~/.local` 的 pip 包会遮蔽 conda 环境，实测）。
- **Check 命令不依赖 PATH 上的 `python` 是哪一个**：`framework-smoke` 运行在 `devflow` 环境里（devflow 本身就在该环境中启动）；产品 Check 一律写成 `["conda", "run", "-n", "<产品环境>", "--no-capture-output", "python", ...]`。
- 旧的 `up3d`、`up3d-dev` 环境中有 editable 安装的 `pyfp3d`，指向 `/home/lrz/code/UP3D/pyfp3d`；目录切换后，它们会静默导入**新**代码。第 0 步处理（见 0.8）。

## 3. 统一的 devflow 生命周期

### 3.1 普通 Change 的进入方式

在主 checkout `/home/lrz/code/UP3D` 中同步干净 `main`，再 propose。`start --create` 前必须满足：本地 `HEAD` 与 `origin/main` 一致；除新 Change 目录 `openspec/changes/<change-id>/` 的未跟踪文件外没有其他改动；`openspec/project.json` 使用 `"integrationBranch": "main"` 和 `"integrationRemote": "origin"`。不得设为 `null` 绕过远程同步检查，不得从 linked worktree 运行 `start`，也不得清理他人工作。

```bash
conda activate devflow
cd /home/lrz/code/UP3D
devflow start --change <change-id> --create
cd <start 报告的 worktree>
devflow next --change <change-id>
```

第 1 步只用第 5 节明确的一次性初始化替代上述进入方式，之后的任务和验收没有豁免。

### 3.2 每个任务和整个 Change 的顺序

字面顺序固定为：`next` 先于实现，提交先于正式验证：

```text
devflow next --change <id> / devflow run --change <id> --task <task-id>
→ 实现当前任务并做快速检查
→ 提交当前任务及全部声明输入
→ devflow verify --change <id> --task <task-id>
→ devflow complete --change <id> --task <task-id>
→ 提交 tasks.md、execution.jsonl 和本次 Run 元数据（openspec/changes/<id>/runs/）形成的状态变更
→ 下一任务，直至全部完成
→ devflow verify --change <id> --profile affected
→ 独立复核（3.3），修正并重新验证，关闭全部发现
→ 复核方用 devflow review 记录当前内容的批准结论
→ 用户明确批准 archive
→ devflow archive --change <id> --approval "<用户批准>"
→ 提交 spec 同步（如有）、acceptance 和归档结果
→ 分别取得用户对 push、PR、merge 的明确授权
```

约束如下：

- **每个任务的必需 Check 写在 `execution.json`**（`tasks.<task-id>.checks`）。没有该文件或某任务未声明，`verify --task` 会以 `Task declares no checks` 拒绝（`execution.py:560-562`），而 `complete` 本身不检查这一点；所以“每个任务至少一个必需 Check”只有写进 `execution.json` 才成立。
- 快速检查不能代替正式 `verify --task`。Check 的 `sourcePaths`、环境文件和 Case 输入必须真实存在并被 Git 跟踪；未提交、未跟踪或被忽略的声明输入不能作为有效 Run 的依据。
- 失败或过期时修正当前内容，不借助旧基线豁免。`complete` 更新任务复选框和执行记录后必须再次提交，不能手工勾选。
- **pytest 类 Check 必须声明 `junit` 证据**：只声明 `console` 时只看退出码，测试全部 skip 时 pytest 退出码为 0，Check 会被判 pass。声明 `junit` 后，devflow 对零用例判 fail、对全 skip 判 skip（`verification.py` `_judge`），必需 Check 都过不了门。
- 根目录 `/runs/` 中的本地证据、worktree 和资源状态不提交；Change 内需要持久化的任务状态、执行记录、评审和 Run 元数据按生命周期提交。
- archive 不自动提交，也不等于 push、PR 或 merge；后三者分别授权。

### 3.3 独立复核

- 复核由**未参与实现的其他上下文**完成，在**独立的只读 worktree** 中阅读被复核的提交，例如在主 checkout 运行 `git worktree add --detach runs/review/<change-id> <commit>`，用完以 `git worktree remove` 删除。复核方不修改产品内容。
- **`devflow review` 由复核方自己执行**（在 Change worktree 中运行，这是复核方唯一的写操作）；不得由实现方代复核方记录。`--reviewer` 取一个与操作系统用户 `lrz` 不同的名字（如 `claude-reviewer`、`kimi-reviewer`）；devflow 只拒绝与执行者同名，无法验证独立性。
- **两次记录**：
  1. 有发现时，复核方先记录拒绝，发现逐条以 `"status": "open"` 登记，devflow 为其分配或保留 `id`：
     ```bash
     devflow review --change <id> --reviewer <复核者> --rejected --findings <findings.json>
     ```
  2. 实现方修正、重新验证后，复核方针对**当前**内容再次复核，逐条按原 `id` 记录 `"closed"`；确因产品改动而失效的才记 `"superseded"`（devflow 要求该发现先前为 open 且产品指纹已变）。没有开放发现时：
     ```bash
     devflow review --change <id> --reviewer <复核者> --approved --findings <findings.json>
     ```
  无发现时直接记录 `--approved`。不得省略仍开放的发现；最终复核必须针对当前产品内容且没有开放发现，否则不得请求 archive 批准。
- **复核有效性只绑定 `productPaths`**（`fingerprint.py:311`），`openspec/changes/` 不在其中。因此复核后再修改 proposal、design、spec delta 或 `tasks.md` 的任务文字，devflow **不会**使复核失效。本项目规定：**任何此类修改都视为复核失效，必须重新复核。**复核后的产品改动由 devflow 自动判为复核过期。

## 4. 第 0 步：收拢、备份与切换

本步不是 devflow Change，是一次性的仓库操作。其中 push、删除远程分支、推送 tag、移动目录**各需使用者单独授权**；每一步执行前先检查目标，执行后核对结果。

### 4.1 收拢分支（实测状态：2026-11-04）

1. 在 `kimi/openspec-devflow-restart` 上提交本方案的修订，以及 `progress.md` 的第 12 轮台账行。
2. 该分支建在 `e00466b` 上，它与 `main` 上的 `1d8e32c`（第 11 轮）树完全相同、哈希不同。先 `git rebase main` 去掉这个孪生提交，再合入 `main`。
3. `main` 当前比 `origin/main`（`12f6e19`）多 1 个提交；合入后 push `main`（需授权）。
4. 删除已合入的分支：
   - 本地：`kimi/phase7-charter`、`kimi/phase7-ibl-ap-projection`、`kimi/phase7-ibl-floor`、`kimi/phase7-xtr030`、`kimi/openspec-devflow-restart`；
   - 远程（需授权）：`origin/kimi/phase7-charter`、`origin/kimi/phase7-ibl-ap-projection`、`origin/kimi/phase7-xtr030`。
   删除前用 `git branch -a --no-merged main` 确认结果为空。
5. 在 `main` 上打 tag 并推送（需授权）：`git tag -a pre-devflow-restart -m "UP3D state before the openspec-devflow clean restart"`、`git push origin pre-devflow-restart`。这是仓库内的回滚锚点，与目录备份互补。

### 4.2 拆除 linked worktree

linked worktree 的 `.git` 指向 `UP3D/.git/worktrees/…`，目录一重命名就会断开，因此必须在重命名之前拆除：

```bash
git -C /home/lrz/code/UP3D worktree remove /home/lrz/code/UP3D_kimi
git -C /home/lrz/code/UP3D worktree remove /home/lrz/code/up3d-baseline-prov
git -C /home/lrz/code/UP3D worktree prune
git -C /home/lrz/code/UP3D worktree list        # 只剩主 checkout
```

`UP3D_kimi` 在合入后没有独有内容（`git status` 为空）；`up3d-baseline-prov` 只有 `__pycache__`，其测量结果位于 `UP3D/bench/results/`（随备份保存），第 11 轮证据已提交到 `bench/studies/baseline_provenance/`。

### 4.3 切换前核对

在 `/home/lrz/code/UP3D` 中逐项确认：当前分支是 `main`；`HEAD` 等于 `origin/main` 且等于 tag `pre-devflow-restart`；`git stash list` 为空；`git status --short` 为空；`git worktree list` 只有一行。

### 4.4 备份（重命名）

```bash
mv /home/lrz/code/UP3D /home/lrz/code/up3d_backup
git -C /home/lrz/code/up3d_backup remote set-url --push origin DISABLED
```

- 用 `mv` 而不是复制：目录共 23 GB，其中 22 GB 是被忽略的 `tools/cfl3d_work`；`mv` 瞬间完成且保留时间戳，结果就是“包括所有不追踪文件的完整副本”。
- 禁用推送：备份与新项目共用同一个 `origin`，备份中的一次误推可能覆盖新的 `main`。
- 备份即只读参考：不在其中提交，也不从其中 push。

### 4.5 备份的运行环境

1. `conda create -n up3d-backup --clone up3d`，然后在新环境中执行 `pip install --no-deps -e /home/lrz/code/up3d_backup`。
2. 从 `up3d` 和 `up3d-dev` 中 `pip uninstall pyfp3d`，避免它们静默导入新代码。
3. 在仓库目录**之外**核对：`conda run -n up3d-backup python -c "import pyfp3d; print(pyfp3d.__file__)"` 必须打印 `up3d_backup` 下的路径。
4. 备份内的已知注意点：
   - `bench/studies/baseline_provenance/run_arms.sh` 硬编码了 `/home/lrz/code/UP3D/...`，运行前必须改路径，否则会写进新项目；
   - 搬迁后第一次运行，numba 可能重编译：按旧纪律 11，先白跑一遍，再刷新或验证证据；
   - `tools/`（cfl3d、cgrid、xfoil）随备份保存；新项目需要时从备份链接或复制，不提交进仓库。

### 4.6 重新 clone 与 devflow 环境

```bash
git clone git@github.com:swayli94/UP3D.git /home/lrz/code/UP3D
cd /home/lrz/code/UP3D
git log -1 --format=%H                          # 等于 pre-devflow-restart
printf '/runs/\nnode_modules/\n' >> .git/info/exclude
```

`.git/info/exclude` 只在本地生效：第 1 步合入前，`main` 的 `.gitignore` 还不忽略 `/runs/`，而 devflow 的证据、锁和 worktree 都写在主 checkout 的 `runs/` 下，不排除的话主 checkout 会一直不干净。

然后建立 `devflow` 环境：

```bash
conda create -n devflow -c conda-forge python=3.12 "nodejs>=22.12"
conda env config vars set -n devflow PYTHONNOUSERSITE=1
conda activate devflow
python -m pip install "openspec-devflow @ git+https://github.com/swayli94/openspec-devflow@v0.2.0"
python -c "import importlib.metadata as m; print(m.distribution('openspec-devflow').read_text('direct_url.json'))"   # commit_id 必须为 cc0a9f12704f…
node --version                                   # >= v22.12
```

### 4.7 Agent 本机状态

Claude 的项目记忆按路径 `/home/lrz/code/UP3D` 存放，新项目会继承全部旧产品记忆。迁移前复制一份到备份路径对应的记忆目录；新项目只保留仍适用的机器、环境、git 类条目，其余标注“属于备份”或删除。2026-11-04 已核实：旧记忆完整保留在 `/home/lrz/.claude/projects/-home-lrz-code-up3d-backup/memory/`，当前项目只保留线程上限、devflow 环境/拓扑和用户重启裁决，本项完成。

Kimi 没有 Claude 式独立持久 `memory` 目录，无需建立 `up3d_backup` 工作区或复制会话。`~/.kimi-code/sessions/` 下的旧 `UP3D_kimi` 与当前 `UP3D` 会话原样保留，仅作为可显式恢复的历史记录，不自动注入新会话；不移动、重命名或删除 session 文件，以免破坏索引。

## 5. 第 1 步：`bootstrap-devflow` —— 最小工具接入

这一步只让框架可用：不删除旧产品代码，不接入旧测试，也不实现求解器功能。由于没有 `devflow init`，且 `next` 需要先加载项目配置和至少一个 Check（`Registry.load` 对零 Check 报 `No checks registered`），严格采用以下**一次性初始化顺序**，全部在主 checkout `/home/lrz/code/UP3D` 中发起：

1. **确认同步干净 `main`**：`git fetch origin`；`HEAD` 与 `origin/main` 一致；`git status --short` 为空。
2. **手工建立独立 branch/worktree**：本 Change 不运行 `devflow start`，其他 Change 均使用 `start --create`。在主 checkout 根目录执行，使相对路径与 `change.json` 中按主 checkout 解析的 `worktree.path` 一致：

   ```bash
   git worktree add -b kimi/bootstrap-devflow runs/worktrees/bootstrap-devflow HEAD
   cd runs/worktrees/bootstrap-devflow
   ```

3. **安装项目适配器与 OpenSpec CLI**：在 `devflow` 环境中，用项目 `package.json` 与 `package-lock.json` 固定 `@fission-ai/openspec` `1.16.0`（`npm install --save-dev --save-exact`），并在**主 checkout** 中执行 `npm ci`（`devflow archive` 先找 Change worktree，再找主 checkout 的 `node_modules/.bin/openspec`）。然后：

   ```bash
   devflow install --host claude --host kimi
   ```

4. **写最小项目配置和 `framework-smoke`**：新建 `openspec/config.yaml`（选择 `devflow` schema）和最小 `openspec/project.json`：

   ```json
   {
     "formatVersion": 1,
     "integrationBranch": "main",
     "integrationRemote": "origin",
     "evidenceRoot": "runs/evidence",
     "resourceStateRoot": "runs/devflow",
     "productPaths": [
       "pyfp3d",
       "AGENTS.md",
       "CLAUDE.md",
       "package.json",
       "package-lock.json",
       ".claude/skills",
       ".kimi-code/skills",
       "openspec/config.yaml",
       "openspec/project.json",
       "openspec/devflow-install.json",
       "openspec/schemas",
       "openspec/checks",
       "openspec/framework_smoke.py"
     ],
     "taskLayerPaths": [],
     "environmentPackages": [],
     "resources": {}
   }
   ```

   - `productPaths` 是第一步的最小完整交付清单，必须逐项存在并被 Git 跟踪。不得改成整个 `openspec/`，也不得纳入 `openspec/changes/`，否则执行记录会进入 review/archive 的产品指纹。
   - `pyfp3d` 当前真实存在，第 2 步清场后也必须在同一路径重建。
   - **`taskLayerPaths` 为空**：这一步唯一的 Check 不导入 `pyfp3d`。把 `pyfp3d` 列进去就等于虚报覆盖（devflow `REFACTORING.md:89` 明令禁止）。留空时，任何 `pyfp3d` 改动都被判为“未覆盖”，所有 Check 都会被选中，这才符合实际状态。
   - 省略 `packageVersionProvider`，不引入环境 provider；第 2 步再加（见 6.3）。

   新建并 Git 跟踪 `openspec/framework_smoke.py` 和 `openspec/checks/framework-smoke.json`：
   - 脚本只导入 Python 标准库，不导入 `pyfp3d`，不读取 `tests/`。逐项检查以下内容，任一不符即以非零退出码退出：
     - Python ≥ 3.12；
     - `openspec-devflow` 的 `direct_url.json` 中 `commit_id` 以 `cc0a9f12704f` 开头；
     - `node --version` ≥ 22.12；
     - 主 checkout 的 `node_modules/.bin/openspec --version` 为 `1.16.0`；
     - `sys.executable -m devflow install --check` 报告 `current`；
     - `list-checks` 与 `list-cases` 成功执行。
   - Check 命令为 `["python", "openspec/framework_smoke.py"]`，在 `devflow` 环境中由 devflow 启动；`required: true`，`platforms: ["Linux"]`，只声明 `console` 证据，不使用 pytest 或 JUnit。
   - `sourcePaths` 逐项覆盖：脚本、`openspec/checks`、`openspec/config.yaml`、`openspec/project.json`、`openspec/devflow-install.json`、`openspec/schemas`、`.claude/skills`、`.kimi-code/skills`。所有声明路径必须真实存在并被 Git 跟踪。`package.json` 与 `package-lock.json` 属于 devflow 默认的环境文件，已自动纳入指纹。

5. **写 tooling-only Change artifacts**：按已安装 schema 建立 `openspec/changes/bootstrap-devflow/`：
   - `.openspec.yaml` 设置 `skip_specs: true`；
   - proposal 的 New/Modified Capabilities 均为空；
   - `change.json`：`branch` 为 `kimi/bootstrap-devflow`，`worktree.path` 为 `runs/worktrees/bootstrap-devflow`，`worktree.machine` 填写真实机器名；`impact.capabilities`、`impact.cases` 为空，`impact.checks` 只含 `framework-smoke`；`requirements.platforms` 为 `["Linux"]`，`independentReview: true`；
   - `tasks.md` 至少一个任务；
   - **`execution.json`** 为每个任务声明 `checks: ["framework-smoke"]`；
   - 不创建 spec delta，不留模板占位值。

6. **首次选择任务**：上述配置、Check 和 Change artifacts 写好后再运行 `devflow next --change bootstrap-devflow`。

**首个任务的内容**（在快速检查之后、正式 `verify --task` 之前，与全部声明输入一起提交：安装生成的 schema 和 adapters、安装记录、最小配置、冒烟脚本、Check、Change artifacts）：

- **新 `AGENTS.md`**：简短，但必须包含以下与产品无关、经本机实测得出的约束：
  - **环境**：两套环境的分工与激活方式（2.2）；`PYTHONNOUSERSITE=1` 必需；`scipy >= 1.12` 这类“代码要求而非偏好”的依赖下限写进依赖清单；
  - **共享机线程上限**：上限 16，含 `NUMBA_NUM_THREADS`、`OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS`；有重活在跑时先查负载，降到 8；
  - **git 安全**：
    - 不在 `main` 上提交；push、PR、merge、删除分支、`reset --hard`、改写已推送历史，均须使用者授权；
    - `git checkout <rev> -- <paths>` 之后只能用 `git checkout HEAD -- <paths>` 还原，不能用裸的 `git checkout -- <paths>`；
    - 会改动被跟踪文件的脚本必须带 `set -euo pipefail` 和 `trap` 还原；在同一会话里跑过文件替换脚本后，不得直接 `git add -A`；
  - **进程**：启动时就记录 PID，并用 `ps -o args -p $PID` 验证它确实是目标进程；不用 `pgrep`/`pkill -f` 按模式查找或杀进程；
  - **贵重计算**：先缓存再出报告；报告代码出错只能多一列，不能丢掉已完成的计算；
  - **拓扑**：唯一主 checkout 规则（2.1）、分支命名、复核流程（3.3）；
  - **备份**：`/home/lrz/code/up3d_backup` 只读、禁推；旧测试按 1.1 第 4 条取回。
- **`CLAUDE.md`**：只含指向 `AGENTS.md` 的入口指针。
- 删除 `openspec/changes/adopt-devflow/` 旧草案。
- `.gitignore`：增加根级 `/runs/` 和 `/node_modules/`；把现有不带锚点的 `tools/`、`reference/` 改为根级的 `/tools/`、`/reference/`（原写法会匹配任意深度的同名目录）。
- `git mv` 本方案到 `docs/migration-plan.md`，避免它在第 2 步随旧阶段文档被删；迁移完成时删除（第 8 节）。

首个任务之后严格按 3.2 执行。初始化写在 `next` 之前，不构成提交、验证或完成门槛的豁免。

第 1 步不做：旧测试接入、旧 bug 修复、Numba/codegen provider、缓存包装、旧网格或外部网格 Case 适配、skip/xfail 适配、CI、roadmap HTML、历史规格迁移。

完成的前提：`devflow install --check` 为 `current`；`list-checks`/`list-cases` 正常；`framework-smoke` 真实通过；OpenSpec strict validation 和 affected verification 通过。随后记录有效的独立复核（3.3），取得用户 archive 批准，提交归档结果，再经分别授权的 push、PR、merge 合入 `main`。

## 6. 第 2 步：`restart-solver-core` —— 清场并重建最小闭环

第 1 步归档并合入后，按 3.1 从同步干净的 `main` 开始：propose 新行为，写新的 spec delta、任务与 `execution.json`，然后执行：

```bash
devflow start --change restart-solver-core --create
cd <start 报告的 worktree>
devflow next --change restart-solver-core
```

### 6.1 顶层路径去留清单

被删除的内容都可以从 `up3d_backup` 或 tag `pre-devflow-restart` 取回。

| 路径 | 处置 |
|---|---|
| `pyfp3d/` | 删除并在同一路径重建为新的最小包 |
| `tests/` | 删除并重建为新测试 |
| `bench/`、`phases/` | 删除 |
| `cases/gates/`、`cases/demo/`、`cases/meshes/` | 删除（含旧网格生成与 case 运行脚本） |
| `cases/reference_data/` | **保留**，只读资料；不自动进入 Check 输入，不构成数值真值或验收门槛 |
| `docs/dev_phase_seven/`、`docs/agent-rules.md`、`docs/overview.md`、`docs/design.md`、`docs/design_track_v.md` | 删除（理论参考与已测死的路线按 1.1 第 4 条从备份取） |
| `docs/migration-plan.md`（本方案） | 保留到迁移完成（第 8 节） |
| `PROJECT_STRUCTURE.md`、`setup.py` | 删除；依赖与打包只保留新依赖清单实际使用的文件 |
| `README.md`、`environment.yml`、`pyproject.toml` | 按新最小闭环重写；用不到的删除 |
| `LICENSE` | 保留 |
| `AGENTS.md`、`CLAUDE.md`、`openspec/`、`package*.json`、`.claude/skills`、`.kimi-code/skills` | 第 1 步产物，保留；`framework-smoke` 清场后仍须真实通过 |

### 6.2 任务

1. **原子产品重启：以下改动必须在同一个任务、同一个提交中完成，不得先提交纯删除。**
   - 删除旧 `pyfp3d/`、`tests/`、`bench/`、`phases/`，同时建立含 `__init__.py`、可导入的新最小 `pyfp3d/` 包；重写 README 与依赖清单，并按依赖清单新建产品环境（2.2）。
   - 只实现一个最小纵向闭环：从一个新的极小四面体网格入口读取输入，施加均匀无黏来流，从真实入口输出速度场和自由流保持误差。新 spec 明确输出契约和误差判据。
   - 此时才为新的 API/E2E 测试选择并固定测试框架，创建新 `tests/` 和对应的必需 Checks。测试只证明新 spec，不引用旧期望值、旧测试账目或历史性能基线。`framework-smoke` 不改用该框架。
   - 在同一提交中更新 `openspec/project.json`：
     - 把 `tests` 和实际保留的依赖清单加入 `productPaths`；
     - 承担 task-profile 的 Check **确实覆盖**的路径才加入 `taskLayerPaths`；
     - 已删除的文件或目录不得留在 `productPaths`、`taskLayerPaths` 或任何 Check 的 `sourcePaths` 中。
2. 按 6.1 删除其余旧路径。
3. 同步清理 `.gitignore`：删除已删除的 case/gate/bench 的专属忽略规则及过时注释；保留根级 `/reference/`、`/docs/references/`、`/tools/` 的本地忽略规则。原始论文和工具留在 `up3d_backup` 中，不为迁移提交 PDF 或工具，也不作为干净 worktree 的输入或验收依赖。
4. 全部任务、验证、复核记录、用户批准、archive 和归档后的提交，均按第 3 节执行。

### 6.3 产品 Check 的配置要求

- 命令经 `conda run -n <产品环境> --no-capture-output python -m pytest …` 调用，不依赖 PATH 上的 `python`。
- 声明 `junit` 证据。`project.json` 增加 `evidenceDirectoryVariable`，Check 命令把 JUnit 报告写到该变量指向的目录。
- `project.json` 增加 `environmentPackages` 与 `packageVersionProvider`（只用标准库）：
  - provider 在 devflow 进程内被调用，因此产品环境中的包版本（新核心实际导入的 numpy、scipy 等）要用子进程查询，例如 `conda run -n <产品环境> python -c "import importlib.metadata …"`；
  - `openspec-devflow` 自身的版本与 `commit_id` 可在进程内直接读取；
  - 这样 Run 会绑定本项目实测踩过的依赖版本（如 scipy 1.11→1.12 的 `rtol` 不兼容），环境一变，旧 Run 不再被复用。
- 新核心若引入 numba，同一 Change 内恢复第 10 轮的缓存方案 A（1.3），并先实测它的两个前提。

本 Change 不做旧网格或外部网格、缓存适配，也不追求跨声速、升力、Kutta、IBL、AMG、复杂边界条件或完整 CFD 功能。完成只以新最小闭环的真实 API/E2E 验证和完整 devflow 生命周期为准；归档提交后，再分别授权 push、PR、merge。

## 7. 第 3 步：重新设计开发流程与计划，再逐个增加能力

`restart-solver-core` 归档、提交并合入后：

1. **先重新设计开发流程和计划**（使用者裁决 1.1 第 3 条）。产物写进 `AGENTS.md`、`openspec/specs/` 和（如采用）`openspec/plan/`，不复活旧的 phase/round 文书体系。同时把 `openspec/specs` 加入 `productPaths`：它在第 2 步归档、规格同步之后才存在，加入后规格变更会使复核失效。
2. 之后一次只规划一个 Change。每个 Change 用一条新的或修改后的 spec、一个真实入口测试和独立 worktree，说明可观察能力及回滚边界；不批量启动 Change。
3. **取回旧测试**：设计中认为有必要复现旧测试或旧结果时，在 proposal 中引用备份来源，写明 `up3d_backup` 中的路径或 tag `pre-devflow-restart` 加上文件路径；在新 spec 下重新推导期望值，以新 Check 引入。在备份里复现旧结果时，使用 `up3d-backup` 环境（4.5）。
4. 每次都按 **同步干净 main → propose → start --create → next** 进入，再重复 3.2。不要求阅读历史会话、旧 phase 台账或旧审计。

## 8. 回滚与迁移完成判据

回滚不恢复旧流程：

- 第 0 步 4.4 之前没有破坏性操作。
- 4.4 之后，旧状态有两份：目录 `/home/lrz/code/up3d_backup`（含全部未跟踪文件），以及 `origin` 上的 tag `pre-devflow-restart`。
- 第 1 步之后、第 2 步之前，只回滚框架接入提交。
- 第 2 步之后如需读取旧实现，从备份读取，或执行 `git worktree add --detach <路径> pre-devflow-restart` 建一个只读 worktree。
- 把旧代码恢复为当前产品，必须另开带 spec delta 的行为 Change，不直接复制旧目录。

只有同时满足以下条件，才称迁移完成：

- `up3d_backup` 存在且已禁推；`conda run -n up3d-backup` 导入的 `pyfp3d` 来自备份路径；tag `pre-devflow-restart` 已在 `origin` 上。
- 最小框架已合入 `main`；干净 clone 能按 2.2 安装，并运行不依赖 pytest/JUnit 的 `framework-smoke`。
- `restart-solver-core` 已在独立 worktree 中走通完整生命周期并经授权合入，新最小闭环的 API/E2E 验证真实通过。
- 当前产品只由新 specs、新代码和新测试描述；旧红、旧文档、旧 case 和历史证据不再阻塞新开发。
- 开发流程与计划已按第 7 节重新设计，下一个 Change 能从 **同步干净 main → propose → start --create → next** 开始，无需旧基线或复杂适配。
- `docs/migration-plan.md`（本方案）在上述条件满足后的那个 Change 中删除。

迁移完成不表示求解器功能完整或没有 bug，只表示 UP3D 已能在 openspec-devflow 下持续进行可验证的新开发。
