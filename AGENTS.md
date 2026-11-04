# UP3D Agent 规则

开发权威只有本文件、`openspec/specs/`、活跃 Change artifacts 和 `devflow protocol execution`。旧文档、旧测试和旧证据不产生当前任务或验收条件；迁移期间的边界见 `docs/migration-plan.md`。

## 环境

- `devflow` 环境只运行 devflow 与 OpenSpec：Python 3.12，从 tag 非 editable 安装 `openspec-devflow` v0.2.0（commit `cc0a9f12704fdce38ae98d5750c32b338207790f`），Node >= 22.12。先 `conda activate devflow` 再运行命令，或使用 `conda run -n devflow --no-capture-output ...`。
- 新产品环境只运行新 `pyfp3d` 及其测试，在 `restart-solver-core` 中按新依赖清单创建并命名。产品 Check 使用 `conda run -n <产品环境> --no-capture-output ...`，不依赖 PATH 中偶然的 `python`。
- 两套环境都必须设置 `PYTHONNOUSERSITE=1`，避免用户 site-packages 遮蔽环境依赖。`scipy >= 1.12` 这类代码要求而非偏好的下限必须写进实际依赖清单。

## 共享机线程

线程上限为 16，包括 `NUMBA_NUM_THREADS`、`OMP_NUM_THREADS` 和 `OPENBLAS_NUM_THREADS`。启动重活前先看机器负载；已有重活时降到 8，并让所有相关线程变量保持一致。

## Git 安全

- 不在 `main` 上提交。push、PR、merge、删除分支、`reset --hard`、改写已推送历史均须分别取得使用者授权；不得 force-push 或改动他人的分支。
- `git checkout <rev> -- <paths>` 会同时改 index 和 worktree；之后只能用 `git checkout HEAD -- <paths>` 还原，不能用裸的 `git checkout -- <paths>`。
- 会改动被跟踪文件的脚本必须使用 `set -euo pipefail`，并用 `trap` 还原。同一会话跑过文件替换脚本后不得直接 `git add -A`；先核对状态、差异和能识别正确内容的 sentinel。

## 进程

启动时记录 PID，并立即用 `ps -o args -p "$PID"` 验证它确实是目标进程；只按已验证的 PID 管理进程。不得用 `pgrep` 或 `pkill -f` 按命令行模式查找或杀进程。

## 贵重计算

先缓存计算结果，再生成报告。报告代码出错只能多一列，不能丢掉已完成的计算；优先增量写出原始结果。不重算贵重已提交 artifact 来替代读取现有证据。

## Checkout 拓扑与生命周期

- `/home/lrz/code/UP3D` 是唯一主 checkout，始终保持干净、与 `origin/main` 一致的 `main`，不在其中开发或提交。`UP3D_kimi` 退役，不再重建。
- Change 在 `runs/worktrees/<change-id>` 开发；Kimi 分支用 `kimi/<topic>`，Claude 分支用 `claude/<topic>`。除一次性 bootstrap 外，只在主 checkout 运行 `devflow start --create` 创建 worktree，不在 devflow 之外创建开发用 linked worktree。
- `openspec/project.json` 保持 `integrationBranch: main`、`integrationRemote: origin`，不得改为 `null` 绕过同步检查，也不得清理他人工作。
- 顺序固定为：`next`/`run` → 实现与快速检查 → 提交任务和全部声明输入 → `verify --task` → `complete` → 提交任务状态与 Run 元数据 → 下一任务 → affected verification → 独立复核 → 用户批准 archive → archive 与归档提交。快速检查不代替正式 verify；不得手工勾任务。push、PR、merge 之后分别授权。

## 独立复核

- 复核必须由未参与实现的其他上下文完成。在主 checkout 用 `git worktree add --detach runs/review/<change-id> <commit>` 建独立只读 worktree，阅读被复核提交，用完 `git worktree remove`；复核方不得在其中修改产品。
- `devflow review` 由复核方本人在 Change worktree 执行，这是复核方唯一的写操作；`--reviewer` 必须与执行者及操作系统用户 `lrz` 不同。
- 发现逐条记录为 `open`；修正并重新验证后，由复核方针对当前内容按原 id 记录 `closed`，只有发现先前为 open 且产品指纹已变时才可记 `superseded`。最终批准必须针对当前内容且没有开放发现。
- devflow 的复核指纹只绑定 `productPaths`；本项目额外规定，复核后任何 proposal、design、spec delta 或 `tasks.md` 任务文字的修改同样使复核失效，必须重新复核。

## 备份与旧测试

`/home/lrz/code/up3d_backup` 是只读备份，禁止修改、提交或 push；其远端推送保持禁用。只有重新设计确有需要时，才在一个新 Change 中引用备份路径或 `pre-devflow-restart` tag 取回旧测试或旧结果：必须在新 spec 下重新推导期望值并配置新 Check，不照抄旧字面值，不整批恢复旧基线。需要在备份复现时使用 `up3d-backup` 环境。
