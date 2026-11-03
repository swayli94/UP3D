# Tasks

> 前置条件（不在本 Change 内，未满足前不 `start`）：审计阻塞 1 —— 迁移机 12 红逐项分类与处置。

## 1. 工具环境

- [ ] 1.1 新增工具环境文件（Python 3.12、Node ≥ 22.12、devflow 固定 `cc0a9f1`）；验证：新建环境后 `devflow --help` 可用、`python -c "import devflow; print(devflow.__version__)"` 输出 `0.2.0`，且 `conda run -n up3d python -c "import sys; print(sys.version)"` 仍为 3.11（求解环境未动）。
- [ ] 1.2 新增 `package.json` + `package-lock.json`，OpenSpec CLI 精确 `1.16.0`；验证：`npm ci` 干净完成且 `npx openspec --version` 输出 `1.16.0`。

## 2. 安装与项目配置

- [ ] 2.1 `devflow install --host claude --host kimi`，提交 schema、适配器与安装记录；验证：`devflow install --check` 报 `current`。
- [ ] 2.2 新增 `openspec/config.yaml` 与 `openspec/project.json`（键按审计 §7.2：integration 分支/远端、`devflowCommand`、`evidenceRoot`、`resourceStateRoot`、`dataEnvironmentVariable=UP3D_DATA_DIR`、`evidenceDirectoryVariable`、`stepTraceVariable`、`environmentFiles`、`environmentPackages`、`productPaths`、`resources`）；验证：`devflow list-checks`、`devflow list-cases` 无错误，`openspec validate adopt-devflow --strict` 通过。
- [ ] 2.3 `.gitignore` 增加 `runs/`、`node_modules/`；验证：`git check-ignore runs/x node_modules/x` 均命中。

## 3. 出处：环境 provider 与网格

- [ ] 3.1 `tools/devflow/` 环境 provider：以子进程查询 `up3d`，返回 D4 所列包版本与 `numba-codegen` 伪包；验证：仅改 `NUMBA_CPU_FEATURES`（或设 `NUMBA_CPU_NAME=generic`）时伪包值改变，且一个已记录的 Run 被判“Run environment differs from this machine”。
- [ ] 3.2 被忽略网格的 Case raw 定义（外部 `UP3D_DATA_DIR`），wrapper 以已提交 `.msh.manifest.json` 校验 SHA-256 后物化；验证：改 raw 字节 ⇒ 旧 Run stale；删 raw ⇒ Check 判 `fail`（不是 skip/pass）；新建 worktree 在主 checkout 无该网格时仍能取得相同输入。

## 4. pytest Check wrapper 与 Check 定义

- [ ] 4.1 `tools/devflow/` pytest wrapper：进入 `up3d`、`PYTHONNOUSERSITE=1`、线程上限、headless，读库确认 BLAS 线程，写原始 JUnit，校验预期 xfail/skip 名单，零测试/意外 skip/XPASS/名单不符 ⇒ 非零退出；验证：合成测试四种情形各判 `fail`，正常情形判 `pass`。
- [ ] 4.2 各 Check 定义以 `--deselect` / `-m` 排除图导出腿；验证：默认 verify 下每个非 xfail Check 的 JUnit `skipped="0"`。
- [ ] 4.3 含 strict xfail 的门文件（C02、C03；E01 进 milestone）各自以整个文件成不声明 `data` 的 Check；验证：每个文件 JUnit 中 `tests > skipped`，devflow 判 `pass`；故意让其中一条 xfail 变 XPASS ⇒ `fail`。
- [ ] 4.4 敏感证据 Check 的暖态自证（D5）：独立进程白跑 + 测量进程零 `cache_misses` 断言；验证前提 (a)：故意让一个 pyfp3d 内核在测量进程内编译 ⇒ 判 `fail`（G-TEETH）；验证前提 (b)：同一 codegen 下两次暖跑全精度逐位一致。**任一不成立即停止本 Change，回使用者重新裁决。**
- [ ] 4.5 “通过后清缓存，再执行 status/archive”：验证 status 仍为 verified（D5 语义），紧接着的冷态新鲜运行判 `fail`。

## 5. 协作规则与资源

- [ ] 5.1 改写 `AGENTS.md`（D9）并把 `CLAUDE.md` 改为指针；验证：`grep -n "在 \`main\` 上直接修改和提交" AGENTS.md` 无结果，`pytest tests/F/test_F07_repo_invariants.py` 通过。
- [ ] 5.2 `resources: {"solver-cpu": {"slots": 1}}`；验证：两个 Change 同时 verify 同一重型 Check 时第二个等待资源而不是并行运行。

## 6. 集成验证

- [ ] 6.1 plan-only verification 与一次真实 worktree smoke（`start` → 一个轻量 task → verify → complete）；验证：Run 记录含本 Change 的 commit、环境与有效 codegen，worktree 不读主 checkout 的偶然文件。
- [ ] 6.2 adoption 边界：`git diff --quiet main -- pyfp3d tests cases/reference_data cases/gates phases` 为真，且 `pytest tests/A/test_A01_freestream_preservation.py` 通过。

## Workflow follow-up

- 独立 reviewer 只读复核最终内容（`independentReview: true`）。
- 使用者分别批准 archive 与 push/PR/merge；任何一步都不授权另一步。
