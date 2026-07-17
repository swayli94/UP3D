# AGENTS.md — 跨 agent 协作规则（Kimi / Claude 共用）

本文件是仓库级 agent 规则。项目工作规则仍以 [CLAUDE.md](CLAUDE.md)、
[docs/agent-rules.md](docs/agent-rules.md)、[docs/roadmap.md](docs/roadmap.md)
为准；本文件只补充**多 agent 并行开发的 git 协作模式**。

## Git 协作模式（最高优先级）

本仓库由两个 AI agent 并行开发：

- **Claude**：在工作副本 `UP3D/` 中直接于 `main` 上修改和提交。
- **Kimi（Kimi Code CLI）**：在工作副本 `UP3D_kimi/` 中工作，**永远不在 `main` 上提交**。

Kimi 必须遵守：

1. 开工前先 `git fetch origin`，从最新的 `origin/main` 新建分支，
   分支命名 `kimi/<topic>`（如 `kimi/agent-rules-and-inspection`）。
2. 所有修改只提交到自己的 `kimi/*` 分支；绝不 `git commit` 到 `main`，
   绝不 force-push，绝不动别人的分支。
3. 完成后通过 **Pull Request 合入 `main`**（PR 目标分支 = `main`）；
   **push 和开 PR 之前必须先获得用户确认**。
4. 分支保持小步、单一主题；`main` 会被 Claude 持续推进，
   因此工作期间要频繁 `git fetch origin` 并同步（rebase/merge 最新 `origin/main`），
   避免基于过期 main 做大规模改动。
5. 任何 git 破坏性操作（reset --hard、rebase 已推送分支、删除分支等）
   必须先经用户确认。

## 项目关键纪律速记（详细版见 docs/agent-rules.md）

- 线程上限 16，含 BLAS/OMP：`NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16`。
- kernel/assembly 改动后先跑 `pytest tests/test_v0_freestream.py`，再跑更大范围。
- 证据纪律：结论必须有已提交的 artifact（PNG/CSV）；只写在 .md 里的数字不算证据。
- 不重算贵重已提交 artifact（P4 heavy demo ~40 min、P5 medium 45–75 min 等），
  以 committed CSV/PNG 为准。
- `cases/reference_data/` 是 ground truth，永不修改。
- 重型 gate 只在 `PYFP3D_TRANSONIC_GATES=1` 下运行；M6 `.msh` 被 gitignore，
  需要先跑 `cases/meshes/onera_m6/generate_onera_m6.py` 生成。
