# Design

## Context

见 `proposal.md` 的 Why。当前 `development-plan.json` 只含一个 P1 nextChange 和一个未展开 backlog；四个历史 Change 已归档，当前唯一活跃 Change 是 `adopt-roadmap`。固定 devflow v0.3.1 只读取 canonical `roadmap.json`，其启动依赖门读取关联 item 的 `deps`，普通视图与 live 视图分别写入 `internal/roadmap.html` 和 `runs/views/roadmap.html`。

主 checkout 在 `start --create` 前只能有本 Change 规划目录，因此实际 graph、plan brief 和视图不得在 start 前写入 `openspec/plan/` 或 `internal/`。规划 artifacts 评审完成后，先在 main 创建声明的 worktree，再在 worktree 中迁移。

`openspec/plan` 已在产品指纹和 workflow Check 输入中；`internal/roadmap.html` 将与生成它的任务同提交加入产品指纹。活跃 `openspec/changes` 不作为永久 Check 的 `sourcePaths`，因此最终 live 任务状态需要一次 milestone 新鲜验证，不能只依赖可复用的 task/affected Run。

## Goals / Non-Goals

**Goals:**

- 用单一 canonical graph 表达最小路线、真实依赖和 Change 关联，并保留现有 P1 详细计划。
- 普通视图可提交、确定且能检测手工漂移；live 视图能反映本机任务、worktree、锁和 merged 状态。
- 在 graph 自关联本 Change 的情况下，仍按单 Change、任务 Checks、独立复核和 archive 门完成交付。
- 新 Checks 可在没有备份和旧 Run 的干净 clone 中运行。

**Non-Goals:**

- 不维护 `globalMap`、`parallelPlan`、自动 SVG DAG 布局、Git hook 或 CI。
- 不修改 devflow 引擎或 `devflow/templates/roadmap.html`，不另写项目自定义路线图模板。
- 不添加 P1 之后的详细 Change artifacts，不实现 P1 梯度或其他数值能力。
- 不改变现有产品 API/CLI、参考数据、历史 acceptance、备份或科学证据。
- 不把 `change.json.planning.dependencies` 或 task `dependsOn` 复制成第二份 roadmap graph。

## Decisions

### 1. 原子迁移到唯一计划源

在 Change worktree 中先建立 `roadmap.json` 和 P1 brief，更新 workflow Check 后移除 `development-plan.json`。graph、旧文件删除、AGENTS 和检查逻辑在同一任务提交中落地，避免中间状态同时存在两份计划或旧 Check 误读新 graph。

主 specs 只在 archive 时由 delta 同步；实施期间不手工改写主 spec 来伪造已归档状态。新 workflow Check 在实施期同时理解当前主 spec 和本 Change delta，archive 后只要求同步后的主契约。

备选方案是保留 `development-plan.json` 作为兼容投影。它仍需要定义生成方向和防漂移规则，首版没有兼容消费者，因此不采用。

### 2. 最小 graph 的精确形状

`roadmap.json` 使用 `formatVersion: 2`、`version: 1.0.0` 和实际迁移日期填写的 `updated`。首版只定义两个 phase，不包含 `globalMap` 或 `parallelPlan`。

| Phase / lane | Item id | code | kind | stored status | deps | changes |
| --- | --- | --- | --- | --- | --- | --- |
| `foundation` / 已归档 Changes | `bootstrap-devflow` | F0 | proposal | done | — | bootstrap-devflow |
| 同上 | `restart-solver-core` | F1 | proposal | done | bootstrap-devflow | restart-solver-core |
| 同上 | `redesign-development-workflow` | F2 | proposal | done | restart-solver-core | redesign-development-workflow |
| 同上 | `upgrade-devflow-0-3-1` | F3 | proposal | done | redesign-development-workflow | upgrade-devflow-0-3-1 |
| `numerical-core` / 最小路线 | `roadmap-adoption` | R0 | proposal | planned → done | upgrade-devflow-0-3-1 | adopt-roadmap |
| 同上 | `tetrahedron-p1-gradient` | P1 | proposal | planned | roadmap-adoption | 暂无 |
| 同上 | `cell-field-reconstruction` | C1 | discussion | deferred | tetrahedron-p1-gradient | 暂无 |

`cell-field-reconstruction` 另设 `backlog: unscheduled`。其他 deferred 字段不得出现在非 deferred 项上。历史节点的 `end` 首版省略，避免把提交、归档和 merge 日期混为一谈。

`roadmap-adoption` 初始存储为 `planned`，关联活跃 `adopt-roadmap` 后自动显示 active。只有在 graph、视图和检查全部完成后才将其存储状态改为 `done`；archive 前它仍显示 active，archive 后显示 archived，归档记录进入 main 后 live 才显示 merged。P1 在 adoption archived 后自动成为 ready。

### 3. P1 详细计划迁入一个 plan brief

新增 `openspec/plan/briefs/tetrahedron-p1-gradient.md`，由 P1 item 以 `brief: briefs/tetrahedron-p1-gradient.md` 引用。内容使用 devflow brief 支持的 Markdown 子集，包含：

- Scope：从节点势函数计算每个四面体单元的常梯度。
- Contract：公共 `run_p1_gradient(path)` 与 `python -m pyfp3d p1-gradient <input.json>`，以及现有输入、N×3 gradients 输出和单元顺序。
- Validation：仿射场、手工非仿射节点值、非法输入、真实 API/E2E 和原自由流回归。
- Rollback：只移除新势函数输入、梯度入口、对应 spec、测试和 Checks。
- 待决事项：数值输出语义、共享网格验证边界及旧裸文件 CLI 的兼容方式，留待 P1 自己的 proposal/design 决定。

`cell-field-reconstruction` 只保留一行 evaluation summary，不写 brief、proposal、design 或 Change id。

### 4. 视图配置与跟踪策略

`openspec/project.json` 增加：

```json
"roadmapView": {
  "title": "UP3D 开发路线图",
  "locale": "zh-CN"
}
```

`internal/roadmap.html` 由 `devflow roadmap-view` 生成并提交，同时加入 `productPaths` 和 `roadmap-view` Check 的输入指纹。生成它的任务必须同提交包含 graph、project 配置和该文件，避免声明不存在的输入。

`runs/views/roadmap.html` 由 `devflow roadmap-view --live` 生成，继续使用现有 `/runs/` 忽略规则，不加入产品指纹或 Git。普通视图不含任务计数；live 视图不要求确定性，但要求与生成时的 Change/worktree/tasks/locks/Git 状态一致。

视图刷新边界为：graph 或 brief 改动、任务提交、PR、merge、archive，以及需要观察本机锁或 worktree 时。本 Change 不安装自动 hook；调用者或项目适配器显式运行生成器。

### 5. 三类 Check 分工

保留并更新 `workflow-contract`：验证 AGENTS 锚点、唯一计划源、最小 graph 规则、brief 必需内容、无旧 `development-plan.json`、单活跃 Change 和永久 Check 可移植性。

新增 `roadmap-contract`，命令为 `python openspec/roadmap_contract.py`：

- 通过 devflow plan loader 读取当前 graph，校验 item id、deps、changes、stored status、deferred/backlog 和 brief 引用。
- 在 `runs/` 下的真实临时 Git 项目中构造最小注册表、plan 和 Change，使用真实 devflow CLI 验证未完成依赖拒绝 start、完成后允许 start/next；临时项目使用 local-only 配置，不读取 UP3D origin、备份或主 worktree。
- 验证 planned/ready、active、done/archived、reopened dependency 等状态推导，以及 deferred 不会自动 ready。

新增 `roadmap-view`，命令为 `python openspec/roadmap_view_check.py`：

- 连续两次生成普通视图并要求字节一致，且与已提交 `internal/roadmap.html` 相同；检查模板替换完成、无 live overlay。
- 生成 live 视图并解析嵌入 JSON，检查当前 worktree、Change 任务计数、历史归档 merged 和锁结构；不断言本 Change 在未 merge 时已 merged。
- 确认 live 输出只位于 `runs/views/`，普通视图不含机器路径或任务计数。

两个新 Check 均为 Linux、required、fast、console evidence、`profiles: ["task", "affected"]`，`data/resources` 为空。`sourcePaths` 只声明当前仓库可提交输入，不包含 `openspec/changes`、runs、备份或历史 Run。`framework_smoke.py` 的 Check 清单同步扩展到六个 id。

### 6. 两阶段执行元数据与最终新鲜验证

初始 `execution.json` 只引用 main 已注册的 `framework-smoke`、`solver-core-api`、`solver-core-e2e`、`workflow-contract`，保证 start 和首次 next 可解析。`change.json.impact.checks` 可以提前声明两个未来 Check。

任务 1.2 同提交创建两个 Check 定义、脚本、普通视图、产品路径和最终 execution 映射，随后 CLI 重新加载注册表，正式 verify 才能使用新 id。不得在定义文件和声明输入尚未提交时提前让 tasks 引用它们。

本 Change 是 roadmap 集成门，设置 `requirements.milestone: true`。全部任务后除 affected verification 外，还运行 milestone，使所有非 manual Checks 新鲜执行；这样最终 live 任务计数、graph 的 done 状态和视图不是旧 Run 复用产生的假证据。现有产品测试很快，无需新增资源锁。

### 7. 任务与状态提交顺序

1. 任务 1.1 原子建立 graph、brief、唯一计划源、AGENTS 和更新后的 workflow 契约；使用现有 Checks 验证。
2. 任务 1.2 注册 roadmap Checks、生成两类视图并提交普通视图，扩展 execution 映射后验证。
3. 任务 1.3 在实现证据齐全后将 adoption item 标为 done、刷新普通和 live 视图，并进行最终任务验证。
4. 全部任务后运行 affected 与 milestone，再由独立上下文复核当前产品、graph、brief、规划 artifacts 和视图。
5. 用户批准 archive 后，由 archive 同步主 spec 并刷新普通视图；提交归档结果。push、PR、merge 分别授权，merge 后重新生成普通和 live 视图。

## Risks / Trade-offs

- [graph 与旧计划双写漂移] → 在任务 1.1 同提交退役旧文件，并由 workflow-contract 禁止其重新出现。
- [本 Change 自关联造成循环或启动死锁] → adoption 只依赖已归档升级项；初始 execution 只用已注册 Check，启动后才迁移 graph。
- [新 Check id 使 next 无法解析] → 任务 1.2 同提交注册定义、脚本、声明输入和最终 execution 映射。
- [普通视图被手工编辑或过期] → 提交生成文件并用双生成字节比较检测漂移，禁止以手改 HTML 完成任务。
- [live 进度被旧 Run 复用遮蔽] → milestone 在全部任务后新鲜执行；live 不作为可提交验收快照。
- [提前展开 backlog 形成伪授权] → graph Check 只允许一个 deferred evaluation 项，新增节点必须另经规划讨论。
- [done 状态提前或复核后失效] → 仅在任务 1.3 标 done；之后任何 graph、proposal、design、spec delta 或任务文字变更都重新验证并重新复核。
- [首版没有全局 SVG 图] → 接受这一展示限制，避免重复维护坐标；graph 稳定后再另开规划讨论。

## Migration Plan

1. 用户审阅全部规划 artifacts；确认后从同步干净 main 执行 `devflow start --change adopt-roadmap --create`，进入返回的 worktree 并运行首次 next。
2. 按任务 1.1 → 1.2 → 1.3 实施，每项先提交声明输入，再 verify --task、complete，最后提交任务状态和 Run 元数据。
3. 运行 affected 和 milestone verification；独立复核通过后请求用户 archive 批准，不手工勾选任务或跳过证据门。
4. archive 后提交主 spec、归档记录和刷新视图；后续 push、PR、merge 分别取得授权，merge 后刷新普通与 live 视图。
5. 合入前放弃时，移除 graph、brief、内部视图、新 Checks 与治理改动并恢复 `development-plan.json`；产品代码、参考数据、历史 Change 和备份不动。合入后回退另开带新 spec 的 Change。
