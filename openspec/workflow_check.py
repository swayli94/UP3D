import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent.parent
DEVFLOW_VERSION = "0.3.1"
DEVFLOW_COMMIT = "1c7ab6f223d4301aea8ea1deeb4dde2b246de37a"
AGENT_ANCHORS = (
    "开发权威只有本文件、`openspec/specs/`、活跃 Change artifacts 和 `devflow protocol execution`",
    "`openspec/plan/roadmap.json` 是唯一协调 graph",
    "active/ready/archived/merged",
    "`devflow roadmap-view`",
    "`runs/views/roadmap.html`",
    "`pyfp3d-core`",
    "v0.3.1",
    DEVFLOW_COMMIT,
    "engineCommit",
    "一次最多规划或执行一个活跃 Change",
    "同步干净 `main` → propose → `devflow start --create` → `devflow next`/`devflow run`",
    "真实入口测试",
    "独立 worktree",
    "可观察能力",
    "回滚边界",
    "任务级 Checks",
    "独立复核",
    "用户明确批准",
    "push、PR、merge",
    "分别取得使用者授权",
)
AGENT_FORBIDDEN = (
    "docs/migration-plan.md",
    "restart-solver-core",
    "docs/dev_phase_seven",
    "phases/",
)
EXPECTED_ITEM_IDS = (
    "bootstrap-devflow",
    "restart-solver-core",
    "redesign-development-workflow",
    "upgrade-devflow-0-3-1",
    "roadmap-adoption",
    "tetrahedron-p1-gradient",
    "cell-field-reconstruction",
)
EXPECTED_DEPS = {
    "bootstrap-devflow": [],
    "restart-solver-core": ["bootstrap-devflow"],
    "redesign-development-workflow": ["restart-solver-core"],
    "upgrade-devflow-0-3-1": ["redesign-development-workflow"],
    "roadmap-adoption": ["upgrade-devflow-0-3-1"],
    "tetrahedron-p1-gradient": ["roadmap-adoption"],
    "cell-field-reconstruction": ["tetrahedron-p1-gradient"],
}
EXPECTED_CHANGES = {
    "bootstrap-devflow": ["bootstrap-devflow"],
    "restart-solver-core": ["restart-solver-core"],
    "redesign-development-workflow": ["redesign-development-workflow"],
    "upgrade-devflow-0-3-1": ["upgrade-devflow-0-3-1"],
    "roadmap-adoption": ["adopt-roadmap"],
    "tetrahedron-p1-gradient": [],
    "cell-field-reconstruction": [],
}
REQUIRED_PRODUCT_PATHS = {
    "openspec/specs",
    "openspec/plan",
    "openspec/workflow_check.py",
    "openspec/devflow-install.json",
}
REQUIRED_WORKFLOW_SOURCES = {
    "AGENTS.md",
    "openspec/specs",
    "openspec/plan",
    "openspec/workflow_check.py",
    "openspec/project.json",
    "openspec/devflow-install.json",
    "openspec/checks",
}
FORBIDDEN_CHECK_DEPENDENCIES = (
    "/home/lrz/code/up3d_backup",
    "up3d-backup",
    "pre-devflow-restart",
    "openspec/changes",
    "runs/evidence",
    "execution.jsonl",
)
KEBAB_CASE = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*")
PLAN_STATUSES = {"planned", "draft", "deferred", "done"}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_json(path):
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path} must contain a JSON object")
    return value


def change_inventory():
    base = ROOT / "openspec/changes"
    paths = sorted(base.glob("*/change.json")) + sorted(base.glob("archive/*/change.json"))
    inventory = {}
    active = []
    for path in paths:
        metadata = read_json(path)
        change_id = metadata.get("id")
        require(isinstance(change_id, str) and change_id, f"Invalid Change id in {path}")
        require(change_id not in inventory, f"Duplicate Change id: {change_id}")
        inventory[change_id] = path
        if path.parent.parent == base:
            active.append(change_id)
    return inventory, sorted(active)


def check_agents():
    text = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    for anchor in AGENT_ANCHORS:
        require(anchor in text, f"AGENTS.md is missing required anchor: {anchor}")
    for forbidden in AGENT_FORBIDDEN:
        require(forbidden not in text, f"AGENTS.md has obsolete authority reference: {forbidden}")
    return f"{len(AGENT_ANCHORS)} anchors, no obsolete authority reference"


def check_framework_manifest():
    manifest = read_json(ROOT / "openspec/devflow-install.json")
    require(manifest.get("formatVersion") == 1, "Install manifest formatVersion must be 1")
    require(
        manifest.get("engineVersion") == DEVFLOW_VERSION,
        f"Install manifest engineVersion must be {DEVFLOW_VERSION}",
    )
    require(
        manifest.get("engineCommit") == DEVFLOW_COMMIT,
        f"Install manifest engineCommit must be {DEVFLOW_COMMIT}",
    )
    return f"engineVersion={DEVFLOW_VERSION}, engineCommit={DEVFLOW_COMMIT}"


def check_plan():
    old_plan = ROOT / "openspec/plan/development-plan.json"
    require(not old_plan.exists(), "Legacy development-plan.json must not coexist with roadmap.json")
    plan = read_json(ROOT / "openspec/plan/roadmap.json")
    require(
        isinstance(plan.get("formatVersion"), int)
        and not isinstance(plan["formatVersion"], bool)
        and plan["formatVersion"] == 2,
        "Roadmap formatVersion must be the integer 2",
    )
    for field in ("version", "updated"):
        require(isinstance(plan.get(field), str) and plan[field], f"Roadmap {field} must be non-empty text")
    require("globalMap" not in plan and "parallelPlan" not in plan, "Initial roadmap must omit optional display maps")

    items = [
        item
        for phase in plan.get("phases", [])
        for lane in phase.get("lanes", [])
        for item in lane.get("items", [])
    ]
    require(tuple(item.get("id") for item in items) == EXPECTED_ITEM_IDS, "Unexpected minimal roadmap items or order")
    by_id = {item["id"]: item for item in items}
    inventory, active = change_inventory()
    require(len(active) <= 1, f"At most one active Change is allowed: {active}")

    linked = {}
    for item in items:
        item_id = item["id"]
        require(KEBAB_CASE.fullmatch(item_id), f"Roadmap item id must be kebab-case: {item_id}")
        for field in ("code", "title", "kind", "status", "summary"):
            require(isinstance(item.get(field), str) and item[field], f"{item_id}.{field} must be non-empty text")
        require(item["status"] in PLAN_STATUSES, f"{item_id} has derived or unsupported stored status")
        require(item.get("deps", []) == EXPECTED_DEPS[item_id], f"Unexpected dependencies for {item_id}")
        require(item.get("changes", []) == EXPECTED_CHANGES[item_id], f"Unexpected Change links for {item_id}")
        for change_id in item.get("changes", []):
            require(change_id in inventory, f"{item_id} links unknown Change: {change_id}")
            require(change_id not in linked, f"Change {change_id} is linked more than once")
            linked[change_id] = item_id
        if item["status"] == "deferred":
            require(item.get("backlog") == "unscheduled", f"Deferred item {item_id} must be unscheduled")
        else:
            require("backlog" not in item, f"Non-deferred item {item_id} must not have a backlog group")
    require(set(linked) == set(inventory), "Every active and archived Change must be linked exactly once")

    for item_id in EXPECTED_ITEM_IDS[:4]:
        require(by_id[item_id]["status"] == "done", f"Historical item {item_id} must store done")
    require(by_id["roadmap-adoption"]["status"] in {"planned", "done"}, "Adoption must be planned or done")
    require(by_id["tetrahedron-p1-gradient"]["status"] == "planned", "P1 must remain planned")
    require(by_id["cell-field-reconstruction"]["status"] == "deferred", "Reconstruction must remain deferred")
    require("brief" not in by_id["cell-field-reconstruction"], "Deferred reconstruction must not have a brief")

    brief_name = by_id["tetrahedron-p1-gradient"].get("brief")
    require(brief_name == "briefs/tetrahedron-p1-gradient.md", "Unexpected P1 brief reference")
    brief = (ROOT / "openspec/plan" / brief_name).read_text(encoding="utf-8")
    for anchor in (
        "run_p1_gradient(path)",
        "python -m pyfp3d p1-gradient <input.json>",
        "vertices",
        "tetrahedra",
        "nodal_potential",
        "N×3 gradients",
        "仿射",
        "非仿射",
        "API/E2E",
        "Rollback",
        "待决事项",
    ):
        require(anchor in brief, f"P1 brief is missing required content: {anchor}")
    return f"items={len(items)}, linkedChanges={len(linked)}, activeChanges={len(active)}, deferred=1"


def check_project_and_checks():
    project = read_json(ROOT / "openspec/project.json")
    product_paths = project.get("productPaths")
    task_layer_paths = project.get("taskLayerPaths")
    require(isinstance(product_paths, list), "project.productPaths must be an array")
    require(isinstance(task_layer_paths, list), "project.taskLayerPaths must be an array")
    require(REQUIRED_PRODUCT_PATHS <= set(product_paths), "productPaths is missing workflow inputs")
    for field, paths in (("productPaths", product_paths), ("taskLayerPaths", task_layer_paths)):
        for path in paths:
            require(isinstance(path, str) and path, f"{field} contains an invalid path")
            require((ROOT / path).exists(), f"{field} path does not exist: {path}")
    require(list((ROOT / "openspec/specs").glob("**/spec.md")), "No main spec found")

    check_files = sorted((ROOT / "openspec/checks").glob("**/*.json"))
    require(check_files, "No permanent Checks found")
    for path in check_files:
        check = read_json(path)
        for field in ("command", "data", "sourcePaths"):
            serialized = json.dumps(check.get(field, []), ensure_ascii=False)
            for forbidden in FORBIDDEN_CHECK_DEPENDENCIES:
                require(
                    forbidden not in serialized,
                    f"{check.get('id', path.name)} has non-portable {field}: {forbidden}",
                )
        source_paths = check.get("sourcePaths", [])
        require(isinstance(source_paths, list), f"{check.get('id', path.name)} sourcePaths must be an array")
        for source in source_paths:
            require(isinstance(source, str) and source, f"{path} contains an invalid sourcePath")
            require((ROOT / source).exists(), f"{check.get('id', path.name)} sourcePath does not exist: {source}")

    workflow = read_json(ROOT / "openspec/checks/workflow-contract.json")
    require(
        workflow.get("command") == ["python", "openspec/workflow_check.py"],
        "workflow-contract must execute openspec/workflow_check.py",
    )
    require(
        REQUIRED_WORKFLOW_SOURCES <= set(workflow.get("sourcePaths", [])),
        "workflow-contract sourcePaths is incomplete",
    )
    return f"{len(check_files)} portable Checks, all declared paths exist"


def main():
    checks = [
        ("AGENTS workflow anchors", check_agents),
        ("devflow install manifest", check_framework_manifest),
        ("canonical roadmap graph", check_plan),
        ("project and permanent Check portability", check_project_and_checks),
    ]
    failures = []
    for name, check in checks:
        try:
            detail = check()
        except Exception as error:
            failures.append(name)
            print(f"FAIL {name}: {error}")
        else:
            print(f"PASS {name}: {detail}")
    if failures:
        print(f"Workflow contract failed: {', '.join(failures)}")
        return 1
    print("Workflow contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
