import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent.parent
DEVFLOW_VERSION = "0.3.1"
DEVFLOW_COMMIT = "1c7ab6f223d4301aea8ea1deeb4dde2b246de37a"
AGENT_ANCHORS = (
    "开发权威只有本文件、`openspec/specs/`、活跃 Change artifacts 和 `devflow protocol execution`",
    "`openspec/plan/development-plan.json` 只是静态能力队列",
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
NEXT_CHANGE_FIELDS = {
    "id",
    "capability",
    "observable",
    "entry",
    "input",
    "output",
    "tests",
    "rollback",
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


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_json(path):
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path} must contain a JSON object")
    return value


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
    plan_path = ROOT / "openspec/plan/development-plan.json"
    plan = read_json(plan_path)
    require(set(plan) == {"formatVersion", "nextChange", "backlog"}, "Invalid plan top-level fields")
    require(
        isinstance(plan["formatVersion"], int)
        and not isinstance(plan["formatVersion"], bool)
        and plan["formatVersion"] == 1,
        "Plan formatVersion must be the integer 1",
    )

    next_change = plan["nextChange"]
    if next_change is not None:
        require(isinstance(next_change, dict), "nextChange must be null or one object")
        require(set(next_change) == NEXT_CHANGE_FIELDS, "Invalid nextChange fields")
        for field, value in next_change.items():
            require(isinstance(value, str) and value.strip(), f"nextChange.{field} must be non-empty text")
        for field in ("id", "capability"):
            require(KEBAB_CASE.fullmatch(next_change[field]), f"nextChange.{field} must be kebab-case")

    backlog = plan["backlog"]
    require(isinstance(backlog, list), "backlog must be an array")
    for index, item in enumerate(backlog):
        require(isinstance(item, dict), f"backlog[{index}] must be an object")
        require(set(item) == {"title"}, f"backlog[{index}] may only contain title")
        require(
            isinstance(item["title"], str) and item["title"].strip(),
            f"backlog[{index}].title must be non-empty text",
        )
    label = next_change["id"] if next_change is not None else "null"
    return f"nextChange={label}, backlog={len(backlog)}"


def check_project_and_checks():
    project = read_json(ROOT / "openspec/project.json")
    product_paths = project.get("productPaths")
    task_layer_paths = project.get("taskLayerPaths")
    require(isinstance(product_paths, list), "productPaths must be an array")
    require(isinstance(task_layer_paths, list), "taskLayerPaths must be an array")
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
        ("static development plan", check_plan),
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
