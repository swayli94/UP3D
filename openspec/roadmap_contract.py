import json
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

from devflow.registry import change_directories, load_plan, plan_items
from devflow.roadmap_view import derive_statuses


ROOT = Path(__file__).resolve().parent.parent
HISTORICAL = (
    "bootstrap-devflow",
    "restart-solver-core",
    "redesign-development-workflow",
    "upgrade-devflow-0-3-1",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run_result(command, cwd):
    env = dict(os.environ)
    env["PYTHONNOUSERSITE"] = "1"
    return subprocess.run(
        command,
        cwd=cwd,
        env=env,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="backslashreplace",
    )


def run(command, cwd):
    result = run_result(command, cwd)
    if result.returncode:
        output = (result.stdout + result.stderr).strip()
        raise RuntimeError(f"{' '.join(map(str, command))} exited {result.returncode}: {output}")
    return result.stdout.strip()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def check_current_plan():
    plan = load_plan(ROOT)
    items = plan_items(plan)
    by_id = {item["id"]: item for item in items}
    directories = change_directories(ROOT)
    active = {change for change, path in directories.items() if not path.startswith("archive/")}
    archived = set(directories) - active
    statuses = derive_statuses(items, active, archived)
    for item_id in HISTORICAL:
        require(statuses[item_id] == "archived", f"{item_id} must be archived")
    require("globalMap" not in plan and "parallelPlan" not in plan, "Optional display maps are not part of the minimal graph")
    require(statuses["cell-field-reconstruction"] == "deferred", "Reconstruction must remain deferred")
    if "adopt-roadmap" in active:
        require(statuses["roadmap-adoption"] == "active", "Adoption must be active during its Change")
        require(statuses["tetrahedron-p1-gradient"] == "planned", "P1 must wait for adoption")
    else:
        require(statuses["roadmap-adoption"] == "archived", "Adoption must be archived after its Change")
        require(statuses["tetrahedron-p1-gradient"] == "ready", "P1 must become ready after adoption")
    return f"current statuses: adoption={statuses['roadmap-adoption']}, p1={statuses['tetrahedron-p1-gradient']}"


def synthetic_item(item_id, status="planned", **fields):
    return {
        "id": item_id,
        "code": item_id.upper(),
        "title": item_id,
        "kind": "proposal",
        "status": status,
        "summary": f"Synthetic {item_id}.",
        **fields,
    }


def check_status_derivation():
    base = synthetic_item("base")
    dependent = synthetic_item("dependent", deps=["base"])
    deferred = synthetic_item("later", status="deferred", deps=["base"], backlog="unscheduled")
    statuses = derive_statuses([base, dependent, deferred], set(), set())
    require(statuses == {"base": "ready", "dependent": "planned", "later": "deferred"}, "Initial ready/deferred derivation is wrong")

    base["status"] = "done"
    statuses = derive_statuses([base, dependent, deferred], set(), set())
    require(statuses["dependent"] == "ready" and statuses["later"] == "deferred", "Done dependency must unlock only planned items")

    base["changes"] = ["old"]
    statuses = derive_statuses([base, dependent, deferred], set(), {"old"})
    require(statuses["base"] == "archived" and statuses["dependent"] == "ready", "Archived dependency must unlock planned items")

    base["changes"] = ["redo"]
    statuses = derive_statuses([base, dependent, deferred], {"redo"}, set())
    require(statuses["base"] == "active" and statuses["dependent"] == "planned", "Reopened dependency must block dependents")
    return "ready/active/archived/reopened/deferred derivation passed"


def initialize_synthetic(root):
    (root / "src").mkdir(parents=True)
    (root / "src/input.txt").write_text("synthetic\n", encoding="utf-8")
    (root / ".gitignore").write_text("runs/\n__pycache__/\n", encoding="utf-8")
    write_json(root / "openspec/project.json", {
        "formatVersion": 1,
        "integrationBranch": "main",
        "integrationRemote": None,
        "evidenceRoot": "runs/evidence",
        "resourceStateRoot": "runs/devflow",
        "productPaths": ["src"],
        "taskLayerPaths": ["src"],
        "environmentPackages": [],
        "resources": {},
        "roadmapView": {"title": "Synthetic roadmap", "locale": "zh-CN"},
    })
    write_json(root / "openspec/checks/smoke.json", {
        "formatVersion": 1,
        "id": "smoke",
        "executionKey": "smoke",
        "title": "Synthetic smoke",
        "capabilities": [],
        "cases": [],
        "command": [sys.executable, "-c", "print('smoke')"],
        "platforms": [platform.system()],
        "data": [],
        "resources": [],
        "cost": "fast",
        "slow": False,
        "required": True,
        "profiles": ["task", "affected"],
        "evidence": ["console"],
        "manual": False,
        "sourcePaths": ["src"],
    })
    write_json(root / "openspec/changes/work/change.json", {
        "formatVersion": 1,
        "id": "work",
        "branch": "feat/work",
        "worktree": {"id": "work", "machine": platform.system(), "path": "runs/worktrees/work"},
        "planning": {"dependencies": [], "lane": "L1"},
        "impact": {"capabilities": [], "cases": [], "checks": ["smoke"], "risks": [], "resources": []},
        "requirements": {
            "platforms": [platform.system()],
            "manualChecks": [],
            "independentReview": False,
            "milestone": False,
        },
    })
    (root / "openspec/changes/work/tasks.md").write_text("- [ ] 1.1 Run the synthetic smoke Check.\n", encoding="utf-8")
    write_json(root / "openspec/changes/work/execution.json", {
        "formatVersion": 1,
        "tasks": {"1.1": {"checks": ["smoke"]}},
    })
    write_json(root / "openspec/plan/roadmap.json", {
        "formatVersion": 2,
        "phases": [{
            "id": "synthetic",
            "title": "Synthetic",
            "lanes": [{
                "title": "Dependency gate",
                "items": [
                    synthetic_item("base"),
                    synthetic_item("work", deps=["base"], changes=["work"]),
                ],
            }],
        }],
    })
    run(["git", "init", "--quiet", "-b", "main"], root)
    run(["git", "config", "user.name", "Roadmap Check"], root)
    run(["git", "config", "user.email", "roadmap@example.invalid"], root)
    run(["git", "add", "--", ".gitignore", "src", "openspec/project.json", "openspec/checks/smoke.json", "openspec/plan/roadmap.json"], root)
    run(["git", "commit", "--quiet", "-m", "Synthetic roadmap"], root)


def check_start_dependency_gate():
    runs = ROOT / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="roadmap-contract-", dir=runs) as directory:
        root = Path(directory)
        initialize_synthetic(root)
        command = [sys.executable, "-m", "devflow", "--root", str(root), "start", "--change", "work", "--create"]
        blocked = run_result(command, root)
        output = blocked.stdout + blocked.stderr
        require(blocked.returncode == 1 and "Plan dependencies are incomplete: ['base']" in output, f"Blocked start was not rejected correctly: {output}")
        require(not (root / "runs/worktrees/work").exists(), "Rejected start created a worktree")
        require(not run(["git", "branch", "--list", "feat/work"], root), "Rejected start created a branch")

        plan_path = root / "openspec/plan/roadmap.json"
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        plan["phases"][0]["lanes"][0]["items"][0]["status"] = "done"
        write_json(plan_path, plan)
        run(["git", "add", "--", "openspec/plan/roadmap.json"], root)
        run(["git", "commit", "--quiet", "-m", "Complete synthetic dependency"], root)
        allowed = run_result(command, root)
        output = allowed.stdout + allowed.stderr
        require(allowed.returncode == 0, f"Completed dependency did not allow start: {output}")
        report = json.loads(allowed.stdout)
        require(report["state"] == "ready" and report["tasks"][0]["id"] == "1.1", "Allowed start did not return its task")
        worktree = Path(report["worktree"])
        require(worktree.is_dir() and (worktree / "openspec/changes/work/tasks.md").is_file(), "Allowed start did not create the declared worktree")
    return "real CLI start rejected an incomplete dependency and allowed a completed one"


def main():
    checks = [
        ("current canonical graph", check_current_plan),
        ("synthetic status derivation", check_status_derivation),
        ("real start dependency gate", check_start_dependency_gate),
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
        print(f"Roadmap contract failed: {', '.join(failures)}")
        return 1
    print("Roadmap contract passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
