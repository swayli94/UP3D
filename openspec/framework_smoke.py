import importlib.metadata
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DEVFLOW_VERSION = "0.3.1"
DEVFLOW_TAG = "v0.3.1"
DEVFLOW_COMMIT = "1c7ab6f223d4301aea8ea1deeb4dde2b246de37a"
OPENSPEC_VERSION = "1.16.0"
CHANGE_ID = "framework-install-preflight"
MANIFEST_PATH = Path("openspec/devflow-install.json")


def run_result(command, cwd=ROOT):
    return subprocess.run(
        command,
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="backslashreplace",
    )


def run(command, cwd=ROOT):
    result = run_result(command, cwd)
    if result.returncode:
        output = (result.stdout + result.stderr).strip()
        raise RuntimeError(
            f"{' '.join(map(str, command))} exited {result.returncode}"
            + (f": {output}" if output else "")
        )
    return result.stdout.strip()


def git_value(*arguments):
    return run(["git", *arguments])


def main_checkout():
    git_dir = Path(git_value("rev-parse", "--path-format=absolute", "--git-dir"))
    common_dir = Path(git_value("rev-parse", "--path-format=absolute", "--git-common-dir"))
    if git_dir == common_dir:
        return Path(git_value("rev-parse", "--show-toplevel"))
    if common_dir.name == ".git":
        return common_dir.parent
    configured = git_value("config", "--get", "core.worktree")
    if configured:
        path = Path(configured)
        return path if path.is_absolute() else (common_dir / path).resolve()
    raise RuntimeError(f"Cannot locate the main checkout from {common_dir}")


def check_python():
    if sys.version_info < (3, 12):
        raise RuntimeError(f"Python >= 3.12 is required, got {sys.version.split()[0]}")
    return f"{sys.version.split()[0]} at {sys.executable}"


def check_devflow_installation():
    distribution = importlib.metadata.distribution("openspec-devflow")
    if distribution.version != DEVFLOW_VERSION:
        raise RuntimeError(f"Unexpected openspec-devflow version: {distribution.version!r}")
    direct_url = distribution.read_text("direct_url.json")
    if direct_url is None:
        raise RuntimeError("openspec-devflow has no direct_url.json")
    metadata = json.loads(direct_url)
    if metadata.get("dir_info") is not None:
        raise RuntimeError("Local or editable openspec-devflow installations are not allowed")
    vcs_info = metadata.get("vcs_info")
    if not isinstance(vcs_info, dict) or vcs_info.get("vcs") != "git":
        raise RuntimeError("openspec-devflow must be a VCS installation")
    if vcs_info.get("requested_revision") != DEVFLOW_TAG:
        raise RuntimeError(f"Unexpected openspec-devflow revision: {vcs_info.get('requested_revision')!r}")
    commit = vcs_info.get("commit_id")
    if commit != DEVFLOW_COMMIT:
        raise RuntimeError(f"Unexpected openspec-devflow commit: {commit!r}")
    return f"{distribution.version} from {DEVFLOW_TAG} at {commit}"


def check_node():
    version = run(["node", "--version"])
    match = re.fullmatch(r"v(\d+)\.(\d+)\.(\d+)", version)
    if match is None:
        raise RuntimeError(f"Cannot parse node version: {version!r}")
    actual = tuple(map(int, match.groups()))
    if actual < (22, 12, 0):
        raise RuntimeError(f"Node >= 22.12.0 is required, got {version}")
    return version


def check_openspec():
    executable = main_checkout() / "node_modules" / ".bin" / "openspec"
    version = run([str(executable), "--version"])
    if version != OPENSPEC_VERSION:
        raise RuntimeError(f"OpenSpec {OPENSPEC_VERSION} is required, got {version!r}")
    return f"{version} at {executable}"


def check_install():
    manifest = json.loads((ROOT / MANIFEST_PATH).read_text(encoding="utf-8"))
    if manifest.get("engineVersion") != DEVFLOW_VERSION:
        raise RuntimeError(f"Unexpected manifest engineVersion: {manifest.get('engineVersion')!r}")
    if manifest.get("engineCommit") != DEVFLOW_COMMIT:
        raise RuntimeError(f"Unexpected manifest engineCommit: {manifest.get('engineCommit')!r}")
    report = json.loads(run([sys.executable, "-m", "devflow", "install", "--check"]))
    if report.get("status") != "current":
        raise RuntimeError(f"devflow install is not current: {report}")
    return json.dumps(report, ensure_ascii=False, sort_keys=True)


def check_list_checks():
    report = json.loads(run([sys.executable, "-m", "devflow", "list-checks"]))
    expected = ["framework-smoke", "roadmap-contract", "roadmap-view", "solver-core-api", "solver-core-e2e", "workflow-contract"]
    if report != expected:
        raise RuntimeError(f"Unexpected list-checks output: {report}")
    return json.dumps(report, ensure_ascii=False)


def check_list_cases():
    report = json.loads(run([sys.executable, "-m", "devflow", "list-cases"]))
    if not isinstance(report, list):
        raise RuntimeError(f"Unexpected list-cases output: {report}")
    return json.dumps(report, ensure_ascii=False)


def check_execution_protocol():
    text = run([sys.executable, "-m", "devflow", "protocol", "execution"])
    for anchor in ("engineCommit", "executionSemantics", "devflow hold"):
        if anchor not in text:
            raise RuntimeError(f"Execution protocol is missing: {anchor}")
    return "engine identity, execution semantics and hold protocol present"


def check_configuration_protocol():
    text = run([sys.executable, "-m", "devflow", "protocol", "configuration"])
    for anchor in ("engineCommit", "outputEncoding", "strict"):
        if anchor not in text:
            raise RuntimeError(f"Configuration protocol is missing: {anchor}")
    return "strict project, Check, Case and Change configuration reference present"


def temporary_directory():
    runs = main_checkout() / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    return tempfile.TemporaryDirectory(prefix="framework-smoke-", dir=runs)


def write_temporary_manifest(root, **updates):
    manifest = json.loads((ROOT / MANIFEST_PATH).read_text(encoding="utf-8"))
    manifest.update(updates)
    path = root / MANIFEST_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def project_snapshot(root):
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file() and ".git" not in path.relative_to(root).parts
    }


def rejection_output(result, marker):
    output = (result.stdout + result.stderr).strip()
    if result.returncode != 1 or marker not in output:
        raise RuntimeError(f"Expected rejection containing {marker!r}, got {result.returncode}: {output}")
    return output


def check_version_guard():
    with temporary_directory() as directory:
        root = Path(directory)
        write_temporary_manifest(root, engineVersion="0.0.0")
        before = project_snapshot(root)
        result = run_result([
            sys.executable, "-m", "devflow", "--root", str(root), "list-checks",
        ])
        rejection_output(result, "Engine version mismatch")
        if project_snapshot(root) != before:
            raise RuntimeError("Version guard rejection changed the temporary project")
    return "version guard rejected before Registry loading without project writes"


def check_commit_preflight():
    with temporary_directory() as directory:
        root = Path(directory)
        shutil.copytree(
            ROOT / "openspec", root / "openspec",
            ignore=shutil.ignore_patterns("changes", "__pycache__"),
        )
        write_temporary_manifest(root, engineCommit="0" * 40)
        change = root / "openspec" / "changes" / CHANGE_ID
        change.mkdir(parents=True)
        metadata = {
            "formatVersion": 1,
            "id": CHANGE_ID,
            "impact": {
                "capabilities": [], "cases": [], "checks": ["framework-smoke"],
                "risks": [], "resources": [],
            },
            "requirements": {
                "platforms": ["Linux"], "manualChecks": [],
                "independentReview": False, "milestone": False,
            },
        }
        (change / "change.json").write_text(json.dumps(metadata), encoding="utf-8")
        (change / "tasks.md").write_text(
            "# Tasks\n\n## 1. Smoke\n\n- [ ] 1.1 Verify installation preflight.\n",
            encoding="utf-8",
        )
        execution = {"formatVersion": 1, "tasks": {"1.1": {"checks": ["framework-smoke"]}}}
        (change / "execution.json").write_text(json.dumps(execution), encoding="utf-8")
        run(["git", "init", "--quiet"], cwd=root)
        before = project_snapshot(root)
        result = run_result([
            sys.executable, "-m", "devflow", "--root", str(root), "verify",
            "--change", CHANGE_ID, "--task", "1.1",
        ])
        rejection_output(result, "Installation preflight failed")
        if project_snapshot(root) != before:
            raise RuntimeError("Commit preflight rejection changed the temporary project")
    return "commit preflight rejected real verify without Run or execution writes"


def main():
    checks = [
        ("Python >= 3.12", check_python),
        ("openspec-devflow installation", check_devflow_installation),
        ("Node >= 22.12.0", check_node),
        (f"main checkout OpenSpec {OPENSPEC_VERSION}", check_openspec),
        ("devflow manifest and install current", check_install),
        ("devflow strict list-checks", check_list_checks),
        ("devflow strict list-cases", check_list_cases),
        ("devflow execution protocol", check_execution_protocol),
        ("devflow configuration protocol", check_configuration_protocol),
        ("devflow version guard rejection", check_version_guard),
        ("devflow commit preflight rejection", check_commit_preflight),
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
        print(f"Framework smoke failed: {', '.join(failures)}")
        return 1
    print("Framework smoke passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
