import importlib.metadata
import json
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DEVFLOW_COMMIT = "cc0a9f12704f"
OPENSPEC_VERSION = "1.16.0"


def run(command):
    result = subprocess.run(
        command,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="backslashreplace",
    )
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


def check_devflow_commit():
    distribution = importlib.metadata.distribution("openspec-devflow")
    direct_url = distribution.read_text("direct_url.json")
    if direct_url is None:
        raise RuntimeError("openspec-devflow has no direct_url.json")
    metadata = json.loads(direct_url)
    if metadata.get("dir_info", {}).get("editable") is True:
        raise RuntimeError("Editable openspec-devflow installations are not allowed")
    commit = metadata.get("vcs_info", {}).get("commit_id", "")
    if not commit.startswith(DEVFLOW_COMMIT):
        raise RuntimeError(f"Unexpected openspec-devflow commit: {commit!r}")
    return commit


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
    report = json.loads(run([sys.executable, "-m", "devflow", "install", "--check"]))
    if report.get("status") != "current":
        raise RuntimeError(f"devflow install is not current: {report}")
    return json.dumps(report, ensure_ascii=False, sort_keys=True)


def check_list_checks():
    report = json.loads(run([sys.executable, "-m", "devflow", "list-checks"]))
    if not isinstance(report, list):
        raise RuntimeError(f"Unexpected list-checks output: {report}")
    return json.dumps(report, ensure_ascii=False)


def check_list_cases():
    report = json.loads(run([sys.executable, "-m", "devflow", "list-cases"]))
    if not isinstance(report, list):
        raise RuntimeError(f"Unexpected list-cases output: {report}")
    return json.dumps(report, ensure_ascii=False)


def main():
    checks = [
        ("Python >= 3.12", check_python),
        ("openspec-devflow commit", check_devflow_commit),
        ("Node >= 22.12.0", check_node),
        (f"main checkout OpenSpec {OPENSPEC_VERSION}", check_openspec),
        ("devflow install current", check_install),
        ("devflow list-checks", check_list_checks),
        ("devflow list-cases", check_list_cases),
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
