"""Does devflow (v0.2.0, cc0a9f1) see a Git-IGNORED file used through Case `extraInputs` / raw?

usage: conda run -n openspec-devflow-dev python -I -c "import sys; sys.path.insert(0, '<devflow checkout>'); exec(open('exp_extra_inputs.py').read())"
Mirrors the `_demo` fixture of devflow/tests/unit/test_devflow_resources.py.
arm1: extra input + in-checkout raw, both git-ignored -> edit each, read hash change and dirty()
arm2: control, nothing ignored -> the same edit must be reported dirty
arm3: ignored extra input ABSENT (as in a fresh Change worktree)
"""
import hashlib, json, subprocess, sys, tempfile
from pathlib import Path
from devflow.fingerprint import GitFingerprints, content_hash
from devflow.registry import Case, CaseStep, Check, Registry

def git(root, *a):
    subprocess.run(["git", "-c", "user.email=t@e", "-c", "user.name=T", *a], cwd=root, check=True, capture_output=True)

def build(root, ignore_lines, extra):
    raw = root / "data" / "demo" / "raw" / "input.txt"
    raw.parent.mkdir(parents=True)
    raw.write_text("first")
    (root / "data/demo/manifest.json").write_text(json.dumps({"files": [{"localName": "input.txt"}]}))
    (root / "cases/meshes").mkdir(parents=True)
    (root / extra).write_text("mesh-v1")
    (root / ".gitignore").write_text("\n".join(ignore_lines) + "\n")
    git(root, "init", "-q", "-b", "main"); git(root, "add", "-A"); git(root, "commit", "-qm", "init")
    case = Case(id="demo", title="Demo", capabilities=("x",), manifest="data/demo/manifest.json",
                role="test", input="test", entry="test", target="regress", non_goal="none", resources=(),
                steps=(CaseStep("s1", (), ("o",), "runs/demo/s1", ("python", "s1.py")),),
                extra_inputs=(extra,), source_paths=("data/demo/manifest.json",))
    check = Check.from_dict({"id": "c", "executionKey": "c", "title": "c", "capabilities": ["x"],
        "cases": ["demo"], "command": ["python"], "platforms": ["Linux"], "data": ["demo"],
        "resources": [], "cost": "test", "slow": False, "required": True, "profiles": ["affected"],
        "evidence": ["console"], "manual": False})
    reg = Registry(root=root, checks={"c": check}, cases={"demo": case}, project={"dataEnvironmentVariable": "XX_UNSET"})
    return GitFingerprints(root, reg, root), check, raw

extra = "cases/meshes/m6.msh"
# Arm 1: extra input git-ignored; raw (inside checkout at data/) also git-ignored
with tempfile.TemporaryDirectory() as d:
    root = Path(d)
    fp, ck, raw = build(root, ["cases/meshes/*.msh", "data/demo/raw/"], extra)
    h0 = content_hash(fp.content(ck)); print("arm1 dirty baseline:", fp.dirty(ck))
    (root / extra).write_text("mesh-v2")
    h1 = content_hash(fp.content(ck))
    print("arm1 extra-input hash changed:", h0 != h1, "| dirty after ignored edit:", fp.dirty(ck))
    raw.write_text("second")
    h2 = content_hash(fp.content(ck))
    print("arm1 ignored in-checkout raw: hash changed:", h1 != h2, "| dirty:", fp.dirty(ck))
# Arm 2 (control): same, nothing ignored -> extra tracked, edit must be dirty
with tempfile.TemporaryDirectory() as d:
    root = Path(d)
    fp, ck, raw = build(root, ["# none"], extra)
    (root / extra).write_text("mesh-v2")
    print("arm2 (tracked control) dirty after edit:", fp.dirty(ck))
# Arm 3: ignored extra input ABSENT (as in a fresh Change worktree)
with tempfile.TemporaryDirectory() as d:
    root = Path(d)
    fp, ck, raw = build(root, ["cases/meshes/*.msh", "data/demo/raw/"], extra)
    (root / extra).unlink()
    c = fp.content(ck)
    print("arm3 absent ignored extra-input:", c["inputs"]["demo"]["extraInputs"], "| dirty:", fp.dirty(ck))
