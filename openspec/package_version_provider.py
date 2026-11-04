import importlib.metadata
import json
import os
import subprocess


PRODUCT_ENVIRONMENT = "pyfp3d-core"
PRODUCT_PACKAGES = {"python", "pytest"}


def _product_package_version(name: str) -> str:
    if name == "python":
        code = "import platform; print(platform.python_version())"
        arguments = []
    else:
        code = (
            "import importlib.metadata as metadata, sys; "
            "print(metadata.version(sys.argv[1]))"
        )
        arguments = [name]
    env = dict(os.environ)
    env["PYTHONNOUSERSITE"] = "1"
    completed = subprocess.run(
        [
            "conda",
            "run",
            "-n",
            PRODUCT_ENVIRONMENT,
            "--no-capture-output",
            "python",
            "-c",
            code,
            *arguments,
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    lines = completed.stdout.strip().splitlines()
    if len(lines) != 1 or not lines[0]:
        raise RuntimeError(f"Cannot read {name} from {PRODUCT_ENVIRONMENT}")
    return lines[0]


def _devflow_version() -> str:
    distribution = importlib.metadata.distribution("openspec-devflow")
    raw_direct_url = distribution.read_text("direct_url.json")
    if raw_direct_url is None:
        raise RuntimeError("openspec-devflow has no direct_url.json")
    direct_url = json.loads(raw_direct_url)
    vcs_info = direct_url.get("vcs_info")
    if not isinstance(vcs_info, dict):
        raise RuntimeError("openspec-devflow direct_url.json has no vcs_info")
    commit = vcs_info.get("commit_id")
    if not isinstance(commit, str) or not commit:
        raise RuntimeError("openspec-devflow direct_url.json has no vcs_info.commit_id")
    return f"{distribution.version}+commit.{commit}"


def package_version(name: str) -> str:
    if name == "openspec-devflow":
        return _devflow_version()
    if name in PRODUCT_PACKAGES:
        return _product_package_version(name)
    raise ValueError(f"Unsupported environment package: {name}")
