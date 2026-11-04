import os
from pathlib import Path


def pytest_load_initial_conftests(early_config, parser, args):
    evidence_directory = os.environ.get("UP3D_EVIDENCE_DIR")
    if evidence_directory:
        report = Path(evidence_directory).resolve() / "pytest.junit.xml"
        early_config.option.xmlpath = str(report)
