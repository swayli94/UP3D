"""How pytest's JUnit report encodes skip and strict xfail (the counts devflow reads).

usage (up3d env): python exp_junit_xfail.py
Builds a throwaway test file in a temp dir: 1 pass, 1 skip, 1 strict xfail.
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SRC = '''import pytest
def test_ok(): pass
@pytest.mark.skip(reason="figure export off")
def test_fig(): pass
@pytest.mark.xfail(strict=True, reason="known negative result")
def test_neg(): assert False
'''
with tempfile.TemporaryDirectory() as d:
    (Path(d) / "test_x.py").write_text(SRC)
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                        "--junitxml=j.xml", "test_x.py"], cwd=d, capture_output=True, text=True)
    print("console:", r.stdout.strip().splitlines()[-1])
    xml = (Path(d) / "j.xml").read_text()
    suite = re.search(r"<testsuite [^>]*>", xml).group(0)
    print("testsuite counts:", " ".join(re.findall(r'(?:tests|skipped|failures|errors)="\d+"', suite)))
    print("skipped types:", re.findall(r'<skipped type="([^"]*)"', xml))
