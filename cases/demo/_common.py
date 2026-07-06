"""
Shared infrastructure for the per-phase demo scripts (cases/demo/*/run_demo.py).

Provides the light-mode chart style used by all demo artifacts, small
save/CSV helpers, and a CheckList recorder that turns each demo into a
self-checking evidence run: every figure is backed by a quantitative check
against the roadmap acceptance criterion, the results land in
results/checks.csv, and the script exits nonzero on any unexpected FAIL
(XFAIL entries -- documented open gates -- do not fail the run).

Everything is headless: matplotlib Agg, no GUI (roadmap Sec 0.1).
"""

import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

MESH_DIR = REPO_ROOT / "cases" / "meshes"
REFERENCE_DIR = REPO_ROOT / "cases" / "reference_data"

# ---------------------------------------------------------------------------
# Chart style: single light-mode palette for static PNG artifacts.
# Categorical hues are assigned in fixed order S1..S5; CRITICAL is reserved
# for gate thresholds / failure annotations, never used as a series color.
# ---------------------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
S1_BLUE = "#2a78d6"
S2_AQUA = "#1baf7a"
S3_YELLOW = "#eda100"
S4_ROSE = "#c65b8a"
S5_VIOLET = "#4a3aa7"
CRITICAL = "#d03b3b"

SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue",
    ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
)
SEQ_MAGMA = "magma"  # perceptually-uniform fallback for noise heatmaps
DIV_BLUE_RED = LinearSegmentedColormap.from_list(
    "div_blue_red",
    ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f2a9a8", "#e34948", "#8f1f1f"],
)


def apply_style():
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": BASELINE,
            "axes.labelcolor": INK_2,
            "axes.titlecolor": INK,
            "axes.titlesize": 12,
            "axes.titleweight": "semibold",
            "axes.labelsize": 10.5,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "grid.linestyle": "-",
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.frameon": False,
            "legend.fontsize": 9.5,
            "legend.labelcolor": INK_2,
            "font.family": "sans-serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "lines.linewidth": 2.0,
            "lines.markersize": 8,
        }
    )


def finish(fig, out_dir: Path, name: str):
    """Save a figure into the demo's results/ directory and close it."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path.relative_to(REPO_ROOT)}")


def write_csv(out_dir: Path, name: str, header: str, rows):
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(v) for v in row) + "\n")
    print(f"  wrote {path.relative_to(REPO_ROOT)}")


# ---------------------------------------------------------------------------
# Acceptance-check recorder
# ---------------------------------------------------------------------------
@dataclass
class CheckList:
    """Collects gate-criterion checks; a demo passes only if every non-XFAIL
    check passes. XFAIL marks a documented open gate (e.g. G1.6) whose
    failure is the expected, recorded state -- it never fails the run, and
    it *unexpectedly passing* is flagged for attention instead."""

    title: str
    checks: list = field(default_factory=list)

    def add(self, gate, name, value, criterion, passed, xfail=False, note=""):
        status = ("XPASS?" if passed else "XFAIL") if xfail else ("PASS" if passed else "FAIL")
        self.checks.append(
            {"gate": gate, "name": name, "value": value,
             "criterion": criterion, "status": status, "note": note}
        )

    def report(self, out_dir: Path) -> int:
        """Print the PASS/FAIL table, write checks.csv, return exit code."""
        w_gate = max(len(c["gate"]) for c in self.checks)
        w_name = max(len(c["name"]) for c in self.checks)
        w_val = max(len(str(c["value"])) for c in self.checks)
        w_crit = max(len(c["criterion"]) for c in self.checks)
        print(f"\n{self.title} -- acceptance checks")
        print("-" * (w_gate + w_name + w_val + w_crit + 18))
        for c in self.checks:
            print(f"  {c['gate']:<{w_gate}}  {c['name']:<{w_name}}  "
                  f"{str(c['value']):<{w_val}}  {c['criterion']:<{w_crit}}  {c['status']}")
        write_csv(out_dir, "checks.csv", "gate,check,value,criterion,status,note",
                  [(c["gate"], c["name"].replace(",", ";"), c["value"],
                    c["criterion"].replace(",", ";"), c["status"],
                    c["note"].replace(",", ";")) for c in self.checks])
        n_fail = sum(c["status"] == "FAIL" for c in self.checks)
        n_xfail = sum(c["status"] == "XFAIL" for c in self.checks)
        n_xpass = sum(c["status"] == "XPASS?" for c in self.checks)
        msg = f"{len(self.checks) - n_fail - n_xfail - n_xpass} passed"
        if n_xfail:
            msg += f", {n_xfail} xfailed (documented open gate)"
        if n_xpass:
            msg += f", {n_xpass} UNEXPECTEDLY passed an open gate -- investigate"
        if n_fail:
            msg += f", {n_fail} FAILED"
        print(msg)
        return 1 if n_fail else 0
