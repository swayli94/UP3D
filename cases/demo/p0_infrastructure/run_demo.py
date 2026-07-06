"""
P0 demo -- mesh infrastructure evidence (gates G0.1-G0.4, all closed).

What this shows, per docs/roadmap.md P0 and docs/demo_report.md:
  1. G0.1 volume exactness: sum(V_e) equals the analytic volume to machine
     precision on meshes with an exact polyhedral volume.
  2. G0.2 gradient exactness: the P1 element gradient reproduces the
     gradient of any linear field to machine precision, even on random
     (badly shaped) tets -- so every later discretization error is the
     PDE approximation, never the geometry kernels.
  3. G0.3 coloring validity: greedy element coloring is valid (no two
     same-color elements share a node) on every committed mesh family,
     with a small color count and balanced classes -> safe prange assembly.
  4. G0.4 VTK round-trip: write -> read of nodes/elements/fields is
     lossless (bit-exact within float64 formatting).

Standalone + self-checking:  python cases/demo/p0_infrastructure/run_demo.py
Outputs: cases/demo/p0_infrastructure/results/{*.png, summary.csv, checks.csv}
Exit code 0 iff every acceptance check passes. Runtime ~30 s.
"""

import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, GRID, INK_2, MESH_DIR, MUTED, REPO_ROOT,
    S1_BLUE, S2_AQUA, S3_YELLOW, S5_VIOLET, SURFACE,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402  (Agg set in _common)

from pyfp3d.mesh.coloring import greedy_coloring, validate_coloring  # noqa: E402
from pyfp3d.mesh.metrics import compute_tet_volumes  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.post.vtk_out import read_vtu, write_vtu  # noqa: E402
from tests.mesh_utils import (  # noqa: E402
    element_gradients_all, generate_structured_cube_mesh,
)
from tests.test_mesh_volume import (  # noqa: E402
    create_unit_cube_mesh, create_unit_octahedron_mesh,
)

OUT = Path(__file__).resolve().parent / "results"


# ---------------------------------------------------------------------------
# 1. G0.1 volume exactness
# ---------------------------------------------------------------------------
def demo_volume(checks):
    print("[1/4] G0.1 volume exactness")
    cases = []

    nodes, elements = create_unit_cube_mesh()
    cases.append(("unit cube, 5 tets", nodes, elements, 1.0))
    nodes, elements = create_unit_octahedron_mesh()
    cases.append(("octahedron, 8 tets", nodes, elements, 4.0 / 3.0))
    for n in (4, 16):
        nodes, elements = generate_structured_cube_mesh(n=n)
        cases.append((f"Kuhn cube n={n}, {len(elements)} tets", nodes, elements, 1.0))
    nodes, elements = create_unit_cube_mesh()
    cases.append(("unit cube scaled x2 (V=8)", 2.0 * nodes, elements, 8.0))

    rows, errs = [], []
    for name, nds, els, exact in cases:
        vol = compute_tet_volumes(nds, els)
        err = abs(vol.sum() - exact) / exact
        rows.append((name.replace(",", ";"), len(els), f"{vol.sum():.15f}",
                     exact, f"{err:.3e}"))
        errs.append(err)
    errs = np.array(errs)

    fig, ax = plt.subplots(figsize=(7.6, 3.6))
    y = np.arange(len(cases))[::-1]
    x_lo = max(errs.min() * 0.3, 1e-18)
    ax.hlines(y, x_lo, np.maximum(errs, x_lo), color=GRID, linewidth=1.0, zorder=2)
    ax.plot(np.maximum(errs, x_lo), y, "o", color=S1_BLUE, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
    for yi, e in zip(y, errs):
        label = f"{e:.1e}" if e > 0 else "0 (exact)"
        ax.text(max(e, x_lo) * 2.2, yi, label, va="center", color=INK_2, fontsize=9)
    ax.axvline(1e-12, color=CRITICAL, linewidth=1.2, linestyle="--", zorder=2)
    ax.set_ylim(-0.75, len(cases) - 0.4)
    ax.text(1.4e-12, -0.55, "gate: 1e-12", color=CRITICAL, fontsize=9, va="center")
    ax.set_xscale("log")
    ax.set_xlim(x_lo, 3e-11)
    ax.set_yticks(y)
    ax.set_yticklabels([c[0] for c in cases], color=INK_2)
    ax.set_xlabel("relative volume error  |sum V_e - V_exact| / V_exact")
    ax.set_title("G0.1: discrete volume matches the analytic volume to machine zero")
    ax.grid(axis="y", visible=False)
    finish(fig, OUT, "g01_volume_exactness.png")
    write_csv(OUT, "g01_volumes.csv", "case,n_tets,sum_volumes,exact,rel_error", rows)

    checks.add("G0.1", "max relative volume error", f"{errs.max():.2e}",
               "< 1e-12", bool(errs.max() < 1e-12))


# ---------------------------------------------------------------------------
# 2. G0.2 gradient exactness on linear fields
# ---------------------------------------------------------------------------
def demo_gradient(checks):
    print("[2/4] G0.2 gradient exactness")
    a = np.array([1.0, 2.0, 3.0])

    def grad_errors(nodes, elements):
        f = nodes @ a + 0.5
        g = element_gradients_all(nodes, elements, f)
        return np.linalg.norm(g - a, axis=1)

    # random tets: same construction as the G0.2 gate test (arbitrary vertex
    # picks -> includes badly-shaped and inverted tets on purpose)
    rng = np.random.RandomState(42)
    rnodes = rng.rand(30, 3)
    relements = np.array([rng.choice(30, 4, replace=False) for _ in range(20)],
                         dtype=np.int32)
    err_random = grad_errors(rnodes, relements)

    nodes, elements = generate_structured_cube_mesh(n=8)
    err_cube = grad_errors(nodes, elements)

    mesh = read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh")
    err_naca = grad_errors(mesh.nodes, mesh.elements)

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    bins = np.linspace(-17, -11, 49)
    for err, color, label in (
        (err_cube, S1_BLUE, f"Kuhn cube n=8 ({len(err_cube)} tets)"),
        (err_naca, S2_AQUA, f"gmsh NACA0012 coarse ({len(err_naca)} tets)"),
        (err_random, S3_YELLOW, f"random tets ({len(err_random)}, incl. slivers)"),
    ):
        ax.hist(np.log10(np.maximum(err, 1e-17)), bins=bins, histtype="step",
                linewidth=2.0, color=color, label=label)
    ax.set_yscale("log")
    ax.set_ylim(top=ax.get_ylim()[1] * 12)  # headroom for the legend
    ax.axvline(-12, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax.text(-11.95, ax.get_ylim()[1] * 0.35, "gate: 1e-12", color=CRITICAL, fontsize=9)
    ax.set_xlabel("log10 |grad(a.x+b) - a|  per element")
    ax.set_ylabel("element count")
    ax.set_title("G0.2: P1 element gradient is exact for linear fields on any tet")
    ax.legend(loc="upper left")
    finish(fig, OUT, "g02_gradient_exactness.png")

    worst = max(err_random.max(), err_cube.max(), err_naca.max())
    checks.add("G0.2", "max linear-field gradient error", f"{worst:.2e}",
               "< 1e-12", bool(worst < 1e-12))


# ---------------------------------------------------------------------------
# 3. G0.3 element coloring across the committed mesh families
# ---------------------------------------------------------------------------
def demo_coloring(checks):
    print("[3/4] G0.3 coloring validity")
    cases = []
    nodes, elements = generate_structured_cube_mesh(n=8)
    cases.append(("Kuhn cube n=8", elements))
    for family, level in (("naca0012_2.5d", "coarse"), ("cylinder_2.5d", "coarse"),
                          ("sphere_shell", "coarse")):
        mesh = read_mesh(MESH_DIR / family / f"{level}.msh")
        cases.append((f"{family} {level}", mesh.elements))

    rows, all_valid, focus = [], True, None
    for name, elements in cases:
        t0 = time.time()
        colors, n_colors = greedy_coloring(elements)
        valid = validate_coloring(elements, colors)
        all_valid &= valid
        sizes = np.bincount(colors, minlength=n_colors)
        rows.append((name.replace(",", ";"), len(elements), n_colors, valid,
                     sizes.min(), sizes.max(), f"{time.time() - t0:.2f}"))
        if "naca" in name:
            focus = (name, sizes)
        print(f"    {name}: {n_colors} colors, valid={valid}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.2))
    labels = [r[0] for r in rows]
    ncols = [r[2] for r in rows]
    y = np.arange(len(rows))[::-1]
    ax1.barh(y, ncols, height=0.55, color=S1_BLUE)
    for yi, nc, r in zip(y, ncols, rows):
        ax1.text(nc + 0.4, yi, f"{nc} colors, valid", va="center", color=INK_2, fontsize=9.5)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, color=INK_2)
    ax1.set_xlim(0, max(ncols) * 1.45)
    ax1.set_xlabel("number of colors (greedy)")
    ax1.set_title("Small, bounded color count per family")
    ax1.grid(axis="y", visible=False)

    name, sizes = focus
    ax2.bar(np.arange(len(sizes)), sizes, color=S1_BLUE, width=0.7)
    ax2.axhline(sizes.mean(), color=MUTED, linewidth=1.2, linestyle="--")
    ax2.text(len(sizes) - 0.4, sizes.mean() * 1.04, "mean", color=MUTED,
             fontsize=9, ha="right")
    ax2.set_xlabel("color class")
    ax2.set_ylabel("elements in class")
    ax2.set_title("naca0012 coarse: balanced classes -> even prange batches")
    ax2.grid(axis="x", visible=False)
    fig.suptitle("G0.3: no two same-color elements share a node (verified), "
                 "enabling race-free colored assembly",
                 fontsize=12.5, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    finish(fig, OUT, "g03_coloring.png")
    write_csv(OUT, "g03_coloring.csv",
              "mesh,n_elements,n_colors,valid,min_class,max_class,seconds", rows)

    checks.add("G0.3", "coloring valid on all 4 mesh families", all_valid,
               "validate_coloring == True", bool(all_valid))


# ---------------------------------------------------------------------------
# 4. G0.4 VTK round-trip
# ---------------------------------------------------------------------------
def demo_roundtrip(checks):
    print("[4/4] G0.4 VTK round-trip")
    mesh = read_mesh(MESH_DIR / "cylinder_2.5d" / "coarse.msh")
    nodes, elements = mesh.nodes, mesh.elements
    fields = {
        "phi_linear": nodes @ np.array([1.0, 2.0, 3.0]) + 0.5,
        "r": np.linalg.norm(nodes, axis=1),
    }

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "roundtrip.vtu"
        write_vtu(path, nodes, elements, point_data=fields)
        nodes_r, elements_r, fields_r = read_vtu(path)

    diffs = [
        ("node coordinates", float(np.max(np.abs(nodes - nodes_r)))),
        ("element connectivity", float(np.max(np.abs(elements - elements_r)))),
    ]
    for name, data in fields.items():
        diffs.append((f"field '{name}'",
                      float(np.max(np.abs(data - fields_r[name].ravel())))))

    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    y = np.arange(len(diffs))[::-1]
    vals = np.array([max(d[1], 1e-18) for d in diffs])
    ax.hlines(y, 1e-18, vals, color=GRID, linewidth=1.0, zorder=2)
    ax.plot(vals, y, "o", color=S5_VIOLET, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
    for yi, (name, d) in zip(y, diffs):
        ax.text(max(d, 1e-18) * 2.5, yi, "0 (bit-exact)" if d == 0 else f"{d:.1e}",
                va="center", color=INK_2, fontsize=9.5)
    ax.axvline(1e-15, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax.set_ylim(-0.75, len(diffs) - 0.4)
    ax.text(1.35e-15, -0.55, "gate: 1e-15", color=CRITICAL, fontsize=9, va="center")
    ax.set_xscale("log")
    ax.set_xlim(1e-18, 1e-13)
    ax.set_yticks(y)
    ax.set_yticklabels([d[0] for d in diffs], color=INK_2)
    ax.set_xlabel("max |written - read| after .vtu round-trip (17k-tet gmsh mesh)")
    ax.set_title("G0.4: VTK I/O is lossless")
    ax.grid(axis="y", visible=False)
    finish(fig, OUT, "g04_vtk_roundtrip.png")

    worst = max(d[1] for d in diffs)
    checks.add("G0.4", "max round-trip difference", f"{worst:.1e}",
               "< 1e-15", bool(worst < 1e-15))
    return diffs


def main():
    apply_style()
    t0 = time.time()
    checks = CheckList("P0 mesh infrastructure (G0.1-G0.4)")

    demo_volume(checks)
    demo_gradient(checks)
    demo_coloring(checks)
    demo_roundtrip(checks)

    write_csv(OUT, "summary.csv", "metric,value",
              [(c["gate"] + " " + c["name"].replace(",", ";"), c["value"])
               for c in checks.checks] + [("runtime_seconds", f"{time.time()-t0:.1f}")])
    code = checks.report(OUT)
    print(f"done in {time.time() - t0:.1f}s -> {OUT.relative_to(REPO_ROOT)}/")
    sys.exit(code)


if __name__ == "__main__":
    main()
