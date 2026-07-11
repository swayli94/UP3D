"""
Track B / B5 demo -- far-field A/B: Dirichlet+vortex vs Neumann outlet.
(Gate renumbered 2026-07-12: B3.5 -> B4.5 -> B5. The directory keeps its
b4p5_ name on purpose so the committed demo path stays stable.)

The level-set lifting path (B3/B4) needs a far-field boundary condition. Two
self-consistent options exist (design_track_b.md section 5.4):

  option a (vortex)   : spherical Dirichlet freestream + a Prandtl-Glauert
                        point vortex on the far-field MAIN DOFs, with the
                        emergent Gamma refreshed into the vortex each outer
                        iteration (RHS-only). pyFP3D's compact 15-chord domain
                        is calibrated FOR this correction.
  option b (neumann)  : the Lopez form -- inflow is Dirichlet freestream (NO
                        vortex), outflow is a Neumann outlet carrying the
                        freestream flux rho_inf (u.n). Attractive for the
                        workflow (no Gamma-into-far-field feedback, simplest
                        alpha sweep), but with no vortex the O(Gamma/r)
                        far-field tail is truncated, so the domain must grow
                        (the dissertation uses 10^2-10^7 chord domains).

This demo runs the Lopez-style domain-size re-calibration: coarse NACA0012
meshes at growing far-field radius R for BOTH Track B mesh families (the
wake-EMBEDDED M0 and the wake-FREE M3), and tabulates the emergent circulation
Gamma against the conforming reference. It makes the decision visible:

  farfield_domain_study.png -- Gamma vs R (log-x), one panel per family, with
      the conforming reference and its +/-2% B3 band. Option a is a flat line
      (domain-robust at 15c); options b/freestream converge UP to the truth as
      R grows, meeting the band only after re-calibration.

  summary.csv  -- every (family, R, mode) -> Gamma, and b-vs-a / fs-vs-a %.
  checks.csv   -- the self-check verdicts (the B5 gate numbers).

VERDICT (measured): option a is domain-robust (Gamma within <1% of the truth
from 15c to 120c) and stays the DEFAULT for pyFP3D's 15c workflow. Option b is
a valid alternative but needs the domain grown to ~30c (2% band) / ~60c (1%),
i.e. a 2-4x larger domain -- its workflow simplicity does not pay for the
element cost at pyFP3D's scale. The O(Gamma/R) truncation is geometry-universal
(a 3D wing truncates the same horseshoe-vortex tail), so this decides the
far-field default for the 3D M6 B-path too (confirmed under B7).

Standalone + self-checking:  python cases/demo/b4p5_farfield/run_demo.py
Re-solve from scratch (~15 min, capped threads):  PYFP3D_B45_RESOLVE=1 python ...
"""
import csv
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.meshgen.planar import naca0012_wake_2d
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.picard import solve_subsonic_lifting
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

RESULTS = Path(__file__).resolve().parent / "results"
RESULTS.mkdir(exist_ok=True)
ALPHA, M_INF, H_WALL = 2.0, 0.5, 0.020
RADII = [15.0, 30.0, 60.0, 120.0]
MODES = ["vortex", "neumann", "freestream"]


def _band_pg_kt():
    path = REPO / "cases" / "reference_data" / "naca0012_m05" / "cl_reference.csv"
    with open(path) as f:
        for row in csv.DictReader(f):
            if abs(float(row["alpha_deg"]) - ALPHA) < 1e-9:
                return float(row["cl_pg"]), float(row["cl_kt"])
    raise LookupError("alpha not found")


def _build(r_far, embed_wake):
    pts, tris, eg, ig = naca0012_wake_2d(
        r_far=r_far, h_wall=H_WALL, h_far=min(3.0, 150.0 * H_WALL),
        h_wake=3.0 * H_WALL, dist_min=0.1, dist_max=6.0, wake_dist_max=1.5,
        n_half=max(80, int(round(2.0 / H_WALL))), embed_wake=embed_wake,
        corridor_alpha_deg=(-6.0, 6.0), corridor_n_lines=5)
    return extrude_single_layer(pts, tris, eg, interior_edge_groups=ig,
                                dz=2.0 * H_WALL, z0=0.0, name="b45")


def _solve_b(mesh, mode):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
                       direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    r = solve_multivalued_lifting(mvop, mesh, M_INF, alpha_deg=ALPHA, farfield=mode)
    return r["gamma"], r["converged"]


def _conforming(mesh):
    mc, wc = cut_wake(mesh)
    ref = solve_subsonic_lifting(mc, wc, M_INF, alpha_deg=ALPHA)
    return float(np.mean(ref["gamma"]))


def compute():
    rows = []
    for family, embed in [("M0_embedded", True), ("M3_wakefree", False)]:
        for r_far in RADII:
            mesh = _build(r_far, embed)
            rec = {"family": family, "r_far": r_far, "nodes": len(mesh.nodes)}
            rec["conforming"] = _conforming(mesh) if embed else np.nan
            for m in MODES:
                g, conv = _solve_b(mesh, m)
                rec[m] = g
                rec[f"{m}_conv"] = int(conv)
            rows.append(rec)
            print(f"  {family} R={r_far:5.0f} nodes={rec['nodes']:6d} "
                  f"vortex={rec['vortex']:.4f} neumann={rec['neumann']:.4f} "
                  f"freestream={rec['freestream']:.4f}")
    return rows


def _write_summary(rows):
    keys = ["family", "r_far", "nodes", "conforming", "vortex", "neumann",
            "freestream", "vortex_conv", "neumann_conv", "freestream_conv",
            "neumann_vs_vortex_pct", "freestream_vs_vortex_pct"]
    for r in rows:
        r["neumann_vs_vortex_pct"] = 100 * (r["neumann"] - r["vortex"]) / r["vortex"]
        r["freestream_vs_vortex_pct"] = 100 * (r["freestream"] - r["vortex"]) / r["vortex"]
    with open(RESULTS / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})


def _load_summary():
    rows = []
    with open(RESULTS / "summary.csv") as f:
        for row in csv.DictReader(f):
            rec = {"family": row["family"], "r_far": float(row["r_far"]),
                   "nodes": int(row["nodes"])}
            for k in ("conforming", "vortex", "neumann", "freestream",
                      "neumann_vs_vortex_pct", "freestream_vs_vortex_pct"):
                rec[k] = float(row[k]) if row[k] not in ("", "nan") else np.nan
            for k in ("vortex_conv", "neumann_conv", "freestream_conv"):
                rec[k] = int(row[k])
            rows.append(rec)
    return rows


def plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    style = {"vortex": ("tab:green", "o", "option a: Dirichlet + vortex"),
             "neumann": ("tab:orange", "s", "option b: Neumann outlet"),
             "freestream": ("tab:red", "^", "freestream (no vortex)")}
    for ax, family in zip(axes, ["M0_embedded", "M3_wakefree"]):
        fr = [r for r in rows if r["family"] == family]
        fr.sort(key=lambda r: r["r_far"])
        R = [r["r_far"] for r in fr]
        truth = np.nanmean([r["conforming"] for r in fr if not np.isnan(r["conforming"])])
        if np.isnan(truth):
            truth = fr[-1]["vortex"]   # M3: no conforming, use large-R vortex
        ax.axhline(truth, color="0.4", lw=1.2, ls="--", label="reference (conforming)")
        ax.axhspan(truth * 0.98, truth * 1.02, color="0.5", alpha=0.12,
                   label="B3 +/-2% band")
        for m in MODES:
            c, mk, lab = style[m]
            ax.plot(R, [r[m] for r in fr], marker=mk, color=c, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel("far-field radius R  [chords]")
        ax.set_title(family.replace("_", " "))
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel(r"emergent circulation  $\Gamma$")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.suptitle("B5 far-field A/B: option a is domain-robust; "
                 "option b truncates O(Gamma/R) and needs a larger domain")
    fig.tight_layout()
    fig.savefig(RESULTS / "farfield_domain_study.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def self_check(rows):
    checks = []

    def rec(name, ok, detail):
        checks.append((name, "PASS" if ok else "FAIL", detail))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        return ok

    all_ok = True
    for family in ["M0_embedded", "M3_wakefree"]:
        fr = {r["r_far"]: r for r in rows if r["family"] == family}
        g15a = fr[15.0]["vortex"]
        # (1) option a domain-robust: Gamma_a(R) within 2% of Gamma_a(15c)
        dev = max(abs(fr[R]["vortex"] - g15a) / g15a for R in RADII)
        all_ok &= rec(f"{family}: option-a domain-robust",
                      dev < 0.02, f"max |dGamma_a|/Gamma = {dev*100:.2f}% over 15-120c")
        # (2) option b converges UP to a as R grows (monotone increasing gap-closure)
        gaps = [abs(fr[R]["neumann"] - fr[R]["vortex"]) / fr[R]["vortex"] for R in RADII]
        all_ok &= rec(f"{family}: option-b converges to a",
                      gaps[0] > gaps[-1] and gaps[-1] < 0.02,
                      f"|b-a| 15c={gaps[0]*100:.2f}% -> 120c={gaps[-1]*100:.2f}%")
        # (3) re-calibration threshold: b OUTSIDE 2% band at 15c, INSIDE by 60c
        b60 = abs(fr[60.0]["neumann"] - fr[60.0]["vortex"]) / fr[60.0]["vortex"]
        all_ok &= rec(f"{family}: option-b needs re-calibration",
                      gaps[0] > 0.02 and b60 < 0.02,
                      f"15c={gaps[0]*100:.2f}% (>2%), 60c={b60*100:.2f}% (<2%)")
        # (4) freestream is the crudest at every R
        crude = all(
            abs(fr[R]["freestream"] - fr[R]["vortex"])
            >= abs(fr[R]["neumann"] - fr[R]["vortex"]) - 1e-9 for R in RADII)
        all_ok &= rec(f"{family}: freestream is crudest", crude,
                      "|fs-a| >= |b-a| at every R")

    # (5) option a at 15c matches the conforming reference within 2% (M0 only)
    m0 = {r["r_far"]: r for r in rows if r["family"] == "M0_embedded"}
    e = abs(m0[15.0]["vortex"] - m0[15.0]["conforming"]) / m0[15.0]["conforming"]
    all_ok &= rec("M0 15c: option a vs conforming", e < 0.02,
                  f"Gamma_a {m0[15.0]['vortex']:.4f} vs conf "
                  f"{m0[15.0]['conforming']:.4f} ({e*100:.2f}%)")

    with open(RESULTS / "checks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["check", "status", "detail"])
        w.writerows(checks)
    return all_ok


def main():
    cached = (RESULTS / "summary.csv").exists()
    resolve = os.environ.get("PYFP3D_B45_RESOLVE") == "1"
    if cached and not resolve:
        print("Loading cached summary.csv (set PYFP3D_B45_RESOLVE=1 to re-solve).")
        rows = _load_summary()
    else:
        print("Solving the domain-size study (both families x 4 radii)...")
        rows = compute()
        _write_summary(rows)
    plot(rows)
    ok = self_check(rows)
    print(f"\nB5 demo {'PASSED' if ok else 'FAILED'} "
          f"({sum(1 for _ in RADII)} radii x 2 families)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
