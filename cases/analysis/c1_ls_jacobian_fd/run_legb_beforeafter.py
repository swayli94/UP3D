"""
B20 GB20.3/GB20.4 — before/after on committed cases, side vs main density.

The operator carries the mode, the driver is identical, so this is a clean A/B.
Prints CSV rows. Subsonic (M0.5) is cheap; transonic (M6 coarse M0.84 via the
B15 Newton ramp) is the heavier leg.

Usage:
  python run_legb_beforeafter.py subsonic     # GB20.3  (B3/B4 Gamma)
  python run_legb_beforeafter.py transonic    # GB20.4  (M6 coarse gamma/M_max)
"""

# ---------------------------------------------------------------------------
# HISTORICAL (B20 close-out, 2026-07-18). This script performed a side-vs-main
# A/B through the `plain_density` knob on MultivaluedOperator. That knob was
# REMOVED when the main-field density was made permanent and non-optional
# (user-arbitrated): there is no longer a "side" mode to compare against, by
# design. The measurements it produced are committed in results/legb_*.csv and
# are the evidence; to re-run the A/B itself, check out commit 5369a84, where
# the knob still existed. The standing (non-historical) check that the Jacobian
# remains exact under the permanent fix lives in tests/test_b19_jacobian_3d.py.
# ---------------------------------------------------------------------------
raise SystemExit(
    "HISTORICAL: the plain_density knob this A/B toggles was removed when B20 "
    "was made permanent. Committed evidence: results/legb_*.csv. To reproduce "
    "the A/B, check out commit 5369a84.")

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import read_mesh                        # noqa: E402
from pyfp3d.meshgen.wing3d import x_te                          # noqa: E402
from pyfp3d.solve.newton_ls import solve_multivalued_newton_transonic  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting   # noqa: E402
from pyfp3d.wake import (                                       # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

RESULTS = Path(__file__).resolve().parent / "results"
MESHES = REPO_ROOT / "cases" / "meshes"
UI = 1.0


def naca_op(level, wakefree, plain_density):
    d = "naca0012_wakefree_2.5d" if wakefree else "naca0012_2.5d"
    mesh = read_mesh(MESHES / d / f"{level}.msh")
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
                       direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                     levelset=wls, plain_density=plain_density)


def m6_op(level, plain_density, alpha=3.06):
    mesh = read_mesh(MESHES / "onera_m6" / f"{level}.msh")
    a = np.radians(alpha)
    b = 1.1963
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(b), 0.0, b]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                     levelset=wls, plain_density=plain_density)


def subsonic():
    rows = []
    for mode in ("side", "main"):
        m0, op0 = naca_op("medium", wakefree=False, plain_density=mode)
        r0 = solve_multivalued_lifting(op0, m0, m_inf=0.5, alpha_deg=1.25,
                                       farfield="vortex", te_kutta="pressure")
        g_emb = float(r0["gamma"])
        m3, op3 = naca_op("medium", wakefree=True, plain_density=mode)
        r3 = solve_multivalued_lifting(op3, m3, m_inf=0.5, alpha_deg=1.25,
                                       farfield="vortex", te_kutta="pressure")
        g_free = float(r3["gamma"])
        n_mp = int(op0._mixed_plain_mask().sum())
        rows.append((mode, g_emb, g_free, abs(g_free - g_emb) / g_emb * 100,
                     n_mp))
        print(f"[{mode}] NACA medium M0.5: Gamma embedded {g_emb:.6f} "
              f"wakefree {g_free:.6f} ({abs(g_free-g_emb)/g_emb*100:.4f}%), "
              f"mixed_plain={n_mp}")
    gs = {m: (e, f) for m, e, f, _, _ in rows}
    print(f"  embedded Gamma side->main: {gs['side'][0]:.6f} -> "
          f"{gs['main'][0]:.6f} "
          f"({(gs['main'][0]-gs['side'][0])/gs['side'][0]*100:+.4f}%)")
    out = RESULTS / "legb_subsonic_ab.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mode", "gamma_embedded", "gamma_wakefree",
                    "wakefree_vs_embedded_pct", "n_mixed_plain"])
        for r in rows:
            w.writerow(r)
    print(f"wrote {out}")


def transonic():
    rows = []
    for mode in ("side", "main"):
        mesh, op = m6_op("coarse", plain_density=mode)
        r = solve_multivalued_newton_transonic(
            op, mesh, 0.84, m_start=0.5, dm=0.05, alpha_deg=3.06,
            farfield="freestream", n_seed=30, n_newton_max=40,
            tol_residual=1e-6, freeze_tol=1e-5, verbose=False)
        gamma = float(r.get("gamma"))
        mmax = float(r.get("mach2_max", np.nan)) ** 0.5
        conv = r.get("converged")
        mfin = r.get("m_final", r.get("m_last_converged"))
        n_mp = int(op._mixed_plain_mask().sum())
        rows.append((mode, mfin, conv, gamma, mmax, n_mp))
        print(f"[{mode}] M6 coarse ramp->M0.84: m_final={mfin} conv={conv} "
              f"gamma={gamma:.6f} M_max={mmax:.4f} mixed_plain={n_mp}")
    out = RESULTS / "legb_transonic_ab.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mode", "m_final", "converged", "gamma", "m_max",
                    "n_mixed_plain"])
        for r in rows:
            w.writerow(r)
    print(f"wrote {out}")


if __name__ == "__main__":
    RESULTS.mkdir(parents=True, exist_ok=True)
    {"subsonic": subsonic, "transonic": transonic}[sys.argv[1]]()
