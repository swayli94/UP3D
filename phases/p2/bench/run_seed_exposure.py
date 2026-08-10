"""Which cases does n_picard_seed = 0 fail on, and WHY?

Asked directly by the user 2026-08-05 after the M1 re-measurement found the seed
default's 5 -> 0 flip breaking NACA0012 M0.80 medium
(docs/dev_phase_two/20260805-1200-m1-remeasure.md sec 2).

HYPOTHESIS, WRITTEN BEFORE RUNNING. GS3.3b's own comment in tests/test_p8_newton.py
scoped its seed change with a stated safety condition:

    "It works because the ramp's FIRST level is subcritical (m_start 0.70), where
     Newton needs no help; do NOT carry this to a recipe whose first level is
     supercritical without re-measuring."

The M1 gate solves at the TARGET Mach directly -- M0.80, supercritical -- with no
ramp, so it is exactly the recipe that comment warned about. Predictions:

  P1  M0.80 medium, no ramp, seed 0   -> FAILS            (already measured)
  P2  M0.80 medium, WITH a ramp from 0.70, seed 0 -> CONVERGES.  If so the mechanism
      is not "the seed" but "no seed AND a supercritical first solve", and GS3.3b's
      safety condition is exactly right and was simply not carried to the call sites.
  P3  If P2 also fails, the seed matters even from a subcritical start, and GS3.3b's
      stated safety condition is INSUFFICIENT -- a stronger and worse result.

And the instrument gap I admitted in the LE factorial verdict is closed here: the
clamped cells are LOCATED, not just counted. CLAUDE.md's classification discipline
asks for "n_limited / n_floored > 0 -- AND WHERE THOSE CELLS ARE"; the previous two
rounds reported counts only, so "it hits the cap" could not be turned into "it hits
the cap at the leading edge / at the shock".

The mask is reconstructed with the solver's OWN condition (newton.py:378-380):
q2n = |grad phi|^2 / u_inf^2, then limit_q2_field(q2n, m_inf, m_cap) != q2n.

Outputs (TRACKED): bench/gate_results/seed_exposure.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.physics.isentropic import limit_q2_field                # noqa: E402
from pyfp3d.post.surface import _element_gradients_and_centroids     # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from run_le14_common_root import classify_failure                   # noqa: E402

CSV = os.path.join(_GATE, "seed_exposure.csv")
M_CAP = 3.0                        # solve_newton_lifting's default
ALPHA = 1.25
#: the M1 gate's recipe verbatim (bench/run_m1_gate.py). ★ alpha_deg is a DRIVER
#: argument, not a newton_kw one -- the first version of this script put it in the
#: dict and passed the dict as newton_kw, so the ramp legs died with
#: "solve_newton_transonic() missing 1 required positional argument: 'alpha_deg'".
#: Tenth wrong-from-structure call this phase; NK is the newton_kw subset.
RECIPE = dict(alpha_deg=ALPHA, upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6,
              freeze_refresh_max=8, precond="direct", direct_refactor_every=4,
              n_newton_max=80)
NK = {k: v for k, v in RECIPE.items() if k != "alpha_deg"}
#: (tag, level, m_inf, seed, ramp) -- ramp=True goes through m_start 0.70, i.e. a
#: SUBCRITICAL first level, which is the condition GS3.3b's comment names.
LEGS = (
    ("m080_medium_seed0_noramp", "medium", 0.80, 0, False),   # P1: the failure
    ("m080_medium_seed5_noramp", "medium", 0.80, 5, False),   # control
    ("m080_medium_seed0_ramp",   "medium", 0.80, 0, True),    # P2: THE discriminator
    ("m080_medium_seed5_ramp",   "medium", 0.80, 5, True),    # its control
    ("m080_coarse_seed0_noramp", "coarse", 0.80, 0, False),   # why does coarse live?
    ("m072_medium_seed0_noramp", "medium", 0.72, 0, False),   # in-envelope (M1a's M)
)


def clamp_map(mc, phi, m_inf):
    """Reconstruct the limiter mask with the solver's own condition and locate it.

    Returns (n_lim, frac_by_band, x_of_max, m_max) where the bands are chordwise
    thirds of the airfoil plus an 'offbody' bucket for elements whose centroid is
    outside [0, 1] in x -- so "at the leading edge" and "in the shock" and "out in
    the field" are distinguishable rather than pooled.
    """
    grad, cen = _element_gradients_and_centroids(mc.nodes, mc.elements, phi)
    q2 = np.einsum("ij,ij->i", grad, grad)
    lim = limit_q2_field(q2, m_inf, M_CAP) == q2
    bad = ~lim
    m2 = q2 * 0.0
    # local Mach^2 from the isentropic relation is monotone in q2, so for LOCATING
    # the peak the raw q2 ordering is enough; the reported m_max comes from the
    # solver itself, not from here.
    x = cen[:, 0]
    buckets = {"LE_0_15": (0.0, 0.15), "MID_15_85": (0.15, 0.85),
               "TE_85_100": (0.85, 1.0)}
    out = {}
    n_bad = int(bad.sum())
    for name, (lo, hi) in buckets.items():
        m = bad & (x >= lo) & (x < hi)
        out[f"clamp_frac_{name}"] = (round(float(m.sum()) / n_bad, 4)
                                     if n_bad else 0.0)
    m = bad & ((x < 0.0) | (x >= 1.0))
    out["clamp_frac_offbody"] = round(float(m.sum()) / n_bad, 4) if n_bad else 0.0
    i = int(np.argmax(q2))
    out["n_lim_reconstructed"] = n_bad
    out["x_of_peak_q2"] = round(float(x[i]), 4)
    out["peak_q2"] = round(float(q2[i]), 5)
    del m2
    return out


def main():
    rows = []
    print("n_picard_seed exposure: which cases fail with seed 0, and why?")
    print("hypothesis (pre-written): the binding condition is not the seed but "
          "'no seed AND a supercritical first solve' -- GS3.3b's own caveat.\n")
    for tag, level, m_inf, seed, ramp in LEGS:
        path = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d",
                            f"{level}.msh")
        if not os.path.exists(path):
            print(f"  {tag}: mesh missing"); continue
        mc, wc = cut_wake(read_mesh(path))
        t0 = time.perf_counter()
        if ramp:
            r = solve_newton_transonic(
                mc, wc, m_inf=m_inf, alpha_deg=ALPHA, m_start=0.70, dm=0.05,
                dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-4,
                newton_kw=dict(NK, n_picard_seed=seed))
        else:
            r = solve_newton_lifting(mc, wc, m_inf=m_inf,
                                     n_picard_seed=seed, **RECIPE)
        wall = time.perf_counter() - t0
        hist = np.asarray(r.get("residual_history", []), dtype=float)
        ch = np.asarray(r.get("clamp_history", []), dtype=float)
        nlim, nflr = int(r.get("n_limited") or 0), int(r.get("n_floored") or 0)
        conv = bool(r.get("converged"))
        if conv:
            mode, ev = "", f"accept_reason={r.get('accept_reason')!r}"
        else:
            mode, ev, _d, _rv = classify_failure(
                hist, ch, np.asarray(r.get("F_history", []), dtype=float),
                int(r.get("n_gmres_stalled") or 0), str(r.get("accept_reason")),
                nlim, nflr)
        loc = clamp_map(mc, np.asarray(r["phi"]), m_inf)
        # when does clamping FIRST appear? distinguishes "the starting field is
        # already past the cap" from "the iteration wanders there".
        first = None
        if ch.ndim == 2 and len(ch):
            nz = np.nonzero(ch[:, 0] + ch[:, 1] > 0)[0]
            first = int(nz[0]) if len(nz) else None
        row = dict(tag=tag, level=level, m_inf=m_inf, n_picard_seed=seed,
                   ramp=ramp, converged=conv, res_final=float(hist[-1]),
                   n_newton=len(hist), n_limited=nlim, n_floored=nflr,
                   m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                   first_clamped_step=first, failure_mode=mode,
                   mode_evidence=ev, accept_reason=r.get("accept_reason"),
                   wall_s=round(wall, 1), **loc)
        rows.append(row)
        print(f"  {tag:28} conv={str(conv):5} |R|={row['res_final']:.2e} "
              f"nk={len(hist):3d} lim/flr={nlim}/{nflr} M_max={row['m_max']:.4f} "
              f"first_clamp@{first} ({wall:.0f}s)", flush=True)
        if nlim or nflr or not conv:
            print(f"      MODE={mode or '-'}  clamped by band: "
                  f"LE {loc['clamp_frac_LE_0_15']:.2f} / "
                  f"MID {loc['clamp_frac_MID_15_85']:.2f} / "
                  f"TE {loc['clamp_frac_TE_85_100']:.2f} / "
                  f"offbody {loc['clamp_frac_offbody']:.2f}   "
                  f"peak q2 at x = {loc['x_of_peak_q2']}", flush=True)
    keys = sorted({k for r in rows for k in r})
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    by = {r["tag"]: r for r in rows}
    print("\n=== reading the pre-written predictions ===")
    a = by.get("m080_medium_seed0_noramp")
    b = by.get("m080_medium_seed0_ramp")
    if a and b:
        if not a["converged"] and b["converged"]:
            print("  P2 HOLDS: seed 0 fails without a ramp and CONVERGES with one")
            print("  => the binding condition is 'no seed AND a supercritical first")
            print("     solve', not the seed alone. GS3.3b's stated safety condition")
            print("     is correct; it was simply never carried to the call sites.")
        elif not a["converged"] and not b["converged"]:
            print("  P3 HOLDS (the worse one): seed 0 fails even from a SUBCRITICAL")
            print("  start => GS3.3b's safety condition is INSUFFICIENT.")
        else:
            print("  neither prediction as written -- read the rows above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
