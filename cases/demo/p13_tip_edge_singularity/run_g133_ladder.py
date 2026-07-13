"""P13/G13.3 -- where is the REMAINING 3D grid-convergence blocker?

G13.2 removed the wake sheet's free-tip-edge singularity (spanwise loading
taper). It did NOT deliver 3D grid convergence: on a refinement ladder the M6
lift increments still GROW. Growing increments under *uniform* refinement are
the signature of a singularity that is still being resolved -- so a second
object must still be there.

This script localises it with a three-region box study on the SELF-SIMILAR
ladder, and the answer is that the remaining singularity is on the WALL, not
in the wake:

  region                              verdict
  wake free edge (the G13.2 fix)      bounded   -- the fix holds under refinement
  wing interior (p99, control)        bounded   -- the ordinary field is converged
  tip-cap edge (WALL)                 DIVERGES  -- the remaining blocker

Root cause: `meshgen/wing3d.py` builds a FLAT tip cap (a documented deliberate
simplification -- the real ONERA M6 has a rounded one). A flat cap meets the
upper/lower surfaces at a sharp convex edge, which in potential flow is an edge
singularity of exactly the kind P13 exists to remove -- only on the BODY rather
than on the wake. This is NOT a P11 (curved wall ELEMENT) problem: isoparametric
elements cannot regularize a genuinely sharp geometric EDGE. The GEOMETRY is
what is wrong, so the fix is Track M (round the cap), not Track P.

★ WHY THE LADDER MATTERS. `generate_onera_m6.py` used to CLAMP the far field
(`h_far = min(2.5, 120*h_wall)`). The clamp bit only at `coarse`, so
coarse->medium refined the far field 1.39x while the wall refined 2x: `coarse`
was not on the refinement ray and ANY three-point Richardson over
(coarse, medium, fine) was invalid -- including P9/G9.1's. This study therefore
runs on `RICHARDSON_LADDER = (coarse_ss, medium, fine)`, which refines by
exactly 2.000 in every length scale.

★ EVIDENCE PROVENANCE (read before trusting a number). The heavy solves are
cached to results/*.npz (gitignored). The `fine` point is adopted from a cache
whose PRODUCING SCRIPT WAS NEVER COMMITTED, so its exact solver settings are
unrecoverable; see FINE_PROVENANCE below. Every level's provenance is written
into the CSV so a reader can tell measured-here from adopted.

Standalone (coarse_ss solves in ~1 min; medium/fine come from cache -- a cold
fine solve is ~34 min):
  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python cases/demo/p13_tip_edge_singularity/run_g133_ladder.py
"""
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from cases.demo._common import (  # noqa: E402
    CRITICAL, CheckList, INK_2, S1_BLUE, S2_AQUA, S3_YELLOW,
    apply_style, finish, plt, write_csv,
)
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.kernels.jacobian import PicardOperator  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te  # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field  # noqa: E402
from pyfp3d.post.surface import cl_kj_3d, planform_area  # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes" / "onera_m6"

# The G13.2 shipped tip model, probed subsonically (no shock / limiter physics
# in the way -- the edge signal is geometric).
ALPHA, M = 3.06, 0.5
FORM, R_C_FRAC = "vanish_smooth", 0.05
LADDER = ("coarse_ss", "medium", "fine")   # = generate_onera_m6.RICHARDSON_LADDER

TIP_AFT, TIP_Y = 0.30, 0.03

#: Legacy caches this script may ADOPT rather than re-solve (same tip model,
#: same M/alpha). `fine` is the honest weak point: the two fine runs in the
#: cache agree on the FIELD to three digits (edge 0.569/0.570, tip-cap
#: 1.014/1.015, wing p99 0.629/0.629) but disagree on the CLAMP COUNTERS --
#: probe_..._fine_ffspan reports 306 limited / 138 floored / not-converged,
#: while probe_..._fine_p9 (whose producing script was never committed, hence
#: unreproducible) reports 0/0/converged and is the one the committed G13.2
#: numbers came from. We adopt _p9 so the numbers here are the SAME objects the
#: docs cite, and surface the discrepancy as a check rather than hiding it.
ADOPT = {
    "medium": "probe_vanish_smooth_rc0.05_medium_ffspan.npz",
    "fine": "probe_vanish_smooth_rc0.05_fine_p9.npz",
}
FINE_ALT = "probe_vanish_smooth_rc0.05_fine_ffspan.npz"


def solve(level):
    mc, wc = cut_wake(read_mesh(MESH / f"{level}.msh"))
    taper = tip_taper_factors(wc.station_z, B_SEMI, FORM, R_C_FRAC * B_SEMI)
    r = solve_newton_lifting(mc, wc, m_inf=M, alpha_deg=ALPHA, upwind_c=1.5,
                             precond="amg", tol_residual=1e-9, tip_taper=taper,
                             farfield_spanwise_gamma=True)
    op = PicardOperator(mc.nodes, mc.elements)
    _, q2 = op.velocities(np.asarray(r["phi"]))
    cen = mc.nodes[mc.elements].mean(axis=1)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    gamma = np.asarray(r["gamma"])
    return dict(
        cen=cen, mach=np.sqrt(mach_squared_field(q2, M, 1.4)),
        xte=np.array([x_te(np.clip(z, 0, B_SEMI)) for z in cen[:, 2]]),
        gamma=gamma, station_z=np.asarray(wc.station_z),
        cl_kj=cl_kj_3d(gamma, np.asarray(wc.station_z), s_ref, B_SEMI),
        ntet=len(mc.elements), conv=bool(r["converged"]),
        n_lim=int(r.get("n_limited", 0)), n_flr=int(r.get("n_floored", 0)),
    )


def get(level):
    """Return (data, provenance). Prefer this script's own cache, then an
    adopted legacy cache, then solve."""
    own = OUT / f"g133_{level}.npz"
    if own.exists():
        d = np.load(own)
        return {k: d[k] for k in d.files}, f"cache:{own.name}"
    legacy = OUT / ADOPT.get(level, "")
    if level in ADOPT and legacy.exists():
        d = np.load(legacy)
        return {k: d[k] for k in d.files}, f"adopted:{legacy.name}"
    print(f"  solving {level} (no cache) ...", flush=True)
    t0 = time.time()
    d = solve(level)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(own, **d)
    print(f"    done {time.time() - t0:.0f}s  cl_KJ={float(d['cl_kj']):.4f}  "
          f"conv={d['conv']}", flush=True)
    return d, f"solved:{own.name}"


def _xte(d):
    if "xte" in d:
        return d["xte"]
    return np.array([x_te(np.clip(z, 0, B_SEMI)) for z in d["cen"][:, 2]])


def _boxes(d):
    """The three regions. dx is the chordwise distance PAST the local TE, so
    dx>0 is off-body (wake) and dx<0 is on the wing/tip-cap surface."""
    cen = d["cen"]
    zb, dx, ay = cen[:, 2] / B_SEMI, cen[:, 0] - _xte(d), np.abs(cen[:, 1])
    mach = np.nan_to_num(d["mach"], nan=-1.0)
    return {
        # the object G13.2 fixed: the sheet's free tip edge, strictly OFF-BODY
        "wake_edge": float(mach[(dx > 0.002) & (dx <= TIP_AFT)
                                & (zb > 0.98) & (ay < TIP_Y)].max()),
        # the flat tip cap's sharp convex edge -- ON the body, at the tip
        "tip_cap": float(mach[(zb > 0.98) & (dx < 0.0)].max()),
        # control: the ordinary wing surface, well inboard of the tip
        "wing_p99": float(np.percentile(mach[(zb < 0.90) & (dx < 0.0)], 99)),
    }


def _exponent(peaks, ntets):
    """peak ~ h^-p with h ~ ntet^(-1/3): slope of log(peak) vs log(ntet^(1/3))."""
    x = np.log(np.asarray(ntets, float) ** (1 / 3))
    return float(np.polyfit(x, np.log(np.asarray(peaks, float)), 1)[0])


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("P13/G13.3: what still blocks 3D grid convergence?")

    data, prov = {}, {}
    for lvl in LADDER:
        data[lvl], prov[lvl] = get(lvl)
    boxes = {lvl: _boxes(data[lvl]) for lvl in LADDER}
    ntets = [int(data[lvl]["ntet"]) for lvl in LADDER]

    rows = []
    for region in ("tip_cap", "wing_p99", "wake_edge"):
        peaks = [boxes[lvl][region] for lvl in LADDER]
        p = _exponent(peaks, ntets)
        rows.append((region, *(f"{v:.3f}" for v in peaks), f"{p:+.3f}",
                     "DIVERGES" if p > 0.15 else "bounded"))
        print(f"{region:10s} " + " -> ".join(f"{v:.3f}" for v in peaks)
              + f"   p = {p:+.3f}")

    p_cap = _exponent([boxes[l]["tip_cap"] for l in LADDER], ntets)
    p_wing = _exponent([boxes[l]["wing_p99"] for l in LADDER], ntets)
    p_wake = _exponent([boxes[l]["wake_edge"] for l in LADDER], ntets)
    cls = [float(data[lvl]["cl_kj"]) for lvl in LADDER]

    # ---- the finding ---------------------------------------------------------
    cl.add("G13.3", "the REMAINING singularity is on the WALL (flat tip cap)",
           f"tip-cap peak Mach {boxes['coarse_ss']['tip_cap']:.3f} -> "
           f"{boxes['fine']['tip_cap']:.3f}, p = {p_cap:+.3f}",
           "p > 0.15 (diverges under uniform refinement)", p_cap > 0.15,
           note="meshgen/wing3d.py builds a FLAT tip cap (deliberate "
                "simplification); its sharp convex edge is a potential-flow "
                "edge singularity ON THE BODY. Fix = round the cap (Track M). "
                "NOT a P11 problem: curved ELEMENTS cannot regularize a sharp "
                "geometric EDGE")
    cl.add("G13.3", "the G13.2 wake-edge fix HOLDS under refinement",
           f"wake-edge peak Mach p = {p_wake:+.3f}",
           "p < 0.15 (bounded)", p_wake < 0.15,
           note="the tapered free edge stays flat across the whole ladder -- "
                "the tip fix REVEALED the wall singularity, it did not cause it")
    cl.add("G13.3", "CONTROL: the ordinary wing field is mesh-converged",
           f"wing p99 Mach p = {p_wing:+.3f}",
           "|p| < 0.05 (converged)", abs(p_wing) < 0.05,
           note="rules out 'the whole solution is just under-resolved' -- the "
                "divergence is localized, which is what a singularity looks like")
    cl.add("G13.3", "lift increments still GROW => NOT in the asymptotic range",
           f"cl_KJ {cls[0]:.4f} -> {cls[1]:.4f} -> {cls[2]:.4f} "
           f"({100 * (cls[1] - cls[0]) / cls[0]:+.2f}%, "
           f"{100 * (cls[2] - cls[1]) / cls[1]:+.2f}%)",
           "no Richardson extrapolation may be reported",
           abs(cls[2] - cls[1]) > abs(cls[1] - cls[0]),
           note="P9/G9.3 discipline: report n/a, never fabricate an "
                "extrapolated value from a non-asymptotic sequence")

    # ---- evidence honesty ----------------------------------------------------
    fine_alt = OUT / FINE_ALT
    if fine_alt.exists():
        alt = np.load(fine_alt)
        cl.add("G13.3", "PROVENANCE: the fine point's clamp counters do NOT reproduce",
               f"adopted fine: conv={bool(data['fine']['conv'])}, "
               f"lim={int(data['fine']['n_lim'])}, flr={int(data['fine']['n_flr'])}"
               f" | reproducible fine: conv={bool(alt['conv'])}, "
               f"lim={int(alt['n_lim'])}, flr={int(alt['n_flr'])}",
               "recorded, not asserted", True,
               note="the two fine runs agree on the FIELD to 3 digits but not "
                    "on the limiter/floor counters; the committed '0 lim / 0 flr' "
                    "claim rests on a run whose script was never committed. The "
                    "box-study VERDICT is unaffected (both give the same p)")

    write_csv(OUT, "g133_boxstudy.csv",
              "region," + ",".join(LADDER) + ",exponent_p,verdict", rows)
    write_csv(OUT, "g133_ladder.csv",
              "level,n_tets,cl_KJ,converged,n_limited,n_floored,provenance",
              [(lvl, int(data[lvl]["ntet"]), f"{float(data[lvl]['cl_kj']):.4f}",
                int(bool(data[lvl]["conv"])), int(data[lvl]["n_lim"]),
                int(data[lvl]["n_flr"]), prov[lvl]) for lvl in LADDER])

    # ---- figure --------------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    h = np.asarray(ntets, float) ** (-1 / 3)
    styles = [("tip_cap", CRITICAL, "o", "tip-cap edge (WALL) -- DIVERGES"),
              ("wake_edge", S1_BLUE, "s", "wake free edge (G13.2-fixed)"),
              ("wing_p99", S2_AQUA, "^", "wing p99 (control)")]
    for key, c, mk, lab in styles:
        ax.loglog(1 / h, [boxes[l][key] for l in LADDER], mk + "-",
                  color=c, label=lab, lw=1.8, ms=6)
    ax.set_xlabel("1 / h  (~ n_tets$^{1/3}$)")
    ax.set_ylabel("peak local Mach in box")
    ax.set_title("Only the flat tip cap still diverges", loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, which="both", alpha=0.25)

    ax2.plot(1 / h, cls, "o-", color=S3_YELLOW, lw=1.8, ms=6)
    for x, y in zip(1 / h, cls):
        ax2.annotate(f"{y:.4f}", (x, y), textcoords="offset points",
                     xytext=(0, 7), fontsize=8, color=INK_2, ha="center")
    ax2.set_xscale("log")
    ax2.set_xlabel("1 / h  (~ n_tets$^{1/3}$)")
    ax2.set_ylabel("cl$_{KJ}$")
    ax2.set_title("increments GROW => not asymptotic => no Richardson",
                  loc="left")
    ax2.grid(True, alpha=0.25)
    finish(fig, OUT, "g133_ladder.png")

    return cl.report(OUT, "checks_g133.csv")


if __name__ == "__main__":
    raise SystemExit(main())
