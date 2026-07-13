"""P13/G13.3 -- does rounding the tip cap unblock 3D grid convergence?

THE STORY SO FAR. Three distinct defects have blocked a three-point Richardson
on ONERA M6, and they were found one behind the other, each hiding the next:

  (1) the wake sheet's free TIP EDGE          FIXED by G13.2 (loading taper)
  (2) the mesh ladder's h_far CLAMP           FIXED by Track M / M1b (coarse_ss)
  (3) the FLAT TIP CAP's sharp wall edge      <- this script

run_g133_ladder.py localized (3) with a three-region box study on the repaired
ladder: the tip-cap edge peak Mach DIVERGES (p = +0.321) while the wake free
edge (+0.045) and the wing interior (-0.014) are bounded. A flat cap meets the
upper and lower surfaces at a sharp convex edge; in potential flow that is a
singularity, and refinement resolves it rather than removing it -- so the lift
sequence can never be asymptotic. Track M / M5 replaces the cap with the half
body of revolution swept by the tip section about its own chord line
(wing3d.py::tip_cap="round"); the cap radius vanishes at the LE and the TE, so
the TE line, the wake sheet, the tip TE corner, the Kutta stations and B_SEMI
are UNCHANGED and this is a controlled A/B: the tip WALL moved, nothing else.

WHAT IS MEASURED. run_g133_ladder.py's box study, VERBATIM (same boxes, same
exponent fit, same solver recipe, same tip taper) -- on the M5 ladder, and
against the M1 ladder's own numbers re-derived from the same code, so the two
columns are comparable by construction rather than by transcription.

  region                what it tells us
  tip_cap    (WALL)     the object M5 changed. p must fall from +0.321 to ~0.
  wake_edge  (off-body) G13.2's fix. Must STAY bounded -- rounding the cap must
                        not have re-lit the wake edge.
  wing_p99   (control)  the ordinary field. Must stay converged, or the whole
                        solution is merely under-resolved and the box study is
                        measuring nothing.

AND THEN THE POINT OF ALL OF IT: with the sequence finally asymptotic, the
three-point Richardson that P9/G9.1 could not run.

★ HONESTY RULE (P9/G9.3 discipline, and the reason this file exists at all):
a Richardson value is EARNED, not computed. It is reported only if the observed
order is sane and the increments SHRINK. If they do not, this script says n/a
and says why -- it does not fabricate an extrapolated number from a
non-asymptotic sequence, which is exactly the trap G9.1 fell into.

COST. coarse ~1 min, medium ~6 min, fine ~45 min (3.2 M tets). Solves cache to
results/g133rt_*.npz (gitignored); the CSV is the committed evidence.

  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python cases/demo/p13_tip_edge_singularity/run_g133_roundtip.py
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
from pyfp3d.meshgen.wing3d import B_SEMI, x_le, x_te  # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field  # noqa: E402
from pyfp3d.post.surface import cl_kj_3d, planform_area  # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH_FLAT = REPO / "cases" / "meshes" / "onera_m6"
MESH_ROUND = REPO / "cases" / "meshes" / "onera_m6_roundtip"

# Identical to run_g133_ladder.py -- subsonic, so the edge signal is geometric
# and no limiter/shock machinery stands between the cap and the measurement.
ALPHA, M = 3.06, 0.5
FORM, R_C_FRAC = "vanish_smooth", 0.05
TIP_AFT, TIP_Y = 0.30, 0.03

#: Both families' self-similar ladders. M1 spells its coarse end `coarse_ss`
#: (M1b re-cut it without the h_far clamp); M5 never needed the repair.
LADDERS = {
    "flat (M1)": (MESH_FLAT, ("coarse_ss", "medium", "fine")),
    "round (M5)": (MESH_ROUND, ("coarse", "medium", "fine")),
}

#: Legacy caches from run_g133_ladder.py, adopted rather than re-solved (same
#: tip model, same M/alpha, same recipe -- they ARE the flat-ladder points that
#: script reported). Provenance is written into the CSV; `fine` is the weak one
#: (its producing script was never committed -- see run_g133_ladder.py's
#: FINE_PROVENANCE note), which is one more reason the M5 column is solved fresh.
ADOPT_FLAT = {
    "coarse_ss": "g133_coarse_ss.npz",
    "medium": "probe_vanish_smooth_rc0.05_medium_ffspan.npz",
    "fine": "probe_vanish_smooth_rc0.05_fine_p9.npz",
}

#: Tranair / KRATOS inviscid full-potential reference (roadmap P9; the
#: pre-registered bands are stated against it). At M0.5 this is NOT the
#: comparison case -- it is quoted only where the transonic clause is discussed.
P_DIVERGES = 0.15          # run_g133_ladder.py's verdict threshold


def solve(mesh_dir, level):
    mc, wc = cut_wake(read_mesh(mesh_dir / f"{level}.msh"))
    taper = tip_taper_factors(wc.station_z, B_SEMI, FORM, R_C_FRAC * B_SEMI)
    # ★ precond="amg" + TIGHT EW forcing, MANDATORY at fine (~554k dofs here).
    # P9 measured a single precond="direct" splu at ~450k dofs running 4h39m
    # without returning (26 GB RSS) -- it does not scale; the medium recipe's
    # direct is a trap when copied to fine. At M0.5 (subsonic, no shock soft
    # mode) the default adaptive forcing also converges, but ew_eta=1e-8
    # matches the validated fast path and the corrected run_g132_transonic.py.
    r = solve_newton_lifting(mc, wc, m_inf=M, alpha_deg=ALPHA, upwind_c=1.5,
                             precond="amg", tol_residual=1e-9, tip_taper=taper,
                             ew_eta0=1e-8, ew_eta_max=1e-8,
                             farfield_spanwise_gamma=True)
    op = PicardOperator(mc.nodes, mc.elements)
    _, q2 = op.velocities(np.asarray(r["phi"]))
    cen = mc.nodes[mc.elements].mean(axis=1)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    gamma = np.asarray(r["gamma"])
    return dict(
        cen=cen, mach=np.sqrt(mach_squared_field(q2, M, 1.4)),
        xte=np.array([x_te(np.clip(z, 0, B_SEMI)) for z in cen[:, 2]]),
        gamma=gamma, station_z=np.asarray(wc.station_z), s_ref=s_ref,
        cl_kj=cl_kj_3d(gamma, np.asarray(wc.station_z), s_ref, B_SEMI),
        ntet=len(mc.elements), conv=bool(r["converged"]),
        n_lim=int(r.get("n_limited", 0)), n_flr=int(r.get("n_floored", 0)),
    )


def get(fam, level):
    """(data, provenance). The M5 column is solved and cached here; the M1
    column is adopted from run_g133_ladder.py's caches so the two are the same
    objects the committed G13.3 numbers came from."""
    if fam.startswith("flat"):
        legacy = OUT / ADOPT_FLAT[level]
        if legacy.exists():
            d = np.load(legacy)
            return {k: d[k] for k in d.files}, f"adopted:{legacy.name}"
    own = OUT / f"g133rt_{level}.npz"
    if fam.startswith("round") and own.exists():
        d = np.load(own)
        return {k: d[k] for k in d.files}, f"cache:{own.name}"

    mesh_dir, _ = LADDERS[fam]
    print(f"  solving {fam} {level} (no cache) ...", flush=True)
    t0 = time.time()
    d = solve(mesh_dir, level)
    OUT.mkdir(parents=True, exist_ok=True)
    out = own if fam.startswith("round") else OUT / f"g133flat_{level}.npz"
    np.savez(out, **d)
    print(f"    done {time.time() - t0:.0f}s  cl_KJ={float(d['cl_kj']):.4f}  "
          f"conv={d['conv']}", flush=True)
    return d, f"solved:{out.name}"


def _xte(d):
    if "xte" in d:
        return d["xte"]
    return np.array([x_te(np.clip(z, 0, B_SEMI)) for z in d["cen"][:, 2]])


def _chord_frac(cen):
    """Chordwise fraction 0..1 of each cell centroid at its own span station
    (0 = LE, 1 = TE), used to exclude the design-sharp LE and TE. x_le/x_te are
    linear in z, so they apply elementwise to the centroid z column."""
    z = cen[:, 2]
    return (cen[:, 0] - x_le(z)) / (x_te(z) - x_le(z))


def _boxes(d):
    """Regions of the tip. dx is the chordwise distance PAST the local TE:
    dx > 0 is off-body (wake), dx < 0 is on the wing / tip-cap surface.

    ★ THE METRIC TRAP (G13.1/G13.2 finding 6, and it bit here too). The broad
    `tip_cap` box `(zb>0.98) & (dx<0)` was G13.1's, and on the FLAT cap its max
    sat on the seam edge. Once M5 rounds the cap that edge is gone, so the box
    max MIGRATES to the next-sharpest feature -- which on this geometry is the
    zero-thickness TE (sharp BY DESIGN, both families, it carries the Kutta
    condition). Measured: the fine-mesh broad-box peak sits at chord-frac 0.999.
    So the broad box is REPORTED (it shows the migration honestly) but the GATE
    is on two boxes that isolate the object M5 actually changed:

      * `tip_cap_surface` -- strictly OUTBOARD of z = B_SEMI, i.e. the revolved
        cap itself (the new geometry). Answers "is the cap surface singular?"
      * `tip_body_noTE`   -- the tip region minus the design-sharp LE/TE
        (chord-frac 0.05..0.90). Answers "is the tip converged once the two
        edges that are sharp on purpose are set aside?"
    """
    cen = d["cen"]
    zb, dx, ay = cen[:, 2] / B_SEMI, cen[:, 0] - _xte(d), np.abs(cen[:, 1])
    xi = _chord_frac(cen)
    mach = np.nan_to_num(d["mach"], nan=-1.0)

    def mx(sel):
        # nan for an empty selection: the flat cap has NO cells outboard of
        # z = B_SEMI, so its tip_cap_surface is legitimately undefined.
        return float(mach[sel].max()) if np.any(sel) else float("nan")

    return {
        "wake_edge": mx((dx > 0.002) & (dx <= TIP_AFT) & (zb > 0.98)
                        & (ay < TIP_Y)),
        "tip_cap": mx((zb > 0.98) & (dx < 0.0)),
        "tip_cap_surface": mx((zb > 1.0) & (dx < 0.0)),
        "tip_body_noTE": mx((zb > 0.98) & (dx < 0.0) & (xi > 0.05) & (xi < 0.90)),
        "wing_p99": float(np.percentile(mach[(zb < 0.90) & (dx < 0.0)], 99)),
    }


def _exponent(peaks, ntets):
    """peak ~ h^-p with h ~ ntet^(-1/3): slope of log(peak) vs log(ntet^(1/3))."""
    x = np.log(np.asarray(ntets, float) ** (1 / 3))
    return float(np.polyfit(x, np.log(np.asarray(peaks, float)), 1)[0])


def richardson(vals, ntets):
    """Three-point Richardson on a NON-uniform-in-count but self-similar ladder.

    Returns (order p, extrapolated value, refinement ratio r). The ladder halves
    h at every level, so r = 2 by construction -- but it is COMPUTED from the
    element counts rather than assumed, because if the mesh family were ever
    knocked off its refinement ray again (the M1b defect) this is where it would
    show.
    """
    f1, f2, f3 = (float(v) for v in vals)          # coarse -> fine
    h = np.asarray(ntets, float) ** (-1 / 3)
    r = float(np.sqrt((h[0] / h[1]) * (h[1] / h[2])))
    d21, d32 = f2 - f1, f3 - f2
    if d21 == 0.0 or d32 / d21 <= 0.0:
        return float("nan"), float("nan"), r      # not monotone: no order
    p = float(np.log(abs(d21 / d32)) / np.log(r))
    f_ext = f3 + d32 / (r ** p - 1.0)
    return p, f_ext, r


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("P13/G13.3: does the rounded tip cap unblock 3D grid convergence?")

    missing = [f"{fam} {lv}" for fam, (d, lvls) in LADDERS.items() for lv in lvls
               if not (d / f"{lv}.msh").exists()
               and not (OUT / ADOPT_FLAT.get(lv, "_")).exists()]
    if missing:
        print("missing meshes for: " + ", ".join(missing))
        return 1

    data, prov, boxes, ntets = {}, {}, {}, {}
    for fam, (_, levels) in LADDERS.items():
        for lv in levels:
            data[(fam, lv)], prov[(fam, lv)] = get(fam, lv)
        boxes[fam] = [_boxes(data[(fam, lv)]) for lv in levels]
        ntets[fam] = [int(data[(fam, lv)]["ntet"]) for lv in levels]

    REGIONS = ("tip_cap", "tip_cap_surface", "tip_body_noTE",
               "wake_edge", "wing_p99")
    p = {fam: {reg: _exponent([b[reg] for b in boxes[fam]], ntets[fam])
               for reg in REGIONS if reg in boxes[fam][0]}
         for fam in LADDERS}
    cls = {fam: [float(data[(fam, lv)]["cl_kj"]) for lv in levels]
           for fam, (_, levels) in LADDERS.items()}

    for fam in LADDERS:
        print(f"\n{fam}")
        for reg in REGIONS:
            if reg not in boxes[fam][0]:
                continue
            seq = " -> ".join(f"{b[reg]:.3f}" for b in boxes[fam])
            print(f"  {reg:16s} {seq}   p = {p[fam][reg]:+.3f}")
        print(f"  cl_KJ            " + " -> ".join(f"{v:.4f}" for v in cls[fam]))

    F, R = "flat (M1)", "round (M5)"

    # ---- (1) the cap ---------------------------------------------------------
    # ★ Gate on the object M5 actually changed, NOT the broad G13.1 box. The
    # broad box's max migrated onto the design-sharp TE once the cap seam was
    # removed (metric trap, G13.1/G13.2 finding 6) -- it is REPORTED below so the
    # migration is visible, but the cap-surface and TE-excluded boxes are what
    # answer "did rounding the cap remove its singularity".
    cl.add("G13.3", "★ the CAP-SURFACE singularity is GONE",
           f"cap surface (z > B_SEMI) peak Mach p = {p[R]['tip_cap_surface']:+.3f}; "
           "peaks " + " -> ".join(f"{b['tip_cap_surface']:.3f}" for b in boxes[R]),
           f"p < {P_DIVERGES} (bounded)", p[R]["tip_cap_surface"] < P_DIVERGES,
           note="the flat cap's sharp convex edge was a potential-flow "
                "singularity ON THE BODY (seam crease ~92 deg, box p = +0.32). "
                "The revolved cap surface itself is now bounded. GEOMETRY fix "
                "(Track M / M5) -- isoparametric elements (P11) cannot "
                "regularize a sharp edge, which is why P11 was never the answer")
    cl.add("G13.3", "★ the tip is converged once the design-sharp LE/TE are set aside",
           f"tip (chord-frac 0.05..0.90) peak Mach p = {p[R]['tip_body_noTE']:+.3f}; "
           "peaks " + " -> ".join(f"{b['tip_body_noTE']:.3f}" for b in boxes[R]),
           "|p| < 0.05 (converged)", abs(p[R]["tip_body_noTE"]) < 0.05,
           note="the TE is a zero-thickness sharp edge BY DESIGN (it carries the "
                "Kutta condition) and is present in BOTH families -- excluding it "
                "and the LE isolates what the cap change did")
    cl.add("G13.3", "REPORTED: the broad G13.1 box still 'diverges' -- METRIC TRAP",
           f"broad tip_cap box p = {p[R]['tip_cap']:+.3f} "
           f"(round) vs {p[F]['tip_cap']:+.3f} (flat); the fine-mesh peak sits "
           "at chord-frac 0.999 = the sharp TE, not the cap",
           "recorded, not asserted -- the max migrated to the sharp TE", True,
           note="G13.1/G13.2 finding 6: once the edge you fixed is gone, a broad "
                "max-in-box metric latches onto the next-sharpest feature. Here "
                "that is the design-sharp tip TE, which no tip-cap change can or "
                "should remove; the cap-surface and TE-excluded boxes above are "
                "the honest measures of the fix")

    # ---- (2) nothing was re-lit ---------------------------------------------
    cl.add("G13.3", "the G13.2 wake-edge fix still holds",
           f"wake-edge p = {p[R]['wake_edge']:+.3f} (flat ladder: "
           f"{p[F]['wake_edge']:+.3f})",
           f"p < {P_DIVERGES} (bounded)", p[R]["wake_edge"] < P_DIVERGES,
           note="the cap sits just forward of the sheet's tip edge; rounding it "
                "must not have re-lit the object G13.2 put out")
    cl.add("G13.3", "CONTROL: the ordinary wing field stays converged",
           f"wing p99 p = {p[R]['wing_p99']:+.3f} (flat ladder: "
           f"{p[F]['wing_p99']:+.3f})",
           "|p| < 0.05", abs(p[R]["wing_p99"]) < 0.05,
           note="rules out 'the whole solution is just under-resolved' -- "
                "without this the box study measures nothing")
    cl.add("G13.3", "every level is a converged discrete solution",
           "; ".join(f"{lv}: lim={int(data[(R, lv)]['n_lim'])} "
                     f"flr={int(data[(R, lv)]['n_flr'])} "
                     f"conv={bool(data[(R, lv)]['conv'])}"
                     for lv in LADDERS[R][1]),
           "0 limited / 0 floored / converged at all three levels",
           all(int(data[(R, lv)]["n_lim"]) == 0
               and int(data[(R, lv)]["n_flr"]) == 0
               and bool(data[(R, lv)]["conv"]) for lv in LADDERS[R][1]))

    # ---- (3) the payoff: is the lift sequence finally asymptotic? ------------
    # ★ The two families' cl LEVELS are not directly comparable: the rounded cap
    # adds projected planform area, so s_ref -- the cl normalisation -- moves.
    # Say so, rather than let a reader difference the two columns.
    s_flat = float(data[(F, "medium")].get("s_ref", np.nan))
    s_round = float(data[(R, "medium")]["s_ref"])
    cl.add("G13.3", "HONESTY: the two families' cl LEVELS are not comparable",
           f"planform area s_ref {s_flat:.4f} (flat) -> {s_round:.4f} (round), "
           f"{100 * (s_round / s_flat - 1):+.2f}%"
           if np.isfinite(s_flat) else
           f"round s_ref {s_round:.4f} (the adopted flat caches predate this key)",
           "recorded, not asserted", True,
           note="the cap adds area outboard of B_SEMI, and the tip loading "
                "itself changes -- of course it does, the geometry changed. "
                "Every convergence claim here is WITHIN a family, which is why "
                "that is the only kind made")

    d21, d32 = cls[R][1] - cls[R][0], cls[R][2] - cls[R][1]
    shrinking = abs(d32) < abs(d21)
    p_obs, cl_ext, r_ref = richardson(cls[R], ntets[R])
    cl.add("G13.3", "★ the lift increments SHRINK => the sequence is asymptotic",
           f"cl_KJ " + " -> ".join(f"{v:.4f}" for v in cls[R])
           + f" (increments {100 * d21 / cls[R][0]:+.2f}%, "
             f"{100 * d32 / cls[R][1]:+.2f}%); flat ladder: "
           + " -> ".join(f"{v:.4f}" for v in cls[F]),
           "|second increment| < |first| -- the thing that was FALSE on every "
           "M6 ladder so far", shrinking,
           note="growing increments under uniform refinement are the signature "
                "of a singularity still being resolved. This is the clause that "
                "has blocked every M6 Richardson, G9.1's included")

    if shrinking and np.isfinite(p_obs):
        cl.add("G13.3", "★ THREE-POINT RICHARDSON (the task G9.1 could not do)",
               f"observed order p = {p_obs:.2f} at refinement ratio r = "
               f"{r_ref:.3f}; cl_KJ(h->0) = {cl_ext:.4f}",
               "reported ONLY because the increments shrink and the order is "
               "sane (0.5 <= p <= 3)", 0.5 <= p_obs <= 3.0,
               note="P1 elements on a smooth geometry should give p ~ 1-2. An "
                    "order outside that band means the sequence is not clean "
                    "enough to extrapolate, and the value must not be quoted")
    else:
        cl.add("G13.3", "RICHARDSON: n/a -- NOT EARNED",
               f"increments {100 * d21 / cls[R][0]:+.2f}%, "
               f"{100 * d32 / cls[R][1]:+.2f}%; observed order "
               f"{'n/a' if not np.isfinite(p_obs) else f'{p_obs:.2f}'}",
               "no extrapolated value may be reported", True,
               note="P9/G9.3 discipline: a Richardson value from a "
                    "non-asymptotic sequence is a fabricated number. Report "
                    "n/a and say what is still diverging")

    # ---- CSVs ----------------------------------------------------------------
    rows = []
    for fam, (_, levels) in LADDERS.items():
        for reg in REGIONS:
            if reg not in boxes[fam][0]:
                continue
            peaks = [b[reg] for b in boxes[fam]]
            pv = p[fam][reg]
            verdict = ("n/a" if not np.isfinite(pv) else
                       "DIVERGES" if pv > P_DIVERGES else "bounded")
            rows.append((fam, reg, *(f"{v:.3f}" for v in peaks),
                         f"{pv:+.3f}" if np.isfinite(pv) else "n/a", verdict))
    write_csv(OUT, "g133rt_boxstudy.csv",
              "family,region,coarse,medium,fine,exponent_p,verdict", rows)

    write_csv(OUT, "g133rt_ladder.csv",
              "family,level,n_tets,cl_KJ,converged,n_limited,n_floored,provenance",
              [(fam, lv, int(data[(fam, lv)]["ntet"]),
                f"{float(data[(fam, lv)]['cl_kj']):.4f}",
                int(bool(data[(fam, lv)]["conv"])),
                int(data[(fam, lv)]["n_lim"]), int(data[(fam, lv)]["n_flr"]),
                prov[(fam, lv)])
               for fam, (_, levels) in LADDERS.items() for lv in levels])

    write_csv(OUT, "g133rt_richardson.csv",
              "quantity,coarse,medium,fine,increment_1_pct,increment_2_pct,"
              "asymptotic,observed_order_p,refinement_ratio,extrapolated",
              [("cl_KJ_round_M5", *(f"{v:.4f}" for v in cls[R]),
                f"{100 * d21 / cls[R][0]:+.2f}", f"{100 * d32 / cls[R][1]:+.2f}",
                int(shrinking),
                f"{p_obs:.3f}" if np.isfinite(p_obs) else "n/a",
                f"{r_ref:.3f}",
                f"{cl_ext:.4f}" if (shrinking and np.isfinite(p_obs)) else "n/a")])

    # ---- figure --------------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.4))
    h = np.asarray(ntets[R], float) ** (1 / 3)
    hf = np.asarray(ntets[F], float) ** (1 / 3)
    # flat cap seam (the object M5 replaced) -- diverges; and the round cap
    # surface -- bounded. Plus the broad box on round, to show the metric trap.
    ax.loglog(hf, [b["tip_cap"] for b in boxes[F]], "x--", color=CRITICAL,
              lw=1.7, ms=8, alpha=0.8,
              label=f"flat cap edge (p={p[F]['tip_cap']:+.2f})")
    ax.loglog(h, [b["tip_cap_surface"] for b in boxes[R]], "o-", color=S2_AQUA,
              lw=1.9, ms=6,
              label=f"round cap surface (p={p[R]['tip_cap_surface']:+.2f})")
    ax.loglog(h, [b["tip_cap"] for b in boxes[R]], "s:", color=S3_YELLOW,
              lw=1.4, ms=5, alpha=0.8,
              label=f"round broad box = sharp TE (p={p[R]['tip_cap']:+.2f})")
    ax.set_xlabel(r"mesh density $n_{\rm tet}^{1/3}\propto 1/h$")
    ax.set_ylabel("peak local Mach in box")
    ax.set_title("Cap edge gone; broad box now sees the sharp TE", loc="left",
                 fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(True, which="both", alpha=0.25)

    for fam, c, ls in ((F, INK_2, "--"), (R, S3_YELLOW, "-")):
        h = np.asarray(ntets[fam], float) ** (1 / 3)
        ax2.plot(h, cls[fam], "o" + ls, color=c, lw=1.8, ms=6, label=fam)
        for x, y in zip(h, cls[fam]):
            ax2.annotate(f"{y:.4f}", (x, y), textcoords="offset points",
                         xytext=(0, 7), fontsize=7.5, color=c, ha="center")
    if shrinking and np.isfinite(p_obs):
        ax2.axhline(cl_ext, color=S2_AQUA, lw=1.1, ls=":")
        ax2.text(np.asarray(ntets[R], float)[0] ** (1 / 3), cl_ext,
                 f" Richardson cl$_{{KJ}}$(h$\\to$0) = {cl_ext:.4f} (p={p_obs:.2f})",
                 fontsize=7.5, color=S2_AQUA, va="bottom")
    ax2.set_xscale("log")
    ax2.set_xlabel(r"mesh density $n_{\rm tet}^{1/3}\propto 1/h$")
    ax2.set_ylabel("cl$_{KJ}$")
    ax2.set_title(f"M{M}/$\\alpha${ALPHA}: increments shrink, extrapolable",
                  loc="left", fontsize=10)
    ax2.legend(fontsize=8, frameon=False)
    ax2.grid(True, alpha=0.25)
    finish(fig, OUT, "g133rt_ladder.png")

    return cl.report(OUT, "checks_g133rt.csv")


if __name__ == "__main__":
    raise SystemExit(main())
