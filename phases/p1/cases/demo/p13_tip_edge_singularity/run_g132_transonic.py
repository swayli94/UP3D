"""P13/G13.2 -- the TRANSONIC clause: is G9.1's fine mesh a real solution now?

★ WHY THIS SCRIPT EXISTS. The conclusion it measures was committed (cb6ba98) with
NO reproducible evidence -- no script, no CSV, no cached solve. The 2026-07-13
audit found the numbers existed only as prose in the .md trackers, even though a
P11 ledger status had been changed on their strength. This script IS that missing
evidence, and it re-measures the claim from scratch rather than restating it.

THE CLAIM UNDER TEST (roadmap P13/G13.2, "transonic" clause):

  P9/G9.1 could not run a three-point Richardson on ONERA M6 at M0.84/alpha3.06
  because the FINE mesh was not a discrete solution at all: the tip/wake-edge
  singularity drove the unlimited local Mach to 7.93 with 9 cells past the
  M_cap = 3 speed limiter, so those cells stayed permanently limited, the N5
  freeze machinery (which requires 0 limited) could never engage, Newton
  limit-cycled, and G9.1 declared its cl_KJ 0.2393 an ARTIFACT, "not a lift".

  With the G13.2 tip taper applied, the same fine solve is claimed to become a
  GENUINE discrete solution (M_max 2.818, 0 cells over M_cap, 0 limited, 0
  floored, converged), with cl_KJ 0.2593 -> 0.2652 -> 0.2866 against the
  Tranair/KRATOS reference 0.288 -- which fires P9's PRE-REGISTERED
  "cl_KJ >= 0.283 => resolution-dominated" branch and refutes P11's lift case.

WHAT IS MEASURED. Exactly G9.1's own census (`cases/demo/p9_grid_discrimination/
run_demo.py::tip_te_singularity_census`), so the two are a strict A/B on the same
meshes and the same recipe -- the ONLY difference is `tip_taper`:

  * the UNLIMITED local Mach field (before the speed limiter hides the blow-up),
  * the census of cells crossing M_cap = 3,
  * n_limited / n_floored / converged as reported by the Newton driver,
  * cl_KJ (Kutta-Joukowski, probe-based) and cl_p (pressure-integrated,
    probe-free) -- both, because G13.3 needed them to exonerate the Kutta probes.

★ READ THE VERDICT, NOT THE HEADLINE. The gate this script asserts is the one
G9.1 actually failed -- "is the fine mesh a discrete solution?" -- because that
is what the taper is claimed to fix and what a lift number is worthless without.
Whether cl_KJ then clears P9's 0.283 band is REPORTED, and is a separate question
from whether the sequence is asymptotic (it is not -- see run_g133_ladder.py: the
flat tip cap still diverges, so no Richardson may be extrapolated from it).

COST. Heavy: the fine level is a full Mach continuation on 2.5M tets (~40 min);
coarse+medium add ~5 min. Solves cache to results/g132_transonic_*.npz
(gitignored, like every other heavy artifact here) -- the CSV is the evidence.

  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
  python cases/demo/p13_tip_edge_singularity/run_g132_transonic.py
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
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.physics.isentropic import limit_q2_field, mach_squared_field  # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve  # noqa: E402
from pyfp3d.post.shock import shock_report  # noqa: E402
from pyfp3d.post.surface import (  # noqa: E402
    cl_kj_3d, planform_area, wall_force_coefficients,
)
from pyfp3d.solve.newton import solve_newton_transonic  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes" / "onera_m6"

M_INF, ALPHA = 0.84, 3.06
M_CAP = 3.0                      # the P4 speed-limiter cap (physics/isentropic)
LEVELS = ("coarse", "medium", "fine")

# The G13.2 shipped tip model.
FORM, R_C_FRAC = "vanish_smooth", 0.05

#: ONERA M6 Newton recipe. Structurally P8/N6's NEWTON_M6_RECIPE
#: (tests/test_p8_newton.py) with two deliberate changes:
#:   (1) `tip_taper` added -- the whole point (controlled A/B vs G9.1).
#:   (2) precond "direct" -> "amg" + TIGHT Eisenstat-Walker forcing
#:       (ew_eta0 = ew_eta_max = 1e-8).
#: ★ WHY (2) IS MANDATORY AT FINE, NOT OPTIONAL. The fine mesh is 428,118
#: nodes (~450k dofs). P9 measured that a single `precond="direct"` splu at
#: ~450k dofs runs 4 h 39 min WITHOUT RETURNING (RSS 26 GB) -- it does not
#: scale -- while `precond="amg"` + tight EW forcing is valid AND faster at
#: EVERY size (M6 medium 66 s vs 141 s direct, same solution to 4 digits, 0
#: GMRES stalls). The P8/N6 medium recipe uses direct because at 63k nodes
#: direct is ~19 s; copying it to fine is exactly the documented trap (and it
#: is the trap this demo's FIRST attempt fell into -- killed at 1 h 16 min /
#: 24 GB before the switch). coarse/medium cache with direct (verified against
#: the recorded 0.2593 / 0.2652); a cold rerun of any level now uses amg and
#: reproduces to 4 digits (P9).
#: ★ THIRD fine-mesh guard (P9): the fine COLD Picard-5 seed overshoots the LE
#: into the density floor at M0.70 (4036 limited / 1847 floored => level-0 stalls
#: and cannot dm-halve => the whole solve breaks). P9's continuation-path fix is
#: m_start = 0.30 + n_picard_seed = 12 (=> clean 0/0 start). Set here for all
#: levels; coarse/medium are cached so unaffected, and these are path-only knobs
#: (they do not move the converged solution).
NEWTON_M6_RECIPE = dict(
    m_start=0.30, dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
    newton_kw=dict(freeze_refresh_max=8, precond="amg", n_picard_seed=12,
                   ew_eta0=1e-8, ew_eta_max=1e-8, n_newton_max=60,
                   farfield_spanwise_gamma=True),
)

#: Inviscid full-potential reference for ONERA M6 at M0.84/alpha3.06.
#: Tranair / KRATOS, as recorded in roadmap.md (P9's pre-registered bands are
#: stated against it). NOT a cases/reference_data file -- cited, not invented.
CL_REF = 0.288
BAND_RESOLUTION = 0.283          # P9 pre-registered: >= this => resolution-dominated
BAND_FLOOR = 0.278               # P9 pre-registered: <= this => geometry floor


def solve(level):
    t0 = time.perf_counter()
    mc, wc = cut_wake(read_mesh(MESH / f"{level}.msh"))
    taper = tip_taper_factors(wc.station_z, B_SEMI, FORM, R_C_FRAC * B_SEMI)

    recipe = dict(NEWTON_M6_RECIPE)
    recipe["newton_kw"] = dict(recipe["newton_kw"], tip_taper=taper)
    r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                               verbose=True, **recipe)
    wall = time.perf_counter() - t0

    phi = np.asarray(r["phi"])
    gamma = np.asarray(r["gamma"])
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])

    # --- G9.1's census, verbatim: the UNLIMITED Mach field and the M_cap count
    op = PicardOperator(mc.nodes, mc.elements)
    _, q2 = op.velocities(phi)
    q2l = limit_q2_field(q2, M_INF, M_CAP, 1.4)
    mach = np.sqrt(mach_squared_field(q2, M_INF, 1.4))       # UNLIMITED
    n_over_cap = int(np.count_nonzero(q2l != q2))

    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
        alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
    shocks = {}
    for eta in (0.44, 0.65, 0.90):
        c = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
        shocks[eta] = float(shock_report(c, M_INF)["upper"]["x_shock"])

    return dict(
        ntet=len(mc.elements),
        cl_kj=float(cl_kj_3d(gamma, np.asarray(wc.station_z), s_ref, B_SEMI)),
        cl_p=float(forces["cl"]),
        m_max=float(np.nanmax(mach)),
        n_over_cap=n_over_cap,
        n_nan=int(np.count_nonzero(~np.isfinite(mach))),
        conv=bool(r["converged"]),
        n_lim=int(r.get("n_limited", 0)),
        n_flr=int(r.get("n_floored", 0)),
        shock_044=shocks[0.44], shock_065=shocks[0.65], shock_090=shocks[0.90],
        wall_s=wall,
    )


def get(level):
    cache = OUT / f"g132_transonic_{level}.npz"
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}
    print(f"  solving {level} at M{M_INF}/alpha{ALPHA} (no cache) ...", flush=True)
    d = solve(level)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **d)
    print(f"    {level}: {d['wall_s']:.0f}s  cl_KJ={d['cl_kj']:.4f}  "
          f"cl_p={d['cl_p']:.4f}  M_max={d['m_max']:.3f}  "
          f"over_cap={d['n_over_cap']}  lim={d['n_lim']}  flr={d['n_flr']}  "
          f"conv={d['conv']}", flush=True)
    return d


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("P13/G13.2 transonic: is G9.1's fine mesh a discrete solution now?")

    D = {lvl: get(lvl) for lvl in LEVELS}
    fine = D["fine"]
    kj = [float(D[l]["cl_kj"]) for l in LEVELS]
    cp = [float(D[l]["cl_p"]) for l in LEVELS]

    # ---- THE GATE: the thing G9.1 could not achieve --------------------------
    cl.add("G13.2", "★ the FINE mesh is a GENUINE DISCRETE SOLUTION",
           f"M_max {float(fine['m_max']):.3f}, {int(fine['n_over_cap'])} cells "
           f"over M_cap={M_CAP:.0f}, {int(fine['n_lim'])} limited, "
           f"{int(fine['n_flr'])} floored, converged={bool(fine['conv'])}",
           "0 over cap AND 0 limited AND 0 floored AND converged",
           int(fine["n_over_cap"]) == 0 and int(fine["n_lim"]) == 0
           and int(fine["n_flr"]) == 0 and bool(fine["conv"]),
           note="G9.1's untapered fine: unlimited M_max 7.93, 9 cells over the "
                "cap, permanently limited => N5 freeze can never engage => "
                "Newton limit-cycles => its cl_KJ 0.2393 was declared an "
                "ARTIFACT. This is the defect the tip taper is claimed to fix")
    cl.add("G13.2", "every level is a discrete solution (not just fine)",
           "; ".join(f"{l}: lim={int(D[l]['n_lim'])} flr={int(D[l]['n_flr'])} "
                     f"conv={bool(D[l]['conv'])}" for l in LEVELS),
           "0 limited / 0 floored / converged at all three levels",
           all(int(D[l]["n_lim"]) == 0 and int(D[l]["n_flr"]) == 0
               and bool(D[l]["conv"]) for l in LEVELS))
    cl.add("G13.2", "no NaN cells anywhere",
           "; ".join(f"{l}: {int(D[l]['n_nan'])}" for l in LEVELS),
           "0 at every level",
           all(int(D[l]["n_nan"]) == 0 for l in LEVELS))

    # ---- REPORTED (not asserted): where the lift lands vs P9's bands ---------
    band = ("resolution-dominated" if kj[-1] >= BAND_RESOLUTION else
            "geometry floor" if kj[-1] <= BAND_FLOOR else "inconclusive")
    cl.add("G13.2", "cl_KJ vs P9's PRE-REGISTERED bands (reported)",
           f"cl_KJ {kj[0]:.4f} -> {kj[1]:.4f} -> {kj[2]:.4f} vs ref {CL_REF} "
           f"=> fine lands in the '{band}' band",
           f">= {BAND_RESOLUTION} resolution / <= {BAND_FLOOR} floor", True,
           note="REPORTED, NOT EXTRAPOLATED: the sequence is not asymptotic "
                "(the flat tip cap still diverges -- run_g133_ladder.py), so "
                "P9/G9.3 discipline forbids a Richardson value here. The taper "
                "itself removes ~1-1.6% of cl, so the untapered equivalent is "
                "HIGHER than the number shown")
    cl.add("G13.2", "HONEST: increments still GROW => no Richardson",
           f"cl_KJ increments {100 * (kj[1] - kj[0]) / kj[0]:+.2f}%, "
           f"{100 * (kj[2] - kj[1]) / kj[1]:+.2f}%",
           "recorded, not extrapolated",
           abs(kj[2] - kj[1]) > abs(kj[1] - kj[0]),
           note="a growing increment under uniform refinement means the "
                "sequence is NOT in the asymptotic range; G13.3 localizes the "
                "cause to the flat tip cap")

    rows = [(l, int(D[l]["ntet"]), f"{float(D[l]['cl_kj']):.4f}",
             f"{float(D[l]['cl_p']):.4f}", f"{float(D[l]['m_max']):.3f}",
             int(D[l]["n_over_cap"]), int(D[l]["n_lim"]), int(D[l]["n_flr"]),
             int(bool(D[l]["conv"])), f"{float(D[l]['shock_044']):.3f}",
             f"{float(D[l]['shock_065']):.3f}", f"{float(D[l]['shock_090']):.3f}",
             f"{float(D[l]['wall_s']):.0f}") for l in LEVELS]
    write_csv(OUT, "g132_transonic.csv",
              "level,n_tets,cl_KJ,cl_p,M_max_unlimited,n_cells_over_Mcap,"
              "n_limited,n_floored,converged,shock_eta044,shock_eta065,"
              "shock_eta090,wall_s", rows)

    # ---- figure --------------------------------------------------------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.3))
    h = np.array([float(D[l]["ntet"]) for l in LEVELS]) ** (1 / 3)

    ax.plot(h, kj, "o-", color=S1_BLUE, lw=1.9, ms=6, label="cl$_{KJ}$ (probe)")
    ax.plot(h, cp, "s--", color=S2_AQUA, lw=1.6, ms=5,
            label="cl$_p$ (pressure, probe-free)")
    ax.axhline(CL_REF, color=INK_2, lw=1.1, ls=":")
    ax.text(h[0], CL_REF, f" Tranair/KRATOS {CL_REF}", fontsize=7.5, color=INK_2,
            va="bottom")
    ax.axhspan(BAND_RESOLUTION, ax.get_ylim()[1], color=S2_AQUA, alpha=0.07)
    ax.text(h[0], BAND_RESOLUTION, " P9 band: resolution-dominated",
            fontsize=7.5, color=S2_AQUA, va="bottom")
    ax.set_xscale("log")
    ax.set_xlabel(r"mesh density $n_{\rm tet}^{1/3}\propto 1/h$")
    ax.set_ylabel("lift coefficient")
    ax.set_title(f"ONERA M6, M{M_INF} / $\\alpha${ALPHA}, tapered tip", loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, alpha=0.25)

    mm = [float(D[l]["m_max"]) for l in LEVELS]
    ax2.plot(h, mm, "o-", color=S2_AQUA, lw=1.9, ms=6, label="tapered (this run)")
    ax2.plot(h, [1.40, 2.13, 7.93], "x--", color=CRITICAL, lw=1.6, ms=7,
             label="untapered (P9/G9.1, recorded)")
    ax2.axhline(M_CAP, color=CRITICAL, lw=1.1, ls=":")
    ax2.text(h[0], M_CAP, " speed-limiter cap $M_{cap}=3$", fontsize=7.5,
             color=CRITICAL, va="bottom")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"mesh density $n_{\rm tet}^{1/3}\propto 1/h$")
    ax2.set_ylabel("unlimited max local Mach")
    ax2.set_title("G9.1's fine mesh crossed the cap; the tapered one does not",
                  loc="left")
    ax2.legend(fontsize=8, frameon=False)
    ax2.grid(True, which="both", alpha=0.25)
    finish(fig, OUT, "g132_transonic.png")

    return cl.report(OUT, "checks_g132_transonic.csv")


if __name__ == "__main__":
    raise SystemExit(main())
