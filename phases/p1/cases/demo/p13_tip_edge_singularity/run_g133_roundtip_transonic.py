"""P13/G13.3 (TRANSONIC) -- the round cap CANNOT complete the M0.84 Mach ramp.

★ RESULT (negative, and sharper than it first looked). Rounding the tip cap
fixes 3D grid convergence SUBSONICALLY (run_g133_roundtip.py: Richardson
p = 2.31). At M0.84 it does the opposite: coarse and medium reach the target and
are clean discrete solutions, but the FINE mesh's Mach continuation BREAKS DOWN
AT M = 0.75 -- it dm-halves until dm < dm_min and gives up, never reaching M0.80,
let alone M0.84. There is no M0.84 fine state on the round mesh at all.

★ THE SITE IS THE SHARP TIP TE -- AND THE ROUND CAP DOES NOT CREATE IT, IT HEATS
IT. The committed FLAT-cap run (run_g132_transonic.py, same taper) completes the
ramp on the SAME refinement level and converges (M_max 2.818, 0 over cap, cl_KJ
0.2866), which at first seems to exonerate the shared trailing edge. It does not:
run_g133_roundtip_transonic_locate.py puts **20 of the 20 fastest cells ON the
sharp tip TE** (z/b 0.97-0.99, chord-frac ~0.998) and **none on the new cap
surface** (z/b > 1). The rounded tip lets flow WRAP AROUND the tip and
accelerate, so the velocity at that same, pre-existing sharp TE is higher at
every level (M_max 1.51 / 2.00 round vs 1.39 / 1.73 flat) -- and at fine that
pushes it past the limiter, where the flat cap's cooler tip flow stays under.
⇒ the cap did not add a singularity; it AMPLIFIED the one that was already there.

⚠ METHOD BUG THIS DEMO ONCE HAD (fixed; see the guard in solve()).
`solve_newton_transonic` returns the FAILED level's state with converged=False --
at THAT level's Mach, not at m_inf. An earlier version of this script censused
that state at M_INF = 0.84 anyway, i.e. applied the WRONG freestream Mach to a
M ~ 0.75 velocity field, and so reported a spurious "M_max 3.97, 5 cells over
M_cap, cl_KJ 0.2415 at M0.84". Those numbers were retracted. solve() now records
the level sequence and REFUSES to census a state whose ramp did not reach the
target; every census column reads "n/a" for such a level.

P9's pre-registered bands therefore CANNOT be fired on the round ladder (only two
of three points exist). They cannot be fired on the flat ladder either: those
solves converge, but their sequence 0.2593 -> 0.2652 -> 0.2866 RISES (not
asymptotic, and on the clamped ladder), so 0.2866 is a single reported point, not
a Richardson extrapolation. ⇒ the "0.019 gap = resolution ⇒ P11 refuted"
conclusion has no clean asymptotic discrete-solution basis on EITHER geometry.

--- original intent, kept for context ---

The subsonic leg (run_g133_roundtip.py) proved the MECHANISM: rounding the tip
cap (Track M / M5) removes the last wall singularity and the three-point
Richardson G9.1 could never run becomes definable (M0.5: p = 2.31,
cl_KJ(h->0) = 0.2050). This script does the leg that fires the VERDICT: the same
M0.84/alpha3.06 transonic three-point Richardson on the rounded ladder, read
against P9's PRE-REGISTERED decision bands.

  cl_KJ(h->0) >= 0.283  => resolution-dominated (the gap to Tranair 0.288 is mesh)
  cl_KJ(h->0) <= 0.278  => a geometry floor is confirmed
  in between            => inconclusive

WHY IT HAS TO BE THE ROUND CAP. The audit-flagged M0.84 claim (cl_KJ 0.2866 =>
"the 0.019 gap is resolution" => "P11's lift case is refuted") was (a) prose-only
/ never reproduced, AND (b) computed on the FLAT cap -- whose sharp tip-cap edge
G13.3 showed diverges under refinement, so its fine mesh was never a clean
discrete solution to extrapolate from. The rounded cap is BOTH reproducible here
AND closer to the real ONERA M6 geometry the 0.288 reference was computed for
(the real wing has a rounded tip; the flat cap was the approximation).

METHOD. run_g132_transonic.py's census + recipe VERBATIM (P8/N6's M6 recipe with
tip taper, precond="amg" + tight EW forcing ew_eta=1e-8 -- NEVER direct, which
does not scale to fine, P9), only the mesh family and the cache names change.
Each level reports whether it is a GENUINE discrete solution (unlimited M_max
bounded, 0 cells over M_cap, 0 limited, 0 floored, converged) before any lift
number is trusted -- the G9.1 discipline.

★ s_ref CAVEAT. The rounded cap adds ~0.9% planform area, so the round-cap cl is
normalised on a slightly larger s_ref than the flat-cap 0.288 reference. Both the
round-normalised cl_KJ and a flat-s_ref-renormalised value are reported, and the
band is read against the flat-s_ref one so the comparison to 0.288 is apples to
apples.

COST. Heavy: fine is a full Mach continuation on 3.26 M tets (~45-60 min);
coarse+medium ~10 min. Caches to results/g133rt_transonic_*.npz (gitignored);
the CSV is the evidence.

  NUMBA_NUM_THREADS=12 OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12 \
  python cases/demo/p13_tip_edge_singularity/run_g133_roundtip_transonic.py
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
MESH = REPO / "cases" / "meshes" / "onera_m6_roundtip"     # the M5 ROUND family
FLAT_MESH = REPO / "cases" / "meshes" / "onera_m6"         # for the s_ref A/B

M_INF, ALPHA = 0.84, 3.06
M_CAP = 3.0
LEVELS = ("coarse", "medium", "fine")
FORM, R_C_FRAC = "vanish_smooth", 0.05

CL_REF = 0.288                    # Tranair / KRATOS (roadmap P9)
BAND_RESOLUTION = 0.283          # P9 pre-registered: >= this => resolution-dominated
BAND_FLOOR = 0.278               # P9 pre-registered: <= this => geometry floor

#: run_g132_transonic.py's corrected recipe, VERBATIM (P8/N6 M6 recipe + the
#: three fine-mesh fixes P9 established): precond="amg" NOT direct (direct splu
#: at ~450k dofs runs 4h39m/26GB, P9), tight EW forcing ew_eta=1e-8 (the loose
#: default stalls the shock soft mode), and m_start=0.30 + n_picard_seed=12 (the
#: cold Picard-5 seed overshoots the LE into the density floor at fine).
NEWTON_M6_RECIPE = dict(
    m_start=0.30, dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5,
    newton_kw=dict(freeze_refresh_max=8, precond="amg", n_picard_seed=12,
                   ew_eta0=1e-8, ew_eta_max=1e-8, n_newton_max=60,
                   farfield_spanwise_gamma=True),
)

#: Flat-cap s_ref (from the M1 family), used to renormalise the round-cap cl so
#: the band read against the 0.288 flat reference is apples to apples. Filled at
#: runtime from the flat mesh if present; falls back to the round s_ref.
_FLAT_SREF = {}


def flat_sref(level):
    if level in _FLAT_SREF:
        return _FLAT_SREF[level]
    # coarse of the round ladder pairs with the flat coarse_ss (self-similar).
    flat_level = "coarse_ss" if level == "coarse" else level
    p = FLAT_MESH / f"{flat_level}.msh"
    if p.exists():
        m = read_mesh(p)
        _FLAT_SREF[level] = float(planform_area(m.nodes, m.boundary_faces["wall"]))
    else:
        _FLAT_SREF[level] = None
    return _FLAT_SREF[level]


def solve(level):
    t0 = time.perf_counter()
    mc, wc = cut_wake(read_mesh(MESH / f"{level}.msh"))
    taper = tip_taper_factors(wc.station_z, B_SEMI, FORM, R_C_FRAC * B_SEMI)
    recipe = dict(NEWTON_M6_RECIPE)
    recipe["newton_kw"] = dict(recipe["newton_kw"], tip_taper=taper)
    r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                               verbose=True, **recipe)
    wall = time.perf_counter() - t0

    # ★★ DID THE MACH RAMP ACTUALLY REACH THE TARGET? This is not a formality.
    # `solve_newton_transonic` dm-halves a level that will not converge and,
    # once dm < dm_min, BREAKS and returns the FAILED level's state with
    # converged=False -- at whatever Mach that level was, NOT at m_inf. Reading
    # a census at M_INF off such a state applies the WRONG freestream Mach to
    # the velocity field and manufactures nonsense (it produced a spurious
    # "M_max 3.97 / 5 cells over cap at M0.84" here before this guard existed:
    # the round fine never got past M = 0.75). So: record the level sequence,
    # and report NO census unless the target was reached.
    lr = r.get("level_results", [])
    m_levels = [float(d["m"]) for d in lr]
    m_final = m_levels[-1] if m_levels else float("nan")
    m_last_ok = max([float(d["m"]) for d in lr if d["converged"]],
                    default=float("nan"))
    reached = bool(r["converged"]) and abs(m_final - M_INF) < 1e-9

    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    nan = float("nan")
    out = dict(
        ntet=len(mc.elements), s_ref=s_ref, wall_s=wall,
        m_target=M_INF, m_final=m_final, m_last_converged=m_last_ok,
        n_levels=len(m_levels), target_reached=reached,
        conv=bool(r["converged"]),
        n_lim=int(r.get("n_limited", 0)), n_flr=int(r.get("n_floored", 0)),
        cl_kj=nan, cl_kj_flatnorm=nan, cl_p=nan, m_max=nan, n_over_cap=-1,
        n_nan=-1, shock_044=nan, shock_065=nan, shock_090=nan,
    )
    if not reached:
        print(f"    ★ {level}: ramp did NOT reach M{M_INF} -- broke at "
              f"M={m_final:.4f} (last converged M={m_last_ok:.4f}). "
              f"No census is reported: there is no M{M_INF} state to census.",
              flush=True)
        return out

    # --- target reached: the state IS an M_INF solution, census it ------------
    phi = np.asarray(r["phi"])
    gamma = np.asarray(r["gamma"])
    st_z = np.asarray(wc.station_z)

    op = PicardOperator(mc.nodes, mc.elements)
    _, q2 = op.velocities(phi)
    q2l = limit_q2_field(q2, M_INF, M_CAP, 1.4)
    mach = np.sqrt(mach_squared_field(q2, M_INF, 1.4))       # UNLIMITED
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
        alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
    shocks = {}
    for eta in (0.44, 0.65, 0.90):
        c = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
        shocks[eta] = float(shock_report(c, M_INF)["upper"]["x_shock"])

    cl_kj = float(cl_kj_3d(gamma, st_z, s_ref, B_SEMI))
    fs = flat_sref(level)   # renormalise onto flat s_ref for the 0.288 band
    out.update(
        cl_kj=cl_kj, cl_kj_flatnorm=cl_kj * (s_ref / fs) if fs else nan,
        cl_p=float(forces["cl"]), m_max=float(np.nanmax(mach)),
        n_over_cap=int(np.count_nonzero(q2l != q2)),
        n_nan=int(np.count_nonzero(~np.isfinite(mach))),
        shock_044=shocks[0.44], shock_065=shocks[0.65], shock_090=shocks[0.90],
    )
    return out


def get(level):
    cache = OUT / f"g133rt_transonic_{level}.npz"
    if cache.exists():
        d = np.load(cache)
        d = {k: d[k] for k in d.files}
        if "target_reached" not in d:
            raise RuntimeError(
                f"{cache.name} predates the ramp-reached guard and its census "
                "may have been computed at the WRONG freestream Mach (see "
                "solve()). Delete it and re-solve.")
        return d
    print(f"  solving {level} at M{M_INF}/alpha{ALPHA} (no cache) ...", flush=True)
    d = solve(level)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **d)
    if bool(d["target_reached"]):
        print(f"    {level}: {d['wall_s']:.0f}s  reached M{M_INF}  "
              f"cl_KJ={d['cl_kj']:.4f} (flatnorm {d['cl_kj_flatnorm']:.4f})  "
              f"M_max={d['m_max']:.3f}  over_cap={d['n_over_cap']}  "
              f"lim={d['n_lim']}  flr={d['n_flr']}", flush=True)
    else:
        print(f"    {level}: {d['wall_s']:.0f}s  ★ RAMP FAILED -- broke at "
              f"M={float(d['m_final']):.4f}, last converged "
              f"M={float(d['m_last_converged']):.4f}; no M{M_INF} state exists",
              flush=True)
    return d


def richardson(vals, ntets):
    f1, f2, f3 = (float(v) for v in vals)
    h = np.asarray(ntets, float) ** (-1 / 3)
    r = float(np.sqrt((h[0] / h[1]) * (h[1] / h[2])))
    d21, d32 = f2 - f1, f3 - f2
    if d21 == 0.0 or d32 / d21 <= 0.0:
        return float("nan"), float("nan"), r
    p = float(np.log(abs(d21 / d32)) / np.log(r))
    return p, f3 + d32 / (r ** p - 1.0), r


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("P13/G13.3 transonic: fire P9's bands on the rounded ladder")

    missing = [lv for lv in LEVELS if not (MESH / f"{lv}.msh").exists()]
    if missing:
        print("missing round meshes: " + ", ".join(missing) + "\n  "
              "python cases/meshes/onera_m6_roundtip/generate_onera_m6_roundtip.py"
              " --all")
        return 1

    D = {lv: get(lv) for lv in LEVELS}
    ntets = [int(D[l]["ntet"]) for l in LEVELS]
    reached = {l: bool(D[l]["target_reached"]) for l in LEVELS}
    fine = D["fine"]

    # ---- (1) did the Mach ramp even reach M0.84? ----------------------------
    cl.add("G13.3", "coarse + medium REACH M0.84 and are clean discrete solutions",
           "; ".join(f"{l}: reached M{float(D[l]['m_final']):.2f}, M_max "
                     f"{float(D[l]['m_max']):.2f}, {int(D[l]['n_over_cap'])} over "
                     f"cap, conv {int(bool(D[l]['conv']))}"
                     for l in ("coarse", "medium")),
           f"ramp reaches M{M_INF}; 0 over cap / 0 lim / 0 flr / converged",
           all(reached[l] and int(D[l]["n_over_cap"]) == 0
               and int(D[l]["n_lim"]) == 0 and int(D[l]["n_flr"]) == 0
               and bool(D[l]["conv"]) for l in ("coarse", "medium")),
           note="the round cap is transonically fine at these two levels "
                "(M_max 1.51 / 2.00, well under the cap)")

    cl.add("G13.3", "★★ NEGATIVE RESULT: the FINE mesh cannot REACH M0.84 at all",
           f"the Mach ramp BREAKS at M={float(fine['m_final']):.4f} "
           f"(last converged M={float(fine['m_last_converged']):.4f}, "
           f"{int(fine['n_levels'])} levels attempted) — it never gets to M0.80, "
           f"let alone M{M_INF}",
           f"target M{M_INF} NOT reached on the fine round mesh",
           not reached["fine"],
           note="★ THIS IS THE FINDING, and it is sharper than 'the fine solve "
                "did not converge'. There is NO M0.84 fine state on the round "
                "mesh at all: the continuation dm-halves at M=0.75 until dm < "
                "dm_min and gives up (1 cell in the density floor, |R| stalling "
                "~8e-6). Reading any census at M_INF off that state applies the "
                "WRONG freestream Mach — an earlier version of this demo did "
                "exactly that and manufactured a spurious 'M_max 3.97 / 5 over "
                "cap at M0.84'. The guard in solve() now refuses to census a "
                "state whose ramp did not reach the target")

    cl.add("G13.3", "★ SITE = the sharp tip TE; the round cap does not CREATE it, "
           "it HEATS it",
           "the 20 fastest cells on the failed fine field: 20/20 on the SHARP "
           "TIP TE (z/b 0.97-0.99, chord-frac ~0.998), 0/20 on the rounded cap "
           "(z/b > 1)  |  and the round M_max is higher at EVERY completed level: "
           f"{float(D['coarse']['m_max']):.2f} / {float(D['medium']['m_max']):.2f} "
           "(round) vs 1.39 / 1.73 (flat)",
           "the singular site is the TE, shared by both families; what the round "
           "cap changes is the SPEED there", not reached["fine"],
           note="★ THE RECONCILIATION. The flat fine DOES complete the ramp "
                "(M_max 2.818, 0 over cap, cl_KJ 0.2866) on the SAME refinement "
                "level, so at first the shared tip TE looks exonerated. It is "
                "not: `run_g133_roundtip_transonic_locate.py` puts every one of "
                "the 20 fastest cells ON the tip TE and NONE on the new cap "
                "surface. The rounded tip lets flow WRAP AROUND the tip and "
                "accelerate — which raises the velocity at that same, "
                "pre-existing sharp TE at every level — and at fine that pushes "
                "it past the limiter, where the flat cap's cooler tip flow stays "
                "under. So the cap did not add a singularity; it amplified the "
                "one that was already there. (Positions are valid on a failed-"
                "level field because argsort(Mach) = argsort(|q|) for any "
                "m_inf; the Mach VALUES on that state are not physical.)")

    # ---- (2) no Richardson, no band --------------------------------------
    cl.add("G13.3", "P9 BAND VERDICT: NOT EARNED — there is no third point",
           f"M0.84 lift exists only at coarse ({float(D['coarse']['cl_kj_flatnorm']):.4f}) "
           f"and medium ({float(D['medium']['cl_kj_flatnorm']):.4f}) "
           f"(flat-norm, {100 * (float(D['medium']['cl_kj_flatnorm']) / float(D['coarse']['cl_kj_flatnorm']) - 1):+.2f}%); "
           "the fine mesh has no M0.84 state at all",
           "two points are not a Richardson; no extrapolation, no band read",
           not reached["fine"],
           note="P9/G9.3 discipline. Coarse and medium agree to 0.2%, which is "
                "suggestive but is NOT a three-point convergence study. ⇒ the P9 "
                "band verdict is not earned on the round cap. On the FLAT cap it "
                "is not earned either: those solves DO converge but the sequence "
                "0.2593 → 0.2652 → 0.2866 RISES (non-asymptotic, and on the "
                "clamped ladder), so 0.2866 is a single reported point, not a "
                "Richardson. ⇒ the 'the 0.019 gap is resolution ⇒ P11 refuted' "
                "conclusion has no clean asymptotic basis on EITHER geometry")

    cl.add("G13.3", "HONESTY: s_ref differs between the families",
           "; ".join(f"{l}: round {float(D[l]['s_ref']):.4f} vs flat "
                     f"{flat_sref(l) if flat_sref(l) else float('nan'):.4f}"
                     for l in LEVELS),
           "recorded -- the band is read on the flat-renormalised cl", True,
           note="the rounded cap adds ~0.9% planform area; cl_kj is round-norm, "
                "cl_kj_flatnorm renormalises onto the flat s_ref for the 0.288 "
                "comparison")

    def fmt(v, n=4):
        v = float(v)
        return "n/a" if not np.isfinite(v) else f"{v:.{n}f}"

    write_csv(OUT, "g133rt_transonic.csv",
              "level,n_tets,m_target,m_final_reached,m_last_converged,"
              "target_reached,s_ref,cl_KJ_roundnorm,cl_KJ_flatnorm,cl_p,"
              "M_max_unlimited,n_over_Mcap,n_limited,n_floored,converged,"
              "shock_044,shock_065,shock_090,wall_s",
              [(l, int(D[l]["ntet"]), f"{float(D[l]['m_target']):.2f}",
                f"{float(D[l]['m_final']):.4f}",
                f"{float(D[l]['m_last_converged']):.4f}",
                int(bool(D[l]["target_reached"])), f"{float(D[l]['s_ref']):.5f}",
                fmt(D[l]["cl_kj"]), fmt(D[l]["cl_kj_flatnorm"]), fmt(D[l]["cl_p"]),
                fmt(D[l]["m_max"], 3),
                int(D[l]["n_over_cap"]) if int(D[l]["n_over_cap"]) >= 0 else "n/a",
                int(D[l]["n_lim"]), int(D[l]["n_flr"]), int(bool(D[l]["conv"])),
                fmt(D[l]["shock_044"], 3), fmt(D[l]["shock_065"], 3),
                fmt(D[l]["shock_090"], 3), f"{float(D[l]['wall_s']):.0f}")
               for l in LEVELS])
    # NOTE: every census column is "n/a" wherever target_reached = 0. That is the
    # point: there is no M0.84 state on that mesh, so there is nothing to census.

    # ---- figure -------------------------------------------------------------
    ok = [l for l in LEVELS if reached[l]]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.4))
    h_ok = np.array([D[l]["ntet"] for l in ok], float) ** (1 / 3)
    h_all = np.asarray(ntets, float) ** (1 / 3)

    ax.plot(h_ok, [float(D[l]["cl_kj_flatnorm"]) for l in ok], "o-",
            color=S1_BLUE, lw=1.9, ms=6, label="round cl$_{KJ}$ (flat-norm)")
    ax.plot(h_ok, [float(D[l]["cl_p"]) for l in ok], "s--", color=S2_AQUA,
            lw=1.5, ms=5, label="round cl$_p$")
    # the flat ladder DOES complete the ramp at all three levels
    ax.plot(h_all, [0.2593, 0.2652, 0.2866], "x:", color=INK_2, lw=1.4, ms=7,
            label="flat cap (g132, all 3 reach M0.84)")
    ax.axhline(CL_REF, color=INK_2, lw=1.0, ls=":")
    ax.text(h_all[0], CL_REF, f" Tranair {CL_REF}", fontsize=7.5, color=INK_2,
            va="bottom")
    ax.axvline(h_all[-1], color=CRITICAL, lw=1.2, ls="--")
    ax.text(h_all[-1], ax.get_ylim()[0], " round FINE:\n ramp dies at M"
            f"{float(fine['m_final']):.2f}", fontsize=7.5, color=CRITICAL,
            va="bottom", ha="right")
    ax.set_xscale("log")
    ax.set_xlabel(r"mesh density $n_{\rm tet}^{1/3}\propto 1/h$")
    ax.set_ylabel("lift coefficient")
    ax.set_title(f"Round cap M{M_INF}: only 2 of 3 levels exist", loc="left",
                 fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(True, alpha=0.25)

    # right: how far up the Mach ramp each mesh got
    ax2.bar(range(3), [float(D[l]["m_final"]) for l in LEVELS], 0.5,
            color=[S2_AQUA if reached[l] else CRITICAL for l in LEVELS])
    ax2.axhline(M_INF, color=INK_2, lw=1.2, ls="--")
    ax2.text(2.45, M_INF, f" target M{M_INF}", fontsize=8, color=INK_2,
             va="bottom", ha="right")
    ax2.set_xticks(range(3))
    ax2.set_xticklabels([f"{l}\n{int(D[l]['ntet']):,} tets" for l in LEVELS],
                        fontsize=8)
    ax2.set_ylim(0.6, 0.9)
    ax2.set_ylabel("highest Mach level the ramp reached")
    ax2.set_title("The fine mesh never gets there", loc="left", fontsize=10)
    ax2.grid(True, axis="y", alpha=0.25)
    finish(fig, OUT, "g133rt_transonic.png")

    return cl.report(OUT, "checks_g133rt_transonic.csv")


if __name__ == "__main__":
    raise SystemExit(main())
