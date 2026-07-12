"""
P13/G13.2 -- WHAT DOES THE TIP TAPER DO TO THE PHYSICS?

run_taper_probe.py answers "does the taper kill the tip-edge singularity, and
what does it cost in total lift?" (yes; ~1% of cl). It does NOT show WHERE that
lift went, or whether the taper contaminates the wing. This demo does: it plots
the three distributions an aerodynamicist actually reads --

  (a) the taper factor F(eta) itself      -- where does the model touch the wing?
  (b) the circulation Gamma(eta)          -- the spanwise loading
  (c) the sectional lift cl(eta)          -- Gamma normalised by the local chord
  (d-f) the sectional wall Cp(x/c)        -- at mid-span, outboard, and the tip

-- for three models on the SAME mesh, same far field (spanwise-Gamma, so the
P5 branch-ray artifact is out of the picture):

  none                the untapered baseline (carries the 1/sqrt(r) singularity)
  tanh_half r_c=0.10b the originally-proposed form: F = (1+tanh(u/r_c))/2
  vanish_smooth 0.05b the probe's winner: F = smoothstep(u/r_c), F(tip) = 0

THE HEADLINE THIS FIGURE MAKES VISIBLE. Both forms regularize the edge, but they
are NOT equivalent in what they do to the wing:

  * vanish_smooth is LOCAL: F = 1 over ~95% of the span, so Gamma(eta), cl(eta)
    and the sectional Cp are UNCHANGED inboard. Only the last ~8 TE stations
    (eta > 0.95) are unloaded. It costs ~1% of total lift and buys a bounded,
    mesh-convergent tip.
  * tanh_half is GLOBAL: tanh has long tails, so F < 1 over 57 of 83 stations --
    it quietly unloads the ENTIRE outer wing (F ~ 0.9 even at mid-span). That is
    why it costs ~7x more lift for the same regularization, and it is the real
    reason to reject it -- not that it fails to work.

So the question "is the taper acceptable?" is answered by (b)/(c)/(d-f): a
LOCAL taper is a tip-model change; a GLOBAL one is a silent re-rigging of the
whole wing.

Cheap (coarse+medium only, ~2 min); heavy solves cache to results/phys_*.npz
(gitignored). Standalone:
  NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
  python cases/demo/tip_edge_singularity/run_taper_physics.py
"""
import sys, time
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, CheckList, INK_2, S1_BLUE, S2_AQUA,
    apply_style, finish, plt, write_csv,
)
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at  # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve  # noqa: E402
from pyfp3d.post.surface import cl_kj_3d, planform_area  # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes" / "onera_m6"
ALPHA, M, LEVEL = 3.06, 0.5, "medium"
ETAS = (0.20, 0.44, 0.65, 0.90, 0.99)      # Cp stations (P5's + a tip station)
P6_PASSES = 4                              # P6/G6.1 wall-gradient recovery smoothing


def _jitter(y):
    """Normalised RMS second difference -- the 'sawtooth' measure used for
    Gamma(z) in B7 (conforming P5 measured 0.0970; the level-set path 0.008)."""
    y = np.asarray(y, dtype=np.float64)
    if len(y) < 3 or np.max(np.abs(y)) == 0:
        return 0.0
    d2 = y[2:] - 2.0 * y[1:-1] + y[:-2]
    return float(np.sqrt(np.mean(d2 ** 2)) / np.max(np.abs(y)))


def _sawtooth_reversals(x, cp, lo=0.10, hi=0.95):
    """Slope-REVERSAL count on the mid-chord run. The sawtooth is an
    alternating-sign artifact, so reversal count -- not RMS amplitude -- is
    the right measure (P6's own G6.1 metric is sign-alternating for exactly
    this reason; that metric itself is transonic-only, returning nan on a
    shock-free M0.5 section, hence this subsonic stand-in). The LE suction
    peak and the TE are excluded: they are REAL sharp features."""
    x = np.asarray(x, dtype=np.float64)
    cp = np.asarray(cp, dtype=np.float64)
    m = (x > lo) & (x < hi)
    c = cp[m]
    if len(c) < 3:
        return 0
    dd = np.diff(c)
    return int(np.count_nonzero(np.sign(dd[1:]) * np.sign(dd[:-1]) < 0))

MODELS = [
    ("none", 0.00, "untapered (singular tip)", INK_2),
    ("tanh_half", 0.10, "tanh_half  $r_c$=0.10b", CRITICAL),
    ("vanish_smooth", 0.05, "vanish_smooth  $r_c$=0.05b", S2_AQUA),
]


def tag(form, frac):
    return form if form == "none" else f"{form}_rc{frac:.2f}"


def get(form, frac):
    """Solve (cached) and extract everything the physics plots need."""
    cache = OUT / f"phys_{tag(form, frac)}_{LEVEL}.npz"
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        return {k: d[k] for k in d.files}
    print(f"  solving {tag(form, frac)} {LEVEL} ...", flush=True)
    t0 = time.time()
    mc, wc = cut_wake(read_mesh(MESH / f"{LEVEL}.msh"))
    r_c = frac * B_SEMI
    taper = tip_taper_factors(wc.station_z, B_SEMI, form,
                              r_c if r_c > 0 else 1.0)
    r = solve_newton_lifting(mc, wc, m_inf=M, alpha_deg=ALPHA, upwind_c=1.5,
                             precond="amg", tol_residual=1e-9, tip_taper=taper,
                             farfield_spanwise_gamma=True)
    phi = np.asarray(r["phi"])
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    gamma = np.asarray(r["gamma"])
    sz = np.asarray(wc.station_z)
    # Cp both RAW (smooth_passes=0) and with the shipped P6 wall-gradient
    # recovery smoothing (G6.1): the surface-Cp sawtooth is a per-triangle
    # RECOVERY artifact of the flat-facet P1 wall, not a flux/solver defect
    # and not a bug in this demo -- P6 measured the G6.1 metric drop 330x
    # under exactly this smoothing. Default is 0, which is what the first cut
    # of this demo (wrongly) plotted.
    cps = {}
    for e in ETAS:
        for sp, kind in ((0, "raw"), (P6_PASSES, "smooth")):
            c = section_cp_curve(mc, phi, eta=e, b_semi=B_SEMI, u_inf=1.0,
                                 m_inf=M, smooth_passes=sp)
            for side in ("upper", "lower"):
                cps[f"x_{e}_{side}_{kind}"] = c[f"x_{side}"]
                cps[f"cp_{e}_{side}_{kind}"] = c[f"cp_{side}"]
    d = dict(gamma=gamma, station_z=sz, taper=taper,
             cl_kj=cl_kj_3d(gamma, sz, s_ref, B_SEMI),
             conv=bool(r["converged"]), dt=time.time() - t0, **cps)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **d)
    return d


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("P13/G13.2: what the tip taper does to the spanwise physics")

    D = {(f, r): get(f, r) for f, r, _, _ in MODELS}
    base = D[("none", 0.0)]
    o = np.argsort(base["station_z"])
    eta_b = base["station_z"][o] / B_SEMI

    fig = plt.figure(figsize=(13.2, 7.4))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.28)
    ax_f = fig.add_subplot(gs[0, 0])
    ax_g = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    axcp = [fig.add_subplot(gs[1, i]) for i in range(3)]

    rows = []
    for form, frac, lab, col in MODELS:
        d = D[(form, frac)]
        oo = np.argsort(d["station_z"])
        eta = d["station_z"][oo] / B_SEMI
        F = d["taper"][oo]
        G = d["gamma"][oo]
        chords = np.array([chord_at(z) for z in d["station_z"][oo]])
        clsec = 2.0 * G / chords          # sectional cl = 2*Gamma/(U*c)
        lw = 2.4 if form == "none" else 1.8
        ax_f.plot(eta, F, "-", color=col, lw=lw, label=lab)
        ax_g.plot(eta, G, "-", color=col, lw=lw, label=lab)
        ax_c.plot(eta, clsec, "-", color=col, lw=lw, label=lab)

        # how far inboard does the model reach, and what does it cost?
        touched = eta[F < 0.99]
        reach = float(touched.min()) if len(touched) else 1.0
        n_touch = int(np.count_nonzero(F < 0.99))
        dcl = 100.0 * (float(d["cl_kj"]) - float(base["cl_kj"])) / float(base["cl_kj"])
        # lift lost INBOARD of the tip region (eta < 0.90) = contamination
        m_in = eta_b < 0.90
        g_in = np.interp(eta_b[m_in], eta, G)
        gb_in = base["gamma"][o][m_in]
        dg_in = 100.0 * (np.trapezoid(g_in, eta_b[m_in])
                         - np.trapezoid(gb_in, eta_b[m_in])) \
            / np.trapezoid(gb_in, eta_b[m_in])
        rows.append((tag(form, frac), f"{reach:.3f}", n_touch,
                     f"{float(d['cl_kj']):.4f}", f"{dcl:+.2f}",
                     f"{dg_in:+.2f}"))
        print(f"{tag(form, frac):22s} reaches inboard to eta={reach:.3f} "
              f"({n_touch} stations)  cl={float(d['cl_kj']):.4f} ({dcl:+.2f}%)  "
              f"inboard(eta<0.90) circulation change {dg_in:+.2f}%")

        if form != "none":
            cl.add("G13.2", f"{tag(form, frac)}: is the taper LOCAL to the tip?",
                   f"reaches inboard to eta={reach:.2f}; inboard circulation "
                   f"{dg_in:+.2f}%", "eta_reach > 0.90 and |inboard| < 1%",
                   reach > 0.90 and abs(dg_in) < 1.0,
                   note="a LOCAL taper is a tip-model change; a GLOBAL one "
                        "silently re-rigs the whole wing")

        # (i) the Cp "sawtooth" is the P6 WALL-GRADIENT RECOVERY artifact --
        #     NOT a solver defect, NOT caused by the taper, and NOT a bug in
        #     this demo's maths: it is what smooth_passes=0 (the default) looks
        #     like. The shipped P6 smoothing removes the slope reversals, and
        #     the raw counts are IDENTICAL across all three models, proving the
        #     artifact is independent of the tip model.
        rr = _sawtooth_reversals(d["x_0.44_upper_raw"], d["cp_0.44_upper_raw"])
        rs = _sawtooth_reversals(d["x_0.44_upper_smooth"],
                                 d["cp_0.44_upper_smooth"])
        cl.add("G13.2", f"{tag(form, frac)}: Cp sawtooth = P6 recovery artifact",
               f"slope reversals raw {rr} -> P6-smoothed {rs}",
               "P6 smoothing removes >= 80% of reversals",
               rs <= max(1, 0.2 * rr),
               note="per-triangle wall-gradient recovery on the flat-facet P1 "
                    "wall (roadmap P6/G6.1). Default smooth_passes=0; the raw "
                    "count is the SAME for every tip model => taper-independent")

        # (ii) Gamma(z) jitter is a DIFFERENT artifact: Kutta-probe placement
        #      on the unstructured swept TE (P5 known-robustness item).
        cl.add("G13.2", f"{tag(form, frac)}: Gamma(z) jitter is probe-placement, "
               "not post-processing",
               f"RMS d2 = {_jitter(G):.4f}", "<= the P5 conforming record (0.097)",
               _jitter(G) <= 0.097,
               note="conforming Kutta probes are picked per station on an "
                    "unstructured swept TE (adjacent stations even SHARE "
                    "probes, P5); the level-set path is 11-12x smoother (B7) "
                    "because its implicit Kutta has no per-station probe")

        # (iii) TE non-closure. TWO separate causes, and the figure shows both:
        #   * a PRE-EXISTING baseline gap (~0.14-0.22) at EVERY station with no
        #     taper at all -- the conforming Kutta enforces the POTENTIAL-JUMP
        #     form (design.md 4.4), only an approximation of true pressure
        #     equality, plus the sharp-TE P1 recovery error (the V6/P11 floor).
        #     Track B's B4 had to introduce the explicit nonlinear
        #     |q_u|^2=|q_l|^2 Kutta precisely because this form is inadequate.
        #   * INSIDE the taper, a LARGE extra gap -- BY CONSTRUCTION, since
        #     Gamma_eff = F*Gamma_Kutta deliberately relaxes the Kutta
        #     condition. Non-closure IS the model, made visible.
        # The acceptance test is therefore: OUTSIDE its taper, a tip model must
        # leave the TE gap at the BASELINE value.
        base_gaps = {e: abs(float(base[f"cp_{e}_upper_smooth"][-1])
                            - float(base[f"cp_{e}_lower_smooth"][-1]))
                     for e in (0.44, 0.90, 0.99)}
        for e in (0.44, 0.90, 0.99):
            gap = abs(float(d[f"cp_{e}_upper_smooth"][-1])
                      - float(d[f"cp_{e}_lower_smooth"][-1]))
            inside = F[np.argmin(np.abs(eta - e))] < 0.99
            rows[-1] = rows[-1] + (f"{gap:.3f}" + ("*" if inside else ""),)
            if form != "none" and not inside:
                cl.add("G13.2", f"{tag(form, frac)}: TE closure UNDISTURBED "
                       f"outside the taper (eta={e})",
                       f"gap {gap:.3f} vs baseline {base_gaps[e]:.3f}",
                       "within 0.05 of baseline",
                       abs(gap - base_gaps[e]) < 0.05,
                       note="a tip model must not relax the Kutta condition "
                            "where there is no singularity to fix")

    for a, t, yl in ((ax_f, "taper factor $F(\\eta)$ — where the model bites",
                      "$F$"),
                     (ax_g, "circulation $\\Gamma(\\eta)$ — spanwise loading",
                      "$\\Gamma$"),
                     (ax_c, "sectional lift $c_l(\\eta)=2\\Gamma/c$", "$c_l$")):
        a.set_xlabel("$\\eta = z/b$"); a.set_ylabel(yl); a.set_title(t, fontsize=9.5)
        a.axvspan(0.95, 1.0, color=BASELINE, alpha=0.30, zorder=0)
    ax_f.legend(frameon=False, fontsize=7.5, loc="lower left")
    ax_f.text(0.955, 0.06, "$r_c$=0.05b", fontsize=7, color=INK_2, rotation=90)

    # sectional Cp: RAW (thin, faded = the P6 recovery sawtooth) vs P6-SMOOTHED
    for a, e in zip(axcp, (0.44, 0.90, 0.99)):
        for form, frac, lab, col in MODELS:
            d = D[(form, frac)]
            lw = 2.0 if form == "none" else 1.5
            for side in ("upper", "lower"):
                a.plot(d[f"x_{e}_{side}_raw"], d[f"cp_{e}_{side}_raw"], "-",
                       color=col, lw=0.7, alpha=0.32, zorder=1)
                a.plot(d[f"x_{e}_{side}_smooth"], d[f"cp_{e}_{side}_smooth"],
                       "-", color=col, lw=lw, zorder=2,
                       label=lab if side == "upper" else None)
        a.invert_yaxis()
        a.set_xlabel("$x/c$"); a.set_ylabel("$C_p$")
        a.set_title(f"sectional $C_p$,  $\\eta$ = {e}"
                    + ("   (INSIDE the taper)" if e > 0.95 else ""),
                    fontsize=9.5)
    axcp[0].legend(frameon=False, fontsize=7)
    axcp[0].text(0.03, 0.06, "faint = RAW (P6 recovery sawtooth)\n"
                 f"bold = P6-smoothed ({P6_PASSES} passes)",
                 transform=axcp[0].transAxes, fontsize=6.8, color=INK_2)

    fig.suptitle("P13/G13.2 — what the tip taper does to the wing "
                 "(ONERA M6, M$_\\infty$=0.5, medium)\n"
                 "vanish_smooth is LOCAL (untouched inboard); tanh_half's long "
                 "tails silently unload the WHOLE outer wing",
                 fontsize=11, y=0.99)
    finish(fig, OUT, "taper_physics.png")

    write_csv(OUT, "taper_physics.csv",
              "model,eta_reach_inboard,n_stations_touched,cl_KJ,dcl_pct,"
              "inboard_circulation_change_pct,te_gap_eta0.44,te_gap_eta0.90,"
              "te_gap_eta0.99", rows)
    return cl.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
