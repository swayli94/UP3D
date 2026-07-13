"""
Track B / B7 demo -- the ONERA M6 3D gate for the level-set solver.

B1-B6 all ran on quasi-2D meshes. A quasi-2D wake sheet is a flat strip: no
sweep, no tip, no spanwise direction. So three pieces of the level-set path had
never been exercised at all before B7:

1. The TE-POLYLINE ruled level set (D9). The sheet is ruled between a swept
   root->tip TE polyline, so its per-segment frame (v, d_hat, n_hat) is OBLIQUE
   -- the span axis is NOT perpendicular to the wake direction. B1 already found
   a real defect here (an orthogonal projection leaked the downstream distance
   into the spanwise coordinate and wrongly clipped ~60% of the true M6 cut set).

2. The SPANWISE CLIP (0 <= q <= span_length), which is what makes the
   circulation vanish at the tip DISCRETELY: elements that the infinite wake
   PLANE crosses but which lie outboard of the sheet are rejected
   (`CutElementMap.beyond_tip_elems`), so no jump is carried there. It is the
   level-set analogue of the conforming path's free-edge rule (Gamma(tip) = 0).

3. The g2 wake-BC component with a FREE spanwise jump gradient -- the trailing
   vortex DOF (D1). On a quasi-2D mesh there is no spanwise direction to be free
   in, so this term was structurally untested.

★ THE FAR FIELD IS THE 3D DECISION (figure `farfield_decision.png`). The B-path
vortex far field (`picard_ls._farfield_main`) is a SPAN-UNIFORM 2D point vortex
of strength mean(Gamma), whose branch cut is the ray y = 0, x > 0 -- the SAME cut
at every span station. On a 3D wing that is wrong in two independent ways, and
this demo measures both on M6 coarse at M0.5:

  (a) NON-COPLANARITY. The sheet is aimed at the incidence, direction
      (cos a, sin a, 0), so by the outlet (x ~ 10c) it has climbed to
      y ~ x tan(a) ~ 0.5 -- far from the vortex's y = 0 cut. The outlet then
      carries a prescribed Gamma jump that NO cut supports. This is exactly B3's
      recorded load-bearing rule ("the wake must be coplanar with the vortex
      branch cut, or the outlet carries an unsupported Dirichlet jump"), now in
      3D. Measured: a near-sonic spurious spot at the outlet, M_max 0.96 vs the
      neumann run's 0.77, with the excess concentrated at x ~ 10, y ~ 0.

  (b) SPAN-UNIFORMITY (P5's branch-ray artifact). Even re-aimed coplanar
      (direction (1,0,0)), one scalar Gamma cannot match Gamma(z), which decays
      to 0 at the tip: the prescribed jump is too small inboard and nonzero
      OUTBOARD OF THE TIP, where there is no cut at all. Measured: the outlet
      artifact shrinks (M_max 0.96 -> 0.83) but does not vanish -- vs 0.52 for
      neumann. The conforming path's fix for this was the Gamma(z) taper
      (`farfield_dirichlet(spanwise_gamma=True)`, P5).

  => B7 uses farfield="neumann" (the Lopez outlet). It carries NO vortex, so
  there is nothing to misapply: neither defect can exist, and no Gamma(z) taper
  is needed on the B path. It is also what let B6 converge near the fold (the
  live Gamma->vortex loop has gain > 1 there). The price is B5's measured
  O(Gamma/R) outlet truncation -- a few % of lift on a compact domain, which is
  why the gate bands are A/B bands, not <1% bands.

Figures:
  gamma_of_z.png         -- Gamma(z) root->tip on BOTH mesh families vs the P5
                            baseline; tip -> 0 (the spanwise clip, discretely).
  section_cp.png         -- upper/lower Cp(x/c) at eta = 0.44 / 0.65 / 0.90, both
                            families, vs P5, with the upper shock marked.
  shock_planform.png     -- the swept shock line on the planform: forward
                            migration toward the tip.
  farfield_decision.png  -- the (a)/(b) far-field study above.

Self-checking: results/checks.csv + summary.csv; exits nonzero on FAIL.
Heavy solves cache to results/*.npz (gitignored; delete to re-solve).

Standalone:  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
             python cases/demo/b7_onera_m6/run_demo.py
"""
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, CheckList, INK, INK_2, MESH_DIR, MUTED, S1_BLUE,
    S2_AQUA, S3_YELLOW, S4_ROSE, apply_style, finish, plt, write_csv,
)
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_le as x_le_fn, x_te  # noqa: E402
from pyfp3d.post.shock import shock_report  # noqa: E402
from pyfp3d.post.surface import cl_kj_3d, planform_area  # noqa: E402
from pyfp3d.post.surface_ls import (  # noqa: E402
    cl_pressure_3d_levelset, section_cp_curve_levelset, wall_cp_levelset,
)
from pyfp3d.solve.picard_ls import (  # noqa: E402
    solve_multivalued_lifting, solve_multivalued_transonic,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = Path(__file__).parent / "results"
M_INF, ALPHA = 0.84, 3.06
ETAS = (0.44, 0.65, 0.90)

# Committed baselines (roadmap P5 / P8 ledgers; cases/demo/p5_onera_m6/results/).
P5 = {"cl_kj": 0.24788, "shocks": (0.596, 0.570, 0.425),
      "gamma_root": 0.097, "gamma_tip": 0.0206, "m_max": 1.398}
P8 = {"cl_kj": 0.2692, "cl_p": 0.2560, "shocks": (0.596, 0.541, 0.362),
      "m_max": 2.13}

FAMILIES = {
    "M1 (wake-embedded)": MESH_DIR / "onera_m6" / "coarse.msh",
    "M4 (wake-free)": MESH_DIR / "onera_m6_wakefree" / "coarse.msh",
}


def setup(path, alpha_deg=ALPHA):
    mesh = read_mesh(path)
    a = np.radians(alpha_deg)
    wls = WakeLevelSet(np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
                       direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return mesh, cm, mvop


def solve_transonic(tag, path):
    """M0.84 / alpha 3.06 Mach-continuation solve, cached."""
    cache = OUT / f"{tag}.npz"
    mesh, cm, mvop = setup(path)
    if cache.exists():
        d = np.load(cache)
        print(f"  [{tag}] cached ({cache.name})")
        return mesh, cm, mvop, {k: d[k] for k in d.files}

    print(f"  [{tag}] solving M{M_INF} alpha {ALPHA} (minutes)...", flush=True)
    t0 = time.time()
    r = solve_multivalued_transonic(
        mvop, mesh, M_INF, alpha_deg=ALPHA, farfield="neumann",
        m_start=0.60, dm=0.04,
        n_outer_seed=120, n_outer_level=600, tol_residual=1e-7,
    )
    wall_s = time.time() - t0
    z = mesh.nodes[cm.te_nodes, 2]
    o = np.argsort(z)
    d = {
        "phi_ext": r["phi_ext"], "z": z[o], "gamma": r["te_jump"][o],
        # NOTE 2026-07-14: the solver's mach2_max now reads element_mach2
        # with the honest mixed_plain="main" default (flipped from "side").
        # The committed summary.csv m_max 1.453 (M1) / 1.368 (M4) came from
        # the cached side-based npz; the main re-reads are in
        # cases/demo/b8_tip_taper_ls/results/mmax_reread.csv.
        "m_max": float(np.sqrt(r["mach2_max"])), "wall_s": wall_s,
        "n_limited": r["n_limited"], "n_floored": r["n_floored"],
        "residual": r["residual_norm"],
        "level_m": np.array([L["m_inf"] for L in r["levels"]]),
        "level_gamma": np.array([L["gamma"] for L in r["levels"]]),
        "level_mmax": np.array([L["mach_max"] for L in r["levels"]]),
        "level_conv": np.array([L["converged"] for L in r["levels"]]),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **d)
    print(f"  [{tag}] {wall_s/60:.1f} min", flush=True)
    return mesh, cm, mvop, d


def farfield_study():
    """The M0.5 far-field decision study (cheap): neumann vs the span-uniform
    vortex, with the sheet aimed at alpha (non-coplanar) and re-aimed to y=0
    (coplanar) -- separating defect (a) from defect (b)."""
    cache = OUT / "farfield_study.npz"
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}

    path = MESH_DIR / "onera_m6" / "coarse.msh"
    mesh = read_mesh(path)
    c = mesh.nodes[mesh.elements].mean(axis=1)
    # the outlet region where the wake sheet leaves the domain
    outlet = (c[:, 0] > 6.0) & (np.abs(c[:, 1]) < 0.3) & (c[:, 2] < B_SEMI)
    rows = {}
    for aim, aim_tag in ((ALPHA, "alpha_aimed"), (0.0, "coplanar")):
        _, cm, mvop = setup(path, alpha_deg=aim)
        for ff in ("neumann", "vortex"):
            r = solve_multivalued_lifting(mvop, mesh, 0.5, alpha_deg=ALPHA,
                                          farfield=ff, n_outer_max=60,
                                          tol_residual=1e-7)
            # mixed_plain="side" pinned 2026-07-14 (default flipped to
            # "main"): the committed farfield_study.csv was measured through
            # the historical side reading.
            m = np.sqrt(mvop.element_mach2(r["phi_ext"], 0.5,
                                           mixed_plain="side"))
            rows[f"{aim_tag}_{ff}"] = np.array(
                [m.max(), m[outlet].max(), r["gamma"]])
            print(f"  farfield study {aim_tag:11s} {ff:9s} "
                  f"M_max {m.max():.3f}  outlet {m[outlet].max():.3f}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **rows)
    return rows


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("B7 -- ONERA M6 3D level-set gate")

    # ---------------------------------------------------------------- solves
    fams, summary = {}, []
    for tag, path in FAMILIES.items():
        key = tag.split()[0]
        if not path.exists():
            print(f"!! {path} missing (gitignored) -- regenerate it; see the "
                  "module header")
            return 2
        mesh, cm, mvop, d = solve_transonic(key, path)
        s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
        phi = d["phi_ext"]
        clkj = cl_kj_3d(d["gamma"], d["z"], s_ref=s_ref, b_semi=B_SEMI)
        cp = wall_cp_levelset(mesh, mvop, phi, m_inf=M_INF)
        clp = cl_pressure_3d_levelset(mesh, cp["cp"], cp["area"], cp["n_out"],
                                      ALPHA, s_ref)
        curves, shocks = {}, {}
        for eta in ETAS:
            curves[eta] = section_cp_curve_levelset(mesh, mvop, phi, eta=eta,
                                                    b_semi=B_SEMI, m_inf=M_INF)
            shocks[eta] = shock_report(curves[eta], M_INF)["upper"]
        fams[key] = dict(tag=tag, mesh=mesh, cm=cm, mvop=mvop, d=d, s_ref=s_ref,
                         clkj=clkj, clp=clp, curves=curves, shocks=shocks)
        summary.append((
            key, f"{clkj:.4f}", f"{clp:.4f}",
            f"{abs(clp - clkj) / abs(clkj) * 100:.2f}",
            *[f"{shocks[e]['x_shock']:.3f}" for e in ETAS],
            f"{d['gamma'][0]:.4f}", f"{d['gamma'][-1]:.4f}",
            f"{float(d['m_max']):.3f}", f"{int(d['n_limited'])}",
            f"{int(d['n_floored'])}", f"{float(d['wall_s']) / 60:.1f}",
        ))

    ff = farfield_study()

    # ---------------------------------------------------------------- checks
    for key, F in fams.items():
        d, sh = F["d"], F["shocks"]
        g = d["gamma"]

        cl.add("B7", f"{key}: field physical (0 limited / 0 floored)",
               f"{int(d['n_limited'])}/{int(d['n_floored'])}", "== 0/0",
               int(d["n_limited"]) == 0 and int(d["n_floored"]) == 0)
        cl.add("B7", f"{key}: M_max physical",
               f"{float(d['m_max']):.3f}", "< 2.5",
               float(d["m_max"]) < 2.5,
               note="P5 1.398 (Picard) / P8 2.13 (Newton)")
        # Honest convergence semantics: the top Mach levels park on the recorded
        # transonic Picard residual TAIL (~1e-6) instead of reaching
        # tol_residual. The field is bounded and physical; the LS Newton is the
        # cure and is deferred on M6 (needs lagged-LU -- P8/N6). We gate on
        # BOUNDED, not on converged.
        n_conv = int(np.sum(d["level_conv"]))
        cl.add("B7", f"{key}: residual bounded (Picard tail, not divergence)",
               f"{float(d['residual']):.1e}", "< 1e-4",
               float(d["residual"]) < 1e-4,
               note=f"{n_conv}/{len(d['level_conv'])} levels reached "
                    "tol_residual=1e-7; the rest park on the tail -> LS Newton "
                    "(deferred on M6)")
        # -- Gamma(z): the 3D-only machinery (spanwise clip => tip -> 0)
        cl.add("B7", f"{key}: Gamma(z) tip -> 0 (spanwise clip)",
               f"{g[-1]:.4f}", "< 0.02 (P5: 0.0206)", abs(g[-1]) < 0.02)
        cl.add("B7", f"{key}: Gamma(z) positive everywhere",
               f"min {g.min():.4f}", "> -1e-3", g.min() > -1e-3)
        half = len(g) // 2
        viol = int((np.diff(g[half:]) > 1e-3).sum())
        cl.add("B7", f"{key}: Gamma(z) decays monotonically on the outer half",
               f"{viol} violations", "== 0", viol == 0)
        # -- lift: gated against the conforming NEWTON truth (P8), not the
        # conforming Picard (P5) -- the B6 user arbitration, which B7 reproduces
        # in 3D: the LS Picard tracks the Newton truth while the conforming
        # Picard under-circulates below it.
        rel_p8 = (F["clkj"] - P8["cl_kj"]) / P8["cl_kj"] * 100
        rel_p5 = (F["clkj"] - P5["cl_kj"]) / P5["cl_kj"] * 100
        cl.add("B7", f"{key}: cl_KJ vs the conforming NEWTON truth",
               f"{F['clkj']:.4f} ({rel_p8:+.1f}% of P8)", "within +-10% of P8",
               abs(rel_p8) < 10.0,
               note=f"P8 Newton {P8['cl_kj']:.4f}; conforming Picard P5 "
                    f"{P5['cl_kj']:.4f} ({rel_p5:+.1f}% away) under-circulates")
        v6 = abs(F["clp"] - F["clkj"]) / abs(F["clkj"]) * 100
        cl.add("B7", f"{key}: V6 circulation/pressure consistency",
               f"{v6:.2f}%", "< 5%", v6 < 5.0,
               note="P5 coarse 2.40%")
        # -- shocks: present, monotone, no expansion shock, forward migration
        for eta in ETAS:
            s = sh[eta]
            ok = (s["has_shock"] and s["monotone"] and not s["expansion_shock"])
            cl.add("B7", f"{key}: eta={eta} upper shock present & monotone",
                   f"x/c {s['x_shock']:.3f}", "has & monotone & no expansion", ok)
        for i, eta in enumerate(ETAS):
            band = abs(sh[eta]["x_shock"] - P5["shocks"][i])
            cl.add("B7", f"{key}: eta={eta} shock within +-0.06 of P5",
                   f"{sh[eta]['x_shock']:.3f} (P5 {P5['shocks'][i]:.3f})",
                   "|dx/c| <= 0.06", band <= 0.06)
        fwd = sh[0.90]["x_shock"] <= sh[0.65]["x_shock"] + 0.05
        cl.add("B7", f"{key}: shock migrates FORWARD toward the tip",
               f"x90 {sh[0.90]['x_shock']:.3f} vs x65 {sh[0.65]['x_shock']:.3f}",
               "x90 <= x65 + 0.05", fwd)

    # ★ Spanwise SMOOTHNESS of Gamma(z) -- an unplanned B7 finding, visible the moment
    # the real P5 curve is overlaid (gamma_of_z.png). The conforming path solves a
    # SEPARATE secant per TE station, so its Gamma(z) carries station-to-station
    # jitter (P5's INVESTIGATION_gamma_smoothing.md chased exactly this, and recorded
    # that smoothing it moved Gamma AWAY from the self-consistent value). The implicit
    # Kutta has no per-station loop to be noisy in: Gamma is a single solution mode of
    # the coupled system, so it comes out spanwise-smooth with NO smoothing applied.
    # Metric: RMS 2nd difference of Gamma(z), normalised by the curve's own range.
    def roughness(z, g):
        d2 = g[:-2] - 2.0 * g[1:-1] + g[2:]
        return float(np.sqrt(np.mean(d2 ** 2)) / (g.max() - g.min()))

    p5_csv = Path(__file__).parent / "p5_gamma_baseline.csv"
    n_skip = sum(1 for ln in p5_csv.read_text().splitlines()
                 if ln.startswith("#") or ln.startswith("z,"))
    p5_curve = np.loadtxt(p5_csv, delimiter=",", comments="#", skiprows=n_skip)
    r_p5 = roughness(p5_curve[:, 0], p5_curve[:, 1])
    for key in ("M1", "M4"):
        r_ls = roughness(fams[key]["d"]["z"], fams[key]["d"]["gamma"])
        cl.add("B7", f"{key}: Gamma(z) spanwise-smooth WITHOUT any smoothing",
               f"{r_ls:.4f} ({r_p5 / r_ls:.0f}x smoother than P5)",
               "< 1/3 of the conforming P5 jitter", r_ls < r_p5 / 3.0,
               note=f"P5 per-station secant jitter {r_p5:.4f}; the implicit Kutta has "
                    "no per-station loop (Gamma is one solution mode)")

    ab = abs(fams["M1"]["clkj"] - fams["M4"]["clkj"]) / abs(fams["M1"]["clkj"]) * 100
    cl.add("B7", "dual-mesh A/B: cl_KJ (M1 embedded vs M4 wake-free)",
           f"{ab:.2f}%", "< 10%", ab < 10.0,
           note="M4 is within 6-9% of M1's tet count at equal h_wall "
                "(the controlled comparison)")

    # -- the far-field decision, measured
    a_neu = float(ff["alpha_aimed_neumann"][1])
    a_vor = float(ff["alpha_aimed_vortex"][1])
    c_vor = float(ff["coplanar_vortex"][1])
    cl.add("B7", "far field (a): span-uniform vortex + alpha-aimed sheet is "
           "NOT coplanar with the branch cut -> outlet artifact",
           f"outlet M {a_vor:.3f} vs neumann {a_neu:.3f}",
           "vortex >> neumann", a_vor > a_neu + 0.2)
    cl.add("B7", "far field (b): re-aiming coplanar SHRINKS but does not remove "
           "it (span-uniform Gamma != Gamma(z))",
           f"outlet M {c_vor:.3f} (alpha-aimed {a_vor:.3f}, neumann {a_neu:.3f})",
           "neumann < coplanar < alpha-aimed",
           a_neu < c_vor < a_vor)

    # --------------------------------------------------------------- figures
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for key, color in (("M1", S1_BLUE), ("M4", S2_AQUA)):
        F = fams[key]
        ax.plot(F["d"]["z"] / B_SEMI, F["d"]["gamma"], "-o", ms=3, color=color,
                label=f"level set, {F['tag']}")
    # committed P5 baseline curve (its own solution .npz is a gitignored cache)
    p5_csv = Path(__file__).parent / "p5_gamma_baseline.csv"
    n_skip = sum(1 for ln in p5_csv.read_text().splitlines()
                 if ln.startswith("#") or ln.startswith("z,"))
    p5g = np.loadtxt(p5_csv, delimiter=",", comments="#", skiprows=n_skip)
    ax.plot(p5g[:, 0] / B_SEMI, p5g[:, 1], "--", color=BASELINE, lw=1.5,
            label="P5 conforming Picard")
    ax.axhline(0.0, color=MUTED, lw=0.8)
    ax.set_xlabel("span fraction  z / b")
    ax.set_ylabel(r"$\Gamma(z)$")
    ax.set_title(r"B7: spanwise circulation, M6 coarse  M$_\infty$=0.84, "
                 r"$\alpha$=3.06$^\circ$" "\n"
                 r"tip $\Gamma\to0$ is the spanwise clip, discretely")
    ax.legend(frameon=False, fontsize=9)
    finish(fig, OUT, "gamma_of_z.png")

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), sharey=True)
    for ax, eta in zip(axes, ETAS):
        for key, color in (("M1", S1_BLUE), ("M4", S2_AQUA)):
            c = fams[key]["curves"][eta]
            ax.plot(c["x_upper"], c["cp_upper"], "-", color=color, lw=1.3,
                    label=f"{key} upper")
            ax.plot(c["x_lower"], c["cp_lower"], "--", color=color, lw=1.0,
                    alpha=0.65, label=f"{key} lower")
            ax.axvline(fams[key]["shocks"][eta]["x_shock"], color=color, lw=0.8,
                       ls=":", alpha=0.8)
        ax.axvline(P5["shocks"][ETAS.index(eta)], color=CRITICAL, lw=1.0,
                   ls="-.", label="P5 shock")
        ax.axhline(shock_report(fams["M1"]["curves"][eta], M_INF)["cp_critical"],
                   color=MUTED, lw=0.8, ls="--")
        ax.set_title(rf"$\eta$ = {eta}")
        ax.set_xlabel("x / c")
    axes[0].set_ylabel(r"$C_p$")
    axes[0].invert_yaxis()
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle(r"B7: section $C_p$ on the level-set path, M6 coarse "
                 r"(dashed grey = $C_p^*$)", y=1.02)
    finish(fig, OUT, "section_cp.png")

    # Two panels, because the two true statements point opposite ways and one plot
    # alone reads as a contradiction: on a SWEPT wing the shock moves AFT in absolute
    # x toward the tip (the whole section is further aft), while it moves FORWARD in
    # x/c -- and x/c is what the gate asserts (x90 <= x65 + 0.05, the P5 G5.1 clause).
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.8))
    zz = np.linspace(0.0, B_SEMI, 60)
    axL.plot([x_le_fn(z) for z in zz], zz, "-", color=INK_2, lw=1.2, label="LE / TE")
    axL.plot([x_te(z) for z in zz], zz, "-", color=INK_2, lw=1.2)
    for key, color, mk in (("M1", S1_BLUE, "o"), ("M4", S2_AQUA, "s")):
        F = fams[key]
        xs = [F["curves"][e]["x_le"] + F["shocks"][e]["x_shock"] * F["curves"][e]["chord"]
              for e in ETAS]
        axL.plot(xs, [e * B_SEMI for e in ETAS], "-", marker=mk, ms=6, color=color,
                 mfc="none", label=f"{key} shock")
    xs5 = [fams["M1"]["curves"][e]["x_le"] + P5["shocks"][ETAS.index(e)]
           * fams["M1"]["curves"][e]["chord"] for e in ETAS]
    axL.plot(xs5, [e * B_SEMI for e in ETAS], "--^", color=BASELINE, ms=6,
             label="P5 shock")
    axL.set_xlabel("x")
    axL.set_ylabel("z (span)")
    axL.set_title("swept planform: the shock line in absolute x\n"
                  "(it tracks the sweep — aft toward the tip)")
    axL.legend(frameon=False, fontsize=8.5, loc="upper left")

    for key, color, mk in (("M1", S1_BLUE, "o"), ("M4", S2_AQUA, "s")):
        F = fams[key]
        axR.plot([F["shocks"][e]["x_shock"] for e in ETAS], [e for e in ETAS],
                 "-", marker=mk, ms=6, color=color, mfc="none", label=f"{key} shock")
    axR.plot(list(P5["shocks"]), list(ETAS), "--^", color=BASELINE, ms=6,
             label="P5 shock")
    axR.set_xlabel("x / c  (fraction of the LOCAL chord)")
    axR.set_ylabel(r"span fraction  $\eta$ = z / b")
    axR.set_title("the gate clause: FORWARD migration in x/c\n"
                  r"($x_{90} \leq x_{65} + 0.05$)")
    axR.legend(frameon=False, fontsize=8.5)
    fig.suptitle("B7: shock position across the span, M6 coarse", y=1.03)
    finish(fig, OUT, "shock_planform.png")

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    labels = ["neumann\n(no vortex)", "vortex,\ncoplanar sheet",
              "vortex,\n" r"$\alpha$-aimed sheet"]
    vals = [a_neu, c_vor, a_vor]
    bars = ax.bar(labels, vals, color=[S2_AQUA, S3_YELLOW, CRITICAL], width=0.55)
    ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=9)
    ax.axhline(1.0, color=INK, lw=1.0, ls="--")
    ax.text(0.02, 1.01, "sonic", transform=ax.get_yaxis_transform(), fontsize=8,
            color=INK)
    ax.set_ylabel("max local Mach at the OUTLET")
    ax.set_ylim(0, 1.15)
    ax.set_title(r"B7 far-field decision (M$_\infty$=0.5, M6 coarse): the "
                 "span-uniform\nvortex prescribes a jump no cut supports; "
                 "neumann has no vortex")
    finish(fig, OUT, "farfield_decision.png")

    # ----------------------------------------------------------------- CSVs
    write_csv(OUT, "summary.csv",
              "mesh,cl_kj,cl_p_3d,v6_pct,shock_eta044,shock_eta065,shock_eta090,"
              "gamma_root,gamma_tip,m_max,n_limited,n_floored,solve_min",
              summary)
    write_csv(OUT, "farfield_study.csv",
              "case,m_max_global,m_max_outlet,gamma",
              [(k, f"{v[0]:.4f}", f"{v[1]:.4f}", f"{v[2]:.4f}")
               for k, v in ff.items()])
    return cl.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
