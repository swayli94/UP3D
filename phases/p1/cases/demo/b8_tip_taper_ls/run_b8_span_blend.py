"""
Track B / B8 re-spec — SPAN BLEND of the wake LS rows (sheet-termination
desingularization), measured on the HONEST metric.

★★ RESULT (2026-07-14): NEGATIVE, with the mechanism pinned. The span blend
does exactly what it was designed to do to its target -- it WELDS the
termination-ring jump (max|delta| 0.026 -> 0.001/0.0003 at r_blend 0.05/0.08)
-- and the price is a ~20% GLOBAL lift loss that disqualifies it as a model:

  (1) the loss is UNIFORM in z (the whole Gamma(z) profile scales down ~0.8x,
      root included, where the blended rows are BITWISE identical to
      baseline) and essentially r_blend-INSENSITIVE (-19.8/-20.2/-21.8% at
      r_blend 0.03/0.05/0.08 medium);
  (2) component isolation (coarse): the beyond-tip straddler weld ALONE costs
      -13.3%, the inboard smooth blend ALONE costs -10.8% -- BOTH components
      of a sheet-side delta-pin are amplified, so this is not a bug in either
      piece. ⇒ the implicit Kutta has NO per-station target: Gamma is ONE
      global solution mode, and ANY constraint intervention on the sheet near
      the tip re-levels it (the G13.2 finding-(2) fixed-point amplification,
      appearing GLOBALLY on the LS path where the conforming secant keeps it
      per-station -- which is why the conforming taper costs -1.6% and this
      costs -20%);
  (3) the loss GROWS under refinement (rb0.08: -16.9% coarse -> -21.8%
      medium), so it CONTAMINATES the refinement exponent: the h-growing
      unloading contributes p_unload ~ -0.10 to the raw p = +0.05, i.e. a
      THIRD of the apparent +0.37 -> +0.05 reduction; the
      confounder-corrected exponent is ~ +0.15 -- suggestive of a real
      partial regularization, but NOT certifiable from a 2-point ladder
      under a 20% global flow distortion, and moot at this cost;
  (4) r_blend 0.03 (~2 elements at medium) is ACTIVELY harmful: it steepens
      the termination (ring jump x2.3, honest p +1.31) -- a blend narrower
      than a few elements re-creates the Heaviside it was meant to remove.

WHAT SURVIVES (from run_b8_termination_diagnosis.py, independent of the
blend): the committed LS exponent p = +1.34 was a x5 METRIC ARTIFACT
(element_mach2 read mixed-side plain/beyond-tip elements from the
aux-substituted side field; the assembled field gives +0.62, no-sliver
+0.37 -- the SAME object as the conforming +0.52); the honest residual
singularity is the finite termination-ring jump (|delta| ~ 0.026, h- and
TE-taper-independent). The cure, if one is needed, must change the FUNCTION
SPACE at the termination (how the clip ends the multivalued region), not
add sheet-side constraints -- every constraint-side route is now measured:
TE rows (old B8: no effect), wake-LS rows (this: global 10x-amplified
damage).

MODEL/BLEND MACHINERY: `MultivaluedOperator(span_blend=(form, r_blend))`,
per non-TE cut node j:
    w_j * [wake LS row]  +  (1 - w_j) * s_j * [phi_aux - phi_main] = 0
w_j = tip_taper_factors(q_j, span_length, form, r_blend); s_j = the row's
own LS magnitude (h-invariance of the blend itself); beyond-tip straddler
nodes (q > span_length) get w = 0 at any r_blend. Default None is
bit-identical (tests/test_b8_span_blend.py, B suite 116 green).

Baseline states are READ from the diagnosis caches (results/
b8_diag_none_*.npz); blend solves cache to results/b8_blend_*.npz;
isolation solves to results/b8_iso_*.npz. Standalone:
  NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
  python cases/demo/b8_tip_taper_ls/run_b8_span_blend.py
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
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te  # noqa: E402
from pyfp3d.post.surface import cl_kj_3d, planform_area  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes" / "onera_m6"
ALPHA, M = 3.06, 0.5
LEVELS = ("coarse", "medium")
TIP_AFT, TIP_Y = 0.30, 0.03            # the G13.1/G13.2/B8 tip-edge box
V_REL_MIN = 0.1                        # sliver filter (committed B8 CSV:
                                       # the rc>=0.05 TE-taper peak sat on a
                                       # V/median = 0.03 tet)

# The diagnosis numbers this demo must reproduce (metric identity).
DIAG_BASE = {"e_ns": (0.551, 0.691)}

VARIANTS = [("none", 0.00),
            ("vanish_smooth", 0.03),
            ("vanish_smooth", 0.05),
            ("vanish_smooth", 0.08)]

# Component-isolation experiments (coarse only, the discriminator for
# finding (2)): custom per-node weights instead of a (form, r_blend) recipe.
ISOLATIONS = ("straddler_only", "inboard_only")


def tag(form, frac):
    return form if form == "none" else f"{form}_rb{frac:.2f}"


_GEOM = {}


def geom(level, span_blend=None):
    key = (level, span_blend)
    if key not in _GEOM:
        mesh = read_mesh(MESH / f"{level}.msh")
        a = np.radians(ALPHA)
        wls = WakeLevelSet(
            np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
            direction=(np.cos(a), np.sin(a), 0.0))
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]))
        mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                   levelset=wls, span_blend=span_blend)
        _GEOM[key] = (mesh, wls, cm, mvop)
    return _GEOM[key]


def _nonte_and_q(cm):
    cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
    is_te = np.zeros(cm.n_main, dtype=bool)
    is_te[cm.te_nodes] = True
    nonte = cut_nodes[~is_te[cut_nodes]]
    return nonte, cm.q[nonte]


def _apply_custom_w(mvop, cm, w):
    """The span-blend _ls_coo transform with an arbitrary per-nonte-node
    weight vector (the isolation experiments; mirrors MultivaluedOperator)."""
    ls_row, ls_col, ls_data = mvop._ls_coo
    nonte, _ = _nonte_and_q(cm)
    aux = cm.ext_dof_of_node[nonte]
    s_of_row = np.zeros(mvop.n_total)
    np.add.at(s_of_row, ls_row, np.abs(ls_data))
    s_of_row *= 0.5
    w_of_row = np.ones(mvop.n_total)
    w_of_row[aux] = w
    ls_data = ls_data * w_of_row[ls_row]
    sel = w < 1.0
    aux_b, nodes_b = aux[sel], nonte[sel]
    coef = (1.0 - w[sel]) * s_of_row[aux_b]
    mvop._ls_coo = (np.concatenate([ls_row, aux_b, aux_b]),
                    np.concatenate([ls_col, aux_b, nodes_b]),
                    np.concatenate([ls_data, coef, -coef]))
    return int(sel.sum())


def _solve(mvop, mesh):
    return solve_multivalued_lifting(
        mvop, mesh, M, alpha_deg=ALPHA, farfield="neumann", upwind_c=0.0,
        n_outer_max=120, tol_residual=1e-7)


def solve_phi(form, frac, level):
    """phi_ext, cached. `none` reads the diagnosis cache (same recipe)."""
    if form == "none":
        cache = OUT / f"b8_diag_none_{level}.npz"
        if not cache.exists():
            raise SystemExit(f"{cache} missing -- run "
                             "run_b8_termination_diagnosis.py first")
        return np.load(cache)["phi_ext"]
    cache = OUT / f"b8_blend_{tag(form, frac)}_{level}.npz"
    if cache.exists():
        return np.load(cache)["phi_ext"]
    sb = (form, frac * B_SEMI)
    mesh, wls, cm, mvop = geom(level, sb)
    print(f"  solving span-blend {tag(form, frac)} {level} "
          f"(n_blended={mvop.n_span_blended}) ...", flush=True)
    t0 = time.time()
    r = _solve(mvop, mesh)
    dt = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, phi_ext=r["phi_ext"], converged=r["converged"],
             n_limited=r.get("n_limited", 0), n_floored=r.get("n_floored", 0),
             dt=dt)
    print(f"    done {dt:.0f}s  conv={r['converged']}", flush=True)
    return r["phi_ext"]


def solve_isolation(which):
    """Coarse solve with only ONE component of the blend active."""
    cache = OUT / f"b8_iso_{which}_coarse.npz"
    if cache.exists():
        return np.load(cache)["phi_ext"]
    mesh = read_mesh(MESH / "coarse.msh")
    a = np.radians(ALPHA)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    nonte, q = _nonte_and_q(cm)
    straddler = q > cm.span_length
    if which == "straddler_only":
        w = np.ones(len(nonte))
        w[straddler] = 0.0
    else:  # inboard_only: the rb0.08 smooth blend, straddlers left FREE
        w = tip_taper_factors(q, cm.span_length, "vanish_smooth",
                              0.08 * B_SEMI)
        w[straddler] = 1.0
    n = _apply_custom_w(mvop, cm, w)
    print(f"  solving isolation {which} coarse (n_blended={n}) ...",
          flush=True)
    r = _solve(mvop, mesh)
    np.savez(cache, phi_ext=r["phi_ext"], converged=r["converged"],
             n_blended=n)
    return r["phi_ext"]


def measure(level, phi_ext, span_blend=None):
    """Gate metrics for one solved state, on the HONEST metric."""
    mesh, wls, cm, mvop = geom(level, span_blend)
    mach_h = np.sqrt(mvop.element_mach2(phi_ext, M, mixed_plain="main"))
    cen = mesh.nodes[np.asarray(mesh.elements)].mean(axis=1)
    xte = np.array([x_te(np.clip(z, 0, B_SEMI)) for z in cen[:, 2]])
    zb, dx = cen[:, 2] / B_SEMI, cen[:, 0] - xte
    box = ((dx > 0.002) & (dx <= TIP_AFT) & (zb > 0.98)
           & (np.abs(cen[:, 1]) < TIP_Y))
    v_rel = mvop.op.V / np.median(mvop.op.V)

    def _peak(mask):
        idx = np.flatnonzero(mask)
        v = np.nan_to_num(mach_h, nan=-1.0)[idx]
        return float(v.max())

    aux_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
    xyz = mesh.nodes[aux_nodes]
    xte_n = np.array([x_te(np.clip(z, 0, B_SEMI)) for z in xyz[:, 2]])
    zb_n, dx_n = xyz[:, 2] / B_SEMI, xyz[:, 0] - xte_n
    sel = (dx_n > 0.002) & (dx_n <= TIP_AFT) & (zb_n > 0.98)
    dlt = np.abs(mvop.node_jump(phi_ext, aux_nodes[sel]))

    te_z = mesh.nodes[cm.te_nodes, 2]
    g = mvop.te_jump(phi_ext)
    o = np.argsort(te_z)
    z_s, g_s = te_z[o], g[o]
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    eta = z_s / B_SEMI
    return dict(
        e_h=_peak(box), e_ns=_peak(box & (v_rel >= V_REL_MIN)),
        d_max=float(dlt.max()),
        gamma_in=float(g_s[eta < 0.90].mean()),
        cl_kj=float(cl_kj_3d(g_s, z_s, s_ref, B_SEMI)),
        eta=eta, gamma=g_s, ntet=len(mesh.elements),
        n_blended=int(getattr(mvop, "n_span_blended", 0)),
    )


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("Track B / B8 re-spec: SPAN BLEND of the wake-LS "
                   "termination.  ANSWER: NEGATIVE -- it welds the ring "
                   "jump but any sheet-side pin is globally amplified")

    res = {}
    for form, frac in VARIANTS:
        sb = None if form == "none" else (form, frac * B_SEMI)
        for level in LEVELS:
            res[(form, frac, level)] = measure(
                level, solve_phi(form, frac, level), sb)
    iso = {w: measure("coarse", solve_isolation(w)) for w in ISOLATIONS}

    def expo(key, field):
        a, b = res[(key[0], key[1], "coarse")], res[(key[0], key[1], "medium")]
        ia, ib = float(a["ntet"]) ** (1 / 3), float(b["ntet"]) ** (1 / 3)
        return float(np.log(b[field] / a[field]) / np.log(ib / ia))

    base_c = res[("none", 0.0, "coarse")]
    base_m = res[("none", 0.0, "medium")]
    p_base_ns = expo(("none", 0.0), "e_ns")

    def dcl(key, level):
        base = base_c if level == "coarse" else base_m
        return 100.0 * (res[(key[0], key[1], level)]["cl_kj"]
                        / base["cl_kj"] - 1.0)

    # ---- identity + health -------------------------------------------------
    cl.add("B8r", "metric identity: baseline honest no-sliver reproduces "
           "the diagnosis",
           f"{base_c['e_ns']:.3f}->{base_m['e_ns']:.3f}, p={p_base_ns:+.3f}",
           f"peaks within 5e-3 of {DIAG_BASE['e_ns']}",
           (abs(base_c["e_ns"] - DIAG_BASE["e_ns"][0]) < 5e-3
            and abs(base_m["e_ns"] - DIAG_BASE["e_ns"][1]) < 5e-3),
           note="honest = element_mach2(mixed_plain='main'); the committed "
                "committee p=+1.341 was the x5 beyond-tip metric artifact "
                "on top of this +0.367")

    # ---- the blend DOES weld its target ------------------------------------
    for form, frac in VARIANTS[2:]:
        r = res[(form, frac, "medium")]
        cl.add("B8r", f"the blend WELDS the termination-ring jump "
               f"({tag(form, frac)})",
               f"max|delta| {base_m['d_max']:.5f} -> {r['d_max']:.5f} "
               f"({r['d_max'] / base_m['d_max']:.2f}x)",
               "< 0.1x baseline (medium)",
               r["d_max"] < 0.1 * base_m["d_max"],
               note="the h- and TE-taper-independent ring jump the TE row "
                    "blend could not reach -- the blend's own target, hit")

    # ---- FINDING 4: a narrow blend is actively harmful ---------------------
    r03 = res[("vanish_smooth", 0.03, "medium")]
    cl.add("B8r", "* FINDING: a ~2-element blend is ACTIVELY harmful "
           "(rb0.03)",
           f"ring jump {base_m['d_max']:.5f} -> {r03['d_max']:.5f} "
           f"({r03['d_max'] / base_m['d_max']:.2f}x), honest p "
           f"{p_base_ns:+.2f} -> {expo(('vanish_smooth', 0.03), 'e_ns'):+.2f}",
           "ring jump > 1.5x AND p worse than baseline",
           (r03["d_max"] > 1.5 * base_m["d_max"]
            and expo(("vanish_smooth", 0.03), "e_ns") > p_base_ns),
           note="a blend narrower than a few elements re-creates the "
                "Heaviside it was meant to remove, STEEPER")

    # ---- FINDING 1: the cost is global and huge ----------------------------
    d_m = [dcl(v, "medium") for v in VARIANTS[1:]]
    cl.add("B8r", "* FINDING: the lift cost is GLOBAL and ~10x the "
           "conforming taper's",
           "  ".join(f"{tag(*v)}:{d:+.1f}%"
                     for v, d in zip(VARIANTS[1:], d_m))
           + "  (conforming vanish_linear_rc0.05: -1.74%)",
           "all medium dcl < -15%", all(d < -15.0 for d in d_m),
           note="Gamma(z) scales down UNIFORMLY root-to-tip -- including "
                "where the blended rows are BITWISE identical to baseline "
                "(test-locked) -- so this is the global circulation MODE "
                "re-leveling, not local tip unloading")

    spread = max(d_m) - min(d_m)
    cl.add("B8r", "* FINDING: the cost is r_blend-INSENSITIVE "
           "(structural, not dose-dependent)",
           f"spread {spread:.1f} points over r_blend 0.03->0.08 "
           f"(n_blended 246->421)",
           "spread < 3 points", spread < 3.0)

    # ---- FINDING 2: component isolation ------------------------------------
    d_str = 100.0 * (iso["straddler_only"]["cl_kj"] / base_c["cl_kj"] - 1.0)
    d_inb = 100.0 * (iso["inboard_only"]["cl_kj"] / base_c["cl_kj"] - 1.0)
    cl.add("B8r", "* FINDING: BOTH components are amplified independently "
           "(coarse isolation)",
           f"straddler-only weld {d_str:+.1f}% | inboard-only blend "
           f"{d_inb:+.1f}%",
           "each alone < -8%", d_str < -8.0 and d_inb < -8.0,
           note="not a bug in either piece: ANY delta-pin on the sheet near "
                "the tip re-levels the whole mode. The implicit Kutta has "
                "no per-station target (Gamma is ONE solution mode), so the "
                "G13.2 finding-(2) fixed-point amplification acts GLOBALLY "
                "-- the conforming secant keeps it per-station, which is "
                "why the same F(z) costs -1.6% there and -20% here")

    # ---- FINDING 3: h-growing cost => the flat p is confounded -------------
    ia = float(base_c["ntet"]) ** (1 / 3)
    ib = float(base_m["ntet"]) ** (1 / 3)
    rows_conf = []
    for form, frac in VARIANTS[1:]:
        rc = 1.0 + dcl((form, frac), "coarse") / 100.0
        rm = 1.0 + dcl((form, frac), "medium") / 100.0
        p_unload = float(np.log(rm / rc) / np.log(ib / ia))
        rows_conf.append(((form, frac), p_unload))
    (f8, r8), p_unl8 = rows_conf[-1][0], rows_conf[-1][1]
    p_raw8 = expo(("vanish_smooth", 0.08), "e_ns")
    cl.add("B8r", "* FINDING: the loss GROWS under refinement => rb0.08's "
           "flat p is CONFOUNDED",
           f"rb0.08 dcl {dcl(('vanish_smooth', 0.08), 'coarse'):+.1f}% -> "
           f"{dcl(('vanish_smooth', 0.08), 'medium'):+.1f}%; unloading "
           f"contributes p_unload={p_unl8:+.2f} to the raw p={p_raw8:+.2f} "
           f"(baseline {p_base_ns:+.2f}; corrected ~{p_raw8 - p_unl8:+.2f})",
           "p_unload <= -0.08 (a material confounder)", p_unl8 <= -0.08,
           note="the h-growing unloading accounts for ~a THIRD of the "
                "apparent reduction; the corrected ~+0.15 hints at a real "
                "partial regularization but is NOT certifiable from a "
                "2-point ladder under a 20% global flow distortion -- and "
                "is moot at this cost")

    # ---- artifacts ---------------------------------------------------------
    rows = []
    for form, frac in VARIANTS:
        a, b = res[(form, frac, "coarse")], res[(form, frac, "medium")]
        rows.append((
            form, f"{frac:.2f}",
            f"{a['e_ns']:.4f}", f"{b['e_ns']:.4f}",
            f"{expo((form, frac), 'e_ns'):+.3f}",
            f"{expo((form, frac), 'e_h'):+.3f}",
            f"{a['d_max']:.6f}", f"{b['d_max']:.6f}",
            f"{dcl((form, frac), 'coarse'):+.2f}" if form != "none" else "+0.00",
            f"{dcl((form, frac), 'medium'):+.2f}" if form != "none" else "+0.00",
            f"{b['cl_kj']:.4f}", f"{b['gamma_in']:.6f}",
            b["n_blended"],
        ))
    for which, d_i in (("straddler_only", d_str), ("inboard_only", d_inb)):
        n_i = int(np.load(OUT / f"b8_iso_{which}_coarse.npz")["n_blended"])
        rows.append((f"iso_{which}", "", f"{iso[which]['e_ns']:.4f}",
                     "", "", "", f"{iso[which]['d_max']:.6f}", "",
                     f"{d_i:+.2f}", "", f"{iso[which]['cl_kj']:.4f}",
                     f"{iso[which]['gamma_in']:.6f}", n_i))
    write_csv(OUT, "b8_span_blend.csv",
              "form,r_blend_over_b,ens_coarse,ens_medium,p_honest_nosliver,"
              "p_honest,delta_max_coarse,delta_max_medium,dcl_pct_coarse,"
              "dcl_pct_medium,cl_KJ_medium,gamma_inboard_medium,"
              "n_blended_medium",
              rows)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.0))
    ib_ax = [float(res[("none", 0.0, lv)]["ntet"]) ** (1 / 3)
             for lv in LEVELS]
    for (form, frac), col in zip(VARIANTS,
                                 [CRITICAL, S3_YELLOW, S1_BLUE, S2_AQUA]):
        p = expo((form, frac), "e_ns")
        ax.plot(ib_ax, [res[(form, frac, "coarse")]["e_ns"],
                        res[(form, frac, "medium")]["e_ns"]],
                "o-", lw=2.2, ms=7, color=col,
                label=f"{tag(form, frac)}  (p={p:+.2f})")
    ax.set_xscale("log")
    ax.set_xlabel("mesh density  $N_{tet}^{1/3}$  (coarse $\\to$ medium)")
    ax.set_ylabel("HONEST tip-edge peak Mach (no-sliver, off-body box)")
    ax.set_title("B8 re-spec / span blend: the flat rb0.08 curve is "
                 "CONFOUNDED\n(the -20% h-growing global unloading "
                 "contributes $p_{unload}\\approx-0.2$)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    for (form, frac), col in zip(VARIANTS,
                                 [CRITICAL, S3_YELLOW, S1_BLUE, S2_AQUA]):
        r = res[(form, frac, "medium")]
        ax2.plot(r["eta"], r["gamma"], "-", lw=2.0, color=col,
                 label=f"{tag(form, frac)}  "
                       f"$\\Delta_{{max}}$={r['d_max']:.4f}")
    ax2.set_xlabel("$\\eta = z/b$")
    ax2.set_ylabel("$\\Gamma(z)$  (TE jump, medium)")
    ax2.set_title("THE finding: $\\Gamma(z)$ drops UNIFORMLY root-to-tip\n"
                  "(a tip-only constraint re-levels the whole implicit-"
                  "Kutta mode)")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    finish(fig, OUT, "b8_span_blend.png")

    return cl.report(OUT, "checks_b8_span_blend.csv")


if __name__ == "__main__":
    sys.exit(main())
