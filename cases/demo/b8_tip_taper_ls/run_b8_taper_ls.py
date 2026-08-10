"""
Track B / B8 — level-set tip-edge desingularization (row-blend tip taper).

THE QUESTION. P13/G13.2 killed the tip/wake-edge singularity on the CONFORMING
path with a spanwise loading taper Gamma_eff(z) = F(z) * Gamma_Kutta(z). The
level-set path cannot use it: it has NO Gamma DOF (the implicit Kutta makes Gamma
a solution mode) and its TE Kutta row

    s . (q_u - q_l) = 0                     (kernels/cut_assembly.py::te_kutta_coo)

is HOMOGENEOUS (RHS = 0), so scaling that row by F is an algebraic NO-OP
(G13.2 finding (8)). B8's model is a convex BLEND of that row with B2's
continuity weld (`te_weld_coo`), per TE node:

    F_i * [ s . (q_u - q_l) ]  +  (1 - F_i) * [ phi_aux - phi_main ]  = 0

  F = 1 (inboard)  -> pure pressure Kutta   (bit-identical to the untapered path)
  F = 0 (tip)      -> pure weld -> jump = 0 -> the tip is UNLOADED

It is not a no-op because the weld row is not proportional to the pressure row
(q_u and q_l come from DIFFERENT element sets, so q_u - q_l is not a jump
gradient). It is a DIFFERENT MODEL from Gamma = F*Gamma_Kutta, so r_c is
calibrated here independently, and the conforming comparison is a PHYSICS A/B,
not a model-identity check.

★★ THE RESULT (2026-07-13): THE ROW BLEND DOES NOT CLOSE B8, AND THE REASON IS
THAT ITS PREMISE IS WRONG. The blend is correctly implemented and behaves exactly
as its model predicts -- it converges cleanly (0 lim / 0 flr at every r_c), it
UNLOADS the tip circulation far past the conforming criterion (Gamma_last ~ h^q
with q = 4.73, criterion q >= 1), and it is perfectly LOCAL (inboard Gamma
+0.01%, cl_KJ +0.03%). And yet the tip-edge peak STILL diverges under every
taper (p = +1.34 untapered -> +1.35...+1.37 tapered; larger r_c is WORSE).

  * G13.2's DISCRETE mechanism does not transfer. There, p ~ 1 - q: the
    outermost TE station sheds its retained Gamma_last as a concentrated vortex
    over the last cell, so the edge velocity ~ Gamma_last/h, and q >= 1 kills
    it. Here q = 4.73 yet p = +1.37 -- nowhere near 1 - q = -3.73. Killing
    Gamma_last does NOT kill the peak.
  * The lift cost is ~0 (+0.03%, vs the conforming taper's -1.74%) because there
    is NOTHING TO UNLOAD: the level-set implicit Kutta ALREADY drives
    Gamma(tip) -> 0 emergently (B7 measured +-3e-4). The conforming path NEEDS
    the taper because its free-edge rule leaves Gamma_last ~ sqrt(h) (q = 0.44).
    The level-set path never had that disease.
  * WHERE the peak actually lives (the mechanism): the peak cell is OUTBOARD of
    the geometric tip (z/b = 1.012), it is a `beyond_tip` element -- one the
    SPANWISE CLIP refuses to cut (a crossing needs q <= span_length) -- it is
    the SAME element tapered or not, and it is NOT a small-cut sliver (volume
    0.71x the median, and not even a cut element). So the level-set tip
    singularity lives in how the embedded sheet TERMINATES, not in the
    circulation it sheds.

  ==> The two paths' tip singularities are DIFFERENT OBJECTS. G13.2 finding
      (8)'s "clean analogue" is a faithful analogue of the conforming MODEL, but
      the level-set path does not have the conforming path's disease, so the
      analogue treats a patient that is not ill. B8 needs a re-spec aimed at the
      sheet TERMINATION (the spanwise clip / beyond-tip zone), not at Gamma.
      The machinery shipped here (tip_taper on the LS path) is correct, tested,
      and bit-identical by default -- it is simply not the cure for this.

WHAT THIS DEMO MEASURES (roadmap Track B B8 gates)
  1. MECHANISM. Peak Mach on the wake sheet's free TIP EDGE, measured strictly
     OFF-BODY (dx > 0, z/b > 0.98 -- G13.2 finding (6)'s metric trap: a box that
     admits wing cells has its max migrate onto the wing suction peak once the
     edge is regularized, making a WORKING fix look worse). Refinement exponent
     p = dlog(edge peak)/dlog(1/h) over M6 coarse -> medium at M0.5.
       untapered  -> p > 0  (DIVERGES: the singularity is being resolved)
       tapered    -> p ~ 0  (BOUNDED: the singularity is gone)
  2. LOCALITY. Inboard Gamma(z) and cl_KJ must be essentially untouched -- a tip
     model that unloads the whole wing (the conforming `tanh_half` failure mode)
     is disqualified however well it regularizes.
  3. r_c CALIBRATION (independent of the conforming 0.05 b_semi).
  4. TWO-PATH PHYSICS A/B against the committed conforming numbers
     (cases/demo/p13_tip_edge_singularity/results/taper_probe.csv -- SAME meshes,
     SAME M0.5/alpha3.06 condition). Read, not recomputed (cost discipline).

Subsonic only (M0.5, upwind_c = 0 -> no limiter, no shock): the tip singularity
is a GEOMETRIC/discrete object, and probing it with the transonic machinery off
is what made G13.1 clean. The M0.84 P9-band question stays with P13/G13.3.

Heavy solves cache to results/b8_*.npz (gitignored). Standalone:
  NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
  python cases/demo/b8_tip_taper_ls/run_b8_taper_ls.py
"""
import csv
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
CONF_CSV = (REPO / "cases" / "demo" / "p13_tip_edge_singularity" / "results"
            / "taper_probe.csv")
ALPHA, M = 3.06, 0.5
LEVELS = ("coarse", "medium")
TIP_AFT, TIP_Y = 0.30, 0.03          # same tip-edge box as G13.1/G13.2

# (form, r_c as a FRACTION of B_SEMI). r_c = 0 is the untapered baseline.
# The r_c sweep is the B8 calibration: the blend is a DIFFERENT model, so the
# conforming 0.05 is a starting point, not an answer.
VARIANTS = [
    ("none", 0.00),
    ("vanish_smooth", 0.03),
    ("vanish_smooth", 0.05),
    ("vanish_smooth", 0.08),
    ("vanish_linear", 0.05),
]


def tag(form, frac):
    return form if form == "none" else f"{form}_rc{frac:.2f}"


def solve(form, frac, level):
    mesh = read_mesh(MESH / f"{level}.msh")
    a = np.radians(ALPHA)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)

    # F(z) on the LS path comes from the per-TE-NODE spanwise arclength cm.q
    # (there is no "station" abstraction here -- Gamma is a solution mode), with
    # the sheet's free edge at cm.span_length. Same tip_taper_factors as the
    # conforming path: only the model that CONSUMES F differs.
    r_c = frac * B_SEMI
    taper = (None if form == "none" else
             tip_taper_factors(cm.q[cm.te_nodes], cm.span_length, form, r_c))

    r = solve_multivalued_lifting(
        mvop, mesh, M, alpha_deg=ALPHA, farfield="neumann", upwind_c=0.0,
        n_outer_max=120, tol_residual=1e-7, tip_taper=taper)

    z = mesh.nodes[cm.te_nodes, 2]
    gamma = mvop.te_jump(r["phi_ext"])
    o = np.argsort(z)
    z, gamma = z[o], gamma[o]
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    cen = mesh.nodes[mesh.elements].mean(axis=1)
    return dict(
        cen=cen,
        # mixed_plain="side" pinned 2026-07-14 (default flipped to "main"):
        # this demo's committed verdict numbers (p=+1.341 untapered etc.)
        # were measured through the historical side reading; the honest
        # exponent lives in run_b8_termination_diagnosis.py.
        mach=np.sqrt(mvop.element_mach2(r["phi_ext"], M, mixed_plain="side")),
        xte=np.array([x_te(np.clip(zz, 0, B_SEMI)) for zz in cen[:, 2]]),
        gamma=gamma, station_z=z,
        taper=(np.ones(len(cm.te_nodes)) if taper is None else taper)[o],
        cl_kj=cl_kj_3d(gamma, z, s_ref, B_SEMI),
        n_st_in_taper=int(0 if taper is None
                          else np.count_nonzero(taper < 0.999)),
        ntet=len(mesh.elements), conv=bool(r["converged"]),
        n_lim=int(r.get("n_limited", 0)), n_flr=int(r.get("n_floored", 0)),
    )


def get(form, frac, level):
    cache = OUT / f"b8_{tag(form, frac)}_{level}.npz"
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}
    print(f"  solving LS {tag(form, frac)} {level} ...", flush=True)
    t0 = time.time()
    d = solve(form, frac, level)
    d["dt"] = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **d)
    print(f"    done {d['dt']:.0f}s  edge={_edge(d):.3f}  "
          f"cl_KJ={float(d['cl_kj']):.4f}  conv={d['conv']}  "
          f"lim/flr={d['n_lim']}/{d['n_flr']}", flush=True)
    return d


def _edge(d):
    """Peak Mach on the wake sheet's FREE TIP EDGE: strictly AFT of the TE
    (OFF-BODY, dx > 0) and at the tip (z/b > 0.98), in the chord plane. This is
    G13.2 finding (6)'s metric -- a box that admits wing cells measures the wing
    suction peak once the edge is fixed, not the singularity."""
    cen, xte = d["cen"], d["xte"]
    zb, dx = cen[:, 2] / B_SEMI, cen[:, 0] - xte
    m = ((dx > 0.002) & (dx <= TIP_AFT) & (zb > 0.98)
         & (np.abs(cen[:, 1]) < TIP_Y))
    return float(np.nan_to_num(d["mach"], nan=-1.0)[m].max())


def _gamma_last(d):
    """Circulation retained by the OUTERMOST TE node -- the discrete free edge
    sheds it as a concentrated vortex over the last cell (edge velocity
    ~ Gamma_last/h), so with Gamma_last ~ h^q the edge peak grows as h^(q-1):
    p ~ 1 - q, and the discrete regularization criterion is q >= 1."""
    return float(abs(np.asarray(d["gamma"])[-1]))


_CUTINFO = {}


def _cut_info(level):
    """(is_cut, is_beyond_tip, V/median) per element on `level` -- the cut map
    the level set built. Cached; no solve."""
    if level in _CUTINFO:
        return _CUTINFO[level]
    from pyfp3d.kernels.jacobian import PicardOperator
    mesh = read_mesh(MESH / f"{level}.msh")
    a = np.radians(ALPHA)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    n = len(mesh.elements)
    is_cut = np.zeros(n, dtype=bool)
    is_cut[cm.cut_elems] = True
    btip = np.zeros(n, dtype=bool)
    if len(getattr(cm, "beyond_tip_elems", [])):
        btip[cm.beyond_tip_elems] = True
    V = PicardOperator(mesh.nodes, mesh.elements).V
    _CUTINFO[level] = (is_cut, btip, V / np.median(V))
    return _CUTINFO[level]


def _peak_cell(d, level="medium"):
    """WHERE the tip-edge peak actually is -- the B8 mechanism evidence."""
    cen, xte = d["cen"], d["xte"]
    zb, dx = cen[:, 2] / B_SEMI, cen[:, 0] - xte
    m = ((dx > 0.002) & (dx <= TIP_AFT) & (zb > 0.98)
         & (np.abs(cen[:, 1]) < TIP_Y))
    idx = np.flatnonzero(m)
    j = int(idx[np.argmax(np.nan_to_num(d["mach"], nan=-1.0)[idx])])
    is_cut, btip, v_rel = _cut_info(level)
    return dict(elem=j, zb=float(zb[j]), dx=float(dx[j]),
                mach=float(d["mach"][j]), cut=bool(is_cut[j]),
                beyond_tip=bool(btip[j]), v_rel=float(v_rel[j]))


def _inboard(d, eta_max=0.90):
    """Mean Gamma inboard of the taper -- the LOCALITY probe."""
    eta = np.asarray(d["station_z"]) / B_SEMI
    return float(np.asarray(d["gamma"])[eta < eta_max].mean())


def _read_conforming():
    """The committed conforming taper numbers (P13/G13.2, SAME meshes and
    condition). Read, never recomputed -- cost discipline (CLAUDE.md)."""
    if not CONF_CSV.exists():
        return {}
    out = {}
    with open(CONF_CSV) as f:
        for row in csv.DictReader(f):
            out[(row["form"], float(row["r_c_over_b"]))] = {
                "p": float(row["p_edge"]), "q": float(row["q"]),
                "e_c": float(row["edge_coarse"]), "e_m": float(row["edge_medium"]),
                "cl": float(row["cl_KJ"]), "dcl": float(row["dcl_pct"]),
            }
    return out


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    cl = CheckList("Track B / B8: does the row-blend taper kill the tip "
                   "singularity on the LEVEL-SET path?  ANSWER: NO -- and why")

    data = {}
    for form, frac in VARIANTS:
        for lvl in LEVELS:
            data[(form, frac, lvl)] = get(form, frac, lvl)

    base_cl = float(data[("none", 0.0, "medium")]["cl_kj"])
    base_in = _inboard(data[("none", 0.0, "medium")])
    res = {}
    for form, frac in VARIANTS:
        a, b = data[(form, frac, "coarse")], data[(form, frac, "medium")]
        ia, ib = float(a["ntet"]) ** (1 / 3), float(b["ntet"]) ** (1 / 3)
        ea, eb = _edge(a), _edge(b)
        ga, gb = _gamma_last(a), _gamma_last(b)
        p = float(np.log(eb / ea) / np.log(ib / ia))
        q = (float(np.log(gb / ga) / np.log(ia / ib))
             if ga > 0 and gb > 0 else float("nan"))
        res[(form, frac)] = dict(
            p=p, q=q, e_c=ea, e_m=eb, g_c=ga, g_m=gb,
            cl=float(b["cl_kj"]),
            dcl=100.0 * (float(b["cl_kj"]) - base_cl) / base_cl,
            din=100.0 * (_inboard(b) - base_in) / abs(base_in),
            n_st=int(b["n_st_in_taper"]),
            conv=bool(a["conv"]) and bool(b["conv"]),
            clean=(int(a["n_lim"]) + int(a["n_flr"])
                   + int(b["n_lim"]) + int(b["n_flr"])) == 0,
        )

    # ---- what the run established (every check asserts a MEASURED fact) ---
    base = res[("none", 0.00)]
    ref = ("vanish_linear", 0.05)          # representative tapered variant
    tr = res[ref]

    cl.add("B8", "every LS solve converged, 0 limited / 0 floored",
           all(r["conv"] and r["clean"] for r in res.values()),
           "all True", all(r["conv"] and r["clean"] for r in res.values()),
           note="the row blend is numerically well behaved at every r_c")

    cl.add("B8", "UNTAPERED LS tip edge DIVERGES (off-body box)",
           f"{base['e_c']:.3f}->{base['e_m']:.3f}, p={base['p']:+.3f}",
           "p > 0.25", base["p"] > 0.25,
           note="reproduces G13.1's level-set exponent (1.34) on the same "
                "meshes -- the metric is measuring the right object")

    # -- the blend does exactly what its MODEL says --------------------------
    cl.add("B8", "the blend DOES unload the tip circulation (q >= 1)",
           f"q={tr['q']:.2f} at {tag(*ref)}", "q >= 1", tr["q"] >= 1.0,
           note="Gamma_last ~ h^q. The conforming criterion (G13.2) is q >= 1; "
                "the blend clears it by a wide margin")
    cl.add("B8", "the blend is LOCAL (inboard Gamma untouched)",
           f"{tr['din']:+.2f}% at {tag(*ref)}", "|change| < 2%",
           abs(tr["din"]) < 2.0)

    # -- ...and yet the edge still diverges. THE FINDING. --------------------
    still = all(r["p"] > 0.25 for k, r in res.items() if k[0] != "none")
    cl.add("B8", "* FINDING: the tip edge STILL diverges under every taper",
           "  ".join(f"{tag(*k)}:p={r['p']:+.2f}"
                     for k, r in res.items() if k[0] != "none"),
           "p > 0.25 for all (the taper does NOT bound it)", still,
           note="THE GATE'S MECHANISM CLAUSE IS NOT MET, and not because the "
                "blend misbehaves -- it clears its own q>=1 criterion. "
                "Killing Gamma_last does not kill the peak")

    cl.add("B8", "* FINDING: p ~ 1-q does NOT hold on the level-set path",
           f"q={tr['q']:.2f} => 1-q={1 - tr['q']:+.2f}, but p={tr['p']:+.2f}",
           "p differs from 1-q by > 1", abs(tr["p"] - (1 - tr["q"])) > 1.0,
           note="G13.2's DISCRETE mechanism (the outermost TE station sheds "
                "Gamma_last as a concentrated vortex over the last cell, so "
                "edge ~ Gamma_last/h and p = 1-q) is what the taper is built "
                "on. It does not transfer: the LS tip peak is NOT driven by "
                "the shed circulation")

    cl.add("B8", "* FINDING: lift cost ~0 because there is nothing to unload",
           f"level-set {tr['dcl']:+.2f}% vs conforming -1.74%",
           "|dcl| < 0.5% (vs conforming's -1.74%)", abs(tr["dcl"]) < 0.5,
           note="the LS implicit Kutta ALREADY drives Gamma(tip)->0 emergently "
                "(B7: +-3e-4). The conforming path NEEDS the taper because its "
                "free-edge rule leaves Gamma_last ~ sqrt(h) (q=0.44). The "
                "level-set path never had that disease")

    # -- where the peak actually lives: the mechanism ------------------------
    pk = _peak_cell(data[("none", 0.0, "medium")])
    pk_t = _peak_cell(data[(ref[0], ref[1], "medium")])
    cl.add("B8", "* MECHANISM: the peak cell is OUTBOARD of the geometric tip",
           f"z/b={pk['zb']:.4f}, dx={pk['dx']:+.4f}", "z/b > 1.0",
           pk["zb"] > 1.0,
           note="it is a `beyond_tip` element -- one the SPANWISE CLIP refuses "
                "to cut (cut_elements.py: crossings need q <= span_length), "
                "i.e. where the embedded sheet TERMINATES")
    cl.add("B8", "* MECHANISM: the taper does not even move the peak cell",
           f"untapered elem {pk['elem']} == tapered elem {pk_t['elem']}",
           "same element", pk["elem"] == pk_t["elem"])
    cl.add("B8", "* MECHANISM: it is NOT a small-cut sliver",
           f"V/median = {pk['v_rel']:.2f}, cut={pk['cut']}",
           "normal volume (V/median > 0.3), not a cut element",
           pk["v_rel"] > 0.3 and not pk["cut"],
           note="rules out the CutFEM small-cut instability -- the element is "
                "ordinary; what is singular is the sheet's TERMINATION")

    # ---- two-path PHYSICS A/B -------------------------------------------
    conf = _read_conforming()
    if conf:
        cb = conf.get(("none", 0.0))
        cl.add("B8", "A/B: conforming baseline also diverges (committed CSV)",
               f"conforming p={cb['p']:+.3f} | level-set p={base['p']:+.3f}",
               "both p > 0.25", cb["p"] > 0.25 and base["p"] > 0.25,
               note="same M6 meshes, same M0.5/alpha3.06 -- READ from "
                    "p13_tip_edge_singularity/results/taper_probe.csv, not "
                    "recomputed (cost discipline)")
        if ref in conf:
            cr = conf[ref]
            cl.add("B8",
                   f"* A/B VERDICT: the taper works on ONE path only ({tag(*ref)})",
                   f"conforming p={cr['p']:+.3f} (bounded) | "
                   f"level-set p={tr['p']:+.3f} (diverges)",
                   "conforming |p| < 0.25 AND level-set p > 0.25",
                   abs(cr["p"]) < 0.25 and tr["p"] > 0.25,
                   note="THE two-path A/B RESULT. The SAME F(z), the SAME "
                        "meshes, the SAME condition: the conforming taper "
                        "bounds its edge, the level-set row blend does not. "
                        "The two paths' tip singularities are DIFFERENT "
                        "OBJECTS, so finding (8)'s 'clean analogue' targets "
                        "the wrong one")

    # ---- artifacts -------------------------------------------------------
    rows = []
    for form, frac in VARIANTS:
        r = res[(form, frac)]
        c = conf.get((form, frac), {})
        pc = _peak_cell(data[(form, frac, "medium")])
        rows.append((
            form, f"{frac:.2f}", f"{r['e_c']:.3f}", f"{r['e_m']:.3f}",
            f"{r['p']:+.3f}", f"{r['g_c']:.5f}", f"{r['g_m']:.5f}",
            f"{r['q']:.2f}", f"{r['cl']:.4f}", f"{r['dcl']:+.2f}",
            f"{r['din']:+.2f}", r["n_st"], int(r["conv"] and r["clean"]),
            pc["elem"], f"{pc['zb']:.4f}", f"{pc['dx']:+.4f}",
            int(pc["cut"]), int(pc["beyond_tip"]), f"{pc['v_rel']:.2f}",
            f"{c.get('p', float('nan')):+.3f}" if c else "",
            f"{c.get('dcl', float('nan')):+.2f}" if c else "",
        ))
    write_csv(OUT, "b8_taper_ls.csv",
              "form,r_c_over_b,edge_coarse,edge_medium,p_edge,"
              "gamma_last_coarse,gamma_last_medium,q,cl_KJ,dcl_pct,"
              "d_inboard_pct,n_te_in_taper,clean,"
              "peak_elem,peak_z_over_b,peak_dx,peak_is_cut,peak_beyond_tip,"
              "peak_V_over_median,conforming_p,conforming_dcl",
              rows)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.0))
    ib = [float(data[(f, r, l)]["ntet"]) ** (1 / 3) for f, r in VARIANTS[:1]
          for l in LEVELS]
    for (form, frac), col in zip(
            VARIANTS, [CRITICAL, S3_YELLOW, S1_BLUE, S2_AQUA, INK_2]):
        r = res[(form, frac)]
        ax.plot(ib, [r["e_c"], r["e_m"]], "o-", color=col, lw=2.2, ms=7,
                label=f"{tag(form, frac)}  (p={r['p']:+.2f})")
    ax.set_xscale("log")
    ax.set_xlabel("mesh density  $N_{tet}^{1/3}$  (coarse $\\to$ medium)")
    ax.set_ylabel("peak Mach on the wake FREE TIP EDGE (off-body box)")
    ax.set_title("B8 / level-set: the row blend does NOT bound the tip edge\n"
                 "(every taper still diverges; M6, $M_\\infty$0.5, "
                 "$\\alpha$3.06 — no limiter, no shock)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    for (form, frac), col in zip(
            VARIANTS, [CRITICAL, S3_YELLOW, S1_BLUE, S2_AQUA, INK_2]):
        d = data[(form, frac, "medium")]
        ax2.plot(np.asarray(d["station_z"]) / B_SEMI, d["gamma"], "-",
                 color=col, lw=2.0, label=tag(form, frac))
    ax2.set_xlabel("$\\eta = z/b$")
    ax2.set_ylabel("$\\Gamma(z)$  (TE jump)")
    ax2.set_title("...because there is nothing to unload: the implicit Kutta\n"
                  "already drives $\\Gamma(tip)\\to 0$ (medium mesh)")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    finish(fig, OUT, "b8_taper_ls.png")

    return cl.report(OUT, "checks_b8.csv")


if __name__ == "__main__":
    sys.exit(main())
