"""GB31.3 candidate C1 -- outboard faded fringe: coarse-mesh mechanics probe.

Pre-reg cases/analysis/b31_tip_termination/PRE_REGISTRATION.md section 3
(GB31.3, candidate 1). This is the CHEAP mechanics rung, ahead of the
medium 0.7875 gate (run separately by the phase): three variants of the
B29 production LS config (B25 inboard fragment clip + B28 flat sheet) on
the wing-body COARSE mesh at M0.5 alpha = 3.06 -- the subsonic F1 regime
(|delta| ~ 0.026 at the last cut ring):

  baseline     production build, untouched
  fringe_only  clip extended w = 0.05*B_SEMI OUTBOARD, span_blend=None --
               isolates the clip extension alone (expected: the Heaviside
               is DISPLACED to the new edge, not removed)
  fringe_fade  same extension + span_blend=("vanish_smooth", w) -- the C1
               treatment: fade 1->0 across the fringe, pure-weld edge

All variants are seeded from scratch with the SAME n_seed (the fringe
changes the DOF count, so cross-variant warm starts are impossible) and
solved with the production strict recipe (wb30.LS_LEVEL_KW).

Measured per variant:
  * ring delta(q) profile across the termination on a SHARED bin grid
    covering the physical tip AND the fringe (wb31.ring_profile_ls
    machinery, adapted to an explicit bin range);
  * honest main-field tip-box peak Mach (wb30.mach2_ls; B8 x5 erratum
    discipline: main-dof readings only);
  * cl_p (wb30.cl_p_ls);
  * inboard (q <= 0.95*span) drift vs baseline: TE-jump Gamma(z) profile
    bins + cl_p (pre-reg guardrail <= 1%).

Decision hint (pre-locked for this probe): fringe_only merely displacing
the Heaviside AND fringe_fade producing a graded delta decay with inboard
drift <= 1% => C1 mechanics CONFIRMED, proceed to the medium gate.
Report either way, with numbers.

C1 RESULT (2026-07-21): graded decay achieved, inboard drift FAILS the
guardrail by 20-30x (cl -19.5%, Gamma -33%) -> C1 recorded X. Next
pre-registered candidate:

  c3_noclip   C3 sheet extension past the tip: outboard_fringe=np.inf
              DROPS the tip clip -- the flat y=0 sheet runs outboard to
              the far field ("no free edge in the fluid", the B25
              analogue; NOT the conforming free-edge model). Measured:
              convergence; honest tip peak M and WHERE the global peak
              sits; ring delta(q) to the extension's end; cl_p/cl_kj/
              gamma drift; a P5 far-field audit (aux DOF census, the
              pin_gamma scalar-gamma vs local-delta mismatch on the
              outboard outflow ring, max |R| over far-field rows);
              plus a MEDIUM classification sanity (no solve).

Run:  python cases/analysis/b31_tip_termination/run_g3.py
Artifacts: results/g3_probe_coarse.csv; solve caches results/g3_*.npz
"""

import csv
import time

import numpy as np

from wb31 import B_SEMI, OUT, z_of_q
import wb30  # noqa: E402  (importable after wb31's sys.path hook)
from pyfp3d.meshgen.fuselage import make_inboard_clip  # noqa: E402
from pyfp3d.meshgen.wingbody import te_polyline  # noqa: E402
from pyfp3d.solve.newton_ls import solve_multivalued_newton  # noqa: E402
from pyfp3d.wake import (CutElementMap, MultivaluedOperator,  # noqa: E402
                         WakeLevelSet)

M_PROBE = 0.5
W_FRAC = 0.05                  # fringe width w = 0.05 * B_SEMI (pre-reg)
FORM = "vanish_smooth"
N_SEED = 40                    # from-scratch Picard-LS seed (all variants)
TAIL_LO = 0.85                 # shared ring-profile grid: [TAIL_LO*span,
N_BINS = 14                    #   span + w], same edges for all variants
VARIANTS = ("baseline", "fringe_only", "fringe_fade", "c3_noclip")


# ------------------------------------------------------------------ builders
def build_variant(mesh, variant):
    """The wb30.build_ls_flat production config plus the C1/C3 knobs."""
    a = np.radians(wb30.ALPHA)
    wls = WakeLevelSet(te_polyline(wb30.FUS),
                       direction=(np.cos(a), np.sin(a), 0.0),
                       sheet_direction=(1.0, 0.0, 0.0))
    if variant == "baseline":
        fringe_w = 0.0
    elif variant == "c3_noclip":
        fringe_w = np.inf          # C3: the tip clip is dropped
    else:
        fringe_w = W_FRAC * B_SEMI
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]),
                       inboard_clip=make_inboard_clip(wb30.FUS),
                       outboard_fringe=fringe_w)
    span_blend = (FORM, fringe_w) if variant == "fringe_fade" else None
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls,
                               span_blend=span_blend)
    return wls, cm, mvop


def run_variant(mesh, mvop, variant):
    """One strict M0.5 level from scratch (production recipe), cached to
    results/g3_{variant}.npz (wb30 cache pattern)."""
    cache = OUT / f"g3_{variant}.npz"
    if cache.exists():
        return wb30._cache_load(cache)
    t0 = time.perf_counter()
    r = solve_multivalued_newton(
        mvop, mesh, M_PROBE, alpha_deg=wb30.ALPHA, phi_init=None,
        gamma_init=0.0, n_seed=N_SEED, upwind_c=1.5, verbose=True,
        **wb30.LS_LEVEL_KW)
    wall_s = time.perf_counter() - t0
    rec = dict(variant=variant, m_inf=float(M_PROBE),
               converged=bool(r["converged"]),
               accept_reason=str(r["accept_reason"]),
               residual_norm=float(r["residual_history"][-1]),
               n_newton=int(r["n_newton"]),
               n_limited=int(r["n_limited"]), n_floored=int(r["n_floored"]),
               mach_max=float(np.sqrt(r["mach2_max"])),
               gamma=float(r["gamma"]), cl_kj=float(r["cl_kj"]),
               wall_s=wall_s)
    wb30._cache_save(cache, rec, r["phi_ext"], r["gamma"])
    rec["phi"], rec["gamma"] = r["phi_ext"], r["gamma"]
    rec["cached"] = False
    print(f"  [g3 {variant}] converged={rec['converged']} "
          f"res={rec['residual_norm']:.2e} n_newton={rec['n_newton']} "
          f"({wall_s:.0f}s)", flush=True)
    return rec


# ------------------------------------------------------------------ metrics
def ring_profile(mesh, wls, cm, mvop, phi_ext, m_inf, q_lo, q_hi, n_bins):
    """wb31.ring_profile_ls with an explicit bin range (covers the fringe);
    per-element |delta| = mean |aux - main| over its aux-bearing nodes,
    M^2 from the honest main-field reading only."""
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    _s, _d, q = wls.evaluate(cents)
    aux_idx = cm.ext_dof_of_node
    has = aux_idx >= 0
    node_jump = np.zeros(mesh.nodes.shape[0])
    node_jump[has] = phi_ext[aux_idx[has]] - phi_ext[np.where(has)[0]]
    m2_main = wb30.mach2_ls(mvop, phi_ext, m_inf)
    is_cut = np.zeros(mesh.elements.shape[0], dtype=bool)
    is_cut[cm.cut_elems] = True
    node_has = has[mesh.elements]
    elem_jump = np.where(node_has,
                         np.abs(node_jump[mesh.elements]), 0.0).sum(axis=1)
    n_aux = node_has.sum(axis=1)
    elem_jump = np.divide(elem_jump, n_aux,
                          out=np.zeros_like(elem_jump), where=n_aux > 0)
    rows = []
    for lo, hi in zip(np.linspace(q_lo, q_hi, n_bins + 1)[:-1],
                      np.linspace(q_lo, q_hi, n_bins + 1)[1:]):
        sel = is_cut & (q >= lo) & (q < hi)
        if not np.count_nonzero(sel):
            continue
        rows.append(dict(section="ring", q_lo=float(lo), q_hi=float(hi),
                         q_mid=float(0.5 * (lo + hi)),
                         z_mid=float(z_of_q(0.5 * (lo + hi))),
                         n_elem=int(np.count_nonzero(sel)),
                         delta_mean=float(elem_jump[sel].mean()),
                         delta_max=float(elem_jump[sel].max()),
                         m2_main_max=float(m2_main[sel].max())))
    return rows


def te_gamma_profile(cm, phi_ext, n_bins=20):
    """TE-jump Gamma(q): phi_ext[aux] - phi_ext[main] at each TE node,
    binned over the PHYSICAL span [0, span_length] (fringe-invariant:
    TE nodes are untouched by the clip)."""
    te = cm.te_nodes
    jump = phi_ext[cm.ext_dof_of_node[te]] - phi_ext[te]
    q = cm.q[te]
    rows = []
    for lo, hi in zip(np.linspace(0.0, cm.span_length, n_bins + 1)[:-1],
                      np.linspace(0.0, cm.span_length, n_bins + 1)[1:]):
        sel = (q >= lo) & (q < hi)
        if not np.count_nonzero(sel):
            continue
        rows.append(dict(section="gamma", q_lo=float(lo), q_hi=float(hi),
                         q_mid=float(0.5 * (lo + hi)),
                         z_mid=float(z_of_q(0.5 * (lo + hi))),
                         n_elem=int(np.count_nonzero(sel)),
                         gamma=float(jump[sel].mean())))
    return rows


def inboard_drift(base_rows, var_rows, key, guard_q):
    """Max / mean |relative drift| vs baseline over the inboard bins
    (q_mid <= guard_q), matched by q_lo."""
    base = {r["q_lo"]: r[key] for r in base_rows
            if r["q_mid"] <= guard_q}
    drifts = []
    for r in var_rows:
        if r["q_mid"] > guard_q or r["q_lo"] not in base:
            continue
        b = base[r["q_lo"]]
        if abs(b) > 1e-14:
            drifts.append(abs(r[key] - b) / abs(b))
    if not drifts:
        return 0.0, 0.0
    return float(np.max(drifts)), float(np.mean(drifts))


# ------------------------------------------------------------------ driver
def fringe_mechanism_rows(mesh, states):
    """GB31.3/C1 mechanism census on the SOLVED states (cache-only): where
    the fringe sheet lives, what seeds its jump, and how the fade's effect
    distributes in (q, x).

    Measured on this coarse mesh (2026-07-21): the fringe band holds cut
    elements only out to x ~ 5 and NONE touch the far-field boundary --
    the band is a NEAR-FIELD STUB, so pin_gamma never sees it and the
    fringe jump is inherited from the physical sheet through the aux DOFs
    shared across the artificial q = span_length line (per-node
    granularity), not pinned at the outflow.
    """
    from pyfp3d.solve.picard_ls import farfield_aux_dofs

    rows = []
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    ej = {}
    for v, s in states.items():
        wls_cm = s["cm"]
        phi = s["rec"]["phi"]
        aux_idx = wls_cm.ext_dof_of_node
        has = aux_idx >= 0
        nh = has[mesh.elements]
        nj = np.zeros(mesh.nodes.shape[0])
        nj[has] = phi[aux_idx[has]] - phi[np.where(has)[0]]
        ejv = np.where(nh, np.abs(nj[mesh.elements]), 0.0).sum(axis=1)
        ejv = ejv / np.maximum(nh.sum(axis=1), 1)
        is_cut = np.zeros(mesh.elements.shape[0], dtype=bool)
        is_cut[wls_cm.cut_elems] = True
        ej[v] = (ejv, is_cut)
        gam = float(np.atleast_1d(s["rec"]["gamma"])[0])
        rows.append(dict(section="fringe_mech", variant=v,
                         metric="gamma_scalar", value=gam))
    cm = states["fringe_fade"]["cm"]
    span, w = cm.span_length, cm.outboard_fringe
    _s, _d, qc = states["fringe_fade"]["wls"].evaluate(cents)
    band = (qc > span) & (qc <= span + w)
    is_cut_f = ej["fringe_fade"][1]
    hosts, _auxs = farfield_aux_dofs(mesh, cm)
    fb = np.isin(mesh.elements,
                 np.unique(mesh.boundary_faces["farfield"])).any(axis=1)
    rows.append(dict(section="fringe_mech", variant="fringe_fade",
                     metric="ff_aux_total", value=len(hosts)))
    rows.append(dict(section="fringe_mech", variant="fringe_fade",
                     metric="ff_aux_in_fringe_band",
                     value=int(np.count_nonzero(
                         (cm.q[hosts] > span) & (cm.q[hosts] <= span + w)))))
    rows.append(dict(section="fringe_mech", variant="fringe_fade",
                     metric="fringe_band_cut_elems_touching_farfield",
                     value=int(np.count_nonzero(is_cut_f & band & fb))))
    for xlo, xhi in ((0.0, 2.0), (2.0, 5.0), (5.0, 30.0)):
        sel = is_cut_f & band & (cents[:, 0] >= xlo) & (cents[:, 0] < xhi)
        rows.append(dict(section="fringe_mech", variant="fringe_fade",
                         metric="fringe_band_cut_elems_by_x",
                         value=int(np.count_nonzero(sel)),
                         x_lo=xlo, x_hi=xhi, n_elem=int(np.count_nonzero(sel))))
        for v in ("fringe_only", "fringe_fade"):
            if np.count_nonzero(sel):
                rows.append(dict(section="fringe_mech", variant=v,
                                 metric="fringe_band_delta_mean",
                                 value=float(ej[v][0][sel].mean()),
                                 x_lo=xlo, x_hi=xhi,
                                 n_elem=int(np.count_nonzero(sel))))
    for v in ("fringe_only", "fringe_fade", "c3_noclip"):
        for xlo, xhi in ((0.0, 2.0), (2.0, 30.0)):
            inb = (qc >= 0.85 * span) & (qc < span) \
                & (cents[:, 0] >= xlo) & (cents[:, 0] < xhi)
            sb = ej["baseline"][1] & inb
            sv = ej[v][1] & inb
            if np.count_nonzero(sb) and np.count_nonzero(sv):
                ratio = float(ej[v][0][sv].mean() / ej["baseline"][0][sb].mean())
                rows.append(dict(section="fringe_mech", variant=v,
                                 metric="inboard_delta_ratio_vs_baseline",
                                 value=ratio, x_lo=xlo, x_hi=xhi,
                                 q_lo=0.85 * span, q_hi=span))
    for r in rows:
        if "value" in r and isinstance(r["value"], float):
            print(f"  [mech] {r['variant']:<12} {r['metric']:<38} "
                  f"x=[{r.get('x_lo', '')},{r.get('x_hi', '')}) "
                  f"{r['value']:.5f}")
    return rows


def c3_farfield_audit(mesh, states):
    """P5 far-field regression audit for the C3 no-clip extension.

    pin_gamma pins EVERY far-field-boundary aux DOF at host phi_inf -/+
    the SCALAR gamma (B17: the wake carries its circulation out). On the
    outboard extension the local sheet jump is NOT gamma (it is inherited
    near the tip and decays outboard), so the pin can impose a MISMATCHED
    jump on the outboard outflow ring -- the P5 branch-ray question in LS
    form. Per variant: pinned-DOF census (count / q range / outboard
    share), pinned jump vs the local element delta on the hosts' cut
    elements (mismatch), extension reach, and max |R| over the main rows
    of far-field-touching elements (LSNewtonSystem.residual at the
    converged state, live selection, recipe constants; free interior rows
    are <= tol by construction, the Dirichlet boundary rows are the
    reactions the boundary absorbs -- both reported).
    """
    from pyfp3d.solve.newton_ls import LSNewtonSystem
    from pyfp3d.solve.picard_ls import farfield_aux_dofs

    UP = wb30.UPWIND_DEFAULT
    rows = []
    ff_nodes = np.unique(mesh.boundary_faces["farfield"])
    ff_touch = np.isin(mesh.elements, ff_nodes).any(axis=1)
    ff_elem_nodes = np.unique(mesh.elements[ff_touch])
    ff_adj_free = ff_elem_nodes[~np.isin(ff_elem_nodes, ff_nodes)]
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    for v in ("baseline", "c3_noclip"):
        cm, mvop, rec = states[v]["cm"], states[v]["mvop"], states[v]["rec"]
        phi = rec["phi"]
        span = cm.span_length
        hosts, auxs = farfield_aux_dofs(mesh, cm)
        jump = phi[auxs] - phi[hosts]
        outb = cm.q[hosts] > span
        # local sheet delta on the hosts' own cut elements
        aux_idx = cm.ext_dof_of_node
        has = aux_idx >= 0
        nh = has[mesh.elements]
        nj = np.zeros(mesh.nodes.shape[0])
        nj[has] = phi[aux_idx[has]] - phi[np.where(has)[0]]
        ej = np.where(nh, np.abs(nj[mesh.elements]), 0.0).sum(axis=1) \
            / np.maximum(nh.sum(axis=1), 1)
        is_cut = np.zeros(mesh.elements.shape[0], dtype=bool)
        is_cut[cm.cut_elems] = True
        local = np.full(len(hosts), np.nan)
        for k, h in enumerate(hosts):
            elems = np.flatnonzero(is_cut
                                   & (mesh.elements == h).any(axis=1))
            if len(elems):
                local[k] = ej[elems].mean()
        mismatch = np.abs(np.abs(jump) - local)
        qc = states[v]["wls"].evaluate(cents)[2]
        sysm = LSNewtonSystem(mvop, M_PROBE, upwind_c=1.5,
                              m_crit=UP["m_crit"], gamma_air=1.4, u_inf=1.0,
                              m_cap=UP["m_cap"], rho_floor=UP["rho_floor"])
        _A, R, _up, _lo = sysm.residual(phi)
        R_main = np.abs(R[:cm.n_main])
        metrics = [
            ("n_ff_aux", len(hosts)),
            ("ff_aux_q_min", float(cm.q[hosts].min()) if len(hosts)
             else np.nan),
            ("ff_aux_q_max", float(cm.q[hosts].max()) if len(hosts)
             else np.nan),
            ("n_ff_aux_outboard", int(np.count_nonzero(outb))),
            ("pinned_jump_abs_mean", float(np.mean(np.abs(jump)))
             if len(hosts) else np.nan),
            ("pinned_jump_abs_mean_outboard",
             float(np.mean(np.abs(jump[outb]))) if np.any(outb) else np.nan),
            ("local_delta_on_pinned_hosts_mean",
             float(np.nanmean(local)) if len(hosts) else np.nan),
            ("local_delta_on_outboard_hosts_mean",
             float(np.nanmean(local[outb])) if np.any(outb) else np.nan),
            ("pin_vs_local_mismatch_max",
             float(np.nanmax(mismatch)) if len(hosts) else np.nan),
            ("maxR_ff_touching_main_rows",
             float(R_main[ff_elem_nodes].max())),
            ("maxR_ff_adjacent_free_rows",
             float(R_main[ff_adj_free].max()) if len(ff_adj_free)
             else np.nan),
            ("n_cut_beyond_span",
             int(np.count_nonzero(is_cut & (qc > span)))),
            ("n_cut_touching_farfield",
             int(np.count_nonzero(is_cut & ff_touch))),
            ("cut_centroid_q_max", float(qc[is_cut].max())),
            ("cut_centroid_z_max", float(cents[is_cut, 2].max())),
            ("cut_centroid_x_max", float(cents[is_cut, 0].max())),
        ]
        for metric, value in metrics:
            rows.append(dict(section="c3_audit", variant=v, metric=metric,
                             value=value))
            print(f"  [c3 audit] {v:<10} {metric:<38} {value:.5g}",
                  flush=True)
    return rows


def medium_classification():
    """GB31.3/C3 MEDIUM classification sanity (no solve): the production
    LS config on the wing-body medium mesh with the no-clip map vs the
    baseline build -- census, extension reach, TE/aux invariant."""
    rows = []
    mesh = wb30.load_mesh(wb30.LS_MESH_DIR / "medium.msh")
    ff_touch = np.isin(mesh.elements,
                       np.unique(mesh.boundary_faces["farfield"])).any(axis=1)
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    a = np.radians(wb30.ALPHA)
    for tag, fringe in (("baseline", 0.0), ("c3_noclip", np.inf)):
        wls = WakeLevelSet(te_polyline(wb30.FUS),
                           direction=(np.cos(a), np.sin(a), 0.0),
                           sheet_direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]),
                           inboard_clip=make_inboard_clip(wb30.FUS),
                           outboard_fringe=fringe)
        is_cut = np.zeros(mesh.elements.shape[0], dtype=bool)
        is_cut[cm.cut_elems] = True
        qc = wls.evaluate(cents)[2]
        invariant = bool((cm.ext_dof_of_node[cm.te_nodes] >= 0).all())
        metrics = [
            ("n_cut", len(cm.cut_elems)),
            ("n_ext", cm.n_ext_dofs),
            ("n_beyond_tip", len(cm.beyond_tip_elems)),
            ("n_te_nodes", len(cm.te_nodes)),
            ("n_cut_beyond_span",
             int(np.count_nonzero(is_cut & (qc > cm.span_length)))),
            ("n_cut_touching_farfield",
             int(np.count_nonzero(is_cut & ff_touch))),
            ("cut_centroid_q_max", float(qc[is_cut].max())),
            ("cut_centroid_z_max", float(cents[is_cut, 2].max())),
            ("cut_centroid_x_max", float(cents[is_cut, 0].max())),
            ("te_aux_invariant", int(invariant)),
        ]
        for metric, value in metrics:
            rows.append(dict(section="medium_class", variant=tag,
                             metric=metric, value=value))
        print(f"  [medium {tag}] cut={len(cm.cut_elems)} "
              f"ext={cm.n_ext_dofs} beyond={len(cm.beyond_tip_elems)} "
              f"te={len(cm.te_nodes)} invariant={invariant} "
              f"reach: z<={cents[is_cut, 2].max():.2f} "
              f"x<={cents[is_cut, 0].max():.2f} "
              f"ff={int(np.count_nonzero(is_cut & ff_touch))}", flush=True)
    return rows


def main():
    mesh = wb30.load_mesh(wb30.LS_MESH_DIR / "coarse.msh")
    w = W_FRAC * B_SEMI
    print(f"=== GB31.3/C1 coarse probe: M{M_PROBE} alpha={wb30.ALPHA}, "
          f"w = {W_FRAC}*B_SEMI = {w:.4f} ===", flush=True)
    states, rings, gammas, summaries = {}, {}, {}, []
    for variant in VARIANTS:
        print(f"--- variant: {variant} ---", flush=True)
        wls, cm, mvop = build_variant(mesh, variant)
        rec = run_variant(mesh, mvop, variant)
        phi = rec["phi"]
        cents = mesh.nodes[mesh.elements].mean(axis=1)
        m2 = wb30.mach2_ls(mvop, phi, M_PROBE)
        tip_box = cents[:, 2] > 0.95 * B_SEMI
        cl_p = wb30.cl_p_ls(mesh, mvop, phi, M_PROBE)
        states[variant] = dict(wls=wls, cm=cm, mvop=mvop, rec=rec, cl_p=cl_p,
                               mmax=float(np.sqrt(m2.max())),
                               tip_mmax=float(np.sqrt(m2[tip_box].max())))
        q_lo = TAIL_LO * cm.span_length
        if variant == "c3_noclip":
            # the no-clip sheet runs to the far field: extend the shared
            # grid to the extension's own reach (finer bins, same start)
            q_hi = float(wls.evaluate(cents)[2][cm.cut_elems].max())
            n_bins = 3 * N_BINS
        else:
            q_hi = cm.span_length + w
            n_bins = N_BINS
        rings[variant] = ring_profile(mesh, wls, cm, mvop, phi, M_PROBE,
                                      q_lo, q_hi, n_bins)
        gammas[variant] = te_gamma_profile(cm, phi)
        i_pk = int(np.argmax(m2))
        summaries.append(dict(
            section="variant", variant=variant,
            n_cut=len(cm.cut_elems), n_ext=cm.n_ext_dofs,
            n_beyond_tip=len(cm.beyond_tip_elems),
            n_span_blended=mvop.n_span_blended,
            converged=rec["converged"], res=rec["residual_norm"],
            n_newton=rec["n_newton"],
            n_limited=rec["n_limited"], n_floored=rec["n_floored"],
            cl_p=cl_p, cl_kj=rec["cl_kj"],
            mmax=states[variant]["mmax"],
            tip_mmax=states[variant]["tip_mmax"],
            peak_x=float(cents[i_pk, 0]), peak_y=float(cents[i_pk, 1]),
            peak_z=float(cents[i_pk, 2]),
            wall_s=round(rec["wall_s"], 1)))
        for r in rings[variant]:
            r["variant"] = variant
        for r in gammas[variant]:
            r["variant"] = variant

    # ---- drift vs baseline (pre-reg guardrail: inboard <= 1%)
    guard_q = 0.95 * states["baseline"]["cm"].span_length
    cl0 = states["baseline"]["cl_p"]
    for variant in ("fringe_only", "fringe_fade", "c3_noclip"):
        s = states[variant]
        cl_drift = abs(s["cl_p"] - cl0) / abs(cl0)
        g_max, g_mean = inboard_drift(gammas["baseline"], gammas[variant],
                                      "gamma", guard_q)
        d_max, d_mean = inboard_drift(rings["baseline"], rings[variant],
                                      "delta_mean", guard_q)
        for row in summaries:
            if row["variant"] == variant:
                row.update(cl_drift_pct=100 * cl_drift,
                           gamma_drift_max_pct=100 * g_max,
                           gamma_drift_mean_pct=100 * g_mean,
                           delta_drift_max_pct=100 * d_max,
                           delta_drift_mean_pct=100 * d_mean)

    mech = fringe_mechanism_rows(mesh, states)
    audit = c3_farfield_audit(mesh, states)
    med = medium_classification()

    # ---- CSV
    keys = ["section", "variant", "n_cut", "n_ext", "n_beyond_tip",
            "n_span_blended", "converged", "res", "n_newton", "n_limited",
            "n_floored", "cl_p", "cl_kj", "mmax", "tip_mmax", "wall_s",
            "peak_x", "peak_y", "peak_z",
            "cl_drift_pct", "gamma_drift_max_pct", "gamma_drift_mean_pct",
            "delta_drift_max_pct", "delta_drift_mean_pct",
            "q_lo", "q_hi", "q_mid", "z_mid", "n_elem", "delta_mean",
            "delta_max", "m2_main_max", "gamma", "metric", "value", "x_lo",
            "x_hi"]
    rows = summaries + [r for v in VARIANTS for r in rings[v]] \
        + [r for v in VARIANTS for r in gammas[v]] + mech + audit + med
    with open(OUT / "g3_probe_coarse.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, restval="")
        wr.writeheader()
        wr.writerows(rows)
    print(f"wrote {OUT / 'g3_probe_coarse.csv'}")

    # ---- table + verdict hint
    print("\n--- variant summary ---")
    hdr = (f"{'variant':<12} {'cut':>5} {'ext':>5} {'beyond':>6} "
           f"{'blend':>5} {'conv':>4} {'res':>9} {'cl_p':>7} "
           f"{'mmax':>5} {'tipM':>5} {'cl drift%':>9} {'Gam drift%':>10}")
    print(hdr)
    for row in summaries:
        print(f"{row['variant']:<12} {row['n_cut']:>5} {row['n_ext']:>5} "
              f"{row['n_beyond_tip']:>6} {row['n_span_blended']:>5} "
              f"{str(row['converged']):>4} {row['res']:>9.1e} "
              f"{row['cl_p']:>7.4f} {row['mmax']:>5.2f} "
              f"{row['tip_mmax']:>5.2f} "
              f"{row.get('cl_drift_pct', 0.0):>9.3f} "
              f"{row.get('gamma_drift_max_pct', 0.0):>10.3f}")
    print("\n--- ring delta(q) across the termination ---")
    print(f"{'variant':<12} {'q_mid':>6} {'z_mid':>6} {'n':>3} "
          f"{'d_mean':>7} {'d_max':>7} {'m2max':>6}")
    for v in VARIANTS:
        for r in rings[v]:
            print(f"{v:<12} {r['q_mid']:>6.3f} {r['z_mid']:>6.3f} "
                  f"{r['n_elem']:>3} {r['delta_mean']:>7.4f} "
                  f"{r['delta_max']:>7.4f} {r['m2_main_max']:>6.3f}")
    fade = next(r for r in summaries if r["variant"] == "fringe_fade")
    graded = rings["fringe_fade"][-1]["delta_mean"] \
        < 0.5 * rings["baseline"][-1]["delta_mean"]
    drift_ok = fade.get("gamma_drift_max_pct", 1e9) <= 1.0 \
        and fade.get("cl_drift_pct", 1e9) <= 1.0
    print("\n--- C1 mechanics hint (pre-reg GB31.3 probe) ---")
    print(f"  graded outboard decay: {graded} "
          f"(last fringe bin {rings['fringe_fade'][-1]['delta_mean']:.4f} "
          f"vs baseline last ring "
          f"{rings['baseline'][-1]['delta_mean']:.4f}); "
          f"inboard drift <= 1%: {drift_ok} "
          f"(cl {fade.get('cl_drift_pct', float('nan')):.3f}%, "
          f"Gamma {fade.get('gamma_drift_max_pct', float('nan')):.3f}%)")
    print(f"  -> C1 mechanics {'CONFIRMED' if graded and drift_ok else 'NOT CONFIRMED -- see numbers'}")

    c3 = next(r for r in summaries if r["variant"] == "c3_noclip")
    base = next(r for r in summaries if r["variant"] == "baseline")
    audit_of = {(r["variant"], r["metric"]): r["value"] for r in audit}
    print("\n--- C3 no-clip hint (pre-reg GB31.3 probe) ---")
    print(f"  converged={c3['converged']} res={c3['res']:.2e} "
          f"n_newton={c3['n_newton']} limited={c3['n_limited']} "
          f"floored={c3['n_floored']}")
    print(f"  inboard drift vs baseline: cl "
          f"{c3.get('cl_drift_pct', float('nan')):.3f}%, Gamma "
          f"{c3.get('gamma_drift_max_pct', float('nan')):.3f}% (guardrail 1%)")
    print(f"  tip-box peak M: baseline {base['tip_mmax']:.3f} -> "
          f"c3 {c3['tip_mmax']:.3f}; global peak sits at "
          f"(x,y,z)=({c3['peak_x']:.2f},{c3['peak_y']:.2f},{c3['peak_z']:.2f}) "
          f"vs baseline ({base['peak_x']:.2f},{base['peak_y']:.2f},"
          f"{base['peak_z']:.2f})")
    print(f"  extension: cut elems beyond span "
          f"{audit_of[('c3_noclip', 'n_cut_beyond_span')]:.0f}, touching "
          f"far field {audit_of[('c3_noclip', 'n_cut_touching_farfield')]:.0f},"
          f" reach z<= {audit_of[('c3_noclip', 'cut_centroid_z_max')]:.2f}")
    print(f"  P5 audit: pinned ff aux outboard "
          f"{audit_of[('c3_noclip', 'n_ff_aux_outboard')]:.0f}/"
          f"{audit_of[('c3_noclip', 'n_ff_aux')]:.0f}, pin-vs-local mismatch "
          f"max {audit_of[('c3_noclip', 'pin_vs_local_mismatch_max')]:.4f}, "
          f"max|R| ff rows "
          f"{audit_of[('c3_noclip', 'maxR_ff_touching_main_rows')]:.2e} "
          f"(baseline "
          f"{audit_of[('baseline', 'maxR_ff_touching_main_rows')]:.2e})")


if __name__ == "__main__":
    main()
