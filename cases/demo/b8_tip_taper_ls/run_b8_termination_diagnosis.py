"""
Track B / B8 re-spec STEP 0 -- termination-zone DIAGNOSIS (one-shot).

B8's negative verdict (the row blend does not bound the LS tip edge,
p = +1.34 -> +1.37) localized the peak to a `beyond_tip` element (93977,
z/b = 1.0118) and attributed the singularity to "how the embedded sheet
TERMINATES". The proposed re-spec (span-blend of the wake LS rows) treats
that attribution as established. This script tests it BEFORE any
implementation, because a code read found THREE distinct mechanisms all
consistent with the committed evidence -- and they need different cures:

  H2a (metric artifact). `element_mach2` reports mixed-side PLAIN elements
      (exactly the beyond_tip class) from the phi_up SIDE field, which
      substitutes AUX values at "-" cut nodes -- but the ASSEMBLY of a plain
      element uses MAIN dofs only (mass_conservation_coo scatters el[plain]).
      The diagnostic field is not the assembled field. This is the same
      disease class as B6's own_side_field te_lower fix; the mixed-side
      plain corner was never audited, and G13.1's LS exponent 1.34 used the
      same chain (p13 run_demo.py:126).
  H2b (density-weight pollution). mass_conservation_coo weights plain_plus
      elements with rho_upper, which on a MIXED-side plain element is
      computed from that same aux-mixed phi_up field -- so aux junk CAN
      enter the assembled operator there, subsonic path has no q2 cap, and
      density_field past vacuum yields NaN (which the demo's
      nan_to_num(nan=-1) silently hides).
  H1' (real jump survives the TE weld). The untapered tip circulation GROWS
      under refinement (Gamma_last 0.00011 -> 0.00218, q = -4.87 -- the
      "implicit Kutta pins Gamma(tip)~0" property degrades at medium), and
      the TE-row blend only welds the TE NODES: the jump delta at the
      DOWNSTREAM termination ring may never follow it. Then the peak really
      is delta_ring/h -- the G13.2 concentrated-vortex mechanism at the
      spanwise termination -- and the span-blend (welding ALL near-tip cut
      nodes) is exactly the cure.
  H3  (constraint-structure Heaviside, the re-spec's own hypothesis) is the
      remainder: honest field still diverges while delta_ring ~ 0.

WHAT IT MEASURES (M6 coarse+medium, M0.5/alpha3.06, farfield=neumann,
upwind_c=0 -- byte-identical recipe to run_b8_taper_ls.py; solves are cached
to results/b8_diag_*.npz WITH phi_ext this time):

  D1  How many nodes of the verdict peak element carry an aux DOF (0 kills
      H2a immediately).
  D2  The tip-edge box peak under THREE metrics, each with its own
      coarse->medium exponent p:
        (a) the committee metric `element_mach2` (must REPRODUCE the
            committed 0.672 -> 1.532 / p = +1.341 -- metric identity check);
        (b) HONEST: same, but mixed-side plain elements read the MAIN-field
            Mach (the field their assembly actually uses);
        (c) = (b) restricted to V/median >= 0.1 (the committed CSV shows the
            rc>=0.05 peak migrating to a V/median = 0.03 sliver, elem 79709).
  D3  |delta| (node_jump) over the aux nodes in the termination box,
      untapered vs tapered (vanish_linear rc0.05): does the TE weld actually
      drive the DOWNSTREAM ring jump to zero?
  D4  Hygiene: NaN count in the box Mach fields BEFORE nan_to_num; range of
      the assembled-side rho on mixed-side plain elements (H2b).

DECISION GATE (user arbitration -- pre-registered in the plan):
  - honest p bounded  & committee p diverges          => H2 (fix the metric,
    errata for G13.1-LS + B8; the span-blend may be unnecessary);
  - honest p diverges & delta_ring does NOT follow TE  => H1' (implement the
    span-blend: it welds exactly those nodes);
  - honest p diverges & delta_ring ~ 0 either way      => H3 (implement, but
    pre-register the possibility of another negative).

Standalone:
  NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
  python cases/demo/b8_tip_taper_ls/run_b8_termination_diagnosis.py
"""
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from cases.demo._common import write_csv  # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te  # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = Path(__file__).parent / "results"
MESH = REPO / "cases" / "meshes" / "onera_m6"
ALPHA, M = 3.06, 0.5
LEVELS = ("coarse", "medium")
TIP_AFT, TIP_Y = 0.30, 0.03          # same tip-edge box as G13.1/G13.2/B8

# untapered baseline + the representative tapered variant of the B8 verdict
VARIANTS = [("none", 0.00), ("vanish_linear", 0.05)]

# The committed B8 numbers this diagnosis must reproduce (metric identity).
COMMITTED = {("none", "coarse"): 0.672, ("none", "medium"): 1.532}
COMMITTED_PEAK_ELEM = 93977          # medium, untapered & vanish_linear


def tag(form, frac):
    return form if form == "none" else f"{form}_rc{frac:.2f}"


_GEOM = {}


def geom(level):
    """(mesh, wls, cm, mvop) for `level` -- deterministic rebuild, cached."""
    if level not in _GEOM:
        mesh = read_mesh(MESH / f"{level}.msh")
        a = np.radians(ALPHA)
        wls = WakeLevelSet(
            np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
            direction=(np.cos(a), np.sin(a), 0.0))
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]))
        mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
        _GEOM[level] = (mesh, wls, cm, mvop)
    return _GEOM[level]


def solve_phi(form, frac, level):
    """phi_ext for (form, frac, level); cached WITH phi_ext (the committed
    b8_*.npz caches do not store it, which is why this re-solve exists)."""
    cache = OUT / f"b8_diag_{tag(form, frac)}_{level}.npz"
    if cache.exists():
        return np.load(cache)["phi_ext"]
    mesh, wls, cm, mvop = geom(level)
    taper = (None if form == "none" else
             tip_taper_factors(cm.q[cm.te_nodes], cm.span_length, form,
                               frac * B_SEMI))
    print(f"  solving {tag(form, frac)} {level} ...", flush=True)
    t0 = time.time()
    r = solve_multivalued_lifting(
        mvop, mesh, M, alpha_deg=ALPHA, farfield="neumann", upwind_c=0.0,
        n_outer_max=120, tol_residual=1e-7, tip_taper=taper)
    dt = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(cache, phi_ext=r["phi_ext"], converged=r["converged"], dt=dt)
    print(f"    done {dt:.0f}s  conv={r['converged']}", flush=True)
    return r["phi_ext"]


def element_classes(level):
    """Masks: is_cut, is_tel, plain, mixed-side (per element), V/median."""
    mesh, wls, cm, mvop = geom(level)
    el = np.asarray(mesh.elements, dtype=np.int64)
    n = len(el)
    is_cut = np.zeros(n, dtype=bool)
    is_cut[cm.cut_elems] = True
    is_tel = np.zeros(n, dtype=bool)
    is_tel[cm.te_lower_elems] = True
    plain = ~is_cut & ~is_tel
    side_e = cm.node_side[el]
    mixed = side_e.min(axis=1) != side_e.max(axis=1)
    btip = np.zeros(n, dtype=bool)
    btip[cm.beyond_tip_elems] = True
    V = mvop.op.V
    return is_cut, is_tel, plain, mixed, btip, V / np.median(V)


def box_mask(level):
    """The G13.1/G13.2/B8 off-body tip-edge box, on element centroids."""
    mesh = geom(level)[0]
    cen = mesh.nodes[np.asarray(mesh.elements)].mean(axis=1)
    xte = np.array([x_te(np.clip(z, 0, B_SEMI)) for z in cen[:, 2]])
    zb, dx = cen[:, 2] / B_SEMI, cen[:, 0] - xte
    return ((dx > 0.002) & (dx <= TIP_AFT) & (zb > 0.98)
            & (np.abs(cen[:, 1]) < TIP_Y)), cen, dx, zb


def mach_fields(level, phi_ext):
    """(mach_committee, mach_honest, n_mixed_plain_changed):
    committee = element_mach2 (the committed metric);
    honest    = same, but mixed-side PLAIN elements read the MAIN field --
                the field their (single-valued, main-DOF) assembly uses."""
    mesh, wls, cm, mvop = geom(level)
    m2_c = mvop.element_mach2(phi_ext, M)
    _, q2_main = mvop.op.velocities(mvop.main_potential(phi_ext))
    m2_main = mach_squared_field(q2_main, M)
    is_cut, is_tel, plain, mixed, btip, v_rel = element_classes(level)
    fix = plain & mixed
    m2_h = m2_c.copy()
    m2_h[fix] = m2_main[fix]
    return np.sqrt(m2_c), np.sqrt(m2_h), np.sqrt(m2_main), int(fix.sum())


def peak(mach, mask):
    """(value, elem) of the box max, NaN treated as -1 (the demo's rule)."""
    idx = np.flatnonzero(mask)
    v = np.nan_to_num(mach, nan=-1.0)[idx]
    j = int(idx[np.argmax(v)])
    return float(v.max()), j


def delta_stats(level, phi_ext):
    """|delta| over the aux nodes in the termination box (node coords), plus
    the outermost-TE-node jump (Gamma_last, ties to CSV anomaly F2)."""
    mesh, wls, cm, mvop = geom(level)
    aux_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
    xyz = mesh.nodes[aux_nodes]
    xte = np.array([x_te(np.clip(z, 0, B_SEMI)) for z in xyz[:, 2]])
    zb, dx = xyz[:, 2] / B_SEMI, xyz[:, 0] - xte
    sel = (dx > 0.002) & (dx <= TIP_AFT) & (zb > 0.98)
    d = np.abs(mvop.node_jump(phi_ext, aux_nodes[sel]))
    te_z = mesh.nodes[cm.te_nodes, 2]
    g = mvop.te_jump(phi_ext)
    g_last = float(abs(g[np.argsort(te_z)][-1]))
    return d, aux_nodes[sel], g_last


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("B8 re-spec STEP 0: termination-zone diagnosis "
          "(H1' real-jump / H2 metric-artifact / H3 constraint-structure)")
    print("=" * 78)

    phi = {(f, r, l): solve_phi(f, r, l)
           for f, r in VARIANTS for l in LEVELS}

    rows = []
    report = {}
    for form, frac in VARIANTS:
        for level in LEVELS:
            mesh, wls, cm, mvop = geom(level)
            p_ext = phi[(form, frac, level)]
            mask, cen, dx, zb = box_mask(level)
            m_c, m_h, m_main, n_fix = mach_fields(level, p_ext)
            is_cut, is_tel, plain, mixed, btip, v_rel = element_classes(level)

            e_c, j_c = peak(m_c, mask)
            e_h, j_h = peak(m_h, mask)
            e_h2, j_h2 = peak(m_h, mask & (v_rel >= 0.1))

            # D4 hygiene: NaNs in the box before nan_to_num; assembled-side
            # rho on mixed plain elements (rho_upper is what
            # mass_conservation_coo consumes on plain_plus).
            nan_c = int(np.isnan(m_c[mask]).sum())
            nan_h = int(np.isnan(m_h[mask]).sum())
            rho_up, rho_lo = mvop.element_densities(p_ext, M)
            fixm = plain & mixed
            r_up = rho_up[fixm]
            rho_min = float(np.nanmin(r_up)) if fixm.any() else float("nan")
            rho_nan = int(np.isnan(r_up).sum())

            dlt, dlt_nodes, g_last = delta_stats(level, p_ext)

            key = (form, frac, level)
            report[key] = dict(
                e_c=e_c, j_c=j_c, e_h=e_h, j_h=j_h, e_h2=e_h2, j_h2=j_h2,
                nan_c=nan_c, nan_h=nan_h, rho_min=rho_min, rho_nan=rho_nan,
                d_max=float(dlt.max()) if dlt.size else 0.0,
                d_p95=float(np.percentile(dlt, 95)) if dlt.size else 0.0,
                n_dnodes=int(dlt.size), g_last=g_last, n_fix=n_fix,
                ntet=len(mesh.elements),
                m_c_at_committed=(float(np.nan_to_num(m_c, nan=-1.0)
                                        [COMMITTED_PEAK_ELEM])
                                  if level == "medium" else float("nan")),
                m_h_at_committed=(float(np.nan_to_num(m_h, nan=-1.0)
                                        [COMMITTED_PEAK_ELEM])
                                  if level == "medium" else float("nan")),
                m_main_at_committed=(float(np.nan_to_num(m_main, nan=-1.0)
                                           [COMMITTED_PEAK_ELEM])
                                     if level == "medium" else float("nan")),
            )
            rows.append((form, f"{frac:.2f}", level,
                         f"{e_c:.4f}", j_c, f"{e_h:.4f}", j_h,
                         f"{e_h2:.4f}", j_h2, nan_c, nan_h,
                         f"{rho_min:.4f}", rho_nan,
                         f"{report[key]['d_max']:.6f}",
                         f"{report[key]['d_p95']:.6f}",
                         report[key]["n_dnodes"], f"{g_last:.6f}", n_fix))

    # ---- D1: aux nodes of the committed verdict element (medium) ----------
    mesh, wls, cm, mvop = geom("medium")
    el = np.asarray(mesh.elements, dtype=np.int64)
    pk_nodes = el[COMMITTED_PEAK_ELEM]
    pk_aux = cm.ext_dof_of_node[pk_nodes] >= 0
    pk_sides = cm.node_side[pk_nodes]
    p_ext0 = phi[("none", 0.0, "medium")]
    pk_delta = mvop.node_jump(p_ext0, pk_nodes)
    p_ext_t = phi[("vanish_linear", 0.05, "medium")]
    pk_delta_t = mvop.node_jump(p_ext_t, pk_nodes)

    print()
    print(f"D1  peak element {COMMITTED_PEAK_ELEM} (medium): "
          f"aux nodes {int(pk_aux.sum())}/4, sides {pk_sides.tolist()}")
    print(f"    |delta| at its nodes  untapered: "
          f"{np.abs(pk_delta).max():.6f}   tapered: "
          f"{np.abs(pk_delta_t).max():.6f}")

    # ---- exponents ---------------------------------------------------------
    def expo(a, b, ka, kb):
        ia = float(report[ka]["ntet"]) ** (1 / 3)
        ib = float(report[kb]["ntet"]) ** (1 / 3)
        return float(np.log(b / a) / np.log(ib / ia))

    print("\nD2  box peak, three metrics (coarse -> medium, exponent p):")
    hdr = f"    {'variant':<22}{'metric':<28}{'coarse':>8}{'medium':>8}{'p':>8}"
    print(hdr)
    p_summary = {}
    for form, frac in VARIANTS:
        ka, kb = (form, frac, "coarse"), (form, frac, "medium")
        for label, field in (("(a) committee element_mach2", "e_c"),
                             ("(b) honest (main on mixed)", "e_h"),
                             ("(c) honest, no sliver", "e_h2")):
            a, b = report[ka][field], report[kb][field]
            p = expo(a, b, ka, kb)
            p_summary[(form, frac, field)] = p
            print(f"    {tag(form, frac):<22}{label:<28}"
                  f"{a:>8.3f}{b:>8.3f}{p:>+8.3f}")

    print("\n    metric identity check vs committed b8_taper_ls.csv:")
    for level in LEVELS:
        got = report[("none", 0.0, level)]["e_c"]
        want = COMMITTED[("none", level)]
        print(f"      none/{level}: committee {got:.3f} vs committed "
              f"{want:.3f}  ({'OK' if abs(got - want) < 5e-3 else 'MISMATCH'})")
    rm = report[("none", 0.0, "medium")]
    print(f"      elem {COMMITTED_PEAK_ELEM} (medium, untapered): "
          f"committee {rm['m_c_at_committed']:.3f} | honest "
          f"{rm['m_h_at_committed']:.3f} | pure-main "
          f"{rm['m_main_at_committed']:.3f}")

    print("\nD3  |delta| over termination-box aux nodes (does the TE weld "
          "reach the downstream ring?):")
    for form, frac in VARIANTS:
        for level in LEVELS:
            r = report[(form, frac, level)]
            print(f"    {tag(form, frac):<22}{level:<8}"
                  f"max {r['d_max']:.6f}  p95 {r['d_p95']:.6f}  "
                  f"(n={r['n_dnodes']})  Gamma_last {r['g_last']:.6f}")

    print("\nD4  hygiene:")
    for form, frac in VARIANTS:
        for level in LEVELS:
            r = report[(form, frac, level)]
            print(f"    {tag(form, frac):<22}{level:<8}"
                  f"box NaN committee/honest {r['nan_c']}/{r['nan_h']}  "
                  f"rho_up(mixed plain) min {r['rho_min']:.4f} "
                  f"NaN {r['rho_nan']}  (n mixed plain = {r['n_fix']})")

    write_csv(OUT, "b8_termination_diagnosis.csv",
              "form,r_c_over_b,level,peak_committee,elem_committee,"
              "peak_honest,elem_honest,peak_honest_nosliver,elem_nosliver,"
              "box_nan_committee,box_nan_honest,rho_up_min_mixed_plain,"
              "rho_up_nan_mixed_plain,delta_max_box,delta_p95_box,"
              "n_aux_nodes_box,gamma_last,n_mixed_plain",
              rows)

    # ---- the decision-gate readout ----------------------------------------
    pc = p_summary[("none", 0.0, "e_c")]
    ph = p_summary[("none", 0.0, "e_h")]
    ph2 = p_summary[("none", 0.0, "e_h2")]
    print("\n" + "=" * 78)
    print(f"GATE  committee p = {pc:+.3f} | honest p = {ph:+.3f} | "
          f"honest-no-sliver p = {ph2:+.3f}")
    if ph < 0.25 <= pc:
        print("GATE  => H2 branch: the divergence is a METRIC artifact "
              "(fix element_mach2 for mixed-side plain elements; errata for "
              "G13.1-LS and B8; span-blend may be unnecessary)")
    else:
        d0 = report[("none", 0.0, "medium")]["d_max"]
        dt_ = report[("vanish_linear", 0.05, "medium")]["d_max"]
        if dt_ > 0.5 * d0 and d0 > 0:
            print("GATE  => H1' branch: honest field still diverges AND the "
                  "downstream ring jump survives the TE weld "
                  f"(delta_max {d0:.5f} -> {dt_:.5f}) => implement the "
                  "span-blend (it welds exactly those nodes)")
        else:
            print("GATE  => H3 branch: honest field diverges but the ring "
                  "jump is already ~0 -- the re-spec's constraint-structure "
                  "hypothesis stands; implement the span-blend but "
                  "pre-register the possibility of another negative")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
