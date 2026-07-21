"""Shared machinery for B31 (C-class wing-tip cure) gates.

Reuses the B30 phase's builders/masks (sys.path hook below) and adds the
GB31.1 tip-atlas diagnostics:

  1. clamp-cell OWNERSHIP on both paths' dying levels -- each limited/
     floored element is classified as cap_wall / wing_wall / cut /
     beyond_tip / straddler (LS) or cap_wall / wing_wall / wake_adjacent /
     field (CONF), with distances to the sheet tip edge and the tip TE
     corner.  Answers: do the dying cells live on the round-cap WALL or on
     the sheet's TIP EDGE?  (B30 located them at z ~ 1.198, riding the
     B_SEMI=1.1963 line; G13.3 says the cap amplifies but does not create
     the site.)
  2. LS termination-ring delta(q) profile over the last 15% span (the
     Heaviside F1: |delta| ~ 0.026 at the last cut ring), honest
     main-field metric only (B8 x5 erratum discipline), plus the B20
     straddler side-vs-main M^2 pollution census.
  3. CONF tip Gamma(z) profile (per-station, wc.station_z) + Gamma_last
     (the G13.2 shed-vortex strength) with the vanish_smooth taper
     overlay F(z)*Gamma(z) -- the GB31.2 cure-leverage prediction.

Pre-reg: cases/analysis/b31_tip_termination/PRE_REGISTRATION.md (GB31.1).
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
B30 = REPO_ROOT / "cases/analysis/b30_transonic_ceiling"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(B30))

import wb30  # noqa: E402  (sets its own sys.path for B23/wb_common)
from pyfp3d.constraints.wake import tip_taper_factors  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.meshgen.wingbody import te_polyline  # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

# tip geometry (analytic spec, shared by both mesh families; verified
# numerically in the B31 kickoff probe: corner (1.1437, 0, 1.1963))
_POLY = te_polyline(wb30.FUS)
TIP_TE = _POLY[-1]                        # tip TE corner (x, 0, B_SEMI)
ROOT_TE = _POLY[0]                        # junction TE point
X_HAT = np.array([1.0, 0.0, 0.0])
# arclength q -> z along the straight 2-point production polyline
_DZ_DQ = (TIP_TE[2] - ROOT_TE[2]) / np.linalg.norm(TIP_TE - ROOT_TE)


def z_of_q(q):
    return ROOT_TE[2] + np.asarray(q, dtype=float) * _DZ_DQ


# ------------------------------------------------------------------ geometry
def tip_distances(cents):
    """(dist_to_sheet_tip_edge, dist_to_tip_TE_corner) per centroid.

    The sheet tip edge is the ray from TIP_TE downstream along +x at
    z = B_SEMI (both mesh families end the sheet there).
    """
    rel = cents - TIP_TE
    t = np.maximum(0.0, rel @ X_HAT)
    proj = TIP_TE + np.outer(t, X_HAT)
    d_edge = np.linalg.norm(cents - proj, axis=1)
    d_te = np.linalg.norm(rel, axis=1)
    return d_edge, d_te


# ------------------------------------------------------------------ ownership
def classify_ls(mesh, cm, cents):
    """Per-element ownership class (priority: cap_wall > wing_wall > cut
    > beyond_tip > straddler > field).  Straddler = uncut element sharing
    >= 1 node with a cut element (the B20 mixed-side-plain class)."""
    n_tets = mesh.elements.shape[0]
    wall_nodes = np.unique(mesh.boundary_faces["wall"])
    touches_wall = np.isin(mesh.elements, wall_nodes).any(axis=1)
    is_cut = np.zeros(n_tets, dtype=bool)
    is_cut[cm.cut_elems] = True
    is_beyond = np.zeros(n_tets, dtype=bool)
    is_beyond[cm.beyond_tip_elems] = True
    aux_node = np.where(cm.ext_dof_of_node >= 0)[0]
    shares_cut = np.isin(mesh.elements, aux_node).any(axis=1)
    cls = np.full(n_tets, "field", dtype=object)
    cls[shares_cut] = "straddler"
    cls[is_beyond] = "beyond_tip"
    cls[is_cut] = "cut"
    cls[touches_wall & (cents[:, 2] <= B_SEMI)] = "wing_wall"
    cls[touches_wall & (cents[:, 2] > B_SEMI)] = "cap_wall"
    return cls


def classify_conf(mc, wc, cents):
    """Per-element ownership class for the conforming path (cap_wall >
    wing_wall > wake_adjacent > field).  wake_adjacent = shares a node
    with the (duplicated) wake sheet."""
    wall_nodes = np.unique(mc.boundary_faces["wall"])
    touches_wall = np.isin(mc.elements, wall_nodes).any(axis=1)
    wake_nodes = np.unique(np.concatenate([wc.master_nodes, wc.slave_nodes]))
    touches_wake = np.isin(mc.elements, wake_nodes).any(axis=1)
    cls = np.full(mc.elements.shape[0], "field", dtype=object)
    cls[touches_wake] = "wake_adjacent"
    cls[touches_wall & (cents[:, 2] <= B_SEMI)] = "wing_wall"
    cls[touches_wall & (cents[:, 2] > B_SEMI)] = "cap_wall"
    return cls


def clamp_cell_rows(leg, mask_dict, cls, m2, cents):
    """One row per clamped cell with ownership + distances + local Mach."""
    d_edge, d_te = tip_distances(cents)
    rows = []
    for clamp, mask in mask_dict.items():
        for e in np.where(mask)[0]:
            rows.append(dict(section="clamp_cells", leg=leg, clamp=clamp,
                             elem=int(e), cls=cls[e],
                             x=float(cents[e, 0]), y=float(cents[e, 1]),
                             z=float(cents[e, 2]),
                             dist_tip_edge=float(d_edge[e]),
                             dist_tip_te=float(d_te[e]),
                             mach=float(np.sqrt(m2[e]))))
    return rows


def ownership_summary(rows):
    """(leg, clamp, cls) -> count pivot of clamp_cell_rows."""
    out = {}
    for r in rows:
        key = (r["leg"], r["clamp"], r["cls"])
        out[key] = out.get(key, 0) + 1
    return [dict(section="clamp_owner", leg=k[0], clamp=k[1], cls=k[2],
                 count=v) for k, v in sorted(out.items())]


# ------------------------------------------------------------------ LS ring δ(q)
def ring_profile_ls(mesh, wls, cm, mvop, phi_ext, cents, m_inf,
                    tail_frac=0.15, n_bins=12):
    """Jump-magnitude profile over the last `tail_frac` of the sheet span
    (cut elements only; per-element |delta| = mean |aux - main| over its
    aux-bearing nodes).  M^2 from the honest main-field reading only
    (B8 x5 erratum discipline)."""
    _s, _d, q = wls.evaluate(cents)
    aux_idx = cm.ext_dof_of_node                     # (n_nodes,) ext dof / -1
    has = aux_idx >= 0
    node_jump = np.zeros(mesh.nodes.shape[0])
    node_jump[has] = phi_ext[aux_idx[has]] - phi_ext[np.where(has)[0]]
    m2_main = wb30.mach2_ls(mvop, phi_ext, m_inf)
    n_tets = mesh.elements.shape[0]
    is_cut = np.zeros(n_tets, dtype=bool)
    is_cut[cm.cut_elems] = True
    node_has = has[mesh.elements]                    # (n_tets, 4)
    elem_jump = np.where(node_has,
                         np.abs(node_jump[mesh.elements]), 0.0).sum(axis=1)
    n_aux = node_has.sum(axis=1)
    elem_jump = np.divide(elem_jump, n_aux,
                          out=np.zeros_like(elem_jump), where=n_aux > 0)
    q0 = (1.0 - tail_frac) * cm.span_length
    bins = np.linspace(q0, cm.span_length, n_bins + 1)
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = is_cut & (q >= lo) & (q < hi)
        if not np.count_nonzero(sel):
            continue
        rows.append(dict(section="ring_profile", leg="LS",
                         q_lo=float(lo), q_hi=float(hi),
                         z_mid=float(z_of_q(0.5 * (lo + hi))),
                         n_elem=int(np.count_nonzero(sel)),
                         delta_mean=float(elem_jump[sel].mean()),
                         delta_max=float(elem_jump[sel].max()),
                         m2_main_max=float(m2_main[sel].max())))
    return rows


def straddler_census_ls(mesh, cm, mvop, phi_ext, cents, m_inf, cls):
    """B20-class pollution census in the tip box (z > 0.95 B_SEMI):
    side-reading M^2 (pre-B20 contaminated) vs honest main reading on
    straddler / beyond-tip elements."""
    sel = ((cls == "straddler") | (cls == "beyond_tip")) \
        & (cents[:, 2] > 0.95 * B_SEMI)
    in_upper, in_lower = mvop._side_element_sets()
    phi_up, phi_lo = mvop.side_potentials(phi_ext)
    m2_side = np.full(mesh.elements.shape[0], np.nan)
    for phi_s, keep in ((phi_up, in_upper), (phi_lo, in_lower)):
        _g, q2 = mvop.op.velocities(phi_s)
        m2 = mach_squared_field(q2, m_inf, 1.4)
        m2_side[keep] = m2[keep]
    m2_main = wb30.mach2_ls(mvop, phi_ext, m_inf)
    rows = []
    for klass in ("straddler", "beyond_tip"):
        k = sel & (cls == klass)
        if not np.count_nonzero(k):
            continue
        smax = float(np.nanmax(m2_side[k]))
        mmax = float(m2_main[k].max())
        rows.append(dict(section="straddler", leg="LS", cls=klass,
                         n_elem=int(np.count_nonzero(k)),
                         m2_side_max=smax, m2_main_max=mmax,
                         ratio=(smax / mmax if mmax > 0 else np.nan)))
    return rows


# ------------------------------------------------------------------ CONF Γ(z)
def gamma_profile_conf(wc, gamma, r_c_frac=0.05):
    """Per-station tip Gamma(z) + the vanish_smooth taper overlay
    (GB31.2 cure-leverage prediction; F3 parameters)."""
    z = np.asarray(wc.station_z, dtype=float)
    g = np.asarray(gamma, dtype=float).ravel()
    f = tip_taper_factors(z, B_SEMI, form="vanish_smooth",
                          r_c=r_c_frac * B_SEMI)
    rows = []
    for j in range(len(z)):
        rows.append(dict(section="gamma", leg="CONF", station=j,
                         z=float(z[j]), gamma=float(g[j]),
                         taper_f=float(f[j]), gamma_eff=float(f[j] * g[j])))
    return rows
