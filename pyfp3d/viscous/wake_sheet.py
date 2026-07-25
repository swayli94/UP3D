"""Track V V6 -- the conforming wake-sheet delta* source (binding:
docs/roadmap/track_v.md GV6.1; pre-registration
cases/analysis/v6_1_wake_sheet/PRE_REGISTRATION.md including the
2026-07-25 addendum; GV6.0 ruling: Option A + producer (i),
cases/analysis/v6_0_design_adjudication/DESIGN_ADJUDICATION.md).

The wake displacement thickness enters the FP solve as a mass-transpiration
sheet source on the wake cut faces, riding the V2 channel: b_wake is
assembled in CUT-MESH volume numbering and added to the loose loop's
``body_source_rhs`` (``viscous/coupling.py::run_loose_coupling``), so the
``constraints/wake.py`` reduce RHS (T^T b) folds the two sides' Galerkin
loads into the master rows -- that fold IS the weak sheet flux (GV6.0
survey section 1.1). No solver-path code is touched: the only touch point
is an additive, flag-gated term in the loose loop's RHS assembly
(``CouplingConfig.wake_transpiration``, default OFF = legacy
bit-identical).

Pieces here (all plain NumPy orchestration like coupling.py -- nothing is
hot next to the kernels):

1. ``build_wake_sheet_case`` -- the wake ``SurfaceMesh`` (V1 layout reused
   as a second instance, design_track_v.md section 6 point 1) on
   ``wc.wake_faces_minus`` plus an open-chain station table (TE = min-x
   endpoint -> downstream end) with the per-node arc length s from the TE,
   and the master->slave fold pairing. The W3 sanity asserts of the
   pre-registration run at build time.

2. ``prescribed_delta_star_wake`` / ``wake_transpiration_source`` -- the
   pinned producer (i): TE confluence by construction (delta*_TE =
   delta*_upper(TE) + delta*_lower(TE), theta_TE likewise, read off the
   wall IBL closure packet at the two TE copies split by ``side_node``),
   then the straight-wake mass-transpiration relaxation

       delta*_wake(s) = theta_TE + (delta*_TE - theta_TE) * exp(-s / L_rel),
       L_rel = WAKE_L_REL_CHORDS * chord  (pinned MODEL CHOICE, recorded;
       the sensitivity is a GV6.2 item, not tuned in this gate)

   with the W2 runtime assert: the constructed delta*_wake(0) equals
   delta*_TE to 1e-12 relative, checked at EVERY call (every loose-loop
   outer iteration). With no wall shear the wake momentum thickness is
   approximately conserved (far-field drag theorem), so H_wake = delta*/theta
   relaxes from its TE value toward 1 downstream.

3. ``wake_edge_velocity`` -- u_e,wake: the per-zone quadratic recovery
   (post/surface.py discipline, same function the wall chain uses) run on
   the wake cut faces, averaged over the two sides:
   u_e(node) = 0.5*(ue_vol[minus] + ue_vol[plus slave]).

4. ``assemble_wake_sheet_rhs`` -- the sheet-source RHS: m_wake =
   div_Gamma(rho_e u_e delta*_wake) (``transpiration_from_delta_star``
   verbatim, SurfaceMesh-generic), scattered HALF PER FACE COPY
   (``0.5*m_wake`` on the minus nodes and on their plus slaves), then
   ``assemble_transpiration_rhs`` over wake_minus U wake_plus.

   WHY THE HALF (the 2026-07-25 pre-registration addendum): the T^T fold
   SUMS the two coincident copies' loads into the master row, so two
   copies at the full m_wake would realize [rho v_n] = 2*m_wake. m_wake is
   the divergence of the TOTAL wake defect mass flux = the sheet's total
   ejection = the jump strength itself (the source-sheet identity
   v_n± = ±m_wake/(2 rho)); the half per copy makes the folded jump equal
   m_wake. The sign convention is the V2 one (positive m = blowing; the
   negated Galerkin load, transpiration.py:19-27,114): a thickening wake
   (m_wake > 0) ejects fluid AWAY from the sheet on both sides. Band (b)'s
   pinned MMS lock (uniform m0: antisymmetry + jump = m0/rho0 within 5%)
   empirically pins both the sign and the factor.

delta*_wake == 0 gives m_wake == 0 and hence the EXACT zero RHS vector
(GV2.1(b) assembly discipline), which is the bit-identity basis of band
(a)(i).
"""

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np

from pyfp3d.physics.isentropic import density_field
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.surface_mesh import SurfaceMesh
from pyfp3d.viscous.transpiration import (
    assemble_transpiration_rhs,
    edge_velocity_per_zone,
    transpiration_from_delta_star,
)

__all__ = [
    "WAKE_L_REL_CHORDS",
    "WakeSheetCase",
    "build_wake_sheet_case",
    "prescribed_delta_star_wake",
    "wake_edge_velocity",
    "wake_transpiration_source",
    "assemble_wake_sheet_rhs",
]


# The single free constant of the prescribed producer, pinned a priori at
# 1.0 chords (PRE_REGISTRATION section 2 -- MODEL CHOICE, RECORDED; the
# L_rel sensitivity sweep is a GV6.2 recorded item, NOT tuned in GV6.1).
WAKE_L_REL_CHORDS = 1.0

# W2/band-(c) tolerance: the TE-continuity construction identity holds to
# round-off; 1e-12 relative is the pinned assert (recipe-error raiser).
_TE_CONTINUITY_RTOL = 1.0e-12


@dataclass
class WakeSheetCase:
    """The wake sheet's SurfaceMesh + station/fold tables (frozen at build).

    sm: SurfaceMesh built on wc.wake_faces_minus (compact wake-surface
        numbering; ``volume_node_of`` maps to the cut mesh's minus-side
        node ids).
    slave_of_surf (n_s,) int64: the plus-side SLAVE volume node paired to
        each wake-surface node (the fold pairing; slave coordinates
        coincide with their master's).
    s (n_s,) float64: arc distance from the TE station along the straight
        wake strip (0 at the TE, downstream positive).
    station_of (n_s,) int64: wake-surface node -> station row (stations
        group the strip's spanwise copies by (x, y), the airfoil-case
        discipline).
    faces_both (2F, 3) int: wake_minus U wake_plus in cut-mesh volume
        numbering (the RHS assembly face set; row f of the minus block and
        row f of the plus block are the two coincident copies of sheet
        face f).
    prescribed_ds: optional (n_s,) prescribed delta*_wake override (the
        band-(a)(i) zero-field hook; when set, the producer formula and
        the W2 construction assert are bypassed -- an override is by
        definition not the TE-continuity construction).
    """

    sm: SurfaceMesh
    slave_of_surf: np.ndarray
    s: np.ndarray
    station_of: np.ndarray
    faces_both: np.ndarray
    prescribed_ds: Optional[np.ndarray] = None


def _station_chain_open(sm: SurfaceMesh) -> Tuple[np.ndarray, np.ndarray]:
    """Open-chain station table on the wake strip: stations = distinct
    (x, y) node positions (the strip's spanwise copies collapse), chained
    from the TE endpoint (min-x) to the downstream end.

    Returns (station_of, s_station): station row per surface node and the
    per-station arc distance from the TE. Raises (W3) unless the station
    graph is a single open chain covering every station with strictly
    monotone arc length.
    """
    xy = np.ascontiguousarray(sm.xyz[:, :2])
    uniq, station_of = np.unique(xy, axis=0, return_inverse=True)
    n_st = len(uniq)
    adj = [set() for _ in range(n_st)]
    for t in range(sm.n_tri):
        rows = np.unique(station_of[sm.triangles[t]])
        for a in rows:
            for b in rows:
                if a != b:
                    adj[int(a)].add(int(b))
                    adj[int(b)].add(int(a))
    ends = [i for i in range(n_st) if len(adj[i]) == 1]
    te = int(np.argmin(uniq[:, 0]))  # TE = the strip's min-x station
    if len(ends) != 2 or te not in ends:
        raise RuntimeError(
            "wake sheet station graph is not the expected open strip "
            f"({len(ends)} endpoint stations, TE endpoint={te in ends})"
        )

    order = [te]
    prev = -1
    while True:
        cur = order[-1]
        nxt = [j for j in adj[cur] if j != prev]
        if not nxt:
            break
        if len(nxt) != 1:
            raise RuntimeError(
                "wake sheet station chain branches (not a quasi-2D strip)"
            )
        prev = cur
        order.append(nxt[0])
    if len(order) != n_st:
        raise RuntimeError(
            f"wake sheet station chain covers {len(order)} of {n_st} "
            "stations (W3: station count does not match the wake strip)"
        )

    s_st = np.zeros(n_st, dtype=np.float64)
    for k in range(1, n_st):
        s_st[order[k]] = s_st[order[k - 1]] + float(
            np.linalg.norm(uniq[order[k]] - uniq[order[k - 1]])
        )
    if not np.all(np.diff(s_st[order]) > 0.0):
        raise RuntimeError(
            "wake sheet arc length is not strictly monotone from the TE (W3)"
        )
    return station_of.astype(np.int64), s_st


def build_wake_sheet_case(mesh_cut, wc) -> WakeSheetCase:
    """Build the wake-sheet case on a wake-cut mesh (W3 sanity at build).

    Args:
        mesh_cut: the cut Mesh from mesh/wake_cut.py::cut_wake (its
            boundary_faces carry wake_minus / wake_plus).
        wc: the matching WakeCut (master/slave duplication map).

    Returns:
        WakeSheetCase. Raises RuntimeError on any W3 violation (recipe
        error): the station chain must cover the strip with strictly
        monotone arc length, node areas must be positive, every
        wake-minus node must pair to exactly one plus-slave at coincident
        coordinates, and volume_node_of must map into the cut mesh.
    """
    nodes = np.asarray(mesh_cut.nodes, dtype=np.float64)
    sm = SurfaceMesh.from_wall_faces(
        nodes, wc.wake_faces_minus, mesh_cut.elements, name="wake"
    )
    station_of, s_st = _station_chain_open(sm)
    s = s_st[station_of]

    # The fold pairing: every wake-minus node is a master with exactly one
    # slave (on quasi-2D sheets every wake node is duplicated; M1 free-edge
    # nodes are single-valued and OUT of scope for GV6.1 -- they would
    # fail this assert honestly).
    slave_of_master = dict(
        zip(wc.master_nodes.tolist(), wc.slave_nodes.tolist())
    )
    if len(slave_of_master) != len(wc.master_nodes):
        raise RuntimeError("WakeCut master/slave map is not one-to-one (W3)")
    vol = sm.volume_node_of
    missing = [int(v) for v in vol if int(v) not in slave_of_master]
    if missing:
        raise RuntimeError(
            f"{len(missing)} wake-minus node(s) have no plus-slave "
            "(free-edge/single-valued sheet nodes are out of GV6.1 scope; W3)"
        )
    slave_of_surf = np.array(
        [slave_of_master[int(v)] for v in vol], dtype=np.int64
    )

    if vol.size == 0 or vol.max() >= len(nodes) or slave_of_surf.max() >= len(nodes):
        raise RuntimeError(
            "wake volume_node_of / slave ids do not map into the cut "
            "mesh's node table (W3)"
        )
    if np.any(sm.node_area <= 0.0):
        raise RuntimeError("wake SurfaceMesh has a non-positive node area (W3)")
    if not np.allclose(nodes[slave_of_surf], nodes[vol], rtol=0.0, atol=0.0):
        raise RuntimeError(
            "wake fold pairing is not coordinate-coincident (W3: a slave "
            "must be its master's duplicated copy)"
        )

    faces_both = np.concatenate(
        [np.asarray(wc.wake_faces_minus), np.asarray(wc.wake_faces_plus)],
        axis=0,
    )
    return WakeSheetCase(
        sm=sm,
        slave_of_surf=slave_of_surf,
        s=s,
        station_of=station_of,
        faces_both=faces_both,
    )


def prescribed_delta_star_wake(
    ds_te: float, th_te: float, s: np.ndarray, l_rel: float
) -> np.ndarray:
    """The pinned producer (i): theta_TE + (delta*_TE - theta_TE) exp(-s/L_rel).

    delta*_wake(0) = delta*_TE by construction (the W2 assert) and
    delta*_wake -> theta_TE downstream (H_wake -> 1: with no wall shear the
    wake momentum thickness is approximately conserved, so the profile
    fills out).
    """
    s = np.asarray(s, dtype=np.float64)
    if l_rel <= 0.0:
        raise ValueError(f"L_rel must be > 0, got {l_rel}")
    return th_te + (ds_te - th_te) * np.exp(-s / l_rel)


def wake_edge_velocity(
    nodes: np.ndarray, wake_case: WakeSheetCase, phi: np.ndarray
) -> np.ndarray:
    """u_e on the wake sheet, (n_s, 3): the loose loop's per-zone recovery
    discipline (quadratic; no LE band on the wake) run on the wake cut
    faces, averaged over the two sides --
    u_e(node) = 0.5*(ue_vol[minus] + ue_vol[plus slave])."""
    ue_vol = edge_velocity_per_zone(
        np.asarray(nodes, dtype=np.float64), wake_case.faces_both, phi
    )
    ue = 0.5 * (
        ue_vol[wake_case.sm.volume_node_of] + ue_vol[wake_case.slave_of_surf]
    )
    if not np.all(np.isfinite(ue)):
        raise RuntimeError(
            "u_e gather hit NaN on the wake sheet (off-sheet leak?)"
        )
    return ue


def wake_transpiration_source(
    wake_case: WakeSheetCase,
    stations,
    outs: np.ndarray,
    ue_surf: np.ndarray,
    m_inf: float,
    gamma_air: float = 1.4,
    chord: float = 1.0,
    l_rel_chords: float = WAKE_L_REL_CHORDS,
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """The prescribed delta*_wake producer + the m_wake field, from the
    current wall-IBL state (one call per loose-loop outer iteration).

    Args:
        wake_case: the wake sheet tables.
        stations: the wall case's AirfoilStations (duck-typed: xc,
            station_of, side_node) -- the TE confluence reads the two TE
            copies off it (side_node distinguishes them).
        outs: (n_wall_surf, C.N_OUT) the wall closure packet of the
            current outer iteration (delta* = OUT_DS1, theta = OUT_TH11).
        ue_surf: (n_s, 3) wake_edge_velocity output.
        m_inf, gamma_air: freestream state for the wall chain's isentropic
            rho_e.
        chord: reference chord (L_rel = l_rel_chords * chord).
        l_rel_chords: the pinned MODEL-CHOICE constant (default
            WAKE_L_REL_CHORDS; the sweep is a GV6.2 item).

    Returns:
        (ds_wake, m_wake, info): ds_wake (n_s,), m_wake (n_s,) =
        div_Gamma(rho_e u_e ds_wake), info with ds_te / th_te /
        mdot_wake_max.

    W2 (band (c)): the constructed ds_wake at s == 0 equals
    delta*_upper(TE) + delta*_lower(TE) to 1e-12 relative, asserted at
    EVERY call; a failure is a TE-copy wiring recipe error. Bypassed only
    under a prescribed_ds override (test hook, not the construction).
    """
    sm = wake_case.sm
    ds_wall = np.asarray(outs, dtype=np.float64)[:, C.OUT_DS1]
    th_wall = np.asarray(outs, dtype=np.float64)[:, C.OUT_TH11]
    te_row = int(np.argmax(stations.xc))
    at_te = stations.station_of == te_row
    upper = at_te & (stations.side_node == 1)
    lower = at_te & (stations.side_node == -1)
    if not (np.any(upper) and np.any(lower)):
        raise RuntimeError(
            "wall station table has no two-sided TE copies (the wake-cut "
            "airfoil case is required; GV6.1 scope)"
        )
    ds_te = float(np.mean(ds_wall[upper]) + np.mean(ds_wall[lower]))
    th_te = float(np.mean(th_wall[upper]) + np.mean(th_wall[lower]))

    if wake_case.prescribed_ds is not None:
        ds_wake = np.asarray(wake_case.prescribed_ds, dtype=np.float64)
        if ds_wake.shape != (sm.n_node,):
            raise ValueError(
                f"prescribed_ds must be ({sm.n_node},), got {ds_wake.shape}"
            )
    else:
        ds_wake = prescribed_delta_star_wake(
            ds_te, th_te, wake_case.s, l_rel_chords * chord
        )
        # W2: the TE-continuity construction identity, every call.
        at_s0 = wake_case.s == 0.0
        if not np.any(at_s0):
            raise RuntimeError("wake strip has no s == 0 (TE) nodes (W2)")
        err = float(np.max(np.abs(ds_wake[at_s0] - ds_te)))
        if err > _TE_CONTINUITY_RTOL * max(abs(ds_te), 1.0e-30):
            raise AssertionError(
                f"W2 TE-continuity violated: |ds_wake(0) - ds_TE| = {err:.3e} "
                f"(ds_TE = {ds_te:.6e}) -- TE-copy wiring recipe error"
            )

    ue_surf = np.asarray(ue_surf, dtype=np.float64)
    q2 = np.sum(ue_surf ** 2, axis=1)
    rho_e = density_field(q2, m_inf, gamma_air)
    m_wake = transpiration_from_delta_star(sm, rho_e, ue_surf, ds_wake)
    info = {
        "ds_te": ds_te,
        "th_te": th_te,
        "mdot_wake_max": float(np.max(np.abs(m_wake))),
    }
    return ds_wake, m_wake, info


def assemble_wake_sheet_rhs(
    nodes: np.ndarray, wake_case: WakeSheetCase, m_wake: np.ndarray
) -> np.ndarray:
    """The wake sheet-source RHS b_wake (n_nodes,) in cut-mesh volume
    numbering: m_wake scattered HALF PER FACE COPY (minus nodes and their
    plus slaves each carry 0.5*m_wake -- the 2026-07-25 addendum: the T^T
    fold sums the two copies, so the folded jump equals m_wake, not
    2*m_wake), assembled over wake_minus U wake_plus by
    ``assemble_transpiration_rhs`` (the V2 sign convention: positive
    m_wake ejects fluid AWAY from the sheet on both sides).

    m_wake == 0 returns the exact zero vector (GV2.1(b) discipline -- the
    band-(a)(i) bit-identity basis).
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    m_wake = np.asarray(m_wake, dtype=np.float64)
    if m_wake.shape != (wake_case.sm.n_node,):
        raise ValueError(
            f"m_wake must be ({wake_case.sm.n_node},), got {m_wake.shape}"
        )
    m_full = np.zeros(len(nodes), dtype=np.float64)
    m_full[wake_case.sm.volume_node_of] = 0.5 * m_wake
    m_full[wake_case.slave_of_surf] = 0.5 * m_wake
    return assemble_transpiration_rhs(nodes, wake_case.faces_both, m_full)
