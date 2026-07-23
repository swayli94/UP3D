"""Track V V3 -- loose viscous-inviscid coupling driver (binding:
docs/roadmap/track_v.md "V3 -- Loose coupling (2.5-D ladder + fuselage
smoke)", 2026-07-22 re-spec).

The loose loop (roadmap text): FP solve -> u_e -> IBL3 solve -> m_dot ->
RHS -> FP re-solve, with under-relaxation on delta*. This module is pure
ORCHESTRATION: the heavy kernels live in viscous/ibl3.py (IBL3 surface
Newton), viscous/transpiration.py (delta* -> m_dot operator and the wall
RHS assembly) and post/surface.py (per-zone u_e recovery, A4 discipline).
No Numba here -- nothing in this file is hot next to the kernels.

Two numbering conventions meet here (the V2 discipline, restated):

* the IBL surface mesh (SurfaceMesh) numbers compact surface nodes
  [0, n_s); delta*, m_dot from ``transpiration_from_delta_star`` and the
  IBL3 state live in this numbering;
* the FP drivers take the wall RHS in VOLUME node numbering of the mesh
  the driver runs on (for the lifting legs: the wake-CUT mesh -- its wall
  carries TE-duplicated nodes, which is exactly what the IBL surface wants
  built on it: upper/lower TE are distinct outflow nodes there).

``SurfaceMesh.volume_node_of`` / ``node_map`` (the design-pinned master-map
hook) mediate: gather u_e with ``ue_vol[sm.volume_node_of]``, scatter m_dot
with ``m_vol[sm.volume_node_of] = m_surf``.

The IBL3 solver consumes u_e once at construction (no setter), so the
loose loop REBUILDS ``IBL3Solver`` every outer iteration and warm-starts
the state U -- the intended usage per the V1 design.

Forced transition follows the track_v.md discipline: regime flags are set
from the prescribed x_tr/c on each side; gate comparisons must match Re
and x_tr explicitly on both sides (XFOIL side: XTR x_tr x_tr).
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from pyfp3d.physics.isentropic import density_field, mach_squared_field
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.surface_mesh import SurfaceMesh
from pyfp3d.viscous.transpiration import (
    assemble_transpiration_rhs,
    edge_velocity_per_zone,
    transpiration_from_delta_star,
)

__all__ = [
    "CouplingConfig",
    "CouplingCase",
    "CouplingResult",
    "FpSeed",
    "build_airfoil_case",
    "build_closed_body_case",
    "build_wing_case",
    "make_picard_lifting_driver",
    "make_newton_lifting_driver",
    "make_picard_nonlifting_driver",
    "run_loose_coupling",
    "station_average",
]


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass
class CouplingConfig:
    """Loose-coupling run parameters (GV3.1/GV3.2 conditions are pinned by
    the PRE_REGISTRATION of cases/analysis/v3_loose_coupling/, not here).

    re_chord: Reynolds number per unit chord. Nondimensionalization follows
        the FP solver (u_inf = rho_inf = chord = 1), so the laminar
        viscosity consumed by the closures is mu = 1 / re_chord (constant;
        XFOIL matches Re through its own input, viscosity ratio 1).
    m_inf, alpha_deg: freestream conditions (alpha only labels the case
        here; the drivers get it themselves).
    x_tr_upper / x_tr_lower: forced-transition x/c per side (turbulent
        flags set on stations with x/c >= x_tr; the LE stagnation station
        stays laminar). D-TR discipline of track_v.md.
    le_band_x: x/c below which u_e uses the A4 LE-band recovery
        (linear + crease-gated smoothed; elsewhere quadratic).
    inflow_band_x: x/c below which nodes get Dirichlet-pinned to the
        laminar seed state (per-node, at the k=1 edge data, then frozen).
        The stagnation BAND must be pinned, not just the min-q station:
        on a closed nose the flow splits into two streams there, and a
        single pinned station row leaves the split under-anchored --
        observed to send the Newton solve into the near-singular
        near-separation basin even at alpha = 0 (V3 debug, 2026-07-22;
        the x/c <= 0.02 band restores the V1 x0-line discipline and
        converges cleanly).
    omega: under-relaxation factor on delta* (ds_hat <- omega*ds_new +
        (1-omega)*ds_hat, ds_hat_0 = 0). Recorded honestly per GV3.2.
    tol_ds: loose-loop convergence, max|ds_new - ds_prev| / max|ds_new|
        (GV3.2's 1e-3 criterion).
    eps_diff / eps_diff_s: IBL3 artificial diffusion (load-bearing V1
        calibration, design_track_v.md §9.4 -- do not zero here).
    """

    re_chord: float
    m_inf: float
    alpha_deg: float = 0.0
    x_tr_upper: float = 0.05
    x_tr_lower: float = 0.05
    le_band_x: float = 0.05
    inflow_band_x: float = 0.02
    omega: float = 1.0
    n_outer_max: int = 10
    tol_ds: float = 1.0e-3
    ibl_tol: float = 1.0e-9
    ibl_max_iter: int = 100
    n_smooth_passes: int = 2
    gamma_air: float = 1.4
    eps_diff: float = 0.005
    eps_diff_s: float = 0.02
    chord: float = 1.0


# ---------------------------------------------------------------------------
# Case: the IBL surface + BC wiring on top of an FP mesh wall
# ---------------------------------------------------------------------------


@dataclass
class AirfoilStations:
    """Quasi-2D station table: one row per distinct (x, y) wall position
    (the extruded strip's span lines collapse; TE-duplicated cut nodes
    share coordinates and hence ONE station row -- the station graph of a
    wake-cut airfoil wall is therefore a closed LOOP slit conceptually at
    the shared TE station, not an open chain).

    station_of (n_s,): surface node -> station row.
    xy (n_st, 2): station coordinates. xc: x/chord.
    s (n_st,): arc distance from the stagnation (LE) station along its side.
    side (n_st,): +1 upper / -1 lower / 0 stagnation station.
    side_node (n_s,): per-NODE side, split by incident-triangle centroid y
        relative to the stagnation station -- this is what distinguishes
        the two TE copies (identical coordinates, opposite sides).
    order (n_st,): station rows in loop order starting at the stagnation
        station (LE -> TE via one side -> back toward LE via the other).
    stag_row: station row of the stagnation (min-x) station.
    le_nbrs: the two station rows adjacent to stag_row in the loop.
    """

    station_of: np.ndarray
    xy: np.ndarray
    xc: np.ndarray
    s: np.ndarray
    side: np.ndarray
    side_node: np.ndarray
    order: np.ndarray
    stag_row: int
    le_nbrs: Tuple[int, ...]


@dataclass
class CouplingCase:
    """Everything the loose loop needs that does NOT change between outer
    iterations (geometry, BC masks, seeding tables).

    nodes / elements / wall_faces: the FP driver's mesh (cut mesh for the
        lifting legs; the IBL surface is built ON those wall faces so the
        RHS scatter lands in the driver's own numbering).
    sm: the IBL SurfaceMesh on wall_faces.
    turbulent_flags (n_s,): D-TR regime per surface node.
    inflow_candidates (n_s,) bool: near-LE (or nose) nodes eligible as the
        Dirichlet inflow region; the loop pins the whole stagnation BAND
        (airfoil: stations x/c <= config.inflow_band_x; closed body: these
        candidates) at the FIRST outer iteration and freezes it (a mid-loop
        discrete switch would churn the fixed point; the single-station
        variant leaves the nose split under-anchored -- V3 debug
        2026-07-22).
    outflow_pin_surf (n_s,) bool or None: closed-body TAIL stagnation
        band (last stag_band_frac of the body length), Dirichlet-pinned
        to per-node seed states alongside the inflow. On a closed
        surface the BL characteristics converge to the aft pole -- there
        is no natural outflow boundary and the Newton Jacobian is
        exactly singular there (measured GV3.3 2026-07-22). The pin
        covers the whole tail band: a NARROW (pole+ring) pin leaves the
        tail-cone BL free to separate and hits the Goldstein singularity
        head-on (same diag). The frozen seed delta* of the pin is
        excluded from the transpiration source (run_loose_coupling masks
        it) so the convergence-zone sink cannot run away against the FP
        velocity. None on airfoil cases (the TE outflow is natural).
    seed_fetch (n_s,): per-node effective laminar fetch for U0 seeding
        (arc distance from the inflow region, floored).
    inflow_fetch: scalar fetch retained for reference (the inflow Dirichlet
        states are seeded per node from seed_fetch at the frozen k=1 edge
        data -- the V1 x0-line discipline over the stagnation band).
    le_band_surf (n_s,) bool: surface nodes in the A4 linear+smoothed u_e
        recovery zone (airfoil: x/c < le_band_x; closed body: the nose and
        tail stagnation bands, GV3.3(c) discipline).
    stations: airfoil-only station table (None for closed bodies).
    """

    nodes: np.ndarray
    elements: np.ndarray
    wall_faces: np.ndarray
    sm: SurfaceMesh
    turbulent_flags: np.ndarray
    inflow_candidates: np.ndarray
    seed_fetch: np.ndarray
    inflow_fetch: float
    le_band_surf: np.ndarray
    stations: Optional[AirfoilStations] = None
    outflow_pin_surf: Optional[np.ndarray] = None

    @property
    def n_volume(self) -> int:
        return len(self.nodes)


def _station_chain(sm: SurfaceMesh, chord: float) -> AirfoilStations:
    """Build the quasi-2D station table from the wall strip.

    Stations are grouped by (x, y) coordinates, so the TE-duplicated
    nodes of a wake-cut mesh collapse into one TE station row and the
    station graph is a closed LOOP (degree 2 everywhere): LE -> one side
    -> TE -> other side -> LE. (An open chain with two endpoints is also
    accepted for robustness.) The per-NODE side split (centroid-y rule)
    is what separates the two TE copies downstream.
    """
    xy = np.ascontiguousarray(sm.xyz[:, :2])
    uniq, station_of = np.unique(xy, axis=0, return_inverse=True)
    n_st = len(uniq)
    adj: List[set] = [set() for _ in range(n_st)]
    for t in range(sm.n_tri):
        rows = np.unique(station_of[sm.triangles[t]])
        for a in rows:
            for b in rows:
                if a != b:
                    adj[int(a)].add(int(b))
                    adj[int(b)].add(int(a))
    ends = [i for i in range(n_st) if len(adj[i]) == 1]
    stag_row = int(np.argmin(uniq[:, 0]))  # LE = min-x station

    if len(ends) == 2:
        # open chain (kept for robustness): walk endpoint -> endpoint
        order = [ends[0]]
        prev = -1
        while order[-1] != ends[1]:
            cur = order[-1]
            nxt = [j for j in adj[cur] if j != prev]
            if len(nxt) != 1:
                raise ValueError(
                    "station chain branches (not a quasi-2D strip)"
                )
            prev = cur
            order.append(nxt[0])
        if len(order) != n_st:
            raise ValueError("station chain does not cover all stations")
        # rotate so the chain starts at the stagnation station
        pos = order.index(stag_row)
        order = order[pos:] + order[:pos]
    elif not ends:
        # closed loop: walk from the LE station back to itself
        order = [stag_row]
        prev, cur = -1, stag_row
        nxt = sorted(adj[cur])[0]
        while nxt != stag_row:
            order.append(nxt)
            prev, cur = cur, nxt
            nxts = [j for j in adj[cur] if j != prev]
            if len(nxts) != 1:
                raise ValueError(
                    "station loop branches (not a quasi-2D strip)"
                )
            nxt = nxts[0]
        if len(order) != n_st:
            raise ValueError("station loop does not cover all stations")
    else:
        raise ValueError(
            f"wall strip graph is neither a chain nor a loop "
            f"({len(ends)} endpoint stations)"
        )
    order = np.asarray(order, dtype=np.int64)

    # split the loop at the TE: the two branches from the LE station
    pos_te = int(np.argmax(uniq[order, 0]))  # TE = max-x station
    branch_a = order[1:pos_te]  # LE -> TE (exclusive)
    branch_b = order[pos_te + 1 :][::-1]  # LE -> TE (exclusive), reversed
    side = np.zeros(n_st, dtype=np.int64)
    ma = float(np.mean(uniq[branch_a, 1])) if len(branch_a) else -np.inf
    mb = float(np.mean(uniq[branch_b, 1])) if len(branch_b) else -np.inf
    upper_branch, lower_branch = (
        (branch_a, branch_b) if ma >= mb else (branch_b, branch_a)
    )
    side[upper_branch] = 1
    side[lower_branch] = -1
    side[int(order[pos_te])] = 0  # TE station: split per-node instead

    s = np.zeros(n_st, dtype=np.float64)
    for branch in (branch_a, branch_b):
        acc = 0.0
        prev_xy = uniq[stag_row]
        for row in branch:  # walk stagnation -> TE
            acc += float(np.linalg.norm(uniq[row] - prev_xy))
            s[row] = acc
            prev_xy = uniq[row]
    # TE station gets the mean of its two branch ends (only used as a
    # fetch for seeding; the TE copies themselves carry per-node sides)
    te_row = int(order[pos_te])
    s[te_row] = max(s[branch_a[-1]], s[branch_b[-1]]) if len(branch_a) and len(branch_b) else 0.0

    # per-node side: mean incident-triangle centroid y vs the stagnation y
    cent = sm.xyz[sm.triangles].mean(axis=1)  # (F, 3)
    ysum = np.zeros(sm.n_node, dtype=np.float64)
    ycnt = np.zeros(sm.n_node, dtype=np.float64)
    np.add.at(ysum, sm.triangles.reshape(-1), np.repeat(cent[:, 1], 3))
    np.add.at(ycnt, sm.triangles.reshape(-1), 1.0)
    y_stag = float(uniq[stag_row, 1])
    side_node = np.where(ysum / np.maximum(ycnt, 1.0) >= y_stag, 1, -1)

    le_nbrs = tuple(
        int(order[j]) for j in (1, len(order) - 1) if len(order) > 2
    )
    return AirfoilStations(
        station_of=station_of.astype(np.int64),
        xy=uniq,
        xc=uniq[:, 0] / chord,
        s=s,
        side=side,
        side_node=side_node,
        order=order,
        stag_row=stag_row,
        le_nbrs=le_nbrs,
    )


def station_average(case: CouplingCase, field_surf: np.ndarray) -> np.ndarray:
    """Per-station mean of a surface-nodal field (averages the spanwise
    copies of a quasi-2D strip). Airfoil cases only."""
    if case.stations is None:
        raise ValueError("station table available on airfoil cases only")
    st = case.stations
    out = np.zeros(len(st.xc), dtype=np.float64)
    cnt = np.zeros(len(st.xc), dtype=np.float64)
    np.add.at(out, st.station_of, field_surf)
    np.add.at(cnt, st.station_of, 1.0)
    return out / np.maximum(cnt, 1.0)


def _turb_seed(s: float, q: float, rho: float, mu: float) -> np.ndarray:
    """Power-law turbulent seed (the V1 recipe, v1_ibl3_standalone/run.py):
    delta = 0.37 x Re_x^-0.2, cf = 0.0576 Re_x^-0.2, A = cf Re_d / 2,
    C_tau1 at the closure's own equilibrium."""
    s = max(s, 1.0e-8)
    q = max(q, 1.0e-8)
    re_x = rho * q * s / mu
    delta = 0.37 * s * re_x ** -0.2
    cf = 0.0576 * re_x ** -0.2
    re_d = rho * q * delta / mu
    st = np.array([delta, 0.5 * cf * re_d, 0.0, 0.0, 1.0e-3, 0.0])
    out, _, _ = C.closure_scalar(st, q=q, rho=rho, mu=mu, turbulent=True)
    st[4] = max((C.C_L_DEFAULT * out[C.OUT_SP1] / out[C.OUT_SD]) ** 2, 1.0e-6)
    return st


def _lam_seed(s: float, q: float, rho: float, mu: float) -> np.ndarray:
    """Laminar Blasius-matched seed (closures.blasius_seed)."""
    return C.blasius_seed(x=max(s, 1.0e-8), q=max(q, 1.0e-8), rho=rho, mu=mu)


def build_airfoil_case(
    nodes: np.ndarray,
    elements: np.ndarray,
    wall_faces: np.ndarray,
    config: CouplingConfig,
) -> CouplingCase:
    """IBL case on the wall of a wake-CUT quasi-2D airfoil mesh.

    The TE duplication of the cut mesh gives each side its own TE outflow
    nodes; the strip chain runs TE -> LE -> TE. Transition flags come from
    config.x_tr_upper/lower; the LE band (x/c < le_band_x) takes the A4
    linear+smoothed u_e recovery.
    """
    sm = SurfaceMesh.from_wall_faces(nodes, wall_faces, elements)
    st = _station_chain(sm, config.chord)
    xc_n = st.xc[st.station_of]
    flags = np.zeros(sm.n_node, dtype=np.int64)
    flags[(st.side_node == 1) & (xc_n >= config.x_tr_upper)] = 1
    flags[(st.side_node == -1) & (xc_n >= config.x_tr_lower)] = 1
    inflow_candidates = st.xc[st.station_of] <= config.le_band_x
    seed_fetch = np.maximum(st.s[st.station_of], 0.0)
    # inflow fetch: mean distance from the LE (stagnation) station to its
    # loop neighbours -- "the BL starts one cell downstream".
    d = [
        float(np.linalg.norm(st.xy[r] - st.xy[st.stag_row]))
        for r in st.le_nbrs
    ]
    inflow_fetch = max(float(np.mean(d)), 1.0e-4 * config.chord)
    return CouplingCase(
        nodes=np.asarray(nodes, dtype=np.float64),
        elements=np.asarray(elements),
        wall_faces=np.asarray(wall_faces),
        sm=sm,
        turbulent_flags=flags,
        inflow_candidates=np.asarray(inflow_candidates, dtype=bool),
        seed_fetch=seed_fetch,
        inflow_fetch=inflow_fetch,
        le_band_surf=np.asarray(inflow_candidates, dtype=bool),
        stations=st,
    )


def build_closed_body_case(
    nodes: np.ndarray,
    elements: np.ndarray,
    wall_faces: np.ndarray,
    config: CouplingConfig,
    x_tr_frac: float = 0.05,
    stag_band_frac: float = 0.05,
) -> CouplingCase:
    """IBL case on a CLOSED wall (GV3.3 body of revolution): no natural
    boundary edges, so the Dirichlet inflow is the nose stagnation pole
    (alpha = 0 puts the stagnation exactly there). Transition is forced at
    x_tr_frac of the body length; the nose/tail stagnation bands get the
    A4 linear+smoothed u_e recovery (GV3.3(c) discipline).
    """
    sm = SurfaceMesh.from_wall_faces(nodes, wall_faces, elements)
    x = sm.xyz[:, 0]
    x0, x1 = float(np.min(x)), float(np.max(x))
    body_len = x1 - x0
    flags = np.zeros(sm.n_node, dtype=np.int64)
    flags[x >= x0 + x_tr_frac * body_len] = 1
    nose = int(np.argmin(x))
    inflow_candidates = np.zeros(sm.n_node, dtype=bool)
    # candidates = the nose pole plus its one ring; the loop pins the
    # single min-q candidate (the pole itself at alpha = 0).
    ring = np.unique(sm.triangles[np.any(sm.triangles == nose, axis=1)])
    inflow_candidates[ring] = True
    h_loc = float(
        np.mean(
            np.linalg.norm(sm.xyz[ring] - sm.xyz[nose], axis=1)[ring != nose]
        )
    )
    # per-node seed fetch: distance along x from the nose (a meridian arc
    # distance is finer-grained than a seed needs)
    seed_fetch = np.maximum(x - x0, h_loc)
    # A4 stagnation-band u_e recovery at BOTH ends (GV3.3(c))
    le_band = (x < x0 + stag_band_frac * body_len) | (
        x > x1 - stag_band_frac * body_len
    )
    # The TAIL stagnation band is also a Dirichlet pin: on a closed
    # surface the BL characteristics converge to the aft pole, leaving the
    # marching system without a natural outflow -- the Newton Jacobian is
    # exactly singular there (measured GV3.3 2026-07-22). The pin covers
    # the whole stagnation band: a NARROW (pole+ring) pin leaves the
    # tail-cone BL free to separate, and the separation (Goldstein)
    # singularity crashes the Newton solve outright (same diag). The
    # pin's artificial seed delta* is excluded from the transpiration
    # source (run_loose_coupling masks it) so the tail-cone convergence
    # sink cannot run away against the FP velocity.
    tail_pin = x > x1 - stag_band_frac * body_len
    return CouplingCase(
        nodes=np.asarray(nodes, dtype=np.float64),
        elements=np.asarray(elements),
        wall_faces=np.asarray(wall_faces),
        sm=sm,
        turbulent_flags=flags,
        inflow_candidates=inflow_candidates,
        seed_fetch=seed_fetch,
        inflow_fetch=max(h_loc, 1.0e-6),
        le_band_surf=le_band,
        stations=None,
        outflow_pin_surf=tail_pin,
    )


def build_wing_case(
    nodes: np.ndarray,
    elements: np.ndarray,
    wall_faces: np.ndarray,
    config: CouplingConfig,
    x_le: Callable[[np.ndarray], np.ndarray],
    chord_at: Callable[[np.ndarray], np.ndarray],
    tip_mask_frac: float = 0.05,
) -> CouplingCase:
    """IBL case on a 3-D lifting wing wall (GV5.0 M6 bridge; GV5.3 gate).

    The wall of the wake-CUT mesh carries TE-duplicated nodes, so the IBL
    surface has natural outflow boundary edges along BOTH TE lines (the
    airfoil discipline, now per span station); the root section (z = 0
    symmetry plane) is an open boundary edge and takes the natural
    zero-flux condition, which IS symmetry. Boundary conditions:

    * INFLOW: the LE band (local x/c <= config.inflow_band_x), Dirichlet-
      pinned to per-node laminar Blasius seeds (the airfoil stagnation-
      band discipline, evaluated on the LOCAL section: x/c =
      (x - x_le(z)) / chord_at(z); planform callables are caller-supplied
      so the builder stays mesh-agnostic).
    * TIP MASK: the band z > z_tip * (1 - tip_mask_frac) -- the
      production tip_taper radius r_c = 0.05 * b_semi (B32; the tip-edge
      singularity zone is outside the VII validity envelope, track_v.md
      scope guards). The band is Dirichlet-pinned to per-node regime seed
      states AND masked out of the transpiration source: mechanically it
      reuses ``outflow_pin_surf`` (pin + m_dot mask), i.e. exactly the
      GV3.3 tail-pin semantics -- the frozen seed delta* there is
      boundary data, not solution, and generates no transpiration.
      ``z_tip`` = max wall z (= B_SEMI on the flat-cap M6 family).
    * The A4 u_e recovery zone (le_band_surf) covers the LE band
      (x/c < config.le_band_x) PLUS the tip mask band (the linear+smoothed
      path is the robust one next to the tip-edge singularity; the band
      is masked from every comparison anyway).

    Transition is forced per side at config.x_tr_upper/lower of the local
    chord; the per-node side split is the incident-triangle centroid y
    sign (the M6 section is symmetric about the chord plane y = 0).
    stations is None: run_loose_coupling takes its closed-body branch,
    pinning candidates | outflow_pin_surf with per-node regime seeds --
    which is exactly the wing's pin set.
    """
    sm = SurfaceMesh.from_wall_faces(nodes, wall_faces, elements)
    x, y, z = sm.xyz[:, 0], sm.xyz[:, 1], sm.xyz[:, 2]
    z_tip = float(np.max(z))
    zc = np.clip(z, 0.0, z_tip)
    chord_n = chord_at(zc)
    xc_n = (x - x_le(zc)) / chord_n

    # per-node side: mean incident-triangle centroid y vs the chord plane
    cent_y = sm.xyz[sm.triangles].mean(axis=1)[:, 1]
    ysum = np.zeros(sm.n_node, dtype=np.float64)
    ycnt = np.zeros(sm.n_node, dtype=np.float64)
    np.add.at(ysum, sm.triangles.reshape(-1), np.repeat(cent_y, 3))
    np.add.at(ycnt, sm.triangles.reshape(-1), 1.0)
    side_node = np.where(ysum / np.maximum(ycnt, 1.0) >= 0.0, 1, -1)

    tip_mask = z > z_tip * (1.0 - tip_mask_frac)
    flags = np.zeros(sm.n_node, dtype=np.int64)
    flags[(side_node == 1) & (xc_n >= config.x_tr_upper)] = 1
    flags[(side_node == -1) & (xc_n >= config.x_tr_lower)] = 1

    inflow_candidates = (xc_n <= config.inflow_band_x) & ~tip_mask
    seed_fetch = np.maximum(x - x_le(zc), 1.0e-4 * chord_n)
    inflow_fetch = (
        float(np.mean(seed_fetch[inflow_candidates]))
        if np.any(inflow_candidates)
        else float(np.mean(seed_fetch))
    )
    le_band = (xc_n < config.le_band_x) | tip_mask
    return CouplingCase(
        nodes=np.asarray(nodes, dtype=np.float64),
        elements=np.asarray(elements),
        wall_faces=np.asarray(wall_faces),
        sm=sm,
        turbulent_flags=flags,
        inflow_candidates=np.asarray(inflow_candidates, dtype=bool),
        seed_fetch=seed_fetch,
        inflow_fetch=inflow_fetch,
        le_band_surf=np.asarray(le_band, dtype=bool),
        stations=None,
        outflow_pin_surf=np.asarray(tip_mask, dtype=bool),
    )


# ---------------------------------------------------------------------------
# FP driver adapters: uniform (rhs, seed) -> (phi, gamma, info) protocol
# ---------------------------------------------------------------------------


@dataclass
class FpSeed:
    """Warm start passed between outer iterations."""

    phi: Optional[np.ndarray] = None
    gamma: Optional[np.ndarray] = None


FpSolve = Callable[[Optional[np.ndarray], Optional[FpSeed]], Tuple[np.ndarray, Optional[np.ndarray], Dict]]


def make_picard_lifting_driver(mc, wc, m_inf: float, alpha_deg: float, **kw) -> FpSolve:
    """Compressible-Picard lifting leg (GV3.1): solve_subsonic_lifting
    with the V2-threaded body_source_rhs (cut-mesh numbering)."""
    from pyfp3d.solve.picard import solve_subsonic_lifting

    def solve(rhs, seed):
        r = solve_subsonic_lifting(
            mc,
            wc,
            m_inf=m_inf,
            alpha_deg=alpha_deg,
            body_source_rhs=rhs,
            phi_init=None if seed is None else seed.phi,
            gamma_init=None if seed is None else seed.gamma,
            **kw,
        )
        return r["phi"], r["gamma"], r

    return solve


def make_newton_lifting_driver(mc, wc, m_inf: float, alpha_deg: float, **kw) -> FpSolve:
    """Conforming-Newton lifting leg (GV3.2's transonic recorded point):
    solve_newton_lifting with the V2 external_rhs channel."""
    from pyfp3d.solve.newton import solve_newton_lifting

    def solve(rhs, seed):
        r = solve_newton_lifting(
            mc,
            wc,
            m_inf=m_inf,
            alpha_deg=alpha_deg,
            external_rhs=rhs,
            phi_init=None if seed is None else seed.phi,
            gamma_init=None if seed is None else seed.gamma,
            **kw,
        )
        return r["phi"], r["gamma"], r

    return solve


def make_picard_nonlifting_driver(
    nodes, elements, dirichlet_nodes, dirichlet_values, m_inf: float, **kw
) -> FpSolve:
    """Non-lifting compressible Picard (GV3.3 body of revolution):
    solve_subsonic with body_source_rhs (uncut mesh, no wake/Kutta)."""
    from pyfp3d.solve.picard import solve_subsonic

    def solve(rhs, seed):
        r = solve_subsonic(
            nodes,
            elements,
            dirichlet_nodes,
            dirichlet_values,
            m_inf=m_inf,
            body_source_rhs=rhs,
            phi_init=None if seed is None else seed.phi,
            **kw,
        )
        return r["phi"], None, r

    return solve


# ---------------------------------------------------------------------------
# The loose loop
# ---------------------------------------------------------------------------


@dataclass
class CouplingResult:
    phi: np.ndarray
    gamma: Optional[np.ndarray]
    U: np.ndarray
    outs: np.ndarray
    delta_star: np.ndarray  # relaxed field used for the final m_dot
    m_dot: np.ndarray  # volume numbering
    ue_surf: np.ndarray
    n_outer: int
    converged: bool
    history: List[Dict] = field(default_factory=list)
    ibl_info: Dict = field(default_factory=dict)
    driver_info: Dict = field(default_factory=dict)


def run_loose_coupling(
    fp_solve: FpSolve,
    case: CouplingCase,
    config: CouplingConfig,
    probe: Optional[Callable[[np.ndarray, Optional[np.ndarray], int], Dict]] = None,
) -> CouplingResult:
    """The V3 loose loop.

    fp_solve(rhs, seed): FP driver adapter (None rhs = inviscid).
    probe(phi, gamma, k): optional per-iteration callback (k = 0 is the
        inviscid baseline) returning a dict merged into the history record
        -- cl/Cp logging lives in the caller, not here.

    Convergence (GV3.2): max|ds_new - ds_prev| / max|ds_new| < tol_ds on
    the RAW successive IBL outputs (the fixed-point residual; the relaxed
    ds_hat is what the FP solve sees). The final FP re-solve still runs so
    the returned phi is consistent with the returned delta_star.
    """
    sm = case.sm
    mu = 1.0 / config.re_chord
    n_vol = case.n_volume
    le_mask_vol = np.zeros(n_vol, dtype=bool)
    le_mask_vol[sm.volume_node_of[case.le_band_surf]] = True

    history: List[Dict] = []

    # -- inviscid baseline --------------------------------------------------
    phi, gamma, dinfo = fp_solve(None, None)
    base = {"k": 0, "mdot_max": 0.0}
    if probe is not None:
        base.update(probe(phi, gamma, 0))
    history.append(base)

    # -- loop ---------------------------------------------------------------
    U = None
    ds_prev = None
    ds_hat = None
    m_dot = np.zeros(n_vol, dtype=np.float64)
    inflow_mask = None
    inflow_state = None
    outs = None
    ue_surf = None
    converged = False
    ibl_info: Dict = {}
    k_done = 0

    for k in range(1, config.n_outer_max + 1):
        k_done = k
        ue_vol = edge_velocity_per_zone(
            case.nodes,
            case.wall_faces,
            phi,
            elements=case.elements,
            le_band_mask=le_mask_vol,
            n_smooth_passes=config.n_smooth_passes,
        )
        ue_surf = ue_vol[sm.volume_node_of]
        if not np.all(np.isfinite(ue_surf)):
            raise RuntimeError(
                "u_e gather hit NaN on the IBL surface (off-wall leak?)"
            )
        q2 = np.sum(ue_surf ** 2, axis=1)
        q = np.sqrt(q2)
        rho_e = density_field(q2, config.m_inf, config.gamma_air)
        mach_e = np.sqrt(mach_squared_field(q2, config.m_inf, config.gamma_air))

        if inflow_mask is None:
            # freeze the inflow BC at the first outer iteration: the whole
            # stagnation BAND is Dirichlet-pinned to per-node laminar seed
            # states (airfoil: stations with x/c <= inflow_band_x; closed
            # body: the nose pole + first ring). Pinning only the min-q
            # station row leaves the two streams splitting at the nose
            # under-anchored -- observed to send the Newton solve into
            # the near-singular near-separation basin even at alpha = 0
            # (V3 debug 2026-07-22; the band restores the V1 x0-line
            # discipline).
            if case.stations is not None:
                xc_n = case.stations.xc[case.stations.station_of]
                inflow_mask = xc_n <= config.inflow_band_x
                seed_kind = None  # laminar Blasius at the airfoil nose
            else:
                inflow_mask = case.inflow_candidates.copy()
                if case.outflow_pin_surf is not None:
                    # closed body: the BL characteristics converge to the
                    # aft pole -- no natural outflow exists and the Newton
                    # Jacobian is exactly singular (measured GV3.3
                    # 2026-07-22). Pin the tail stagnation band too, with
                    # TURBULENT seed states (it sits aft of x_tr).
                    inflow_mask |= case.outflow_pin_surf
                seed_kind = case.turbulent_flags  # per-node lam/turb seed
            idx_in = np.where(inflow_mask)[0]
            i_stag = idx_in[int(np.argmin(q[idx_in]))]
            q_ref = float(np.percentile(q[idx_in], 75))
            q_floor = 0.05 * max(q_ref, 1.0e-12)
            inflow_state = np.stack([
                (_turb_seed(case.seed_fetch[i], max(float(q[i]), q_floor),
                            1.0, mu)
                 if seed_kind is not None and seed_kind[i] else
                 _lam_seed(case.seed_fetch[i], max(float(q[i]), 1.0e-8),
                           1.0, mu))
                for i in idx_in
            ])
            stag_note = {
                "inflow_node": int(i_stag),
                "inflow_n_pinned": int(inflow_mask.sum()),
                "inflow_q_ref": q_ref,
            }
        else:
            stag_note = {}

        solver = IBL3Solver(
            sm,
            ue_surf,
            rho_e,
            mu,
            mach_e,
            case.turbulent_flags,
            inflow_mask,
            inflow_state,
            eps_diff=config.eps_diff,
            eps_diff_s=config.eps_diff_s,
        )
        if U is None:
            q_floor = 0.02 * max(float(np.max(q)), 1.0e-12)
            U = np.zeros((sm.n_node, 6), dtype=np.float64)
            for i in range(sm.n_node):
                qq = max(q[i], q_floor)
                if case.turbulent_flags[i]:
                    U[i] = _turb_seed(case.seed_fetch[i], qq, 1.0, mu)
                else:
                    U[i] = _lam_seed(case.seed_fetch[i], qq, 1.0, mu)
        U, ibl_info = solver.solve(
            U, tol=config.ibl_tol, max_iter=config.ibl_max_iter
        )
        if not np.all(np.isfinite(U)):
            raise RuntimeError(
                f"IBL3 produced a non-finite state at outer iter {k}"
            )
        outs = np.empty((sm.n_node, C.N_OUT), dtype=np.float64)
        douts = np.empty((sm.n_node, C.N_OUT, 6), dtype=np.float64)
        douts_e = np.empty((sm.n_node, C.N_OUT, 2), dtype=np.float64)
        C.closure_all(
            U, q, rho_e, np.full(sm.n_node, mu), mach_e,
            case.turbulent_flags, C.C_L_DEFAULT, outs, douts, douts_e,
        )
        ds_new = outs[:, C.OUT_DS1]

        if ds_prev is None:
            change_rel = np.inf
        else:
            change_rel = float(
                np.max(np.abs(ds_new - ds_prev))
                / max(float(np.max(np.abs(ds_new))), 1.0e-30)
            )
        ds_hat = config.omega * ds_new + (1.0 - config.omega) * (
            ds_hat if ds_hat is not None else 0.0
        )
        n_neg = int(np.count_nonzero(ds_hat < 0.0))
        ds_hat = np.maximum(ds_hat, 0.0)  # physical floor; count recorded

        m_surf = transpiration_from_delta_star(sm, rho_e, ue_surf, ds_hat)
        if case.outflow_pin_surf is not None:
            # the tail pin's delta* is a frozen seed (boundary data, not
            # solution): in the tail-cone convergence geometry its
            # transpiration sink feeds back with the FP velocity and runs
            # away (k=3 sink 3.3 -> k=4 FP non-converged -> k=5 blowup,
            # GV3.3 diag 2026-07-22). Boundary data generates no
            # transpiration; the resulting net-source imbalance is tiny
            # and absorbed by the Dirichlet far field.
            m_surf = m_surf.copy()
            m_surf[case.outflow_pin_surf] = 0.0
        m_dot = np.zeros(n_vol, dtype=np.float64)
        m_dot[sm.volume_node_of] = m_surf
        rhs = assemble_transpiration_rhs(case.nodes, case.wall_faces, m_dot)
        phi, gamma, dinfo = fp_solve(rhs, FpSeed(phi=phi, gamma=gamma))
        if dinfo.get("converged") is False:
            # a non-converged FP solve only feeds garbage forward (GV3.3
            # runaway diag 2026-07-22: k=4 phi ~ 2e23, k=5 NaN) -- fail
            # loudly at the iteration that actually broke
            raise RuntimeError(
                f"FP driver did not converge at outer iter {k} "
                f"(mdot_max={float(np.max(np.abs(m_surf))):.3e})"
            )

        rec = {
            "k": k,
            "ds_max": float(np.max(ds_new)),
            "ds_change_rel": change_rel,
            "ds_neg_floored": n_neg,
            "mdot_max": float(np.max(np.abs(m_surf))),
            "ibl_n_iter": ibl_info["n_iter"],
            "ibl_converged": bool(ibl_info["converged"]),
            "ibl_final_residual": ibl_info["final_residual"],
        }
        rec.update(stag_note)
        if probe is not None:
            rec.update(probe(phi, gamma, k))
        history.append(rec)

        if change_rel < config.tol_ds:
            converged = True
            ds_prev = ds_new
            break
        ds_prev = ds_new

    return CouplingResult(
        phi=phi,
        gamma=gamma,
        U=U,
        outs=outs,
        delta_star=ds_hat,
        m_dot=m_dot,
        ue_surf=ue_surf,
        n_outer=k_done,
        converged=converged,
        history=history,
        ibl_info=ibl_info,
        driver_info=dinfo if isinstance(dinfo, dict) else {},
    )
