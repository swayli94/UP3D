"""
Level-set (Track B) Newton solver -- the post-B6 re-derivation
(design_track_b.md section 5.5). Closes the FP-fold cases (e.g. medium
M0.7875) that the B6 Picard reaches only as a bounded stall: a Picard method,
conforming or level-set, does not reach the isolated fold solution (P4-erratum
/ P8 record), so the medium quantitative gate needs Newton.

Structure (simpler than the P8 conforming Newton, section 5.5): NO Gamma DOF
(the implicit Kutta makes the TE jump a solution mode -- no delta-Gamma
elimination, no Woodbury, no far-field vortex column when the far field is the
B6 neumann outlet); the wake-LS rows are LINEAR in phi (constant Jacobian
block); Terms 1-3 (kernels/jacobian.py physics) are reused per side through the
DOF indirection.

The residual is the multivalued fixed-point residual R = A(phi) phi - b, where
A is `MultivaluedOperator.assemble_matrix` at the current phi -- exact for all
rows including the TE Kutta row (its factorization makes row.x = |q_u|^2 -
|q_l|^2 exactly). The Jacobian adds, on top of that Picard matrix:
  - Terms 2+3 on the mass rows (per side, P7 frozen-selection sensitivities);
  - the EXACT quadratic derivative of the TE Kutta row (replacing the
    frozen-mean linearized row -- the frozen row is exact as a value but not
    as a slope, which is what quadratic convergence needs).
The wake-LS rows need no correction (already exact-linear).

Nonsymmetric (wake-LS aux rows + supersonic upwind) -> sparse-direct LU
(scipy spsolve) by default, consistent with B2/B3's choice. B11 (2026-07-14)
adds the `precond` kwarg (None = the bit-identical direct default; "ilu"/"amg"
run preconditioned GMRES on the fused Jacobian): the iterative escape for when
the LU dominates at 3D sizes (P8/N6 measured true-3D LU fill ~100x the 2.5D
cost). **Use "ilu"** -- it factors the real fused Jacobian and converges; "amg"
(the SPD surrogate `picard_ls._amg_surrogate_preconditioner`) STALLS on the
wake_ls-closure operator (measured; convection-like aux rows), so it is not the
Newton escape (design_track_b.md §5.3). B11 measured that ILU-GMRES is a
COARSE-only escape on this fused matrix (it diverges at 2.5D medium lifting and
`factor_failed`s at M6 medium), so at medium/M6 sizes sparse-direct is the only
converging tool -- and the cost driver is then the NUMBER of factorizations.
B12 (2026-07-14) adds the lagged-LU direct-reuse path for exactly that regime:
with `precond=None` and `direct_refactor_every > 1`, refactor the LU only every
k-th Newton step and drive the steps in between with GMRES on the FRESH
Jacobian preconditioned by the stale (exact) LU (the N6 mechanism ported from
solve/newton.py, minus the Woodbury -- the level-set system has no Gamma
coupling). `direct_refactor_every=1` (default) keeps the bit-identical
per-step `spsolve`.

B14 (2026-07-17) adds `precond="schur"`, the STRUCTURAL iterative escape
(solve/schur_ls.py): eliminate the aux thin-strip block exactly per step and
precondition the reduced main-free operator with AMG on the SPD Picard block
-- no springs, no full-size factorization, the conforming-Newton operator
shape. The ramp wrapper forwards it via **newton_kw.
"""

import time
from typing import Dict, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pyfp3d.constraints.dirichlet import freestream_phi, vortex_phi_2d
from pyfp3d.kernels.cut_assembly import (
    newton_terms23_side_coo,
    te_kutta_jacobian_coo,
    te_kutta_residual,
    te_weld_coo,
)
from pyfp3d.solve.linear import build_ilu_preconditioner, solve_gmres
from pyfp3d.solve.schur_ls import SchurReducedSystem, main_block_preconditioner
from pyfp3d.solve.timing import (
    finalize,
    new_timings,
    phase,
    snapshot,
    step_delta,
    sum_timings,
)
from pyfp3d.solve.picard_ls import (
    _amg_surrogate_preconditioner,
    _farfield_main,
    _farfield_split,
    _neumann_outlet_rhs,
    solve_multivalued_lifting,
)


def _seed_from_picard(mvop, mesh, m_inf, alpha_deg, u_inf, gamma_air,
                      upwind_c, m_crit, farfield, phi_init, gamma_init,
                      n_seed, damping_theta, omega_rho, precond=None,
                      direct_refactor_every=1, direct_reuse_rtol=1e-10):
    """Warm start: a short B6 Picard-LS solve to land near the solution before
    Newton (the transonic ramp / a lower-Mach state is the caller's job).

    B15: the seed inherits B13's lagged-LU (`direct_refactor_every`). Without it
    the seed pays a FULL sparse factorization on every one of its outers -- at
    M6 medium that is ~17 s x n_seed, i.e. ~11 min of pure warm-start on a 40-
    outer seed, which dwarfed the Newton solve it was seeding. Note the seed
    keeps the LIFTING default `direct_reuse_rtol=1e-10`, NOT Newton's 1e-8: a
    Picard fixed point is pinned only by its lag tolerances, so an inexact reuse
    step SHIFTS where it stops (B13, picard_ls.py:361-366).
    """
    r = solve_multivalued_lifting(
        mvop, mesh, m_inf, alpha_deg=alpha_deg, u_inf=u_inf,
        gamma_air=gamma_air, upwind_c=upwind_c, m_crit=m_crit,
        farfield=farfield, damping_theta=damping_theta, omega_rho=omega_rho,
        n_outer_max=n_seed, gamma_init=gamma_init, phi_init=phi_init,
        tol_residual=None, precond=precond,
        direct_refactor_every=direct_refactor_every,
        direct_reuse_rtol=direct_reuse_rtol,
    )
    return r["phi_ext"], r["gamma"]


class LSNewtonSystem:
    """Residual + exact Jacobian of the level-set Newton system at a given
    upwind SELECTION (live or frozen).

    Extracted at B15 so the driver and the finite-difference gate share ONE
    assembly path: the Jacobian must be the derivative of exactly the residual
    the solver evaluates, and an FD test that re-implements the assembly would
    only be testing its own copy.

    R = A(phi) phi - b (A = `MultivaluedOperator.assemble_matrix`, exact on the
    mass, wake-LS and frozen-mean TE Kutta rows), with the TE aux rows
    overwritten by the true nonlinear Kutta residual. J adds, on top of that
    Picard matrix: Terms 2+3 per side (P7 sensitivities AT THE CURRENT
    SELECTION -- frozen ones under a freeze) and the exact quadratic derivative
    of the TE Kutta row. The wake-LS rows need no correction (exact-linear).
    """

    def __init__(self, mvop, m_inf, upwind_c=1.5, m_crit=0.95, gamma_air=1.4,
                 u_inf=1.0, m_cap=3.0, rho_floor=0.05, b_base=None,
                 tip_taper=None, timings=None):
        # A1: the driver hands in its own timings dict; None (the FD tests and
        # any other direct user) keeps the class free of accounting.
        self.timings = timings
        self.mvop = mvop
        self.m_inf = float(m_inf)
        self.upwind_c = float(upwind_c)
        self.m_crit = float(m_crit)
        self.gamma_air = float(gamma_air)
        self.u_inf = float(u_inf)
        self.m_cap = float(m_cap)
        self.rho_floor = float(rho_floor)
        n_total = mvop.n_total
        self.n_total = n_total
        self.b_base = (np.zeros(n_total) if b_base is None
                       else np.asarray(b_base, dtype=np.float64))

        self.te_aux = mvop.cm.ext_dof_of_node[
            mvop.cm.te_nodes].astype(np.int64)
        self.te_main = np.asarray(mvop.cm.te_nodes, dtype=np.int64)

        # B8 tip taper (roadmap Track B B8; P13/G13.2 finding (8)). The blended
        # TE Kutta residual is F_i*(|q_u|^2-|q_l|^2) + (1-F_i)*(phi_aux-phi_main)
        # and its Jacobian is F*J_kutta + (1-F)*weld. tip_taper=None (or all
        # ones) -> pure pressure Kutta, bit-identical. Because D_keep zeros the
        # TE aux rows of the Picard matrix and R[te_aux] is overwritten, the
        # blend is applied ONLY here -- assemble_matrix keeps tip_taper=None.
        tt = None
        if tip_taper is not None:
            tt = np.asarray(tip_taper, dtype=np.float64)
            if tt.shape != self.te_main.shape:
                raise ValueError(
                    f"tip_taper must be ({self.te_main.size},), got {tt.shape}")
            if np.all(tt == 1.0):
                tt = None
        self.tt = tt
        self.f_of_row = np.ones(n_total, dtype=np.float64)  # F_i by TE aux dof
        if tt is not None:
            self.f_of_row[self.te_aux] = tt

        # Terms 2/3 are the Jacobian of the MASS-conservation residual, so their
        # rows must follow the SAME row mapping assemble_matrix applies to
        # mass_conservation_coo (design_track_b.md 5.5): the non-TE aux rows are
        # DROPPED (replaced by the phi-independent wake LS) and the TE aux rows
        # are REROUTED to the TE main row (their mass balance moved there when
        # the TE aux row became the Kutta condition). Without this, Terms 2/3
        # leak the dropped mass Jacobian onto the linear wake-LS rows (FD-caught).
        nonte_aux = np.asarray(mvop._nonte_rows, dtype=np.int64)
        row_map = np.arange(n_total, dtype=np.int64)
        row_map[nonte_aux] = -1                     # drop
        row_map[self.te_aux] = mvop.cm.te_nodes     # reroute TE aux -> TE main
        self.row_map = row_map
        # a diagonal that zeros the TE aux rows (replaced by the exact Kutta
        # Jacobian) while keeping every other row of the Picard matrix
        keep_row = np.ones(n_total, dtype=np.float64)
        keep_row[self.te_aux] = 0.0
        self.D_keep = sp.diags(keep_row)

    def te_residual(self, phi):
        r = te_kutta_residual(self.mvop, phi)
        if self.tt is None:
            return r
        return (self.tt * r + (1.0 - self.tt)
                * (phi[self.te_aux] - phi[self.te_main]))

    def freeze(self, phi):
        """Per-side (upstream, branch) selection captured at `phi` (B15/N5)."""
        return self.mvop.freeze_side_state(
            phi, self.m_inf, self.upwind_c, self.m_crit, self.gamma_air,
            self.u_inf, self.m_cap, self.rho_floor)

    def residual(self, phi, frozen=None):
        """(A, R, up, lo) at `phi` under `frozen` (None = the live selection).
        Refreshes mvop's n_limited / n_floored monitors for THAT state.

        A1 cost accounting: this residual CONTAINS an assemble_matrix, so the
        three sub-costs are split here rather than being timed as one block by
        the caller -- timing the whole call as `residual` would silently charge
        the level-set path's assembly to the wrong phase and make the two wake
        models incomparable, which is the whole point of the exercise.
        """
        mv = self.mvop
        tm = self.timings
        t0 = time.perf_counter()
        up, lo = mv.newton_side_data(
            phi, self.m_inf, self.upwind_c, self.m_crit, self.gamma_air,
            self.u_inf, self.m_cap, self.rho_floor, frozen=frozen)
        t1 = time.perf_counter()
        A = mv.assemble_matrix(
            rho_tilde=(up["rho_tilde"], lo["rho_tilde"]), closure="wake_ls",
            te_kutta="pressure", phi_ext=phi)
        t2 = time.perf_counter()
        R = A @ phi - self.b_base
        t3 = time.perf_counter()
        R[self.te_aux] = self.te_residual(phi)
        if tm is not None:
            t4 = time.perf_counter()
            tm["assembly"] += t2 - t1
            tm["residual"] += (t1 - t0) + (t3 - t2)
            tm["kutta"] += t4 - t3
        return A, R, up, lo

    def _remap(self, r, c, d):
        if not len(r):
            return r, c, d
        rr = self.row_map[r]
        keep = rr >= 0
        return rr[keep], c[keep], d[keep]

    def jacobian(self, A, up, lo, phi):
        mv = self.mvop
        Jc = (self.D_keep @ A).tocoo()
        parts_r = [Jc.row]
        parts_c = [Jc.col]
        parts_d = [Jc.data]
        # blended TE Kutta Jacobian: F*J_kutta (+ (1-F)*weld appended below)
        kr, kc, kd = te_kutta_jacobian_coo(mv, phi)
        if self.tt is not None and len(kr):
            kd = kd * self.f_of_row[kr]
        jac_parts = [(kr, kc, kd)]
        if self.tt is not None:
            jac_parts.append(te_weld_coo(mv.cm, self.tt))  # (1-F)*[aux-main]
        jac_parts += [
            self._remap(*newton_terms23_side_coo(mv.op, up, self.u_inf)),
            self._remap(*newton_terms23_side_coo(mv.op, lo, self.u_inf)),
        ]
        for r_, c_, d_ in jac_parts:
            if len(r_):
                parts_r.append(r_)
                parts_c.append(c_)
                parts_d.append(d_)
        return sp.coo_matrix(
            (np.concatenate(parts_d),
             (np.concatenate(parts_r), np.concatenate(parts_c))),
            shape=(self.n_total, self.n_total)).tocsr()


def solve_multivalued_newton(
    mvop,
    mesh,
    m_inf: float,
    alpha_deg: float = 1.25,
    u_inf: float = 1.0,
    gamma_air: float = 1.4,
    upwind_c: float = 1.5,
    m_crit: float = 0.95,
    m_cap: float = 3.0,
    rho_floor: float = 0.05,
    farfield: str = "neumann",
    vortex_center=(0.25, 0.0),
    phi_init: Optional[np.ndarray] = None,
    gamma_init: float = 0.0,
    n_seed: int = 40,
    seed_damping_theta: Optional[float] = 0.2,
    seed_omega_rho: float = 0.5,
    n_newton_max: int = 40,
    tol_residual: float = 1e-10,
    lam_min: float = 0.05,
    freeze: bool = True,
    freeze_tol: Optional[float] = None,
    freeze_refresh_max: int = 8,
    freeze_max_reverts: int = 3,
    freeze_max_clamped: int = 0,
    stall_window: int = 8,
    tol_residual_loose: Optional[float] = None,
    tol_residual_rel: Optional[float] = None,
    accept_on_stall: bool = False,
    verbose: bool = False,
    tip_taper: Optional[np.ndarray] = None,
    precond: Optional[str] = None,
    seed_precond: Optional[str] = None,
    linear_rtol: float = 1e-10,
    gmres_restart: int = 60,
    gmres_maxiter: int = 50,
    amg_rebuild_every: int = 2,
    direct_refactor_every: int = 1,
    direct_reuse_rtol: float = 1e-8,
) -> Dict[str, object]:
    """Newton solve on the level-set path (B6-Newton; design_track_b.md 5.5).

    Args:
        phi_init: warm-start extended state (e.g. the previous Mach level's).
            None -> seed with a short B6 Picard-LS solve.
        n_seed: Picard-LS seed iterations when phi_init is None (or to settle
            a warm start onto the current Mach before Newton).
        farfield: "neumann" (B6 transonic default -- no Gamma feedback, the
            far-field column vanishes), "vortex" (option a; adds the low-rank
            Gamma(z) far-field column -- handled by refreshing the RHS each
            Newton step, frozen within the step) or "freestream".
        freeze / freeze_tol (B15; N5 on the level-set path): freeze the
            PER-SIDE upwind selection between Newton steps. `freeze_tol=None`
            (default) = live re-selection every step, bit-identical to the
            pre-B15 behaviour. Set `freeze_tol` (recipe: 1e-6) to arm the
            machinery: once the live residual falls below it -- or the live
            iteration PLATEAUS (`live_stalled`) -- with 0 limited/floored, the
            per-side (upstream, branch) assignment is captured
            (`mvop.freeze_side_state`) and held. Within a frozen assignment the
            residual is smooth (no max(nu_e,nu_u) near-tie churn) so Newton
            recovers its terminal quadratic rate, and the frozen sweep applies
            NO floor on branches 0-2 => n_floored == 0 by design, which is what
            lets the 0-clamped convergence gate below actually fire at a shock.
            On frozen convergence the LIVE residual is re-evaluated (honesty):
            accept if it too is tight; else accept the assignment-discontinuity
            floor (two-cycle: live residual stopped improving) or refresh the
            assignment and continue, up to `freeze_refresh_max` times. Safety
            net: a frozen path that diverges (10x the freeze residual, or any
            limited/floored) REVERTS to the freeze point and unfreezes.
            `freeze=False` is a hard master-off switch.
        tol_residual: converged when max|R_free| < tol (with 0 limited/floored).
        tol_residual_loose / tol_residual_rel / accept_on_stall (B15, the
            G10.2 loose-acceptance set, for INTERMEDIATE Mach-ramp levels):
            additional accept routes -- an absolute loose tol, a relative drop
            vs the level's first residual, and accepting a measured plateau.
            All require >= 1 Newton step and never fire inside a frozen phase.
            Defaults (None/None/False) => bit-identical strict behaviour.
            `accept_reason` reports which route fired.
        precond (B11): inner linear solver -- None (default) = the bit-identical
            sparse-direct `spsolve`; **use "ilu"** for the iterative escape from
            the true-3D splu wall (factors the real fused Jacobian, converges).
            "amg" is wired but STALLS on this wake_ls operator (measured; §5.3),
            not the Newton escape. Inexact steps (on_fail="return") are absorbed
            by the existing backtracking line search. `seed_precond` (default
            None) is the same knob for the Picard warm-start solve.
            `linear_rtol=1e-10` keeps the Newton terminus unperturbed at the
            tol_residual gate. **"schur" (B14)** is the STRUCTURAL escape
            (solve/schur_ls.py): exact per-step elimination of the aux
            thin-strip block (`lu_aa = splu(J_aa)`, n_ext-sized,
            milliseconds) + GMRES on the reduced main-free operator
            preconditioned by AMG on the SPD Picard block, NO springs -- the
            surrogate's jump==0 bias is structurally absent. A stalled
            reduced GMRES falls back to a full fused spsolve in the same
            step (`n_schur_fallback`), never below per-step-direct
            robustness.
        direct_refactor_every (B12): lagged-LU direct-reuse, active only on the
            `precond=None` path. 1 (default) = refactor every step, the
            bit-identical `spsolve` behaviour. > 1 = refactor the LU every k-th
            step and drive the intermediate steps with GMRES preconditioned by
            the stale (exact) LU, converged to `direct_reuse_rtol` (1e-8). This
            is the medium/M6-scale escape (ILU diverges/factor-fails there),
            amortising the sparse factorization over k Newton steps. Ignored
            when precond is "ilu"/"amg"/"schur".

    Returns: dict with phi_ext, phi, gamma, cl_kj, te_jump, converged,
    residual_history, n_newton, n_limited, n_floored, mach2_max, nu_max,
    precond, n_gmres_total, n_gmres_stalled (B11 monitors), n_refactor (B12),
    n_schur_fallback (B14).
    """
    if farfield not in ("neumann", "vortex", "freestream"):
        raise ValueError(f"farfield={farfield!r} unknown")
    if len(mvop.cm.te_nodes) == 0:
        # A wake level-set that matches NO trailing-edge wall node carries no
        # Kutta condition at all: Gamma is unpinned, the TE aux rows are empty,
        # and the solve silently produces NaN (measured on ONERA M6 medium: a
        # hand-rolled TE polyline off by ~2e-4 in x -> 0 TE nodes -> 340k
        # limited cells, gamma = mean([]) = NaN). Fail loudly instead. Build the
        # polyline from the AUTHORITATIVE geometry (e.g. meshgen.wing3d.x_te),
        # whose endpoints are exact wall nodes.
        raise ValueError(
            "the wake level-set matched 0 TE nodes on this mesh: the TE "
            "polyline does not lie on the wall. Without TE nodes there is no "
            "Kutta condition and the solution is meaningless. Check that the "
            "polyline comes from the mesh's own geometry (e.g. "
            "pyfp3d.meshgen.wing3d.x_te / B_SEMI) and that wall_nodes was "
            "passed to CutElementMap.")
    if precond not in (None, "ilu", "amg", "schur"):
        raise ValueError(
            f"precond={precond!r} unknown (None|'ilu'|'amg'|'schur')")
    beta = float(np.sqrt(1.0 - m_inf**2))

    # --- far-field split (fixed across the solve) --------------------------
    if farfield == "neumann":
        ff_nodes = _farfield_split(mesh, alpha_deg, u_inf)[3]      # inflow
        ff_vals = freestream_phi(mesh.nodes[ff_nodes], alpha_deg, u_inf)
        b_base = _neumann_outlet_rhs(mesh, alpha_deg, u_inf, mvop.n_total)
    elif farfield == "freestream":
        # Dirichlet freestream on the WHOLE far field -- constant across the
        # solve (no vortex), so build the values ONCE. (Previously left None
        # here and only set for the vortex option, so a freestream Newton
        # solve wrote `phi_ext[ff_nodes] = None` -> non-finite; the committed
        # LS Newton runs all use neumann, so this path was never exercised.
        # B9 wing-body needs it: the fuselage blockage makes the Lopez
        # inlet-Dirichlet/outlet-Neumann outlet unbounded, and the 25-MAC
        # domain makes full-freestream Dirichlet accurate.)
        ff_nodes = np.unique(mesh.boundary_faces["farfield"])
        ff_vals = freestream_phi(mesh.nodes[ff_nodes], alpha_deg, u_inf)
        b_base = np.zeros(mvop.n_total)
    else:  # vortex: (re)built per step with the current gamma
        ff_nodes = np.unique(mesh.boundary_faces["farfield"])
        b_base = np.zeros(mvop.n_total)
        ff_vals = None

    # --- warm start (B15: the seed inherits the B13 lagged-LU, else it pays a
    # full factorization per outer -- ~11 min of pure seed at M6 medium) ------
    t_wall0 = time.perf_counter()
    timings = new_timings()
    t_seed0 = time.perf_counter()
    phi_ext, gamma = _seed_from_picard(
        mvop, mesh, m_inf, alpha_deg, u_inf, gamma_air, upwind_c, m_crit,
        farfield, phi_init, gamma_init, n_seed, seed_damping_theta,
        seed_omega_rho, precond=seed_precond,
        direct_refactor_every=direct_refactor_every)
    timings["seed"] = time.perf_counter() - t_seed0

    # the shared residual/Jacobian assembly (B15: one code path, FD-gated)
    sysm = LSNewtonSystem(mvop, m_inf, upwind_c, m_crit, gamma_air, u_inf,
                          m_cap, rho_floor, b_base=b_base, tip_taper=tip_taper,
                          timings=timings)
    te_aux = sysm.te_aux
    n_total = mvop.n_total
    is_dir = np.zeros(n_total, dtype=bool)
    is_dir[ff_nodes] = True
    free = np.flatnonzero(~is_dir)

    residual_history = []
    step_records = []
    converged = False
    M_pre = None
    n_gmres_total = 0
    n_gmres_stalled = 0
    lu_direct = None       # B12: the lagged direct LU (precond is None, k > 1)
    lu_age = 0
    n_refactor = 0
    n_schur_fallback = 0   # B14: full-spsolve fallbacks (precond == "schur")

    # --- B15 freeze state (N5 ported to the per-side level-set operator) ----
    frozen = None            # ((upstream, branch) upper, (.., ..) lower)
    freeze_point = None      # (phi_ext copy, residual at the freeze)
    freeze_clamped_at = 0    # #clamped cells present when the freeze was taken
    freeze_cooldown = 0
    n_freeze_refresh = 0
    n_freeze_reverts = 0
    n_assignment_stale = 0
    residual_unfrozen = None
    prev_live_residual = None
    assignment_cycle = False
    accept_reason = None
    r_best = np.inf
    armed = bool(freeze) and freeze_tol is not None

    _system = sysm.residual        # (A, R, up, lo) at (phi, selection)
    _jacobian = sysm.jacobian      # J at (A, up, lo, phi) -- same selection

    def _stale(fz_a, fz_b):
        """#elements whose (upstream, branch) assignment moved, both sides."""
        n = 0
        for (ua, ba), (ub, bb) in zip(fz_a, fz_b):
            moved = ua != ub
            n += int(np.count_nonzero(moved))
            n += int(np.count_nonzero((ba != bb) & ~moved))
        return n

    # A1: state ON ENTRY to a step + the cost OF that step, so a record is
    # closed once its step has been taken (same shape as solve/newton.py).
    prev_snap = snapshot(timings)
    prev_counts = [0, 0]           # n_gmres_total, n_refactor at last close

    def _close_step(rec):
        nonlocal prev_snap
        rec.update(step_delta(timings, prev_snap))
        rec["n_lin_iters"] = n_gmres_total - prev_counts[0]
        rec["n_refactor"] = n_refactor - prev_counts[1]
        rec["wall_cum_s"] = time.perf_counter() - t_wall0
        prev_snap = snapshot(timings)
        prev_counts[0] = n_gmres_total
        prev_counts[1] = n_refactor

    for it in range(n_newton_max):
        # far-field Dirichlet values (vortex option refreshes with gamma)
        if farfield == "vortex":
            ff_nodes_v, ff_vals = _farfield_main(mesh, alpha_deg, gamma, u_inf,
                                                 vortex_center, beta)
            phi_ext[ff_nodes_v] = ff_vals
            dir_nodes = ff_nodes_v
        else:
            phi_ext[ff_nodes] = ff_vals
            dir_nodes = ff_nodes

        # residual under the CURRENT selection (frozen or live)
        A, R, up, lo = _system(phi_ext, frozen)
        res = float(np.max(np.abs(R[free])))
        n_lim, n_flr = mvop.n_limited, mvop.n_floored
        residual_history.append(res)
        if step_records:
            _close_step(step_records[-1])
        # gamma is recomputed every step at the bottom of the loop but was
        # never stored before A1 -- the LS Newton had no circulation
        # trajectory at all, which is exactly what the bottleneck study needs.
        rec = {
            "i": it,
            "residual": res,
            "gamma_mean": gamma,
            "gamma_root": gamma,
            "n_limited": int(n_lim),
            "n_floored": int(n_flr),
            "frozen": frozen is not None,
            "n_lin_solves": 0,
            "lam": None,
        }
        step_records.append(rec)
        if verbose:
            print(f"    newton {it}: |R|={res:.3e} gamma={gamma:.5f} "
                  f"lim/flr={n_lim}/{n_flr}"
                  f"{' [frozen]' if frozen is not None else ''}")

        # --- freeze-revert safety net: the frozen sweep is floor-free, so a
        # stale assignment can walk the state into clamp land. Revert to the
        # freeze point, unfreeze, and cool down before re-arming.
        if frozen is not None and freeze_point is not None and (
                res > 10.0 * freeze_point[1]
                or n_lim + n_flr > max(freeze_max_clamped, freeze_clamped_at)):
            phi_ext = freeze_point[0].copy()
            frozen = None
            freeze_point = None
            freeze_cooldown = 5
            n_freeze_reverts += 1
            lu_direct = None          # the stale LU belongs to the frozen J
            if precond == "schur":
                M_pre = None          # B14: AMG built on the frozen rho
            r_best = np.inf           # new selection epoch: rescale fail-fast
            # FAIL-SAFE: a freeze that keeps diverging is worse than no freeze.
            # Disarm for the rest of the level and fall back to the live path
            # (measured: without this the solver freeze-reverts-re-arms forever
            # and never converges, while plain live Newton converges). The
            # freeze may only ever HELP; it can never cost convergence.
            if n_freeze_reverts >= freeze_max_reverts:
                armed = False
            continue
        if freeze_cooldown > 0:
            freeze_cooldown -= 1
        r_best = min(r_best, res)

        # Fail-fast: a wandering level is the continuation's cue to halve dm.
        # NOTE r_best is scoped to the CURRENT SELECTION EPOCH (it is reset on
        # every freeze / refresh / revert below). Residuals from a FROZEN phase
        # and from the LIVE phase are not comparable: the frozen system happily
        # descends to ~1e-11, and when the assignment is then refreshed the
        # residual legitimately returns to the live scale (~1e-3). Comparing
        # across that boundary makes the fail-fast read a 1e8x "blow-up" and
        # kill a perfectly healthy freeze-refresh cycle (measured on M6 medium
        # M0.65: frozen phase reached 1.5e-11, the first refresh was then
        # aborted by this very test).
        if res > 100.0 * r_best and res > 1e-3:
            break

        # --- acceptance (strict "tol" always; the loose routes are the G10.2
        # set, opt-in, >= 1 step, and never inside a frozen phase) -----------
        hist = residual_history
        # STALL = "no new best in `stall_window` steps". NOT the conforming
        # band test (within 4x-up / 2x-down of the best of the previous 6):
        # that is tuned for a trace that TRULY limit-cycles, and it MISFIRES on
        # the level-set path, whose live residual bounces +-2x WHILE STILL
        # DESCENDING (measured, medium M0-embedded M0.75: it fires at ~5e-6, the
        # freeze then locks a not-yet-settled assignment, the frozen step
        # diverges and reverts -- 4 reverts, no progress -- while the untouched
        # live path converges to 7.5e-12 in 54 steps). Requiring the ALL-TIME
        # best not to have improved by >10% over the window fires only on a real
        # plateau, never on a noisy descent.
        live_stalled = False
        if len(hist) > stall_window and res < 1e-3:
            recent_best = min(hist[-stall_window:])
            prior_best = min(hist[:-stall_window])
            live_stalled = recent_best > 0.9 * prior_best
        # Clamp precondition. LIVE phase: strictly 0 limited/floored. FROZEN
        # phase: n_floored counts branch==3, i.e. the cells that were clamped AT
        # THE FREEZE POINT -- that count is stale bookkeeping and never falls,
        # so testing it against 0 would refuse a machine-precision solution
        # forever (measured on M6 medium M0.70: the frozen phase drove |R| to
        # 7.8e-14 with the floored cell CLEARED in the live field, and the gate
        # still would not fire). The frozen phase therefore only has to be no
        # worse than it was at the freeze; the definitive clamp check is the
        # LIVE re-evaluation in the honesty branch below, which is strict.
        clamp_ok = (n_lim == 0 and n_flr == 0 if frozen is None
                    else n_lim + n_flr <= freeze_clamped_at)
        accept = None
        if clamp_ok:
            if res < tol_residual:
                accept = "tol"
            elif frozen is None and it > 0:
                if tol_residual_loose is not None and res < tol_residual_loose:
                    accept = "loose_tol"
                elif (tol_residual_rel is not None
                        and res <= tol_residual_rel * hist[0]):
                    accept = "rel_drop"
                elif accept_on_stall and live_stalled:
                    accept = "stall"

        if accept is not None:
            if frozen is None:
                converged = True
                accept_reason = accept
                break
            # The FROZEN system converged. HONESTY: re-evaluate the LIVE
            # residual -- the frozen solution is only a real solution of the
            # discrete equations if the live selection agrees there.
            _, R_l, _, _ = _system(phi_ext, None)
            residual_unfrozen = float(np.max(np.abs(R_l[free])))
            fz_new = sysm.freeze(phi_ext)
            n_assignment_stale = _stale(frozen, fz_new)
            if (residual_unfrozen < tol_residual
                    and mvop.n_limited == 0 and mvop.n_floored == 0):
                converged = True
                accept_reason = "tol"           # live is tight too: honest pass
                break
            if (prev_live_residual is not None
                    and residual_unfrozen > 0.5 * prev_live_residual):
                # the live residual stopped improving across refreshes: this is
                # the intrinsic assignment-discontinuity floor, not a defect.
                assignment_cycle = True
                converged = True
                accept_reason = "assignment_cycle"
                break
            if n_freeze_refresh >= freeze_refresh_max:
                converged = True
                accept_reason = "refresh_budget"
                break
            # REFRESH the active set and carry on from the live state
            prev_live_residual = residual_unfrozen
            n_freeze_refresh += 1
            frozen = fz_new
            freeze_point = (phi_ext.copy(), residual_unfrozen)
            lu_direct = None                    # assignment moved: J changed
            if precond == "schur":
                M_pre = None                    # B14: rho changed with it
            r_best = np.inf                     # new selection epoch
            continue

        # --- freeze trigger: a SETTLED, clean live state ---------------------
        # NOTE the deliberate difference from the conforming driver, which also
        # freezes on `live_stalled` (newton.py:573). MEASURED on this path
        # (medium M0-embedded M0.75, roadmap B15): the level-set live residual
        # bounces +-2x for tens of steps WHILE STILL DESCENDING (gamma travels
        # 0.183 -> 0.243 over that stretch -- slow progress in a stiff
        # direction, not a stall). Freezing there locks a still-MOVING upwind
        # assignment; the frozen step then diverges and reverts (3 reverts, no
        # convergence), while the same solve with the stall trigger removed
        # freezes late (|R| < 1e-6, assignment settled) and converges cleanly:
        # 53 steps, |R| 2.1e-12, 0 reverts, and the LIVE residual at the frozen
        # solution is 2.1e-12 too -- an honest pass landing on exactly the live
        # gamma (0.243305). So: arm on freeze_tol ONLY.
        # `freeze_max_clamped` (default 0 = the conforming N5 rule) relaxes the
        # 0-clamped precondition. Why it exists: at M6 medium M0.70 a SINGLE
        # persistently-floored cell (of 330k) blocks the freeze at ANY
        # freeze_tol, and the live iteration then sits in an exact period-7
        # limit cycle forever. The frozen sweep REPRESENTS a clamped cell
        # exactly (branch 3: nu=0, rho=rho_floor, s_e=s_u=0, a flat clamp with
        # zero derivative), so the 0-clamped precondition is stricter than the
        # machinery needs.
        #
        # !! READ THIS BEFORE QUOTING A CONVERGED STATE (self-caught, B15) !!
        # Setting freeze_max_clamped > 0 RELAXES THE CONVERGENCE SEMANTICS. The
        # `assignment_cycle` / `refresh_budget` accept routes below do NOT
        # re-check the clamp count, so a returned converged=True state MAY CARRY
        # up to that many clamped cells -- the M6 medium M0.84 solution does
        # carry 3 of 330k (which happens to match the Picard's <=3). Only the
        # strict `tol` route still demands live 0-limited/0-floored. And the
        # clamped cells do NOT "clear themselves": in the shipped ramp they
        # persist at every level from M0.70 up.
        #
        # P9/G9.1 is CITED, NOT RE-TESTED: its record is about permanently
        # LIMITED cells on the CONFORMING path (solve/newton.py, which still has
        # the hard 0-clamped rule); ours are mostly FLOORED. Whether relaxing the
        # precondition would unblock G9.1's conforming fine mesh is an UNTESTED
        # hypothesis.
        if (armed and frozen is None and freeze_cooldown == 0 and m_inf > 0.0
                and n_lim + n_flr <= freeze_max_clamped
                and res < freeze_tol):
            frozen = sysm.freeze(phi_ext)
            freeze_point = (phi_ext.copy(), res)
            freeze_clamped_at = n_lim + n_flr
            lu_direct = None
            if precond == "schur":
                M_pre = None                    # B14: new selection epoch
            r_best = np.inf                     # new selection epoch
            # The frozen flux equals the live flux AT the freeze state, so the
            # residual is unchanged -- but the JACOBIAN needs the frozen
            # sensitivities, so re-derive the side data under the freeze.
            A, R, up, lo = _system(phi_ext, frozen)

        with phase(timings, "assembly"):
            J = _jacobian(A, up, lo, phi_ext)

        rec["n_lin_solves"] = 1
        t_lin0 = time.perf_counter()
        # reduce Dirichlet and solve J_free d = -R_free (nonsymmetric)
        if precond is None:
            J_free = J[free][:, free].tocsc()
            rhs = -R[free]
            if direct_refactor_every <= 1:
                # default: sparse-direct spsolve every step (bit-identical)
                d_free = np.atleast_1d(spla.spsolve(J_free, rhs))
            else:
                # B12 lagged-LU (design_track_b.md 5.5 / roadmap B12): refactor
                # the LU only every direct_refactor_every steps; between
                # refactors drive the FRESH Jacobian with GMRES preconditioned
                # by the stale (exact) LU, converged to direct_reuse_rtol. The
                # stale LU is a near-exact preconditioner so a reuse step
                # converges in ~1-2 Krylov iterations. Because the level-set
                # system has NO Gamma coupling, this is a plain solve -- no
                # Woodbury (unlike solve/newton.py). A reuse step whose GMRES
                # fails falls back to a refactor + exact solve in the same
                # iteration, so robustness never drops below every-step-direct.
                d_free = None
                if lu_direct is not None and lu_age < direct_refactor_every:
                    A_op = spla.LinearOperator(
                        (free.size, free.size), matvec=J_free.dot)
                    M_lag = spla.LinearOperator(
                        (free.size, free.size), matvec=lu_direct.solve)
                    d_try, n_it, info = solve_gmres(
                        A_op, rhs, M=M_lag, rtol=direct_reuse_rtol,
                        restart=gmres_restart, maxiter=2, on_fail="return")
                    n_gmres_total += n_it
                    if info == 0:
                        d_free = np.atleast_1d(d_try)
                    else:
                        n_gmres_stalled += 1     # stale LU exhausted: refactor
                if d_free is None:
                    t_lu = time.perf_counter()
                    lu_direct = spla.splu(J_free)
                    dt_lu = time.perf_counter() - t_lu
                    timings["precond"] += dt_lu
                    t_lin0 += dt_lu        # the factorization is not a solve
                    lu_age = 0
                    n_refactor += 1
                    d_free = np.atleast_1d(lu_direct.solve(rhs))
                lu_age += 1
        else:
            # B11: preconditioned GMRES on the fused Jacobian. Inexact steps
            # (on_fail="return") are caught by the backtracking line search.
            # B14 "schur": exact per-step elimination of the aux thin-strip
            # block, GMRES on the reduced main-free operator with AMG on the
            # SPD Picard block (built from the CURRENT -- possibly frozen --
            # side densities, the conforming newton.py analogue).
            J_free = J[free][:, free].tocsr()
            t_pre = time.perf_counter()
            schur = None
            if precond == "ilu":
                M_pre = build_ilu_preconditioner(J_free)
            elif precond == "schur":
                schur = SchurReducedSystem(J_free, free, mvop.n_main,
                                           n_aux_expected=mvop.n_ext)
                if it % amg_rebuild_every == 0 or M_pre is None:
                    M_pre = main_block_preconditioner(
                        mvop, (up["rho_tilde"], lo["rho_tilde"]),
                        schur.main_free)
            elif it % amg_rebuild_every == 0 or M_pre is None:
                M_pre = _amg_surrogate_preconditioner(
                    mvop, J, (up["rho_tilde"], lo["rho_tilde"]), free)
            dt_pre = time.perf_counter() - t_pre
            timings["precond"] += dt_pre
            t_lin0 += dt_pre
            if precond == "schur":
                d_free, n_it, info = schur.solve(
                    -R[free], M=M_pre, rtol=linear_rtol,
                    restart=gmres_restart, maxiter=gmres_maxiter)
                n_gmres_total += n_it
                if info != 0:
                    # Safety net (the lagged-LU pattern): refactor-and-solve
                    # the full fused system in the same step, so robustness
                    # never drops below per-step-direct. The line search
                    # WOULD absorb an inexact step, but a silently degraded
                    # Krylov solve was B11's failure mode -- count it and
                    # solve exactly instead.
                    n_gmres_stalled += 1
                    n_schur_fallback += 1
                    d_free = np.atleast_1d(
                        spla.spsolve(J_free.tocsc(), -R[free]))
                d_free = np.atleast_1d(d_free)
            else:
                d_free, n_it, info = solve_gmres(
                    J_free, -R[free], M=M_pre, rtol=linear_rtol,
                    restart=gmres_restart, maxiter=gmres_maxiter,
                    on_fail="return")
                d_free = np.atleast_1d(d_free)
                n_gmres_total += n_it
                n_gmres_stalled += int(info != 0)
        timings["linsolve"] += time.perf_counter() - t_lin0

        # Backtracking line search on the free-DOF residual (safety only),
        # evaluated under the SAME selection the step was built with. BEST-OF-
        # TRIED, not smallest-lambda: an unconditional smallest-lambda
        # acceptance was measured (P8/N5) to sustain residual GROWTH.
        best = None
        lam = 1.0
        while True:
            phi_try = phi_ext.copy()
            phi_try[free] = phi_ext[free] + lam * d_free
            phi_try[dir_nodes] = phi_ext[dir_nodes]
            _, R_t, _, _ = _system(phi_try, frozen)
            m_t = float(np.max(np.abs(R_t[free])))
            if np.isfinite(m_t) and (best is None or m_t < best[1]):
                best = (phi_try, m_t)
            if np.isfinite(m_t) and m_t < res:
                break
            lam *= 0.5
            if lam < lam_min:
                break
        if best is None:
            raise RuntimeError(
                "LS Newton line search: every trial state was non-finite")
        rec["lam"] = float(lam)
        phi_ext = best[0]
        with phase(timings, "kutta"):
            gamma = float(np.mean(mvop.te_jump(phi_ext)))

    if step_records:
        _close_step(step_records[-1])
    # Final monitors from the LIVE selection: the reported n_limited/n_floored
    # must describe the real flow, never a frozen assignment's (floor-free)
    # bookkeeping -- a frozen finish always shows 0 floored BY DESIGN, so
    # reporting the frozen counts would be self-congratulatory.
    _system(phi_ext, None)
    with phase(timings, "residual"):
        mach2_max = float(np.max(mvop.element_mach2(phi_ext, m_inf, gamma_air,
                                                    u_inf)))
    finalize(timings, time.perf_counter() - t_wall0)
    return {
        "phi_ext": phi_ext,
        "phi": mvop.main_potential(phi_ext),
        "gamma": gamma,
        "cl_kj": 2.0 * gamma / u_inf,
        "te_jump": mvop.te_jump(phi_ext),
        "converged": converged,
        "residual_history": residual_history,
        "step_records": step_records,
        "timings": timings,
        "n_newton": len(residual_history),
        "n_limited": mvop.n_limited,
        "n_floored": mvop.n_floored,
        "nu_max": mvop.nu_max,
        "mach2_max": mach2_max,
        "farfield": farfield,
        "precond": precond,
        "n_gmres_total": n_gmres_total,
        "n_gmres_stalled": n_gmres_stalled,
        "n_refactor": n_refactor,
        "n_schur_fallback": n_schur_fallback,
        # B15 freeze monitors
        "accept_reason": accept_reason,
        "froze": frozen is not None,
        "n_freeze_refresh": n_freeze_refresh,
        "n_freeze_reverts": n_freeze_reverts,
        "residual_unfrozen": residual_unfrozen,
        "n_assignment_stale": n_assignment_stale,
        "assignment_cycle": assignment_cycle,
    }


# Level-set Newton transonic recipes (B15). The analogues of
# NEWTON_TRANSONIC_RECIPE / NEWTON_M6_RECIPE in tests/test_p8_newton.py.
B_NEWTON_TRANSONIC_DEFAULTS = dict(   # 2.5D fold zone: small steps, strict
    m_start=0.60, dm=0.025, dm_min=0.003, freeze_tol=1e-6,
    freeze_refresh_max=8, intermediate_tol=None,
)
# ★ MEASURED on ONERA M6 medium M0.84 (roadmap GB15.4). The two non-obvious
# entries are the whole game -- both were found the hard way:
#   freeze_tol=1e-3       must sit ABOVE the upwind-selection CHURN FLOOR, which
#                         RISES with Mach (<1e-6 at M0.60, 8.6e-6 at M0.65,
#                         2.7e-4 at M0.70). Set below it, the residual is thrown
#                         back by a selection flip before the freeze can ever arm
#                         and the ramp dies. Same law as "tol_residual must sit
#                         above the Picard plateau".
#   freeze_max_clamped=8  a SINGLE persistently-floored cell (of 330k) otherwise
#                         blocks the freeze at ANY freeze_tol -- the P9/G9.1 wall
#                         ("permanently-limited cells block the N5 freeze
#                         machinery" -- CITED, NOT re-tested; that record is about
#                         LIMITED cells on the conforming path, ours are FLOORED).
#                         The frozen sweep represents a clamped cell exactly
#                         (branch 3), so the precondition was stricter than the
#                         machinery needs. NOTE the clamped cells do NOT clear:
#                         the converged M0.84 state CARRIES 3 of 330k (matching
#                         the Picard's <=3), and this recipe therefore RELAXES the
#                         convergence semantics -- see solve_multivalued_newton.
B_NEWTON_M6_DEFAULTS = dict(          # true 3D, far from the fold: lagged-LU
    m_start=0.60, dm=0.05, dm_min=0.01, freeze_tol=1e-3,
    freeze_refresh_max=8, freeze_max_clamped=8, direct_refactor_every=1000,
)


def solve_multivalued_newton_transonic(
    mvop,
    mesh,
    m_target: float,
    alpha_deg: float = 1.25,
    m_start: float = 0.60,
    dm: float = 0.05,
    dm_min: float = 0.01,
    intermediate_tol: Optional[float] = None,
    freeze_tol: Optional[float] = 1e-6,
    upwind_c: float = 1.5,
    upwind_c_post=None,
    n_seed: int = 40,
    n_newton_max: int = 60,
    phi_init: Optional[np.ndarray] = None,
    gamma_init: float = 0.0,
    verbose: bool = False,
    **newton_kw,
) -> Dict[str, object]:
    """Mach-continuation ramp for the LEVEL-SET Newton (B15) -- the analogue of
    `solve/newton.py::solve_newton_transonic`, and the plateau-removing cure for
    the B6/B7 transonic Picard.

    Why it exists: the LS transonic Picard (`solve_multivalued_transonic`) parks
    its top Mach levels on the shock-position residual PLATEAU (P4/B6/N5 soft
    mode) and burns its whole outer budget there -- the measured 24-38 min of the
    M6-medium M0.84 workflow solve. Newton has no such soft mode; with the B15
    per-side freeze it reaches a strict, 0-clamped solution instead of a bounded
    stall.

    Ramp mechanics (mirrored from the conforming N5/N6 driver):
      * upward-only schedule (`continuation.mach_schedule`), warm-started from
        the last CONVERGED level only -- a failed level's state is never a seed;
      * a level that fails to converge HALVES the Mach step and inserts a new
        level BELOW it (down to `dm_min`); inserted retry levels always run
        STRICT (the halving cascade is the robustness fallback -- a loose retry
        cannot repair a bad seed);
      * `upwind_c_post`: an optional DECREASING list of dissipation constants
        (Lopez Table 4.13, 2.0 -> 1.6) re-solved strictly at the target Mach.

    intermediate_tol (the speed knob): loosen the STOPPING tolerance on every
    level except the final one, so intermediate levels take a few warm-started
    Newton steps and move on. None (default) = every level strict.

    ** The freeze stays ARMED on loose intermediate levels **, which is where
    this mask deliberately DIFFERS from the conforming one (`newton.py` sets
    freeze_tol=None there). Reason: the accept gate requires 0 limited/floored
    on EVERY route, loose ones included, and on a 0.60->0.84 ramp the shock forms
    MID-ramp -- those intermediate levels carry limiter cells and can only reach
    a 0-clamped accept THROUGH a freeze. Loosen the tolerance, keep the
    mechanism. (The conforming fold contraindication -- a loose level leaving an
    untracked Gamma seed, G10.2 -- has no analogue here: the level-set path has
    no Gamma DOF, Gamma is a solution mode carried inside phi_ext.)

    Returns the FINAL level's dict plus `levels` (per-level records),
    `mach_schedule`, and the honesty fields `target_reached` / `m_last_converged`
    (P13/G13.3 erratum: NEVER census a state whose ramp did not reach the target
    -- the returned state is then at the LAST level's Mach, not m_target).
    """
    from pyfp3d.solve.continuation import mach_schedule

    if m_target < m_start:
        raise ValueError(f"m_target={m_target} < m_start={m_start} "
                         "(the ramp is upward-only)")
    levels = list(mach_schedule(m_target, m_start=m_start, dm=dm))
    loose = [True] * (len(levels) - 1) + [False]     # final level is strict

    last_good = None            # (phi_ext, gamma) of the last CONVERGED level
    level_results = []
    res = None
    i = 0
    t_wall0 = time.perf_counter()

    def _run(m, lvl_kw):
        t0 = time.perf_counter()
        r = solve_multivalued_newton(
            mvop=mvop, mesh=mesh, m_inf=m, alpha_deg=alpha_deg,
            upwind_c=upwind_c, n_newton_max=n_newton_max, verbose=verbose,
            **lvl_kw)
        r["_wall_s"] = time.perf_counter() - t0
        return r

    def _record(m, r, tag):
        # A1: wall_s, timings, residual_history and step_records were dropped
        # here before -- the ramp reported one residual_norm per level.
        level_results.append({
            "m_inf": float(m), "tag": tag,
            "gamma": float(r["gamma"]), "cl_kj": float(r["cl_kj"]),
            "mach_max": float(np.sqrt(r["mach2_max"])),
            "n_newton": int(r["n_newton"]),
            "converged": bool(r["converged"]),
            "accept_reason": r["accept_reason"],
            "residual_norm": float(r["residual_history"][-1]),
            "n_limited": int(r["n_limited"]),
            "n_floored": int(r["n_floored"]),
            "froze": bool(r["froze"]),
            "n_freeze_refresh": int(r["n_freeze_refresh"]),
            "n_freeze_reverts": int(r["n_freeze_reverts"]),
            "n_refactor": int(r["n_refactor"]),
            "n_schur_fallback": int(r["n_schur_fallback"]),
            "wall_s": float(r["_wall_s"]),
            "timings": r["timings"],
            "residual_history": list(r["residual_history"]),
            "step_records": r["step_records"],
            "n_lin_iters": int(r["n_gmres_total"]),
            "n_lin_solves": int(sum(s.get("n_lin_solves", 0)
                                    for s in r["step_records"])),
        })

    while i < len(levels):
        m = levels[i]
        lvl_kw = dict(newton_kw)
        lvl_kw["freeze_tol"] = freeze_tol
        if intermediate_tol is not None and loose[i]:
            # loosen the stopping tol; the freeze stays armed (see docstring)
            lvl_kw.update(tol_residual_loose=intermediate_tol,
                          tol_residual_rel=1e-3, accept_on_stall=True)
        if last_good is not None:
            lvl_kw.update(phi_init=last_good[0], gamma_init=last_good[1],
                          n_seed=0)
        else:
            lvl_kw.update(phi_init=phi_init, gamma_init=gamma_init,
                          n_seed=n_seed)

        res = _run(m, lvl_kw)
        _record(m, res, "loose" if (intermediate_tol is not None
                                    and loose[i]) else "strict")
        if verbose:
            print(f"  [ramp] M={m:.4f} conv={res['converged']} "
                  f"reason={res['accept_reason']} n={res['n_newton']} "
                  f"|R|={res['residual_history'][-1]:.2e} "
                  f"gamma={res['gamma']:.5f}")

        if not res["converged"]:
            m_prev = m_start if i == 0 else levels[i - 1]
            dm_new = 0.5 * (m - m_prev)
            if i == 0 or dm_new < dm_min:
                break                      # ramp dies: target NOT reached
            levels.insert(i, m_prev + dm_new)
            loose.insert(i, False)         # retry levels run STRICT
            continue
        last_good = (res["phi_ext"], res["gamma"])
        i += 1

    target_reached = bool(last_good is not None and res is not None
                          and res["converged"]
                          and abs(levels[i - 1] - m_target) < 1e-12)

    # --- optional post-ramp dissipation staging (Lopez 2.0 -> 1.6), strict ----
    if target_reached and upwind_c_post:
        for c_post in upwind_c_post:
            lvl_kw = dict(newton_kw)
            lvl_kw.update(freeze_tol=freeze_tol, n_seed=0,
                          phi_init=last_good[0], gamma_init=last_good[1])
            r = solve_multivalued_newton(
                mvop=mvop, mesh=mesh, m_inf=m_target, alpha_deg=alpha_deg,
                upwind_c=c_post, n_newton_max=n_newton_max, verbose=verbose,
                **lvl_kw)
            _record(m_target, r, f"post_c{c_post}")
            if not r["converged"]:
                break
            res, last_good = r, (r["phi_ext"], r["gamma"])

    out = dict(res)
    out.pop("_wall_s", None)
    out["levels"] = level_results
    out["mach_schedule"] = [float(x) for x in levels]
    out["target_reached"] = target_reached
    out["m_last_converged"] = (float(levels[i - 1]) if last_good is not None
                               else None)
    # the Mach the RETURNED state actually lives at (== m_target only when
    # target_reached; otherwise it is the level the ramp died on)
    out["m_final"] = float(levels[min(i, len(levels) - 1)])
    out["timings_total"] = sum_timings([lv["timings"] for lv in level_results])
    out["wall_s"] = time.perf_counter() - t_wall0
    return out
