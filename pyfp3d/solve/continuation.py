"""
Mach continuation + transonic circulation strategy (roadmap P4;
design.md Sec 8 accelerations 3-4 and Sec 12.4 mitigation ladder).

The P4 transonic recipe, every ingredient evidence-driven (see the
roadmap G4.1 entry for the measurement trail):

1. MACH CONTINUATION: converge M_inf - dM first, restart from its
   (phi, Gamma). A cold transonic start diverges outright (measured:
   NaN at M0.80 direct).

2. FROZEN-GAMMA DENSITY SOLVES WITH PSEUDO-TIME at supercritical
   levels: the nested/interleaved Kutta iterations both fail at
   transonic conditions (nested: Gamma runaway 0.115 -> 4.99; damped
   interleaved: limit cycle -- and any relaxed fixed-point update
   provably diverges once |d target/d Gamma| >= 1). The density
   iteration at FIXED Gamma plus a pseudo-transient term A + D is
   stable and physical (measured M_max 1.36 at M0.80 where the
   undamped iteration blows through the M_cap limiter); the residual
   settles to a slowly-decaying, bounded shock-cell tail (see
   solve_subsonic_lifting docstring), with cl drifting < 1e-3 over the
   last hundreds of iterations -- the engineering-converged regime the
   P4 gates measure in. Newton (P6) is the designed cure for the tail.
   D defaults to the LOCAL `damping_theta * diag(A_free)` form
   (solve/picard.py docstring): the original GLOBAL
   diag(m_lumped/dtau) form (`pseudo_dt`) was calibrated on the coarse
   mesh but its damping ratio vs the operator weakens ~4x under
   refinement to medium and the medium G4.1 gate DIVERGED on it
   (roadmap G4.1, 2026-07-07) -- the local form is mesh- and
   shock-strength-independent by construction and is the one measured
   stable stepping M0.75 -> M0.80 (theta = 0.2). The two are mutually
   exclusive (solve_subsonic_lifting raises if both given); pass
   `pseudo_dt` explicitly to fall back to the retired global form.

3. GAMMA AS AN OUTER SCALAR ROOT-FIND (per station): impose Kutta
   OUTSIDE the density iteration -- secant on

       F(Gamma) = kutta_target(density-converged phi at fixed Gamma) - Gamma

   warm-starting each evaluation from the previous one (later
   evaluations cost a fraction of the first). tol_gamma is matched to
   the evaluation noise of the frozen-Gamma solves (~1e-4 = the
   measured cl drift scale), not to the subsonic 1e-8.
"""

from typing import Dict, Optional

import numpy as np

#: One parameter set for the whole G4.3 robustness sweep. The iteration
#: budgets are part of the set: at transonic the frozen-Gamma evals never
#: meet their tol_rho early-exit (the pseudo-time density lag is a bounded
#: tail, see module docstring), so every eval runs its FULL n_picard_eval
#: and total cost is exactly levels x n_evals x n_picard_eval (measured
#: coarse G4.1: 9664 = 64 seed + 12 x 800).
TRANSONIC_DEFAULTS = dict(
    upwind_c=1.5,
    m_crit=0.95,
    damping_theta=0.2,
    n_picard_seed=400,
    n_picard_eval=800,
    max_gamma_evals=12,
    tol_gamma=2e-4,
    omega_seed=0.9,
    forcing_seed=0.05,
)


def mach_schedule(m_target: float, m_start: float = 0.70, dm: float = 0.05):
    """Ascending Mach ramp ending exactly at m_target."""
    if m_target <= m_start:
        return [float(m_target)]
    ms = list(np.arange(m_start, m_target - 1e-12, dm))
    ms.append(m_target)
    return [float(m) for m in ms]


def solve_transonic_lifting(
    mesh_cut,
    wc,
    m_inf: float,
    alpha_deg: float,
    m_start: float = 0.70,
    dm: float = 0.05,
    upwind_c: float = TRANSONIC_DEFAULTS["upwind_c"],
    m_crit: float = TRANSONIC_DEFAULTS["m_crit"],
    upwind_weighted: bool = False,
    upwind_mode: str = "kernel",
    upwind_p_weight: float = 3.0,
    upwind_reach_frac: float = 1.0,
    upwind_sigma_s_frac: float = 0.35,
    upwind_sigma_p_frac: float = 0.35,
    upwind_nbr_depth: int = 3,
    damping_theta: Optional[float] = TRANSONIC_DEFAULTS["damping_theta"],
    pseudo_dt: Optional[float] = None,
    tol_gamma: float = TRANSONIC_DEFAULTS["tol_gamma"],
    n_picard_seed: int = TRANSONIC_DEFAULTS["n_picard_seed"],
    n_picard_eval: int = TRANSONIC_DEFAULTS["n_picard_eval"],
    max_gamma_evals: int = TRANSONIC_DEFAULTS["max_gamma_evals"],
    rtol: float = 1e-10,
    maxiter: int = 3000,
    u_inf: float = 1.0,
    farfield_spanwise_gamma: bool = False,
    omega_rho: float = 1.0,
    n_kutta_polish: int = 0,
    omega_rho_polish: float = 0.5,
    verbose: bool = False,
) -> Dict[str, object]:
    """
    Transonic lifting solve: Mach continuation; the first (subcritical)
    level seeds (phi, Gamma) with the P3 nested solve, every later level
    closes the Kutta condition by a per-station secant around
    frozen-Gamma pseudo-time density solves (module docstring).

    `damping_theta`/`pseudo_dt` select the pseudo-transient stabilizer
    for every supercritical density solve (see solve_subsonic_lifting
    docstring for the two forms); mutually exclusive (passing both
    raises), `damping_theta` is the default (local, mesh/shock-
    independent) -- fall back to the retired global mass-lumped form
    with `damping_theta=None, pseudo_dt=2e-3` (its coarse calibration).

    `rtol`/`maxiter` are the inner CG tolerance/cap forwarded to every
    density solve. The default 1e-10 is bit-identical to the pre-P5 path
    (2.5D gates). On the 3D M6 mesh 1e-10 is ~5x more inner CG work than
    the outer density fixed point (tol_rho 1e-6/1e-8) needs, so P5 runs
    with a looser `rtol=1e-7` (measured M_max identical to 5 digits, ~5.5x
    faster per iter) -- speeding the solver proper is P7's job.

    `farfield_spanwise_gamma` selects the 3D spanwise-tapered vortex far
    field (Gamma(z) per station, 0 at/beyond the sheet tip) instead of the
    span-uniform mean -- see solve_subsonic_lifting/farfield_dirichlet.
    Default False keeps the 2.5D paths bit-identical; P5 3D runs enable it.

    `omega_rho` under-relaxes the density update inside every frozen-Gamma
    eval density solve of the CONTINUATION secant loop (forwarded to
    solve_subsonic_lifting). Default 1.0 is bit-identical to the P4 path.
    NOTE omega_rho<1 inside the active-secant continuation is NOT the P5
    fix -- measured to destabilise the top-Mach secant (M_max blew to ~29);
    the continuation stays at 1.0 and the closure is fixed by the polish
    phase below instead.

    `n_kutta_polish` (P5, default 0 = pre-P5 behaviour) appends a
    fixed-Gamma Kutta-closure polish AFTER the Mach continuation. The
    per-station secant regularises by early stopping: on the 3D medium mesh
    it does not converge Gamma at the top Mach level (the steepest-Gamma
    station diverges if the secant is pushed -- measured M_max 29 at
    max_gamma_evals=16), and stopping early leaves that one station ~32%
    under-circulated, driving a spurious outboard-TE M>2 cluster (roadmap P5
    re-diagnosis 2026-07-08, T1-T4; INVESTIGATION_kutta_closure.md). The
    polish replaces the secant with a damped fixed-point iteration at the
    final Mach level -- apply the measured Kutta target, let the
    under-relaxed (`omega_rho_polish`, default 0.5) density catch up, repeat.
    Secant-free, so no overshoot; contractive toward the self-consistent
    Gamma. Measured on the P5 medium: 3-4 steps take M_max 5.2->2.0, 0
    floored/limited. Exact no-op on single-station (2.5D) meshes only if
    n_kutta_polish stays 0; leave it 0 for every non-3D path.

    Returns:
        dict: the final level's solve_subsonic_lifting result, plus
        gamma (root), kutta_mismatch, n_gamma_evals, mach_levels,
        n_picard_total, converged (physical field AND |F| < tol_gamma:
        M_max below the limiter cap and no floored/limited cells)
    """
    from pyfp3d.constraints.wake import kutta_targets
    from pyfp3d.solve.picard import solve_subsonic_lifting

    levels = mach_schedule(m_inf, m_start, dm)
    n_picard_total = 0

    # Subcritical seed level: P3 nested Kutta (stable there).
    r = solve_subsonic_lifting(
        mesh_cut, wc, m_inf=levels[0], alpha_deg=alpha_deg, u_inf=u_inf,
        omega=TRANSONIC_DEFAULTS["omega_seed"], upwind_c=upwind_c,
        m_crit=m_crit, upwind_weighted=upwind_weighted,
        upwind_mode=upwind_mode, upwind_p_weight=upwind_p_weight,
        upwind_reach_frac=upwind_reach_frac,
        upwind_sigma_s_frac=upwind_sigma_s_frac,
        upwind_sigma_p_frac=upwind_sigma_p_frac,
        upwind_nbr_depth=upwind_nbr_depth,
        tol_rho=1e-6, n_picard_max=n_picard_seed,
        forcing=TRANSONIC_DEFAULTS["forcing_seed"], rtol=rtol, maxiter=maxiter,
        farfield_spanwise_gamma=farfield_spanwise_gamma,
    )
    phi, gamma = r["phi"], r["gamma"].copy()
    n_picard_total += r["n_picard"]
    mismatch = 0.0
    n_evals_last = 1
    if verbose:
        print(f"  M={levels[0]:.3f} seed: n={r['n_picard']} "
              f"gamma={gamma[0]:.5f} Mmax={np.sqrt(r['mach2_max']):.3f}")

    def _density_solve(m, g, phi_seed, omr=omega_rho):
        return solve_subsonic_lifting(
            mesh_cut, wc, m_inf=m, alpha_deg=alpha_deg, u_inf=u_inf,
            omega=1.0, upwind_c=upwind_c, m_crit=m_crit,
            upwind_weighted=upwind_weighted, upwind_mode=upwind_mode,
            upwind_p_weight=upwind_p_weight,
            upwind_reach_frac=upwind_reach_frac,
            upwind_sigma_s_frac=upwind_sigma_s_frac,
            upwind_sigma_p_frac=upwind_sigma_p_frac,
            upwind_nbr_depth=upwind_nbr_depth,
            damping_theta=damping_theta, pseudo_dt=pseudo_dt,
            pseudo_dt_max_ratio=1.0,
            tol_rho=1e-8, n_picard_max=n_picard_eval, forcing=0.0,
            phi_init=phi_seed, gamma_fixed=g, rtol=rtol, maxiter=maxiter,
            farfield_spanwise_gamma=farfield_spanwise_gamma,
            omega_rho=omr,
        )

    for m in levels[1:]:
        g = gamma.copy()
        r = _density_solve(m, g, phi)
        n_picard_total += r["n_picard"]
        F = kutta_targets(r["phi"], wc) - g
        prev_g = prev_F = None
        n_evals = 1
        while float(np.max(np.abs(F))) > tol_gamma and n_evals < max_gamma_evals:
            if prev_g is None:
                g_new = g + 0.7 * F
            else:
                dg = g - prev_g
                slope = np.where(
                    np.abs(dg) > 1e-12,
                    (F - prev_F) / np.where(dg == 0, 1, dg),
                    0.0,
                )
                g_new = np.where(
                    np.abs(slope) > 1e-3,
                    g - F / np.where(slope == 0, 1, slope),
                    g + 0.7 * F,
                )
            prev_g, prev_F = g, F
            g = g_new
            r = _density_solve(m, g, r["phi"])
            n_picard_total += r["n_picard"]
            F = kutta_targets(r["phi"], wc) - g
            n_evals += 1
        gamma = g
        phi = r["phi"]
        mismatch = float(np.max(np.abs(F)))
        n_evals_last = n_evals
        if verbose:
            print(f"  M={m:.3f}: {n_evals} gamma evals, |F|={mismatch:.2e}, "
                  f"gamma={gamma[0]:.5f}, Mmax={np.sqrt(r['mach2_max']):.3f}, "
                  f"res={r['residual_history'][-1]:.1e}")

    # Fixed-Gamma Kutta-closure polish (see docstring `n_kutta_polish`):
    # a secant-free, damped fixed-point closure at the final Mach level that
    # drives the station(s) the continuation left under-converged to their
    # self-consistent Gamma. No-op when n_kutta_polish == 0.
    for i in range(n_kutta_polish):
        gamma = kutta_targets(r["phi"], wc)
        r = _density_solve(levels[-1], gamma, r["phi"], omr=omega_rho_polish)
        n_picard_total += r["n_picard"]
        mismatch = float(np.max(np.abs(kutta_targets(r["phi"], wc) - gamma)))
        n_evals_last = i + 1
        if verbose:
            print(f"  polish {i + 1}/{n_kutta_polish}: |F|={mismatch:.2e}, "
                  f"gamma={gamma[0]:.5f}, Mmax={np.sqrt(r['mach2_max']):.3f}, "
                  f"floored/limited={r['n_floored']}/{r['n_limited']}")

    physical = (
        r["mach2_max"] < 9.0  # below the m_cap=3 limiter ceiling
        and r["n_limited"] == 0
        and r["n_floored"] == 0
    )
    out = dict(r)
    out["gamma"] = gamma
    out["kutta_mismatch"] = mismatch
    out["n_gamma_evals"] = n_evals_last
    out["mach_levels"] = levels
    out["n_picard_total"] = n_picard_total
    out["kutta_converged"] = mismatch < tol_gamma
    out["converged"] = bool(physical and mismatch < tol_gamma)
    return out
