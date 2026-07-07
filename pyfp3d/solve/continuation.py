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
    damping_theta: Optional[float] = TRANSONIC_DEFAULTS["damping_theta"],
    pseudo_dt: Optional[float] = None,
    tol_gamma: float = TRANSONIC_DEFAULTS["tol_gamma"],
    n_picard_seed: int = TRANSONIC_DEFAULTS["n_picard_seed"],
    n_picard_eval: int = TRANSONIC_DEFAULTS["n_picard_eval"],
    max_gamma_evals: int = TRANSONIC_DEFAULTS["max_gamma_evals"],
    u_inf: float = 1.0,
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
        m_crit=m_crit, tol_rho=1e-6, n_picard_max=n_picard_seed,
        forcing=TRANSONIC_DEFAULTS["forcing_seed"],
    )
    phi, gamma = r["phi"], r["gamma"].copy()
    n_picard_total += r["n_picard"]
    mismatch = 0.0
    n_evals_last = 1
    if verbose:
        print(f"  M={levels[0]:.3f} seed: n={r['n_picard']} "
              f"gamma={gamma[0]:.5f} Mmax={np.sqrt(r['mach2_max']):.3f}")

    def _density_solve(m, g, phi_seed):
        return solve_subsonic_lifting(
            mesh_cut, wc, m_inf=m, alpha_deg=alpha_deg, u_inf=u_inf,
            omega=1.0, upwind_c=upwind_c, m_crit=m_crit,
            damping_theta=damping_theta, pseudo_dt=pseudo_dt,
            pseudo_dt_max_ratio=1.0,
            tol_rho=1e-8, n_picard_max=n_picard_eval, forcing=0.0,
            phi_init=phi_seed, gamma_fixed=g,
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
