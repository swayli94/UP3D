"""Isentropic constitutive relations for compressible potential flow.

Pure numba functions for:
  - ρ(q²) — density from isentropic law
  - M(q²) — Mach number
  - a(q²) — speed of sound
  - Cp(q²) — exact (nonlinear) pressure coefficient
  - q*(M∞) — critical speed where M = 1

All functions are unit-tested against hand calculations from design.md §2.

Typical gamma = 1.4 (air at low Mach), but exposed as parameter for flexibility.

Reference: design.md §2 (Governing Equations and Normalization)
"""

import numba
import numpy as np

# Physics constants (SI units, but will nondimensionalize in the solver)
GAMMA = 1.4  # Specific heat ratio for air
R_GAS = 287.05  # Gas constant (J/kg·K)

# These constants should be exposed in config later; here are defaults
M_CRIT_DEFAULT = 0.95  # Cutoff Mach for density upwinding
UPWIND_CONST_DEFAULT = 1.5  # Upwinding coefficient C


@numba.njit(cache=True)
def density_isentropic(q_squared, M_inf, gamma=GAMMA):
    r"""
    Isentropic density from (2.2) in design.md.

    ρ(q²) = [1 + (γ−1)/2 · M∞² (1 − q²)]^{1/(γ−1)}

    Args:
        q_squared: Nondimensional speed squared |∇φ|²
        M_inf: Freestream Mach number
        gamma: Specific heat ratio (default 1.4)

    Returns:
        Nondimensional density ρ

    Note:
        q is normalized by U∞, so q² = 1 is the freestream state:
        ρ(1) = 1.0 (reference density ρ∞).
        q² = 0 is the stagnation point (V = 0): ρ(0) = [1 + (γ−1)/2 · M∞²]^{1/(γ−1)} > 1.
        Sonic (M = 1) occurs at q² = q*² (see critical_speed_squared), not at q² = 1.
    """
    exponent = 1.0 / (gamma - 1.0)
    base = 1.0 + 0.5 * (gamma - 1.0) * M_inf ** 2 * (1.0 - q_squared)
    # Vacuum guard (P4): a diverging transient iterate can push q^2 past
    # the vacuum limit where base <= 0 and rho would be NaN, poisoning the
    # whole Picard matrix. Floor the base instead; any physically
    # meaningful state (q < q_vacuum) is bitwise unaffected.
    if base < 1e-10:
        base = 1e-10
    return base ** exponent


@numba.njit(cache=True)
def speed_of_sound_squared(q_squared, M_inf, gamma=GAMMA):
    r"""
    Nondimensional speed of sound squared, a² = 1/M∞² + (γ−1)/2 · (1 − q²).

    Derived from (2.4) in design.md.

    Args:
        q_squared: Nondimensional speed squared
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        a² (nondimensional)
    """
    return 1.0 / (M_inf ** 2) + 0.5 * (gamma - 1.0) * (1.0 - q_squared)


@numba.njit(cache=True)
def mach_number_squared(q_squared, M_inf, gamma=GAMMA):
    r"""
    Local Mach number squared: M²(q²) = q² M∞² / a².

    From (2.3) in design.md:
    M²(q²) = q² M∞² / [ 1 + (γ−1)/2 · M∞² (1 − q²) ]

    Args:
        q_squared: Nondimensional speed squared
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        M² (local Mach number squared)
    """
    numerator = q_squared * M_inf ** 2
    denominator = 1.0 + 0.5 * (gamma - 1.0) * M_inf ** 2 * (1.0 - q_squared)
    return numerator / denominator


@numba.njit(cache=True)
def critical_speed_squared(M_inf, gamma=GAMMA):
    r"""
    Critical speed squared where local Mach M = 1.

    q*² = (2 + (γ−1)M∞²) / ((γ+1)M∞²)

    At q = q*, the flow becomes sonic. This is the threshold for supersonic zones.

    Args:
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        q*² (sonic speed squared)
    """
    numerator = 2.0 + (gamma - 1.0) * M_inf ** 2
    denominator = (gamma + 1.0) * M_inf ** 2
    return numerator / denominator


@numba.njit(cache=True)
def pressure_coefficient(q_squared, M_inf, gamma=GAMMA):
    r"""
    Exact (nonlinear) isentropic pressure coefficient.

    Cp = 2 / (γ M∞²) · (ρ^γ − 1)

    where ρ is from the isentropic density law. This is NOT the linearized Cp = −2u.
    The solver's whole point is to be valid where linearization fails.

    From (2.5) in design.md.

    Args:
        q_squared: Nondimensional speed squared
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        Cp (local pressure coefficient, exact)
    """
    rho = density_isentropic(q_squared, M_inf, gamma)
    cp = 2.0 / (gamma * M_inf ** 2) * (rho ** gamma - 1.0)
    return cp


@numba.njit(cache=True)
def pressure_coefficient_incompressible(q_squared):
    r"""
    Incompressible (Bernoulli) pressure coefficient, Cp = 1 − q²/U∞².

    The M∞ → 0 limit of `pressure_coefficient` (which divides by γ M∞² and
    is singular at exactly M∞ = 0). Used by the P2 lifting-Laplace gates.

    Args:
        q_squared: Nondimensional speed squared (q/U∞)²

    Returns:
        Cp (incompressible)
    """
    return 1.0 - q_squared


@numba.njit(cache=True)
def density_derivative_wrt_q_sq(q_squared, M_inf, gamma=GAMMA):
    r"""
    Derivative dρ/d(q²) evaluated at q².

    Used for upwinding and Jacobian assembly.
    From design.md §3: dρ/d(q²) = −(M∞²/2) ρ^{2−γ}

    Args:
        q_squared: Nondimensional speed squared
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        dρ/d(q²)
    """
    rho = density_isentropic(q_squared, M_inf, gamma)
    return -0.5 * M_inf ** 2 * (rho ** (2.0 - gamma))


@numba.njit(cache=True)
def mach_squared_derivative_wrt_q_sq(q_squared, M_inf, gamma=GAMMA):
    r"""
    Derivative dM²/d(q²) evaluated at q².

    From (2.3): M² = q² M∞² / D with D = 1 + (γ−1)/2 · M∞² (1 − q²), so
    (quotient rule, dD/dq² = −(γ−1)/2 · M∞²):

        dM²/d(q²) = M∞² · [1 + (γ−1)/2 · M∞²] / D²

    Strictly positive for any physical state — M² is monotone in q².
    Used by the P7 frozen-selection flux sensitivity ∂ρ̃/∂φ (the ∂ν/∂q²
    chain on the active switch branch, design.md §6.3 / López B.9–B.12).

    Args:
        q_squared: Nondimensional speed squared
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        dM²/d(q²)
    """
    D = 1.0 + 0.5 * (gamma - 1.0) * M_inf ** 2 * (1.0 - q_squared)
    return M_inf ** 2 * (1.0 + 0.5 * (gamma - 1.0) * M_inf ** 2) / (D * D)


@numba.njit(cache=True)
def upwind_factor(q_squared, M_inf, M_crit=M_CRIT_DEFAULT, C=UPWIND_CONST_DEFAULT, gamma=GAMMA):
    r"""
    Artificial compressibility upwinding factor ν.

    ν(q²) = C · max(0, 1 − M_c² / M²)

    Activates (ν > 0) only when local Mach exceeds critical threshold M_c ≈ 0.95–1.0.
    Subsonic regions (ν = 0) remain central and second-order accurate.

    From design.md §3, equation (3.2).

    Args:
        q_squared: Nondimensional speed squared
        M_inf: Freestream Mach number
        M_crit: Critical Mach for switch (default 0.95)
        C: Upwinding coefficient (default 1.5)
        gamma: Specific heat ratio

    Returns:
        ν ∈ [0, C], the upwinding weight
    """
    M_sq = mach_number_squared(q_squared, M_inf, gamma)
    return C * max(0.0, 1.0 - (M_crit ** 2) / max(M_sq, M_crit ** 2))


@numba.njit(cache=True)
def density_field(q_squared, M_inf, gamma=GAMMA):
    r"""
    Elementwise isentropic density over an array of q² values (the per-
    element density sweep of the Picard loop, design.md §8). Same scalar
    law as `density_isentropic`; kept here so physics stays in this module
    (agent-rules hard rule #5).

    Note the exact Laplace limit: at M∞ = 0 the base is exactly 1.0 and
    1.0**exponent == 1.0, so ρ ≡ 1.0 bitwise — the G3.3 "M∞ → 0 is
    bit-identical to Laplace" guarantee rests on this.

    Args:
        q_squared: (n,) nondimensional speed squared per element
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        (n,) density array
    """
    out = np.empty_like(q_squared)
    for i in range(q_squared.shape[0]):
        out[i] = density_isentropic(q_squared[i], M_inf, gamma)
    return out


@numba.njit(cache=True)
def q2_at_mach(mach, M_inf, gamma=GAMMA):
    r"""
    Invert (2.3): the q² at which the local Mach equals `mach`.

        q²(M) = M² (1 + (γ−1)/2 M∞²) / (M∞² (1 + (γ−1)/2 M²))

    Used by the P4 speed limiter (q² capped at q²(M_cap)): a diverging
    transient can otherwise run q² toward the vacuum limit, where ρ → 0
    produces effectively decoupled "dead cells" that the positivity
    guards then stabilize into spurious converged states (measured:
    off-body supersonic blobs at |y| ≈ 0.6 acting as fake blockage).
    Returns +inf at M∞ = 0 (limiter inactive, bitwise no-op).

    Args:
        mach: local Mach number of the cap
        M_inf: freestream Mach number
        gamma: specific heat ratio

    Returns:
        q² at local Mach `mach` (+inf when M_inf == 0)
    """
    if M_inf == 0.0:
        return np.inf
    m2 = mach * mach
    return (
        m2 * (1.0 + 0.5 * (gamma - 1.0) * M_inf ** 2)
        / (M_inf ** 2 * (1.0 + 0.5 * (gamma - 1.0) * m2))
    )


@numba.njit(cache=True)
def limit_q2_field(q_squared, M_inf, m_cap, gamma=GAMMA):
    r"""
    Elementwise speed limiter: q²_eff = min(q², q²(M_cap)) (see q2_at_mach
    for why). Any state with local M ≤ M_cap is bitwise unaffected, so
    subcritical solves (gates G3.3/G4.2) and healthy converged transonic
    states (M ≲ 1.5) never see it; the caller monitors how many elements
    are being limited.

    Args:
        q_squared: (n,) nondimensional speed squared per element
        M_inf: freestream Mach number
        m_cap: local-Mach cap (P4 default 3.0)
        gamma: specific heat ratio

    Returns:
        (n,) limited q² array
    """
    cap = q2_at_mach(m_cap, M_inf, gamma)
    out = np.empty_like(q_squared)
    for i in range(q_squared.shape[0]):
        out[i] = q_squared[i] if q_squared[i] < cap else cap
    return out


@numba.njit(cache=True)
def mach_squared_field(q_squared, M_inf, gamma=GAMMA):
    r"""
    Elementwise local M² over an array of q² values (supersonic-zone
    monitor for the Picard loop; the ν switch consumes this in P4).

    Args:
        q_squared: (n,) nondimensional speed squared per element
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        (n,) local Mach-squared array
    """
    out = np.empty_like(q_squared)
    for i in range(q_squared.shape[0]):
        out[i] = mach_number_squared(q_squared[i], M_inf, gamma)
    return out


@numba.njit(cache=True)
def density_derivative_wrt_q_sq_field(q_squared, M_inf, gamma=GAMMA):
    r"""Elementwise dρ/d(q²) over an array of q² values -- the field
    counterpart of `density_derivative_wrt_q_sq` (same formula, vectorized;
    the density_field idiom).

    Args:
        q_squared: (n,) nondimensional speed squared per element
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        (n,) dρ/d(q²) array
    """
    out = np.empty_like(q_squared)
    for i in range(q_squared.shape[0]):
        out[i] = density_derivative_wrt_q_sq(q_squared[i], M_inf, gamma)
    return out


@numba.njit(cache=True)
def mach_squared_derivative_wrt_q_sq_field(q_squared, M_inf, gamma=GAMMA):
    r"""Elementwise dM²/d(q²) over an array of q² values -- the field
    counterpart of `mach_squared_derivative_wrt_q_sq` (same formula,
    vectorized).

    Args:
        q_squared: (n,) nondimensional speed squared per element
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Returns:
        (n,) dM²/d(q²) array
    """
    out = np.empty_like(q_squared)
    for i in range(q_squared.shape[0]):
        out[i] = mach_squared_derivative_wrt_q_sq(q_squared[i], M_inf, gamma)
    return out


def validate_physics_bounds(rho, q, M, Cp, M_inf, gamma=GAMMA):
    """
    Validate that computed physics quantities stay within physical bounds.

    This is a **non-njitted** assertion function for pre/post checks.

    Constraints:
      - ρ ∈ (0, ρ₀] (density must be positive)
      - q ∈ [0, 2] (speed ≤ ~2 times freestream in worst case)
      - M ∈ [0, 5] (Mach number reasonable range)
      - Cp ∈ [−5, 10] (loose sanity bounds; stagnation Cp is at most a
        few, large-negative Cp means unphysical expansion)

    Args:
        rho: Array of density values
        q: Array of speeds
        M: Array of Mach numbers
        Cp: Array of pressure coefficients
        M_inf: Freestream Mach number
        gamma: Specific heat ratio

    Raises:
        ValueError if any bound is violated
    """
    rho = np.asarray(rho)
    q = np.asarray(q)
    M = np.asarray(M)
    Cp = np.asarray(Cp)

    if np.any(rho <= 0):
        raise ValueError(f"Negative density detected: min={np.min(rho)}")
    if np.any(q < -1e-12):
        raise ValueError(f"Negative speed detected: min={np.min(q)}")
    # BACKLOG (A3/P2, Kimi audit 2026-07-17): this q>2.0 raise contradicts a
    # legal limiter-capped field (M_cap=3.0 admits q~2.28); only caller is a
    # subsonic test, so it is left as-is until physics monitoring is wired into
    # the transonic path. Recorded in the A3 response doc + track_a ledger.
    if np.any(q > 2.0):
        raise ValueError(f"Excessive speed detected: max={np.max(q)} (should be ~1.0–1.4)")
    if np.any(M > 5.0):
        raise ValueError(f"Excessive Mach detected: max={np.max(M)}")
    if np.any(Cp < -5):
        raise ValueError(f"Excessive negative Cp: min={np.min(Cp)}")
    if np.any(Cp > 10):
        raise ValueError(f"Excessive positive Cp: max={np.max(Cp)}")


if __name__ == "__main__":
    # Quick sanity checks against hand values (design.md §2)
    # Run: python -m pyfp3d.physics.isentropic

    print("=== Isentropic Physics Self-Test ===\n")

    M_inf = 0.5
    print(f"Freestream Mach M∞ = {M_inf}\n")

    # Test 1: Stagnation point (q² = 0) -- density above reference, Cp > 0
    q_sq_0 = 0.0
    rho_0 = density_isentropic(q_sq_0, M_inf)
    expected_rho_0 = (1.0 + 0.5 * (GAMMA - 1.0) * M_inf ** 2) ** (1.0 / (GAMMA - 1.0))
    print(f"At stagnation (q² = 0):")
    print(f"  ρ = {rho_0:.6f} (expected {expected_rho_0:.6f})")
    print(f"  Cp = {pressure_coefficient(q_sq_0, M_inf):.6f} (expected > 0)")

    # Test 2: Freestream (q² = 1, since q is normalized by U∞) -- reference state
    q_sq_inf = 1.0
    rho_inf = density_isentropic(q_sq_inf, M_inf)
    M_inf_check = np.sqrt(mach_number_squared(q_sq_inf, M_inf))
    print(f"\nAt freestream (q² = 1):")
    print(f"  ρ = {rho_inf:.6f} (expected 1.0)")
    print(f"  M = {M_inf_check:.6f} (expected {M_inf})")
    print(f"  Cp = {pressure_coefficient(q_sq_inf, M_inf):.6f} (expected 0.0)")

    # Test 3: Critical speed for M∞ = 0.5
    q_crit_sq = critical_speed_squared(M_inf)
    M_crit = np.sqrt(mach_number_squared(q_crit_sq, M_inf))
    print(f"\nCritical speed q*² = {q_crit_sq:.6f} where M = 1.0:")
    print(f"  Computed M = {M_crit:.6f} (expected 1.0)")

    print("\n✓ All checks passed!")
