"""GS1.1 Case B: the Laval-nozzle shock bench -- a shock whose position is
UNIQUELY determined by the boundary data, with an exact quasi-1-D reference.

Why Case B replaced Case A (recorded in
phases/p2/docs/dev_phase_two/20260728-1640-s1-shock-bench.md §4): in a CONSTANT-area duct
every shock position is an exact solution AND a smooth (shock-free) solution
satisfies the same two Dirichlet conditions -- measured: the Newton, started
exactly on the shocked solution, converges to residual 1e-16 on the smooth
sonic branch (u = 1.2143 vs u* = 1.2119). The bench could not distinguish
"operator moved the shock" from "another legitimate solution exists", so it
cannot answer GS1.1's question.

A converging-diverging duct removes the ambiguity:

  * the flow can only reach supersonic through a SONIC THROAT, which fixes the
    mass flux mdot = f(u*) H_throat (f(u) = rho(u^2) u);
  * for that mdot, u(x) is one of the two roots of f(u) = mdot / H(x);
  * the solution family is then parameterised by the shock position x_s alone,
    and Delta_phi = integral u dx is STRICTLY MONOTONE in x_s -- so prescribing
    phi at inlet and outlet picks exactly one x_s;
  * the fully-subsonic (unchoked) family covers a DISJOINT, lower range of
    Delta_phi (verified numerically by `verify_uniqueness`), so the shocked
    solution is the only one satisfying the imposed Delta_phi.

Geometry (slender on purpose, so the quasi-1-D reference is accurate to
O((dH/dx)^2) ~ 1e-4):

    H(x) = H_t (1 + a_in  ((x_t - x)/x_t)^2)          x <  x_t
         = H_t (1 + a_out ((x - x_t)/(L - x_t))^2)    x >= x_t

with L = 20, x_t = 6, H_t = 0.5, a_in = 0.6, a_out = 0.361 -- the last chosen
so that the target shock at x_s = 12 sits at M ~ 1.3, i.e. the strength of a
real airfoil shock and inside the full-potential validity envelope.
"""

import numpy as np

from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.physics.isentropic import (GAMMA, critical_speed_squared,
                                       density_isentropic,
                                       mach_number_squared)

# geometry
LENGTH = 20.0
X_T = 6.0
H_T = 0.5
A_IN = 0.6
A_OUT = 0.361
X_S_TARGET = 12.0


def height(x):
    """Duct half-height H(x) (the 2-D channel spans y in [0, H])."""
    x = np.asarray(x, dtype=np.float64)
    conv = H_T * (1.0 + A_IN * ((X_T - x) / X_T) ** 2)
    div = H_T * (1.0 + A_OUT * ((x - X_T) / (LENGTH - X_T)) ** 2)
    return np.where(x < X_T, conv, div)


# ---------------------------------------------------------------------------
# exact quasi-1-D solution
# ---------------------------------------------------------------------------

def _f(u, m_inf, gamma):
    return density_isentropic(u * u, m_inf, gamma) * u


def _roots(target, m_inf, gamma, u_star):
    """Both roots of f(u) = target (target <= f(u_star)), by bisection."""
    fs = _f(u_star, m_inf, gamma)
    t = min(target, fs)
    lo, hi = 1e-12, u_star
    for _ in range(200):                       # subsonic root
        mid = 0.5 * (lo + hi)
        if _f(mid, m_inf, gamma) < t:
            lo = mid
        else:
            hi = mid
    u_sub = 0.5 * (lo + hi)
    lo, hi = u_star, u_star * 4.0
    for _ in range(200):                       # supersonic root
        mid = 0.5 * (lo + hi)
        if _f(mid, m_inf, gamma) > t:
            lo = mid
        else:
            hi = mid
    u_sup = 0.5 * (lo + hi)
    return u_sub, u_sup


def exact_solution(m_inf: float, x_s: float = X_S_TARGET,
                   gamma: float = GAMMA, n_quad: int = 40001):
    """Exact choked quasi-1-D solution with the shock at x_s.

    Returns dict with u_star, mdot, the sampled profile (x, u, mach) and
    phi_of_x (a callable interpolant), plus delta_phi and the inlet/exit Mach.
    """
    u_star = float(np.sqrt(critical_speed_squared(m_inf, gamma)))
    mdot = float(_f(u_star, m_inf, gamma) * H_T)
    xs = np.linspace(0.0, LENGTH, n_quad)
    H = height(xs)
    u = np.empty_like(xs)
    for i, (xi, Hi) in enumerate(zip(xs, H)):
        u_sub, u_sup = _roots(mdot / Hi, m_inf, gamma, u_star)
        u[i] = u_sup if (X_T <= xi < x_s) else u_sub
    phi = np.concatenate(([0.0], np.cumsum(0.5 * (u[1:] + u[:-1]) * np.diff(xs))))
    mach = np.sqrt(mach_number_squared(u * u, m_inf, gamma))
    return dict(u_star=u_star, mdot=mdot, x=xs, u=u, mach=mach, phi=phi,
                delta_phi=float(phi[-1]), x_s=float(x_s),
                m_inlet=float(mach[0]), m_shock_up=float(np.max(mach)),
                phi_of_x=lambda q: np.interp(np.asarray(q, dtype=np.float64),
                                             xs, phi))


def verify_uniqueness(m_inf: float, x_s: float = X_S_TARGET,
                      gamma: float = GAMMA):
    """Check that the imposed Delta_phi cannot be met by an unchoked
    (fully subsonic) solution: the subsonic family's Delta_phi is maximal at
    the critical mass flux, and that maximum must be BELOW the shocked one.

    Returns (delta_phi_shocked, delta_phi_subsonic_max, unique: bool).
    """
    ex = exact_solution(m_inf, x_s, gamma)
    u_star = ex["u_star"]
    xs = np.linspace(0.0, LENGTH, 20001)
    H = height(xs)
    mdot_crit = ex["mdot"]
    u = np.array([_roots(mdot_crit / Hi, m_inf, gamma, u_star)[0]
                  for Hi in H])                      # subsonic branch only
    dphi_sub = float(np.trapezoid(u, xs)) if hasattr(np, "trapezoid") \
        else float(np.trapz(u, xs))
    return ex["delta_phi"], dphi_sub, bool(ex["delta_phi"] > dphi_sub)


# ---------------------------------------------------------------------------
# mesh
# ---------------------------------------------------------------------------

def nozzle_mesh(nx: int, ny: int, dz: float = 0.1, jitter: float = 0.0,
                seed: int = 7):
    """Single-layer 2.5-D tet mesh of the nozzle; `jitter` perturbs interior
    2-D nodes (the irregular-mesh leg). Groups: inlet, outlet, wall."""
    xs = np.linspace(0.0, LENGTH, nx + 1)
    Hs = height(xs)
    pts = np.empty(((nx + 1) * (ny + 1), 2), dtype=np.float64)

    def nid(i, j):
        return i * (ny + 1) + j

    for i in range(nx + 1):
        for j in range(ny + 1):
            pts[nid(i, j)] = (xs[i], Hs[i] * j / ny)

    if jitter > 0.0:
        rng = np.random.default_rng(seed)
        for i in range(1, nx):
            for j in range(1, ny):
                hx = 0.5 * (xs[i + 1] - xs[i - 1])
                hy = Hs[i] / ny
                k = nid(i, j)
                pts[k, 0] += jitter * hx * (rng.random() - 0.5)
                pts[k, 1] += jitter * hy * (rng.random() - 0.5)

    tris = []
    for i in range(nx):
        for j in range(ny):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            if (i + j) % 2 == 0:
                tris += [[a, b, c], [a, c, d]]
            else:
                tris += [[a, b, d], [b, c, d]]
    tris = np.asarray(tris, dtype=np.int64)

    inlet = np.array([[nid(0, j), nid(0, j + 1)] for j in range(ny)],
                     dtype=np.int64)
    outlet = np.array([[nid(nx, j), nid(nx, j + 1)] for j in range(ny)],
                      dtype=np.int64)
    wall = np.array([[nid(i, 0), nid(i + 1, 0)] for i in range(nx)]
                    + [[nid(i, ny), nid(i + 1, ny)] for i in range(nx)],
                    dtype=np.int64)
    return extrude_single_layer(
        pts, tris, {"inlet": inlet, "outlet": outlet, "wall": wall},
        dz=dz, name=f"nozzle_{nx}x{ny}")


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------

def shock_from_profile(xc, ux, u_star, n_bins):
    """Sonic crossing of the bin-averaged u(x). u* is a CONSTANT in this
    normalisation, so the crossing is unambiguous."""
    edges = np.linspace(0.0, LENGTH, n_bins + 1)
    idx = np.clip(np.digitize(xc, edges) - 1, 0, n_bins - 1)
    ub = np.full(n_bins, np.nan)
    for b in range(n_bins):
        m = idx == b
        if m.any():
            ub[b] = ux[m].mean()
    xb = 0.5 * (edges[:-1] + edges[1:])
    ok = ~np.isnan(ub)
    xb, ub = xb[ok], ub[ok]
    above = ub > u_star
    # the LAST supersonic-to-subsonic crossing (the shock; the throat crossing
    # is subsonic-to-supersonic and is skipped by the direction test)
    cross = np.where(above[:-1] & ~above[1:])[0]
    if cross.size == 0:
        return float("nan"), int(np.count_nonzero(above)), xb, ub
    k = int(cross[-1])
    u0, u1 = ub[k], ub[k + 1]
    t = (u0 - u_star) / (u0 - u1) if u0 != u1 else 0.5
    return (float(xb[k] + t * (xb[k + 1] - xb[k])),
            int(np.count_nonzero(above)), xb, ub)
