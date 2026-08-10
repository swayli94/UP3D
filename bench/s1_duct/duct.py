"""GS1.1 shock-operator bench: exact quasi-1-D duct states, meshes, solver,
shock extraction.

The bench case is a straight (constant-area) 2-D channel carrying a single
stationary normal shock, meshed as a single-layer 2.5-D tet mesh -- the same
element idiom as the production airfoil meshes.

Why this case (docs/dev_phase_two/20260728-1640-s1-shock-bench.md §1.0):

  * In full potential a shock is an ISENTROPIC, mass-conserving jump, so in a
    constant-area duct rho(u^2) u = mdot has two roots u_sup > u* > u_sub and

        phi_exact(x) = u_sup x                      x <  x_s
                     = u_sup x_s + u_sub (x - x_s)  x >= x_s

    is an EXACT weak solution for any x_s -- piecewise linear, hence also
    exactly representable by the P1 field away from the shock cell.
  * In both uniform states rho_e == rho_upstream, so the artificial-density
    term nu (rho_e - rho_up) vanishes identically: the background dissipation
    is exactly zero and the operator acts ONLY in the shock cells. Any measured
    shock motion is therefore attributable to the shock operator alone.

Everything numerical here goes through the SHIPPED kernels (PicardOperator,
UpwindOperator) so the bench measures the production discretisation, not a
re-implementation. The only bench-local code is the outer Newton loop, and its
job is only to reach a state with ||R||_inf ~ 1e-12 -- at that point the state
is a solution of the shipped discrete equations regardless of how it was found.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pyfp3d.kernels.entropy import EntropyOperator
from pyfp3d.kernels.jacobian import PicardOperator
from pyfp3d.kernels.upwind import UpwindOperator
from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.physics.isentropic import (GAMMA, critical_speed_squared,
                                       density_field, density_isentropic,
                                       mach_number_squared, q2_at_mach)


# ---------------------------------------------------------------------------
# 1. exact states
# ---------------------------------------------------------------------------

def duct_states(m_inf: float, m_sup: float, gamma: float = GAMMA):
    """The two roots of rho(u^2) u = mdot for a shock with upstream Mach
    `m_sup`, in the solver's normalisation (q^2 = 1 is freestream).

    Returns dict(u_sup, u_sub, mdot, q2_sup, q2_sub, q2_star, m_sub).
    """
    q2_sup = q2_at_mach(m_sup, m_inf, gamma)
    u_sup = float(np.sqrt(q2_sup))
    mdot = float(density_isentropic(q2_sup, m_inf, gamma) * u_sup)
    q2_star = float(critical_speed_squared(m_inf, gamma))

    # subsonic root: solve f(u) = rho(u^2) u = mdot on (0, u*) by bisection
    # (f is strictly increasing there, f(0) = 0 < mdot).
    def f(u):
        return float(density_isentropic(u * u, m_inf, gamma) * u)

    lo, hi = 1e-12, float(np.sqrt(q2_star))
    assert f(hi) > mdot, "mdot above the sonic maximum -- no subsonic root"
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < mdot:
            lo = mid
        else:
            hi = mid
    u_sub = 0.5 * (lo + hi)
    return dict(u_sup=u_sup, u_sub=u_sub, mdot=mdot, q2_sup=q2_sup,
                q2_sub=u_sub ** 2, q2_star=q2_star,
                m_sub=float(np.sqrt(mach_number_squared(u_sub ** 2, m_inf,
                                                        gamma))))


def phi_exact(x, x_s: float, st: dict):
    """Exact piecewise-linear potential of the stationary shock at x_s."""
    x = np.asarray(x, dtype=np.float64)
    return np.where(x < x_s, st["u_sup"] * x,
                    st["u_sup"] * x_s + st["u_sub"] * (x - x_s))


# ---------------------------------------------------------------------------
# 2. meshes
# ---------------------------------------------------------------------------

def duct_mesh(nx: int, ny: int, length: float = 4.0, height: float = 1.0,
              dz: float = 0.1, jitter: float = 0.0, seed: int = 0):
    """Single-layer 2.5-D tet mesh of a straight duct.

    `jitter` > 0 displaces INTERIOR 2-D nodes by up to jitter*h in a
    deterministic pseudo-random pattern -- the "irregular mesh" leg, which is
    how the production meshes actually look (sliver prism-split tets, edge
    lengths spanning decades). Boundary nodes are never moved, so the geometry
    and the boundary groups are identical between the regular and irregular
    legs.

    Boundary groups: "inlet" (x = 0), "outlet" (x = L), "wall" (y = 0, y = H).
    """
    xs = np.linspace(0.0, length, nx + 1)
    ys = np.linspace(0.0, height, ny + 1)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel()])

    def nid(i, j):
        return i * (ny + 1) + j

    if jitter > 0.0:
        rng = np.random.default_rng(seed)
        hx, hy = length / nx, height / ny
        for i in range(1, nx):
            for j in range(1, ny):
                k = nid(i, j)
                pts[k, 0] += jitter * hx * (rng.random() - 0.5)
                pts[k, 1] += jitter * hy * (rng.random() - 0.5)

    tris = []
    for i in range(nx):
        for j in range(ny):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            # alternate the diagonal so the triangulation has no global bias
            if (i + j) % 2 == 0:
                tris += [[a, b, c], [a, c, d]]
            else:
                tris += [[a, b, d], [b, c, d]]
    tris = np.asarray(tris, dtype=np.int64)

    inlet = np.array([[nid(0, j), nid(0, j + 1)] for j in range(ny)], dtype=np.int64)
    outlet = np.array([[nid(nx, j), nid(nx, j + 1)] for j in range(ny)],
                      dtype=np.int64)
    wall = np.array(
        [[nid(i, 0), nid(i + 1, 0)] for i in range(nx)]
        + [[nid(i, ny), nid(i + 1, ny)] for i in range(nx)], dtype=np.int64)

    return extrude_single_layer(
        pts, tris, {"inlet": inlet, "outlet": outlet, "wall": wall},
        dz=dz, name=f"duct_{nx}x{ny}")


# ---------------------------------------------------------------------------
# 3. the shipped discrete operator + a Newton driver
# ---------------------------------------------------------------------------

class DuctSystem:
    """Shipped-kernel residual / Jacobian of the duct problem with Dirichlet
    phi at inlet and outlet (walls take the natural zero-flux BC)."""

    def __init__(self, mesh, m_inf: float, upwind_c: float = 1.5,
                 m_crit: float = 0.95, rho_floor: float = 0.05,
                 gamma: float = GAMMA, entropy: bool = False,
                 entropy_refresh_max: int = 8):
        self.mesh = mesh
        self.nodes = np.ascontiguousarray(mesh.nodes, dtype=np.float64)
        self.elements = np.ascontiguousarray(mesh.elements)
        self.op = PicardOperator(self.nodes, self.elements)
        self.upw = UpwindOperator(self.nodes, self.elements, weighted=False)
        self.m_inf, self.C, self.m_crit = m_inf, upwind_c, m_crit
        self.rho_floor, self.gamma = rho_floor, gamma
        # GS1b.3: the entropy-corrected density, with the SAME frozen-sigma
        # semantics as pyfp3d.solve.newton (sigma held over a step, refreshed
        # between steps, refresh stopped after a cap because the post-shock SET
        # limit-cycles) -- so this bench measures the shipped mechanism, not a
        # bench-local variant of it.
        self.entropy = bool(entropy)
        self.entropy_refresh_max = int(entropy_refresh_max)
        self.ent = EntropyOperator(self.op.n_tets) if entropy else None
        self.sigma = None
        self.n_sigma_refresh = 0

        dir_nodes = np.unique(np.concatenate([
            mesh.boundary_faces["inlet"].ravel(),
            mesh.boundary_faces["outlet"].ravel()]))
        self.dir_nodes = dir_nodes
        mask = np.zeros(len(self.nodes), dtype=bool)
        mask[dir_nodes] = True
        self.free = np.where(~mask)[0]

    def state(self, phi):
        grad, q2 = self.op.velocities(phi)
        grad, q2 = grad.copy(), q2.copy()
        rho = density_field(q2, self.m_inf, self.gamma)
        if self.sigma is not None:
            rho = rho * self.sigma
        rho_t = self.upw.rho_tilde(grad, q2, rho, self.m_inf, self.C,
                                   self.m_crit, self.gamma,
                                   self.rho_floor).copy()
        return grad, q2, rho, rho_t

    def refresh_sigma(self, phi):
        """Rebuild the frozen sigma from `phi` using the donor map of the walk
        that the flux itself uses (upw._upstream is filled by state())."""
        _, q2, _, _ = self.state(phi)
        prev = self.sigma
        sig = self.ent.sigma(q2, self.upw._upstream, self.m_inf, self.gamma)
        self.sigma = sig.copy()
        self.n_sigma_refresh += 1
        return (0.0 if prev is None
                else float(np.max(np.abs(self.sigma - prev))))

    def residual(self, phi):
        _, _, _, rho_t = self.state(phi)
        return self.op.assemble_residual(phi, rho_t), rho_t

    def newton(self, phi0, n_max: int = 60, tol: float = 1e-12,
               verbose: bool = False):
        """Damped Newton at frozen upstream selection (the shipped
        sensitivities). Returns (phi, info)."""
        phi = np.asarray(phi0, dtype=np.float64).copy()
        if self.entropy:
            self.refresh_sigma(phi)
        R, _ = self.residual(phi)
        hist = [float(np.max(np.abs(R[self.free])))]
        n_used = 0
        for it in range(n_max):
            grad, q2, rho, rho_t = self.state(phi)
            s_e, s_u, ups = self.upw.rho_tilde_sensitivities(
                grad, q2, rho, self.m_inf, self.C, self.m_crit, self.gamma,
                self.rho_floor)
            J = self.op.assemble_newton_jacobian(
                phi, rho_t, s_e.copy(), s_u.copy(), ups.copy())
            R = self.op.assemble_residual(phi, rho_t)
            Jff = J[self.free][:, self.free].tocsc()
            try:
                d = spla.spsolve(Jff, -R[self.free])
            except Exception:                                  # noqa: BLE001
                return phi, dict(converged=False, reason="linear solve failed",
                                 residual_history=hist, n_newton=n_used)
            if not np.all(np.isfinite(d)):
                return phi, dict(converged=False, reason="non-finite step",
                                 residual_history=hist, n_newton=n_used)
            lam, best = 1.0, None
            r0 = float(np.max(np.abs(R[self.free])))
            for _ in range(12):
                trial = phi.copy()
                trial[self.free] += lam * d
                Rt, _ = self.residual(trial)
                rt = float(np.max(np.abs(Rt[self.free])))
                if np.isfinite(rt) and rt < r0:
                    break
                if best is None or (np.isfinite(rt) and rt < best[0]):
                    best = (rt, lam, trial)
                lam *= 0.5
            else:
                if best is None:
                    return phi, dict(converged=False, reason="line search",
                                     residual_history=hist, n_newton=n_used)
                rt, lam, trial = best
            phi = trial
            if self.entropy and self.n_sigma_refresh < self.entropy_refresh_max:
                self.refresh_sigma(phi)
                Rt, _ = self.residual(phi)
                rt = float(np.max(np.abs(Rt[self.free])))
            hist.append(rt)
            n_used = it + 1
            if verbose:
                print(f"    newton {it:2d}: |R|={rt:.3e} lam={lam:.3g}")
            if rt < tol:
                return phi, dict(converged=True, reason="tol",
                                 residual_history=hist, n_newton=n_used)
        return phi, dict(converged=False, reason="cap",
                         residual_history=hist, n_newton=n_used)


# ---------------------------------------------------------------------------
# 4. measurements
# ---------------------------------------------------------------------------

def element_u(system, phi):
    """Element-wise x-velocity and centroid x."""
    grad, _ = system.op.velocities(phi)
    xc = system.nodes[system.elements].mean(axis=1)[:, 0]
    return xc, grad[:, 0].copy()


def shock_position(system, phi, st, n_bins: int = None):
    """x where the bin-averaged u crosses the sonic value u* = sqrt(q2_star).

    Bin-average over x so the quasi-1-D profile is recovered from the
    unstructured element cloud, then linearly interpolate the sonic crossing.
    Returns (x_shock, n_cells_transition) or (nan, -1) if u never crosses.
    """
    xc, ux = element_u(system, phi)
    L = float(system.nodes[:, 0].max())
    if n_bins is None:
        # one bin per element column: keep it as fine as the mesh
        n_bins = max(8, int(round(L / (np.ptp(xc) / max(len(np.unique(
            np.round(xc, 6))), 1) + 1e-30))))
        n_bins = min(n_bins, 400)
    edges = np.linspace(0.0, L, n_bins + 1)
    idx = np.clip(np.digitize(xc, edges) - 1, 0, n_bins - 1)
    ub = np.full(n_bins, np.nan)
    for b in range(n_bins):
        m = idx == b
        if m.any():
            ub[b] = ux[m].mean()
    xb = 0.5 * (edges[:-1] + edges[1:])
    ok = ~np.isnan(ub)
    xb, ub = xb[ok], ub[ok]
    ustar = float(np.sqrt(st["q2_star"]))
    above = ub > ustar
    if above.all() or (~above).all():
        return float("nan"), -1
    k = int(np.argmax(above[:-1] & ~above[1:]))     # last supersonic bin
    u0, u1 = ub[k], ub[k + 1]
    t = (u0 - ustar) / (u0 - u1) if u0 != u1 else 0.5
    x_shock = xb[k] + t * (xb[k + 1] - xb[k])
    # transition width: bins between 99 % of u_sup and 101 % of u_sub
    hi = st["u_sup"] - 0.01 * (st["u_sup"] - st["u_sub"])
    lo = st["u_sub"] + 0.01 * (st["u_sup"] - st["u_sub"])
    n_cells = int(np.count_nonzero((ub < hi) & (ub > lo)))
    return float(x_shock), n_cells


def mass_flux_profile(system, phi, st, n_bins: int = 40):
    """Bin-averaged discrete mass flux rho_tilde * u_x, normalised by mdot.

    The exact solution has rho u == mdot everywhere; the artificial density
    replaces rho by rho_tilde in the flux the scheme conserves, so the
    deviation of rho_tilde*u from mdot is exactly the conservation error the
    shock operator introduces. Returns (x_bins, flux/mdot).
    """
    _, _, _, rho_t = system.state(phi)
    grad, _ = system.op.velocities(phi)
    xc = system.nodes[system.elements].mean(axis=1)[:, 0]
    f = rho_t * grad[:, 0]
    L = float(system.nodes[:, 0].max())
    edges = np.linspace(0.0, L, n_bins + 1)
    idx = np.clip(np.digitize(xc, edges) - 1, 0, n_bins - 1)
    fb = np.full(n_bins, np.nan)
    for b in range(n_bins):
        m = idx == b
        if m.any():
            fb[b] = f[m].mean()
    return 0.5 * (edges[:-1] + edges[1:]), fb / st["mdot"]
