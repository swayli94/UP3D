"""C05 的载体：**准一维 Laval 喷管** —— 项目里唯一带激波的精确解。

★★★ **来源与「复活」的性质。** 本模块是 phase 2 的 GS1.1 激波试验台
（`phases/p2/bench/s1_duct/{nozzle,duct}.py`）**按 AST 精确抽取**回活树的结果：
`nozzle.py` 全体 + `duct.py` 的 `DuctSystem` / `element_u`（抽取时验过这两块
**不再引用 duct.py 的任何其它名字**）。Case A（等截面管）的部分**没有带回来** ——
它被 phase 2 自己否掉了（等截面管里每个激波位置都是精确解，且存在满足同样两个
Dirichlet 条件的**光滑解**，实测 Newton 从激波解出发收敛到光滑声速支，
所以那个试验台回答不了它要问的问题）。

★★★ **这道门护住什么、护不住什么 —— 必须说准。**
`DuctSystem` 用的是**生产核**：`PicardOperator.assemble_residual` /
`assemble_newton_jacobian`、`UpwindOperator.rho_tilde_sensitivities`、
`EntropyOperator`。⇒ C05 **护住 `pyfp3d/kernels/` 的人工密度、迎风与熵修正**，
对着一个**含激波的精确解**。
★ 它**护不住** `pyfp3d/solve/newton.py`：全局化那一层（线搜索、EW forcing、
freeze、Kutta 行）是本模块自己的阻尼 Newton —— 因为喷管是**进出口 Dirichlet**，
生产驱动假定尾迹 + Kutta。**把 C05 的绿读成「求解器能正确算激波」是错的读法。**

几何（刻意做得细长，使准一维参照精确到 O((dH/dx)^2) ~ 1e-4）::

    H(x) = H_t (1 + a_in  ((x_t - x)/x_t)^2)          x <  x_t
         = H_t (1 + a_out ((x - x_t)/(L - x_t))^2)    x >= x_t

L = 20, x_t = 6, H_t = 0.5, a_in = 0.6, a_out = 0.361 —— 最后一个是这样选的：
让目标激波落在 x_s = 12 且 M ~ 1.3，**即真实翼型激波的强度**，且在全速势的有效包络内。

★★ **唯一性不是假设，有 `verify_uniqueness` 可查**：喉部把质量流固定，解族只由激波
位置 x_s 参数化，Δφ 对 x_s **严格单调** ⇒ 给定进出口 φ 只挑出一个 x_s；而全亚声速
（未壅塞）族的 Δφ 落在**不相交的、更低的**区间。
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


def element_u(system, phi):
    """Element-wise x-velocity and centroid x."""
    grad, _ = system.op.velocities(phi)
    xc = system.nodes[system.elements].mean(axis=1)[:, 0]
    return xc, grad[:, 0].copy()