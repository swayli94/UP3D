"""GS4.1 strip core: a **two-dimensional** integral boundary layer marched
along the streamwise coordinate, independent of the three-dimensional
`ibl3.py` system.

Scope and provenance (GS4.1 round 1, pre-registration
`docs/dev_phase_four/20260818-2100-gs41-strip-core-prereg.md` + addendum #1):

- the closure family is **reused verbatim** from `closures.py` -- this module
  contains no closure formula and no closure constant of its own, so a
  measurement made here is a statement about the *discretization*, never
  about a re-implemented closure (guard G-REUSE);
- `ibl3.py` is **frozen** (roadmap: kept, not deleted); nothing here imports
  its assembly or its solver;
- ★ the zero-pressure-gradient march is **not new**: `march_2d`/`_rhs_2d` in
  `bench/studies/v1_ibl3_standalone/run.py` is the GV1.1 reference arm and
  predates this module. What this module adds is (i) library status, (ii) the
  pressure-gradient terms, and (iii) the self-similar sweep those enable.

Why this is cheaper than `ibl3.py`, which is the whole point of the gate: the
three-dimensional system assembles and factorizes a coupled 6n system over
the surface, while a strip marches one station at a time and never assembles
a global system. Cost is therefore linear in stations with a tiny constant.

Equations (two-dimensional reduction of D13, incompressible edge; Drela's
momentum + kinetic-energy integral pair, in conservative form):

    dtheta/dxi      = c_f/2 - (H + 2) (theta /u_e) du_e/dxi
    dtheta_star/dxi = 2 c_D - 3      (theta_star/u_e) du_e/dxi
    dk_tau/dxi      = S_tau                       (turbulent only)

with theta, theta_star, c_f, c_D, H supplied by the closure packet and
k_tau = delta * Ctau1 * ku1 its shear-stress conserved density. The unknowns
are the closure's own state components, `y = (delta, A)` laminar and
`(delta, A, Ctau1)` turbulent, so the system is *implicit*: the chain rule
turns it into `M y' = F` with `M` the closure's analytic state Jacobian of
the conserved quantities. The crossflow components (B, Psi, Ctau2) are
identically zero here -- that is what makes this a strip.

★ The `q`-dependence is handled exactly rather than assumed away. Laminar
thicknesses are measurably independent of `re_d` (addendum #1 §4), but the
turbulent ones are not, so a varying edge speed contributes an explicit
`(d theta/d re_d)(rho delta/mu) du_e/dxi` term. It is carried on both
branches; on the laminar branch it is exactly zero and costs nothing.

★★ UNVERIFIED COMBINATION, stated rather than hidden: round 1 verifies the
laminar branch with and without a pressure gradient (Blasius, Falkner-Skan)
and the turbulent branch **only at zero pressure gradient**. The turbulent
branch under a pressure gradient runs, but nothing in this round measures it;
do not quote it as verified.
"""

import numpy as np

from pyfp3d.viscous import closures as C

# Marching-only numerics. These are properties of the ODE integrator, not of
# the boundary-layer model -- every physical constant lives in closures.py.
A_MIN = 0.05           # state floor mirroring the GV1.1 reference march
CT_MIN = 1.0e-10       # stress floor, likewise
_H_BRANCH_LO = 5.0     # bracket of the physical A-branch (H(A) is
_H_BRANCH_HI = 8.4     # non-monotone: the same H has two roots, addendum #1)


class StripState:
    """One marched strip: stations and the closure readings along them.

    Attributes are plain arrays, one entry per recorded station, so a caller
    can compare against an analytic solution without touching the solver.
    """

    __slots__ = ("x", "delta", "A", "ctau", "H", "ds1", "theta", "theta_star",
                 "cf", "cD", "re_theta", "re_x", "ue", "turbulent",
                 "wall_time", "n_substep")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    def as_dict(self):
        return {k: getattr(self, k) for k in self.__slots__}


def _state_vector(y, turbulent):
    """Closure state (6,) from the marched unknowns, with the floors applied."""
    st = np.zeros(6)
    st[0] = max(y[0], C.DELTA_MIN)
    st[1] = max(y[1], A_MIN)
    st[4] = max(y[2], CT_MIN) if turbulent else C.CTAU_LAM
    return st


def _rhs(y, ue, due, rho, mu, turbulent, c_l):
    """`y'` from `M y' = F` at one station.

    `ue`/`due` are the edge speed and its streamwise derivative there. The
    returned `out` is the closure packet, so a caller marching and a caller
    recording share one evaluation.
    """
    st = _state_vector(y, turbulent)
    out, dout, dout_e = C.closure_scalar(st, q=ue, rho=rho, mu=mu,
                                         turbulent=turbulent, c_l=c_l)
    theta = out[C.OUT_TH11]
    theta_s = out[C.OUT_THS1]
    H = out[C.OUT_H1]

    # Explicit edge-speed dependence of the thicknesses, through re_d. Zero on
    # the laminar branch (measured), nonzero on the turbulent one.
    dre_dq = rho * st[0] / mu
    ex_th = dout_e[C.OUT_TH11, 0] * dre_dq * due
    ex_ts = dout_e[C.OUT_THS1, 0] * dre_dq * due

    f1 = 0.5 * out[C.OUT_CF1] - (H + 2.0) * theta / ue * due - ex_th
    f2 = 2.0 * out[C.OUT_CD] - 3.0 * theta_s / ue * due - ex_ts

    if not turbulent:
        M = np.array([[dout[C.OUT_TH11, 0], dout[C.OUT_TH11, 1]],
                      [dout[C.OUT_THS1, 0], dout[C.OUT_THS1, 1]]])
        return np.linalg.solve(M, np.array([f1, f2])), out

    # Shear-stress lag: conserved density k_tau = delta * Ctau1 * ku1, source
    # S_tau from the closure's own stress-source routine (no formula here).
    ct = st[4]
    src = np.empty(3)
    dsrc = np.empty((3, 6))
    C.stress_source(st, ue, rho, c_l, out[C.OUT_SP1], out[C.OUT_SD], 0,
                    src, dsrc)
    f3 = src[2]
    de, ku1 = st[0], out[C.OUT_KU1]
    M = np.array([
        [dout[C.OUT_TH11, 0], dout[C.OUT_TH11, 1], dout[C.OUT_TH11, 4]],
        [dout[C.OUT_THS1, 0], dout[C.OUT_THS1, 1], dout[C.OUT_THS1, 4]],
        [ct * ku1 + de * ct * dout[C.OUT_KU1, 0],
         de * ct * dout[C.OUT_KU1, 1],
         de * ku1 + de * ct * dout[C.OUT_KU1, 4]],
    ])
    return np.linalg.solve(M, np.array([f1, f2, f3])), out


def similar_seed(theta, target_H, ue=1.0, rho=1.0, mu=1.0e-5):
    """Seed `(delta, A)` with a prescribed momentum thickness and shape factor.

    Both targets are **explicit arguments**: an earlier version derived theta
    from the Blasius correlation internally, which silently imposed a
    flat-plate thickness on wedge flows and marched them off the physical
    branch. The caller owns the similarity solution; this only maps it onto
    the closure family.

    `target_H` is resolved on the physical `A`-branch -- the one continuously
    connected to Blasius. `H(A)` folds, so the same `H` has a second root that
    a bare Newton can land on.
    """
    A = _branch_A(target_H, ue, rho, mu)
    out, _, _ = C.closure_scalar((1.0, A, 0.0, 0.0, C.CTAU_LAM, 0.0),
                                 q=ue, rho=rho, mu=mu, turbulent=False)
    return np.array([theta / out[C.OUT_TH11], A])   # theta is delta-linear


def _branch_A(target_H, ue, rho, mu):
    """`A` with `H(A) = target_H`, restricted to the physical branch."""
    def h_of(A):
        out, _, _ = C.closure_scalar((1.0, A, 0.0, 0.0, C.CTAU_LAM, 0.0),
                                     q=ue, rho=rho, mu=mu, turbulent=False)
        return out[C.OUT_H1]

    lo, hi = _H_BRANCH_LO, _H_BRANCH_HI
    h_lo, h_hi = h_of(lo), h_of(hi)
    if not (min(h_lo, h_hi) <= target_H <= max(h_lo, h_hi)):
        raise ValueError(
            f"target H={target_H:.4f} outside the physical A-branch "
            f"[{min(h_lo, h_hi):.4f}, {max(h_lo, h_hi):.4f}] -- the closure "
            "family cannot represent it (report, do not extrapolate)")
    for _ in range(200):                       # bisection: robust on a fold
        mid = 0.5 * (lo + hi)
        if (h_of(mid) - target_H) * (h_lo - target_H) > 0.0:
            lo, h_lo = mid, h_of(mid)
        else:
            hi = mid
        if hi - lo < 1.0e-13:
            break
    return 0.5 * (lo + hi)


def march_correlation(stations, y0, x_start, ue_fn, rho=1.0, mu=1.0e-5,
                      n_substep=2000, x_tr=None, arm="new", lag=False):
    """March the strip with the **correlation** closure (`closures_2d.py`).

    GS4.1 round 3, route (a2). The state is `(theta, H)` and the system is
    explicit -- no implicit `M y' = F`, no quadrature, no state Jacobian -- which
    is the whole cost argument. Laminar only.

    ★ `lag=True` (GS4.1 round 9 leg B) extends the state to `(theta, H, sqrt(Ctau))`
    and evolves the shear-stress lag instead of holding `Ctau = CtauEQ`. The
    third component is seeded at the transition crossing from
    `closures_2d.s_tau_at_transition`, three decades below equilibrium, exactly as
    `xblsys.f:1393/1403` does. `lag=False` is the default and is bit-identical to
    rounds 3-9A -- asserted in `tests/test_gs41_turbulent_closure.py`.
    ★ A zero-pressure-gradient plate cannot discriminate the two arms: they were
    measured 1.3 % apart downstream, against 4.3x just after transition.

    `y0 = (theta, H)`. Kept as a separate entry point rather than a branch inside
    `march` so that the profile path stays byte-for-byte what round 1 measured
    (guard G-LEGACY); the two closures are separate authorities, not two
    implementations of one model.

    `x_tr` forces transition: stations upstream of it use the laminar closure,
    stations at or beyond it the turbulent one at local equilibrium. ★ With
    `x_tr=None` (the default) this is byte-for-byte the laminar-only march of
    round 3 -- guard G-LEGACY -- so adding turbulence moves no existing reading.
    `arm` selects the turbulent `H*` correlation, "new" (what XFOIL 6.99 runs)
    or "old" (the form commented out above it); it is inert when `x_tr is None`.

    ★★ Transition is an INSTANTANEOUS switch at `x_tr` holding theta and H
    continuous, not the e^N ramp of D13 eq (34)(35). That matches the project's
    existing D-TR choice on the 3-D side, and its known consequence is already
    on record: GV3.1 measured cf +44 % at the first post-trip station against
    XFOIL's ramp. A large deviation at that one station is the transition model,
    not the correlations.

    Raises `closures_2d.ClosureRangeError` if the march reaches separation
    (`H >= H_SEPARATION_GUARD`) or leaves the correlation's range -- reported,
    never clamped, because the direct form is singular there by construction.
    """
    import time

    from pyfp3d.viscous import closures_2d as C2

    xs = np.asarray(stations, dtype=float)
    if xs.ndim != 1 or xs.size == 0:
        raise ValueError("stations must be a non-empty 1-D array")
    if np.any(np.diff(xs) <= 0.0):
        raise ValueError("stations must be strictly increasing")
    if xs[0] <= x_start:
        raise ValueError("stations must lie downstream of x_start")
    if x_start <= 0.0:
        raise ValueError("x_start must be positive")

    y = np.array(y0, dtype=float)
    if lag and y.size == 2:
        y = np.append(y, 0.0)          # seeded at the transition crossing below
    seeded = bool(lag and y.size == 3 and y[2] > 0.0)
    rec = {k: [] for k in ("x", "theta", "H", "ds1", "theta_star", "cf", "cD",
                           "re_theta", "re_x", "ue", "delta", "A", "ctau")}

    def _turbulent_at(xx):
        return x_tr is not None and xx >= x_tr

    def _f(yy, xx):
        ue, due = ue_fn(xx)
        if _turbulent_at(xx):
            if lag:
                return np.array(C2.rhs_turb(max(yy[0], C.DELTA_MIN), yy[1], ue,
                                            due, rho=rho, mu=mu, arm=arm,
                                            s_tau=max(yy[2], 1.0e-14)))
            return np.array(C2.rhs_turb(max(yy[0], C.DELTA_MIN), yy[1], ue, due,
                                        rho=rho, mu=mu, arm=arm))
        if yy[1] >= C2.H_SEPARATION_GUARD:
            raise C2.ClosureRangeError(
                f"H = {yy[1]:.4f} reached the separation guard "
                f"{C2.H_SEPARATION_GUARD} -- the direct two-equation form is "
                "singular at H = 4 (GS4.2's motivation); stopping this leg")
        out = C2.rhs(max(yy[0], C.DELTA_MIN), yy[1], ue, due, rho=rho, mu=mu)
        return np.array(out + (0.0,) if lag else out)

    t0 = time.perf_counter()
    x = float(x_start)
    span_log = np.log(xs[-1] / x_start)
    for i in range(xs.size):
        seg_log = np.log(xs[i] / x)
        nsub = max(1, int(round(n_substep * seg_log / span_log)))
        ratio = (xs[i] / x) ** (1.0 / nsub)
        for _ in range(nsub):
            dx = x * (ratio - 1.0)
            if lag and not seeded and _turbulent_at(x):
                # xblsys.f:1393/1403 -- the transition seed, three decades below
                # equilibrium, which is what makes the two arms distinguishable
                ue_s, _ = ue_fn(x)
                p_s = C2.packet_turb(max(y[0], C.DELTA_MIN), y[1], ue_s,
                                     rho=rho, mu=mu, arm=arm)
                y[2] = C2.s_tau_at_transition(y[1], p_s["Ctau_eq"])
                seeded = True
            k1 = _f(y, x)
            k2 = _f(y + 0.5 * dx * k1, x + 0.5 * dx)
            k3 = _f(y + 0.5 * dx * k2, x + 0.5 * dx)
            k4 = _f(y + dx * k3, x + dx)
            y = y + dx * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            y[0] = max(y[0], C.DELTA_MIN)
            if lag:
                y[2] = max(y[2], 0.0)
            x += dx
        x = xs[i]
        ue, _ = ue_fn(x)
        p = (C2.packet_turb(y[0], y[1], ue, rho=rho, mu=mu, arm=arm)
             if _turbulent_at(x) else C2.packet(y[0], y[1], ue, rho=rho, mu=mu))
        rec["x"].append(x)
        rec["theta"].append(y[0])
        rec["H"].append(y[1])
        rec["ds1"].append(y[0] * y[1])
        rec["theta_star"].append(y[0] * p["H_star"])
        rec["cf"].append(p["cf"])
        rec["cD"].append(p["cD"])
        rec["re_theta"].append(p["re_theta"])
        rec["re_x"].append(rho * ue * x / mu)
        rec["ue"].append(ue)
        rec["delta"].append(np.nan)     # the correlation state carries no delta
        rec["A"].append(np.nan)         # nor a wall-slope parameter
        rec["ctau"].append(y[2] ** 2 if (lag and _turbulent_at(x))
                           else p.get("Ctau_eq", np.nan))
    wall = time.perf_counter() - t0

    return StripState(turbulent=False, wall_time=wall, n_substep=n_substep,
                      **{k: np.asarray(v) for k, v in rec.items()})


def march(stations, y0, x_start, ue_fn, rho=1.0, mu=1.0e-5, turbulent=False,
          c_l=C.C_L_DEFAULT, n_substep=2000):
    """March the strip from `x_start` and record the closure at `stations`.

    Uses the **profile** closure (`closures.py`). The correlation-closure
    counterpart is `march_correlation`; see its docstring and `closures_2d.py`
    for which closure is authoritative for what.

    `ue_fn(x)` returns `(u_e, du_e/dx)`. Integration is classical RK4 with the
    substeps distributed so that the march **lands exactly on every recording
    station** -- interpolating instead sets a ~1e-6 comparison noise floor,
    which is what masked a refinement signal in GV1.1's first execution.

    Substeps are geometric in `x` and are budgeted per segment by that
    segment's share of `log(x)`, not of `x`. A boundary layer is self-similar
    in `log x`, so a station layout spanning decades would otherwise starve
    the upstream segments -- where the march is stiffest -- of nearly every
    step. (GV1.1's reference arm used linear station spacing, where the two
    rules coincide.)

    `n_substep` is the total substep budget over the span; it is the knob a
    refinement study varies, and the recording stations never move with it.
    """
    import time

    xs = np.asarray(stations, dtype=float)
    if xs.ndim != 1 or xs.size == 0:
        raise ValueError("stations must be a non-empty 1-D array")
    if np.any(np.diff(xs) <= 0.0):
        raise ValueError("stations must be strictly increasing")
    if xs[0] <= x_start:
        raise ValueError("stations must lie downstream of x_start")
    if x_start <= 0.0:
        raise ValueError("x_start must be positive (the leading edge is "
                         "singular; seed downstream of it)")

    y = np.array(y0, dtype=float)
    if turbulent and y.size == 2:
        y = np.append(y, CT_MIN * 10.0)
    rec = {k: [] for k in ("x", "delta", "A", "ctau", "H", "ds1", "theta",
                           "theta_star", "cf", "cD", "re_theta", "re_x", "ue")}

    def _f(yy, xx):
        ue, due = ue_fn(xx)
        return _rhs(yy, ue, due, rho, mu, turbulent, c_l)[0]

    def _clip(yy):
        yy[0] = max(yy[0], C.DELTA_MIN)
        yy[1] = max(yy[1], A_MIN)
        if turbulent:
            yy[2] = max(yy[2], CT_MIN)
        return yy

    t0 = time.perf_counter()
    x = float(x_start)
    span_log = np.log(xs[-1] / x_start)
    for i in range(xs.size):
        seg_log = np.log(xs[i] / x)
        nsub = max(1, int(round(n_substep * seg_log / span_log)))
        ratio = (xs[i] / x) ** (1.0 / nsub)
        for _ in range(nsub):
            dx = x * (ratio - 1.0)
            k1 = _f(y, x)
            k2 = _f(_clip(y + 0.5 * dx * k1), x + 0.5 * dx)
            k3 = _f(_clip(y + 0.5 * dx * k2), x + 0.5 * dx)
            k4 = _f(_clip(y + dx * k3), x + dx)
            y = _clip(y + dx * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0)
            x += dx
        x = xs[i]                                  # kill accumulated round-off
        ue, due = ue_fn(x)
        _, out = _rhs(y, ue, due, rho, mu, turbulent, c_l)
        rec["x"].append(x)
        rec["delta"].append(y[0])
        rec["A"].append(y[1])
        rec["ctau"].append(y[2] if turbulent else C.CTAU_LAM)
        rec["H"].append(out[C.OUT_H1])
        rec["ds1"].append(out[C.OUT_DS1])
        rec["theta"].append(out[C.OUT_TH11])
        rec["theta_star"].append(out[C.OUT_THS1])
        rec["cf"].append(out[C.OUT_CF1])
        rec["cD"].append(out[C.OUT_CD])
        rec["re_theta"].append(rho * ue * out[C.OUT_TH11] / mu)
        rec["re_x"].append(rho * ue * x / mu)
        rec["ue"].append(ue)
    wall = time.perf_counter() - t0

    return StripState(turbulent=turbulent, wall_time=wall, n_substep=n_substep,
                      **{k: np.asarray(v) for k, v in rec.items()})


def flat_plate_ue(u_inf=1.0):
    """Zero-pressure-gradient edge distribution."""
    return lambda x: (u_inf, 0.0)


def falkner_skan_ue(m, c=1.0):
    """Wedge-flow edge distribution `u_e = c x^m` and its derivative."""
    return lambda x: (c * x ** m, 0.0 if m == 0.0 else c * m * x ** (m - 1.0))


def similarity_fixed_point(m=0.0, rho=1.0, mu=1.0e-5, ue=1.0):
    """The closure family's **own** self-similar state for wedge exponent `m`.

    Solved algebraically from the two integral equations under the similarity
    ansatz (`A` constant, `theta = delta t(A)`, `theta_star = delta s(A)`),
    with no marching and no discretization whatsoever. It exists to
    cross-check the march: a march that converges to something else is a
    marching defect, not a closure property.

    Returns `(A_star, H, cf_sqrt_rex)`.
    """
    from scipy.optimize import brentq

    def pack(A):
        out, _, _ = C.closure_scalar((1.0, A, 0.0, 0.0, C.CTAU_LAM, 0.0),
                                     q=ue, rho=rho, mu=mu, turbulent=False)
        return out, out[C.OUT_TH11], out[C.OUT_THS1]

    # Similarity ansatz delta ~ x^k with k = (1-m)/2, so delta' = k delta/x.
    # With theta = t delta and theta_star = s delta the two equations read
    #   t (k + (H+2) m) delta/x = c_f/2
    #   s (k + 3 m)     delta/x = 2 c_D
    # and c_f, c_D both carry one factor mu/(rho u_e delta), so dividing them
    # removes delta and x entirely and leaves one algebraic condition on A.
    def _coeffs(A):
        out, t, s = pack(A)
        k = 0.5 * (1.0 - m)
        return out, t, s, k + (out[C.OUT_H1] + 2.0) * m, k + 3.0 * m

    def resid(A):
        out, t, s, d1, d2 = _coeffs(A)
        return (t * d1) / (s * d2) - (0.5 * out[C.OUT_CF1]) / (2.0 * out[C.OUT_CD])

    A_star = brentq(resid, _H_BRANCH_LO + 1.0, _H_BRANCH_HI, xtol=1.0e-13)
    out, t, _, d1, _ = _coeffs(A_star)
    # Back-substituting delta(x) into c_f sqrt(Re_x) collapses every
    # dimensional factor: c_f sqrt(Re_x) = 2 sqrt(A t d1).
    return A_star, out[C.OUT_H1], 2.0 * np.sqrt(A_star * t * d1)
