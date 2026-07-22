"""GV1.1 standalone IBL3 verification (Track V1 gate).

Binding text: docs/roadmap/track_v.md GV1.1(a)-(e); pre-registered bands:
cases/analysis/v1_ibl3_standalone/PRE_REGISTRATION.md (written before the
first execution). Regenerates every CSV/PNG in results/ and exits 0 iff all
pre-registered assertions hold (honest FAIL otherwise).

References are the closure's OWN 2-D ODE marches (design doc §3.2): the same
closure packet marched in x with the 2-D reduction of the governing
equations (von Karman + kinetic energy + stress lag), so gates judge the FE
implementation against its own 2-D limit. The 1/7th power-law turbulent
correlation is RECORDED (not pass/fail).

Run:  python cases/analysis/v1_ibl3_standalone/run.py
"""

import os
import sys

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.viscous import closures as C
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.surface_mesh import (
    SurfaceMesh,
    structured_rectangle_surface,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

RHO = 1.0
MU = 1.0e-5
Q = 1.0
EPS_DIFF = 0.005
EPS_DIFF_S = 0.02  # streamwise-tensor diffusion (D-HB follow-up; calibrated
                   # knee value: strict (e) decrease + H order ~1.0 + margin)
X0, X1 = 0.2, 2.2
ZH = 0.2

SUMMARY = []  # (gate, metric, band, measured, verdict)


def _record(gate, metric, band, measured, ok):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append((gate, metric, band, measured, verdict))
    print(f"  [{gate}] {metric}: band={band} measured={measured} "
          f"-> {verdict}")


# ---------------------------------------------------------------------------
# 2-D ODE reference marches (closure's own 2-D limit)
# ---------------------------------------------------------------------------

def _rhs_2d(y, turb):
    """Implicit-ODE right-hand side: solve M (delta', A' [, Ctau']) = F."""
    st = np.zeros(6)
    st[0] = max(y[0], 1.0e-8)
    st[1] = max(y[1], 0.05)
    if turb:
        st[4] = max(y[2], 1.0e-10)
    else:
        st[4] = C.CTAU_LAM
    out, dout = C.closure_scalar(st, q=Q, rho=RHO, mu=MU, turbulent=turb)
    f1 = 0.5 * out[C.OUT_CF1]
    f2 = 2.0 * out[C.OUT_CD]
    if not turb:
        M = np.array([
            [dout[C.OUT_TH11, 0], dout[C.OUT_TH11, 1]],
            [dout[C.OUT_THS1, 0], dout[C.OUT_THS1, 1]],
        ])
        F = np.array([f1, f2])
        return np.linalg.solve(M, F), out
    ct = st[4]
    sp1 = out[C.OUT_SP1]
    sd = out[C.OUT_SD]
    ku1 = out[C.OUT_KU1]
    f3 = 2.0 * C.A1_BRADSHAW * (
        ct * sp1 - (1.0 / C.C_L_DEFAULT) * np.sqrt(ct) * ct * sd
    )
    de = st[0]
    M = np.array([
        [dout[C.OUT_TH11, 0], dout[C.OUT_TH11, 1], 0.0],
        [dout[C.OUT_THS1, 0], dout[C.OUT_THS1, 1], 0.0],
        [ct * ku1 + de * ct * dout[C.OUT_KU1, 0],
         de * ct * dout[C.OUT_KU1, 1],
         de * ku1],
    ])
    F = np.array([f1, f2, f3])
    return np.linalg.solve(M, F), out


def march_2d(x_stations, y0, turb, x_start, nstep=4000):
    """RK4 march of the 2-D reduction; records closure outputs per station.

    The integration starts at the physical inflow station x_start (where the
    seed y0 lives), NOT at the first recording station: the ODE is
    x-autonomous, so teleporting the seed to xs[0] shifts the whole
    trajectory by xs[0] - x_start (first-execution bug, caught by the
    (e) diagnosis: the "march" recorded the un-marched inflow seed as its
    first station and stalled delta* over [x0, xs[0]]). Each segment
    [x_prev, xs[i]] gets an integer number of RK4 substeps so the march
    lands exactly on every station: linear recording interpolation sets a
    ~1e-6 comparison noise floor, which masked the 100x16 -> 200x32 error
    decrease in the (e) measurement.
    """
    xs = np.asarray(x_stations)
    y = np.array(y0, dtype=float)
    rec = {k: [] for k in ("x", "delta", "A", "ct", "H", "ds1", "cf", "th",
                           "re_th")}
    # record helper
    def _record_at(xx, yy):
        st = np.zeros(6)
        st[0] = max(yy[0], 1.0e-8)
        st[1] = max(yy[1], 0.05)
        st[4] = max(yy[2], 1.0e-10) if turb else C.CTAU_LAM
        out, _ = C.closure_scalar(st, q=Q, rho=RHO, mu=MU, turbulent=turb)
        rec["x"].append(xx)
        rec["delta"].append(yy[0])
        rec["A"].append(yy[1])
        rec["ct"].append(yy[2] if turb else C.CTAU_LAM)
        rec["H"].append(out[C.OUT_H1])
        rec["ds1"].append(out[C.OUT_DS1])
        rec["cf"].append(out[C.OUT_CF1])
        rec["th"].append(out[C.OUT_TH11])
        rec["re_th"].append(out[C.OUT_TH11] * Q / MU)

    x = x_start
    span = xs[-1] - x_start
    for i in range(xs.size):
        seg = xs[i] - x
        nsub = max(1, int(round(nstep * seg / span)))
        dx = seg / nsub
        for _ in range(nsub):
            k1, _ = _rhs_2d(y, turb)
            k2, _ = _rhs_2d(y + 0.5 * dx * k1, turb)
            k3, _ = _rhs_2d(y + 0.5 * dx * k2, turb)
            k4, _ = _rhs_2d(y + dx * k3, turb)
            y = y + dx * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            y[0] = max(y[0], 1.0e-8)
            y[1] = max(y[1], 0.05)
            if turb:
                y[2] = max(y[2], 1.0e-10)
        x = xs[i]
        _record_at(x, y)
    return {k: np.array(v) for k, v in rec.items()}


def _turb_seed(x):
    re_x = RHO * Q * x / MU
    delta = 0.37 * x * re_x ** -0.2
    cf = 0.0576 * re_x ** -0.2
    re_d = RHO * Q * delta / MU
    A = 0.5 * cf * re_d  # U_tau^2 = cf/2, U_tau = sqrt(A/Re_delta)
    st = np.array([delta, A, 0.0, 0.0, 1.0e-3, 0.0])
    out, _ = C.closure_scalar(st, q=Q, rho=RHO, mu=MU, turbulent=True)
    st[4] = max((C.C_L_DEFAULT * out[C.OUT_SP1] / out[C.OUT_SD]) ** 2, 1.0e-6)
    return st


# ---------------------------------------------------------------------------
# FE case driver
# ---------------------------------------------------------------------------

def run_fe(nx, nz, turbulent, ue_fn, seed_fn, x0=X0, x1=X1, zh=ZH):
    xyz, tris = structured_rectangle_surface(x0, x1, -zh, zh, nx, nz)
    sm = SurfaceMesh.from_wall_faces(xyz, tris)
    n = xyz.shape[0]
    u_e = np.zeros((n, 3))
    for i in range(n):
        u_e[i] = ue_fn(xyz[i, 0])
    inflow = np.abs(xyz[:, 0] - x0) < 1.0e-12
    flags = np.full(n, 1 if turbulent else 0, dtype=np.int64)
    st_bc = seed_fn(x0)
    solver = IBL3Solver(sm, u_e, RHO, MU, 0.0, flags, inflow, st_bc,
                        eps_diff=EPS_DIFF, eps_diff_s=EPS_DIFF_S)
    U0 = np.zeros((n, 6))
    for i in range(n):
        U0[i] = seed_fn(max(xyz[i, 0], 1.0e-3))
    U, info = solver.solve(U0, tol=1.0e-9, max_iter=100)
    # closure outputs at the solution
    outs = np.empty((n, C.N_OUT))
    douts = np.empty((n, C.N_OUT, 6))
    C.closure_all(U, np.full(n, Q), np.full(n, RHO), np.full(n, MU),
                  np.zeros(n), flags, C.C_L_DEFAULT, outs, douts)
    h = (x1 - x0) / nx
    interior = (
        (xyz[:, 0] > x0 + 2.0 * h)
        & (xyz[:, 0] < x1 - 2.0 * h)
        & (np.abs(xyz[:, 2]) < zh - 1.5 * h)
    )
    return {
        "sm": sm, "U": U, "outs": outs, "info": info, "h": h,
        "interior": interior, "flags": flags,
    }


def _centerline(fe, col):
    """Centerline station profile of closure output `col`."""
    xyz = fe["sm"].xyz
    xs = np.unique(np.round(xyz[:, 0], 12))
    prof = []
    xout = []
    for xx in xs:
        m = (np.abs(xyz[:, 0] - xx) < 1.0e-9) & fe["interior"]
        if np.any(m):
            prof.append(np.mean(fe["outs"][m, col]))
            xout.append(xx)
    return np.array(xout), np.array(prof)


# ---------------------------------------------------------------------------
# Gate (a) + (e): laminar plate, refinement family
# ---------------------------------------------------------------------------

def gate_a_e():
    print("== GV1.1(a)/(e): laminar flat plate ==")
    stations = np.linspace(X0 + 0.1, X1 - 0.1, 40)
    seed0 = C.blasius_seed(X0, q=Q, rho=RHO, mu=MU)
    march = march_2d(stations, (seed0[0], seed0[1]), False, X0)

    fes = []
    for nx, nz in ((50, 8), (100, 16), (200, 32)):
        fe = run_fe(nx, nz, False, lambda x: np.array([Q, 0.0, 0.0]),
                    lambda x: C.blasius_seed(x, q=Q, rho=RHO, mu=MU))
        print(f"  mesh {nx}x{nz}: converged={fe['info']['converged']} "
              f"it={fe['info']['n_iter']} |R|={fe['info']['final_residual']:.2e}")
        assert fe["info"]["converged"], f"laminar FE solve failed on {nx}x{nz}"
        fes.append(fe)

    # --- write profile CSV (all meshes) ---
    rows = ["mesh,x,H,ds1,cf,theta"]
    for (nx, nz), fe in zip(((50, 8), (100, 16), (200, 32)), fes):
        for col_name, col in (("H", C.OUT_H1),):
            pass
        xc, H = _centerline(fe, C.OUT_H1)
        _, ds1 = _centerline(fe, C.OUT_DS1)
        _, cf = _centerline(fe, C.OUT_CF1)
        _, th = _centerline(fe, C.OUT_TH11)
        for i in range(xc.size):
            rows.append(f"{nx}x{nz},{xc[i]:.4f},{H[i]:.8f},{ds1[i]:.8f},"
                        f"{cf[i]:.8f},{th[i]:.8f}")
    with open(os.path.join(RESULTS, "gv1_1a_profiles.csv"), "w") as f:
        f.write("\n".join(rows) + "\n")
    mrows = ["x,H_march,ds1_march,cf_march,theta_march"]
    for i in range(march["x"].size):
        mrows.append(f"{march['x'][i]:.4f},{march['H'][i]:.8f},"
                     f"{march['ds1'][i]:.8f},{march['cf'][i]:.8f},"
                     f"{march['th'][i]:.8f}")
    with open(os.path.join(RESULTS, "gv1_1a_march.csv"), "w") as f:
        f.write("\n".join(mrows) + "\n")

    # --- (a) H band on the finest mesh ---
    fe2 = fes[-1]
    H_int = fe2["outs"][fe2["interior"], C.OUT_H1]
    H_dev = float(np.max(np.abs(H_int - 2.59)) / 2.59)
    _record("a", "max |H-2.59|/2.59 over interior", "<= 2%",
            f"{100.0 * H_dev:.2f}% (H range "
            f"[{H_int.min():.4f}, {H_int.max():.4f}])", H_dev <= 0.02)

    # --- (a) delta* power law on the finest mesh ---
    xc, ds1 = _centerline(fe2, C.OUT_DS1)
    p = np.polyfit(np.log(xc), np.log(ds1), 1)[0]
    _record("a", "delta*(x) power-law exponent", "[0.48, 0.52]",
            f"{p:.4f}", 0.48 <= p <= 0.52)
    # RECORDED diagnostics for the exponent gap: the measured window is
    # dominated by the family fixed-point adjustment transient (H rises
    # 2.59 -> H* ~ 2.71), during which delta* grows faster than sqrt(x);
    # the [0.48, 0.52] band assumes near-similarity. Quantify the cause:
    # same fit on the march reference, and a downstream-half FE fit where
    # H has nearly equilibrated.
    ds_mi = np.interp(xc, march["x"], march["ds1"])
    p_march = float(np.polyfit(np.log(xc), np.log(ds_mi), 1)[0])
    _record("a", "delta* exponent of the march reference (same fit)",
            "diagnostic only", f"{p_march:.4f}", None)
    late = xc > 0.5 * (X0 + X1)
    p_late = float(np.polyfit(np.log(xc[late]), np.log(ds1[late]), 1)[0])
    _record("a", "delta* exponent, FE downstream half (x > 1.2)",
            "diagnostic only", f"{p_late:.4f}", None)

    # --- (e) refinement: error vs closure's own march, evaluated exactly
    # at each mesh's centerline stations (per-mesh march, exact landing;
    # no recording interpolation -> no ~1e-6 comparison noise floor) ---
    errs_H = []
    errs_ds = []
    hs = []
    for (nx, nz), fe in zip(((50, 8), (100, 16), (200, 32)), fes):
        xc, H = _centerline(fe, C.OUT_H1)
        _, ds1 = _centerline(fe, C.OUT_DS1)
        mref = march_2d(xc, (seed0[0], seed0[1]), False, X0)
        errs_H.append(float(np.max(np.abs(H - mref["H"]))))
        errs_ds.append(float(np.max(np.abs(ds1 - mref["ds1"]))))
        hs.append((X1 - X0) / nx)
    ord_H = np.log2(errs_H[0] / errs_H[1]), np.log2(errs_H[1] / errs_H[2])
    ord_ds = np.log2(errs_ds[0] / errs_ds[1]), np.log2(errs_ds[1] / errs_ds[2])
    ok_e = errs_H[2] < errs_H[1] < errs_H[0] and errs_ds[2] < errs_ds[1] < errs_ds[0]
    _record("e", "refinement error decrease (H, delta* vs march)",
            "strictly decreasing",
            f"errH={['%.2e' % e for e in errs_H]} orderH="
            f"[{ord_H[0]:.2f},{ord_H[1]:.2f}] errds={['%.2e' % e for e in errs_ds]}"
            f" orderds=[{ord_ds[0]:.2f},{ord_ds[1]:.2f}]", ok_e)
    with open(os.path.join(RESULTS, "gv1_1e_refinement.csv"), "w") as f:
        f.write("h,err_H_max,err_ds1_max\n")
        for i in range(3):
            f.write(f"{hs[i]:.5f},{errs_H[i]:.8e},{errs_ds[i]:.8e}\n")

    # --- (d) lock on the laminar solution ---
    B = float(np.max(np.abs(fe2["U"][:, 2])))
    Psi = float(np.max(np.abs(fe2["U"][:, 3])))
    ct2_lam = np.unique(np.round(fe2["U"][:, 5], 12))
    ok_d = B < 1.0e-10 and Psi < 1.0e-10 and np.allclose(ct2_lam, C.CTAU_LAM)
    _record("d", "quasi-2-D lock (laminar)", "|B|,|Psi|<1e-10, Ctau2=CTAU_LAM",
            f"|B|={B:.2e} |Psi|={Psi:.2e} Ctau2={ct2_lam}", ok_d)
    return fes, march


# ---------------------------------------------------------------------------
# Gate (b): turbulent plate
# ---------------------------------------------------------------------------

def gate_b():
    print("== GV1.1(b): turbulent flat plate ==")
    stations = np.linspace(X0 + 0.1, X1 - 0.1, 40)
    st0 = _turb_seed(X0)
    march = march_2d(stations, (st0[0], st0[1], st0[4]), True, X0)
    fe = run_fe(100, 16, True, lambda x: np.array([Q, 0.0, 0.0]), _turb_seed)
    print(f"  FE converged={fe['info']['converged']} it={fe['info']['n_iter']} "
          f"|R|={fe['info']['final_residual']:.2e}")
    assert fe["info"]["converged"], "turbulent FE solve failed"

    xc, cf = _centerline(fe, C.OUT_CF1)
    _, th = _centerline(fe, C.OUT_TH11)
    _, H = _centerline(fe, C.OUT_H1)
    re_th_fe = th * Q / MU
    # compare at matched Re_theta (gate: C_f(Re_theta))
    lg_fe = np.log(re_th_fe)
    lg_ma = np.log(march["re_th"])
    lo = max(lg_fe.min(), lg_ma.min())
    hi = min(lg_fe.max(), lg_ma.max())
    grid = np.linspace(lo, hi, 30)
    cf_fe_i = np.interp(grid, lg_fe, cf)
    cf_ma_i = np.interp(grid, lg_ma, march["cf"])
    rel = np.abs(cf_fe_i - cf_ma_i) / cf_ma_i
    re_grid = np.exp(grid)
    # RECORDED external reference: 1/7th power law cf = 0.02507 Re_theta^-0.25
    cf_pw = 0.02507 * re_grid ** -0.25
    with open(os.path.join(RESULTS, "gv1_1b_cf_retheta.csv"), "w") as f:
        f.write("re_theta,cf_FE,cf_march,rel_err,cf_powerlaw_recorded\n")
        for i in range(grid.size):
            f.write(f"{re_grid[i]:.2f},{cf_fe_i[i]:.8f},{cf_ma_i[i]:.8f},"
                    f"{rel[i]:.6f},{cf_pw[i]:.8f}\n")
    _record("b", "max |cf_FE - cf_march|/cf_march at matched Re_theta",
            "<= 5%", f"{100.0 * float(rel.max()):.2f}% "
            f"(Re_theta in [{re_grid[0]:.0f}, {re_grid[-1]:.0f}])",
            float(rel.max()) <= 0.05)

    # --- (d) lock on the turbulent solution ---
    B = float(np.max(np.abs(fe["U"][:, 2])))
    Psi = float(np.max(np.abs(fe["U"][:, 3])))
    Ct2 = float(np.max(np.abs(fe["U"][:, 5])))
    ok_d = B < 1.0e-10 and Psi < 1.0e-10 and Ct2 < 1.0e-10
    _record("d", "quasi-2-D lock (turbulent)", "|B|,|Psi|,|Ctau2|<1e-10",
            f"|B|={B:.2e} |Psi|={Psi:.2e} |Ctau2|={Ct2:.2e}", ok_d)

    # profile CSV for the plot
    with open(os.path.join(RESULTS, "gv1_1b_profiles.csv"), "w") as f:
        f.write("x,cf_FE,H_FE,theta_FE,cf_march,H_march,theta_march\n")
        cf_ma = np.interp(xc, march["x"], march["cf"])
        H_ma = np.interp(xc, march["x"], march["H"])
        th_ma = np.interp(xc, march["x"], march["th"])
        for i in range(xc.size):
            f.write(f"{xc[i]:.4f},{cf[i]:.8f},{H[i]:.6f},{th[i]:.8f},"
                    f"{cf_ma[i]:.8f},{H_ma[i]:.6f},{th_ma[i]:.8f}\n")
    return fe, march, (re_grid, cf_fe_i, cf_ma_i, cf_pw)


# ---------------------------------------------------------------------------
# Gate (c): Falkner-Skan decelerating branch
# ---------------------------------------------------------------------------

def gate_c():
    print("== GV1.1(c): decelerating u_e, FS branch ==")
    x0, x1 = 0.4, 1.4
    results = {}
    for m in (-0.05, -0.0904):
        def ue_fn(x, m=m):
            return np.array([x ** m, 0.0, 0.0])

        def seed_fn(x, m=m):
            return C.blasius_seed(x, q=max(x ** m, 1.0e-3), rho=RHO, mu=MU)

        fe = run_fe(80, 8, False, ue_fn, seed_fn, x0=x0, x1=x1)
        print(f"  m={m}: converged={fe['info']['converged']} "
              f"it={fe['info']['n_iter']} |R|={fe['info']['final_residual']:.2e}")
        assert fe["info"]["converged"], f"FS run failed for m={m}"
        xc, H = _centerline(fe, C.OUT_H1)
        results[m] = (xc, H, fe)
    with open(os.path.join(RESULTS, "gv1_1c_fs.csv"), "w") as f:
        f.write("m,x,H\n")
        for m, (xc, H, _fe) in results.items():
            for i in range(xc.size):
                f.write(f"{m},{xc[i]:.4f},{H[i]:.8f}\n")
    # P1: monotone non-decreasing (scatter < 0.5% of total rise)
    xc, Hm, _ = results[-0.0904]
    rise = Hm[-1] - Hm[0]
    dips = np.diff(Hm)
    ok_p1 = rise > 0.0 and np.all(dips > -0.005 * max(rise, 1.0e-12))
    _record("c/P1", "H(x) monotone non-decreasing (m=-0.0904)",
            "monotone within 0.5% of total rise",
            f"rise={rise:.4f}, worst dip={float(dips.min()):.2e}", ok_p1)
    # P2: rise ratio band
    ratio = float(Hm[-1] / Hm[0])
    _record("c/P2", "H rise ratio (m=-0.0904)", "[1.05, 4.0]",
            f"{ratio:.3f} (H {Hm[0]:.3f} -> {Hm[-1]:.3f})",
            1.05 <= ratio <= 4.0)
    # P3: trend with m
    _, Hw, _ = results[-0.05]
    ratio_w = float(Hw[-1] / Hw[0])
    _record("c/P3", "m=-0.05 rise < m=-0.0904 rise", "strictly smaller",
            f"{ratio_w:.3f} < {ratio:.3f}", ratio_w < ratio)
    return results


# ---------------------------------------------------------------------------
# Figures + summary
# ---------------------------------------------------------------------------

def make_figure(fes_a, march_a, b_pack, c_res):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    for (nx, nz), fe in zip(((50, 8), (100, 16), (200, 32)), fes_a):
        xc, H = _centerline(fe, C.OUT_H1)
        ax.plot(xc, H, ".", ms=2, label=f"FE {nx}x{nz}")
    ax.plot(march_a["x"], march_a["H"], "k-", lw=1, label="2-D march (closure)")
    ax.axhline(2.59, color="r", ls="--", lw=1, label="Blasius 2.59")
    ax.axhspan(2.59 * 0.98, 2.59 * 1.02, color="r", alpha=0.08)
    ax.set_xlabel("x")
    ax.set_ylabel("H")
    ax.set_title("(a) laminar plate H(x)")
    ax.legend(fontsize=7)
    ax = axes[0, 1]
    re_g, cf_fe, cf_ma, cf_pw = b_pack
    ax.loglog(re_g, cf_fe, "b.-", ms=3, label="FE")
    ax.loglog(re_g, cf_ma, "g--", label="2-D march (closure ref)")
    ax.loglog(re_g, cf_pw, "k:", label="1/7th power law (RECORDED)")
    ax.set_xlabel(r"Re_theta")
    ax.set_ylabel("cf")
    ax.set_title("(b) turbulent plate cf(Re_theta)")
    ax.legend(fontsize=7)
    ax = axes[1, 0]
    for m, (xc, H, _fe) in c_res.items():
        ax.plot(xc, H, ".-", ms=3, label=f"m={m}")
    ax.set_xlabel("x")
    ax.set_ylabel("H")
    ax.set_title("(c) FS decelerating branch: H rise indicator")
    ax.legend(fontsize=7)
    ax = axes[1, 1]
    dat = np.loadtxt(os.path.join(RESULTS, "gv1_1e_refinement.csv"),
                     delimiter=",", skiprows=1)
    ax.loglog(dat[:, 0], dat[:, 1], "o-", label="max |H - march|")
    ax.loglog(dat[:, 0], dat[:, 2], "s-", label="max |delta* - march|")
    ax.loglog(dat[:, 0], dat[:, 1][0] * dat[:, 0] / dat[0, 0], "k:",
              label="O(h)")
    ax.set_xlabel("h")
    ax.set_ylabel("error")
    ax.set_title("(e) refinement (vs closure 2-D march)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "gv1_1_panels.png"), dpi=150)
    plt.close(fig)


def main():
    fes_a, march_a = gate_a_e()
    fe_b, march_b, b_pack = gate_b()
    c_res = gate_c()
    make_figure(fes_a, march_a, b_pack, c_res)
    with open(os.path.join(RESULTS, "summary.csv"), "w") as f:
        f.write("gate,metric,band,measured,verdict\n")
        for g, m, b, meas, v in SUMMARY:
            f.write(f"{g},\"{m}\",\"{b}\",\"{meas}\",{v}\n")
    n_fail = sum(1 for s in SUMMARY if s[4] == "FAIL")
    print(f"\nGV1.1 summary: {len(SUMMARY) - n_fail} PASS, {n_fail} FAIL")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
