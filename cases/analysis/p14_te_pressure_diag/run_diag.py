"""
P14 Stage-D diagnostic -- the spec-mandated go/no-go BEFORE wiring the
pressure-equality Kutta residual into any solver (roadmap/track_p.md P14:
"Diagnostic-first at open: build the conforming TE control volumes and
verify the recovered two-sided velocity is non-degenerate in Gamma before
wiring the residual").

Legs (all cheap; the only solve is the seconds-scale NACA coarse Laplace):

  D1  CV construction on NACA coarse + M6 coarse (+ M6 medium when the
      local caches are present): every TEControlVolumes constructor assert
      (two-sided non-empty wall-adjacent fans via exact wall-face
      ownership, probe-membership side identity, zero far-field Dirichlet
      contact) + fan-size stats.
  D2  Gamma non-degeneracy + Jacobian plumbing: dense D = dF/dGamma at
      freestream and at converged states -- sign uniformity, diagonal
      dominance margin, and a central-difference FD check of D at fixed
      phi_red (F is exactly quadratic in Gamma, so FD is exact to
      roundoff; this validates the whole _rows_cut @ G chain).
  D3  Estimator plausibility on converged probe-path states: implied
      Gamma* vs the cached Gamma / probe targets (expected: percent-scale
      SMOOTH offset -- the pressure closure is a *different* closure, the
      offset IS S2; gate-grade agreement comes from the Stage-1/2 solves,
      not here) + roughness_d2 preview of Gamma*(z) vs the probe-target
      curve on the SAME cached field (S1 preview).

Run:  python cases/analysis/p14_te_pressure_diag/run_diag.py
Meshes/caches are gitignored; absent legs are SKIPPED with a message.
Committed evidence: results/diag_checks.csv + results/diag_states.csv.
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "cases/analysis/a2_te_kutta_fidelity"))

import _metrics as M                                        # noqa: E402
from pyfp3d.constraints.dirichlet import farfield_dirichlet  # noqa: E402
from pyfp3d.constraints.te_pressure import TEControlVolumes  # noqa: E402
from pyfp3d.constraints.wake import WakeConstraint, kutta_targets  # noqa: E402
from pyfp3d.kernels.residual import assemble_stiffness_matrix  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                     # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                    # noqa: E402
from pyfp3d.solve.picard import solve_laplace_lifting        # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
MESHES = REPO_ROOT / "cases/meshes"
P5_RES = REPO_ROOT / "cases/demo/p5_onera_m6/results"
A1_RES = REPO_ROOT / "cases/analysis/a1_solver_bottleneck/results"

ALPHA_M6 = 3.06
checks = []          # (check, value, criterion, ok)
state_rows = []      # per-state summary CSV rows


def check(name, value, crit, ok):
    checks.append({"check": name, "value": value, "criterion": crit,
                   "ok": "PASS" if ok else "FAIL"})
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value}  ({crit})")


def build_cvs(mc, wc, alpha_deg):
    """CVs with the Dirichlet-contact guard armed (D1)."""
    dir_nodes, _ = farfield_dirichlet(mc, wc, alpha_deg,
                                      np.zeros(wc.n_stations))
    return TEControlVolumes(mc, wc, dirichlet_nodes=dir_nodes)


def fd_check_D(cvs, con, wc, phi_red, gamma, eps=1e-6):
    """Central-difference dF/dGamma at FIXED phi_red vs the exact dense D.

    F is exactly quadratic in Gamma (density never enters F), so central
    differences are exact to roundoff -- this validates the whole
    _rows_cut @ G plumbing, not just non-degeneracy."""
    phi_cut = con.expand(phi_red, gamma)
    D = cvs.gamma_jacobian(phi_cut, mode="exact")
    n_st = wc.n_stations
    D_fd = np.empty_like(D)
    for j in range(n_st):
        e = np.zeros(n_st)
        e[j] = eps
        Fp = cvs.residual_stations(con.expand(phi_red, gamma + e))
        Fm = cvs.residual_stations(con.expand(phi_red, gamma - e))
        D_fd[:, j] = (Fp - Fm) / (2.0 * eps)
    scale = max(np.abs(D).max(), 1e-30)
    return float(np.abs(D - D_fd).max() / scale), D


def d_stats(D):
    """Sign uniformity + diagonal dominance of the dense Gamma Jacobian."""
    dj = np.diag(D)
    off = np.abs(D).sum(axis=1) - np.abs(dj)
    dom = np.abs(dj) - off
    bw = 0
    nz = np.argwhere(np.abs(D) > 1e-12 * np.abs(D).max())
    if len(nz):
        bw = int(np.abs(nz[:, 0] - nz[:, 1]).max())
    return {"djj_min": float(np.abs(dj).min()),
            "djj_max": float(np.abs(dj).max()),
            "sign_uniform": bool(np.all(np.sign(dj) == np.sign(dj[0]))),
            "dom_margin_min": float(dom.min()),
            "bandwidth": bw,
            "cond": float(np.linalg.cond(D))}


# ---------------------------------------------------------------- D1+D2 NACA
print("== D1/D2/D3: NACA0012 coarse (live Laplace, alpha=2) ==")
naca_path = MESHES / "naca0012_2.5d/coarse.msh"
mc, wc = cut_wake(read_mesh(str(naca_path)))
t0 = time.perf_counter()
cvs = build_cvs(mc, wc, alpha_deg=2.0)
fs = cvs.fan_stats()
print(f"  CV build ok ({time.perf_counter() - t0:.2f}s): {fs}")
check("naca_cv_asserts", "all constructor asserts passed", "no raise", True)

r = solve_laplace_lifting(mc, wc, alpha_deg=2.0)
gamma = r["gamma"]
phi = r["phi"]
con = WakeConstraint(assemble_stiffness_matrix(mc.nodes, mc.elements), wc)
phi_red = phi[: wc.n_nodes_orig]

fd_err, D = fd_check_D(cvs, con, wc, phi_red, gamma)
ds = d_stats(D)
check("naca_fd_D_relerr", f"{fd_err:.2e}", "< 1e-9 (F exactly quadratic)",
      fd_err < 1e-9)
check("naca_D_nondegenerate", f"|D_jj| min {ds['djj_min']:.3g}",
      "> 0 with margin", ds["djj_min"] > 1e-6)

g_star = cvs.implied_targets(con.expand(phi_red, gamma), gamma)
tgt = kutta_targets(phi, wc)
rel = float(np.abs(g_star - gamma).max() / np.abs(gamma).max())
print(f"  gamma={gamma}, probe target={tgt}, implied gamma*={g_star}")
check("naca_gamma_star_band", f"{100 * rel:.2f}%",
      "< 10% of converged probe gamma (plausibility; S2 offset expected)",
      rel < 0.10)
state_rows.append({"state": "naca_coarse_a2_laplace", "n_st": wc.n_stations,
                   **fs, **ds, "fd_relerr": fd_err,
                   "gamma_ref": float(gamma[0]),
                   "gamma_star_reldev": rel,
                   "rough_probe": np.nan, "rough_star": np.nan})

# ------------------------------------------------------- D1+D2 M6 freestream
print("== D1/D2: ONERA M6 coarse (freestream phi = x) ==")
m6_path = MESHES / "onera_m6/coarse.msh"
if m6_path.exists():
    mc6, wc6 = cut_wake(read_mesh(str(m6_path)))
    t0 = time.perf_counter()
    cvs6 = build_cvs(mc6, wc6, alpha_deg=ALPHA_M6)
    fs6 = cvs6.fan_stats()
    print(f"  CV build ok ({time.perf_counter() - t0:.2f}s): {fs6}")
    check("m6_coarse_cv_asserts", f"{fs6}", "no raise, fans non-empty", True)

    con6 = WakeConstraint(
        assemble_stiffness_matrix(mc6.nodes, mc6.elements), wc6)
    phi_fs = mc6.nodes[: wc6.n_nodes_orig, 0].copy()   # phi = x, reduced
    g0 = np.zeros(wc6.n_stations)
    fd_err6, D6 = fd_check_D(cvs6, con6, wc6, phi_fs, g0)
    ds6 = d_stats(D6)
    check("m6_coarse_fd_D_relerr", f"{fd_err6:.2e}", "< 1e-9",
          fd_err6 < 1e-9)
    check("m6_coarse_D_sign_uniform", str(ds6["sign_uniform"]),
          "uniform D_jj sign", ds6["sign_uniform"])
    check("m6_coarse_D_well_conditioned",
          f"cond {ds6['cond']:.1f}, min dominance margin "
          f"{ds6['dom_margin_min']:.2f}, "
          f"|D_jj| in [{ds6['djj_min']:.1f}, {ds6['djj_max']:.1f}], "
          f"bandwidth {ds6['bandwidth']}",
          "cond(D) < 100 (dominance is only sufficient; D is solved dense)",
          ds6["cond"] < 100.0)
    state_rows.append({"state": "m6_coarse_freestream",
                       "n_st": wc6.n_stations, **fs6, **ds6,
                       "fd_relerr": fd_err6, "gamma_ref": 0.0,
                       "gamma_star_reldev": np.nan,
                       "rough_probe": np.nan, "rough_star": np.nan})
else:
    print(f"[skip] {m6_path} missing (generate_onera_m6.py ~30 s)")

# ------------------------------------------------- D2+D3 cached M6 states
CACHED = {
    "conf_picard_coarse": (P5_RES / "coarse_solution.npz",
                           MESHES / "onera_m6/coarse.msh"),
    "conf_picard_medium": (P5_RES / "medium_solution.npz",
                           MESHES / "onera_m6/medium.msh"),
    "conf_newton_medium": (A1_RES / "a1_m6_conf_newton.npz",
                           MESHES / "onera_m6/medium.msh"),
}
_cut_cache = {}
for name, (npz_path, mesh_path) in CACHED.items():
    print(f"== D2/D3: {name} ==")
    if not npz_path.exists() or not mesh_path.exists():
        print(f"[skip] {name}: missing "
              f"{npz_path.name if not npz_path.exists() else mesh_path}")
        continue
    key = str(mesh_path)
    if key not in _cut_cache:
        mcx, wcx = cut_wake(read_mesh(key))
        conx = WakeConstraint(
            assemble_stiffness_matrix(mcx.nodes, mcx.elements), wcx)
        cvx = build_cvs(mcx, wcx, alpha_deg=ALPHA_M6)
        _cut_cache[key] = (mcx, wcx, conx, cvx)
    mcx, wcx, conx, cvx = _cut_cache[key]
    d = np.load(npz_path, allow_pickle=True)
    if "_span_gamma" in d.files:
        gx = np.atleast_1d(np.asarray(d["_span_gamma"], dtype=np.float64))
    else:
        gx = np.asarray(d["gamma"], dtype=np.float64)
    phix = np.asarray(d["phi"], dtype=np.float64)
    if len(gx) != wcx.n_stations:
        print(f"[skip] {name}: station count mismatch")
        continue
    phix_red = phix[: wcx.n_nodes_orig]

    phi_cut = conx.expand(phix_red, gx)
    D_x = cvx.gamma_jacobian(phi_cut, mode="exact")
    dsx = d_stats(D_x)
    check(f"{name}_D_sign_uniform", str(dsx["sign_uniform"]),
          "uniform D_jj sign at converged state", dsx["sign_uniform"])
    check(f"{name}_D_nondegenerate",
          f"|D_jj| in [{dsx['djj_min']:.1f}, {dsx['djj_max']:.1f}], "
          f"min dominance margin {dsx['dom_margin_min']:.2f}",
          "|D_jj| bounded away from 0", dsx["djj_min"] > 1.0)

    g_star = cvx.implied_targets(phi_cut, gx)
    tgtx = kutta_targets(phi_cut, wcx)
    o = np.argsort(wcx.station_z)
    rough_tgt = M.roughness_d2(tgtx[o])
    rough_star = M.roughness_d2(g_star[o])
    reldev = float(np.median(np.abs(g_star - gx))
                   / max(np.abs(gx).max(), 1e-30))
    print(f"  median|gamma*-gamma|/max|gamma| = {100 * reldev:.2f}%  "
          f"roughness probe-target {rough_tgt:.4f} vs gamma* {rough_star:.4f}")
    check(f"{name}_gamma_star_plausible", f"{100 * reldev:.2f}%",
          "< 15% median (S2 offset expected, not a closure yet)",
          reldev < 0.15)
    check(f"{name}_gamma_star_smoother",
          f"{rough_tgt:.4f} -> {rough_star:.4f}",
          "gamma*(z) roughness <= probe-target roughness (S1 preview)",
          rough_star <= rough_tgt)
    state_rows.append({"state": name, "n_st": wcx.n_stations,
                       **cvx.fan_stats(), **dsx, "fd_relerr": np.nan,
                       "gamma_ref": float(np.abs(gx).max()),
                       "gamma_star_reldev": reldev,
                       "rough_probe": rough_tgt, "rough_star": rough_star})

# ---------------------------------------------------------------- write CSVs
import csv                                                   # noqa: E402

with open(OUT / "diag_checks.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["check", "value", "criterion", "ok"])
    w.writeheader()
    w.writerows(checks)
if state_rows:
    with open(OUT / "diag_states.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(state_rows[0].keys()))
        w.writeheader()
        w.writerows(state_rows)

n_fail = sum(1 for c in checks if c["ok"] == "FAIL")
print(f"\n{len(checks) - n_fail}/{len(checks)} checks PASS "
      f"-> {OUT / 'diag_checks.csv'}")
sys.exit(1 if n_fail else 0)
