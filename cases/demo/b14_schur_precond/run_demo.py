"""
Track B / B14 demo -- `precond="schur"`: Schur-eliminated aux block + AMG on
the SPD Picard main block (pyfp3d/solve/schur_ls.py; design_track_b.md §5.3,
roadmap §B14).

What it demonstrates, on the roadmap's pre-registered terms:

  * GB14.1 diagnostic-first: J_aa (the aux thin-strip block: wake-LS g1+g2 +
    nonlinear TE-Kutta rows) FACTORS and its measured 1-norm condition
    estimate is finite -- "measure, don't assume".
  * GB14.2 correctness: schur lands on the spsolve gamma (|dgamma| < 1e-8) on
    the 2.5D coarse lifting + Newton solves, 0 stalls / 0 fallbacks -- on the
    exact operator where the B11 spring surrogate stalled to gamma 0.0033.
  * GB14.3 the discriminating tier: 2.5D MEDIUM lifting, where ILU DIVERGED
    (gamma = -136.99, 77 stalls, B11 solver_ab.csv). "Passing there is what
    'a real escape' means" (roadmap B14).
  * GB14.4 3D capability (gated): ONERA M6 wake-free COARSE + MEDIUM, subsonic
    M0.5 lifting AND transonic M0.84 Newton ramp, A/B against the committed
    lagged-LU recipes in the same session; physics locks vs the committed
    GB15.4 state (gamma 0.088338, M_max 2.4938).
  * Timing analysis (RECORDED, NOT a gate -- the design says the medium-scale
    gain over lagged-LU is uncertain; the unique designed value is the FINE
    memory-bounded path, out of scope here): wall-clock A/B, per-phase
    breakdown (the A1 "LS Newton is 42.6% precond-bound" target), GMRES
    iters/step, fallback counts.

Artifacts: results/checks.csv, results/jaa_diag.csv, results/schur_ab.csv,
results/b14_*.png. Heavy M6-medium states are cached to gitignored
results/*.npz (delete to force a re-solve).

Usage:
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
        python cases/demo/b14_schur_precond/run_demo.py
    # parts 4-5 (M6 coarse+medium, ~40-60 min cold) additionally need
    PYFP3D_TRANSONIC_GATES=1
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE,
    CRITICAL,
    CheckList,
    INK_2,
    MESH_DIR,
    S1_BLUE,
    S2_AQUA,
    S3_YELLOW,
    S4_ROSE,
    S5_VIOLET,
    apply_style,
    finish,
    write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from pyfp3d.constraints.dirichlet import freestream_phi  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_le, x_te  # noqa: E402
from pyfp3d.solve.newton_ls import (  # noqa: E402
    B_NEWTON_M6_DEFAULTS,
    solve_multivalued_newton,
    solve_multivalued_newton_transonic,
)
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.solve.schur_ls import (  # noqa: E402
    SchurReducedSystem,
    jaa_diagnostic,
)
from pyfp3d.wake import (  # noqa: E402
    CutElementMap,
    MultivaluedOperator,
    WakeLevelSet,
)

apply_style()
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"

NACA_DIR = MESH_DIR / "naca0012_2.5d"
M6_DIR = MESH_DIR / "onera_m6_wakefree"
ALPHA_NACA = 2.0
ALPHA_M6 = 3.06
M6_SUB, M6_TRANS = 0.5, 0.84

# Committed context numbers -- quoted, never re-derived here.
ILU_MEDIUM_GAMMA = -136.99      # B11 solver_ab.csv: the divergence B14 escapes
# B21 re-baseline (2026-07-19): the committed GB15.4 state is now the post-N1
# freeze-capture-fix one (n1_freeze_fix_sweep.csv; pre-B20 values were
# gamma 0.088338 / M_max 2.4938 / 657.4 s -- the B20-era stall at M0.6625 was
# B20's own patch gap, overturned by B21).
GB15_GAMMA = 0.088343           # committed M6 medium M0.84 Newton-ramp state
GB15_MMAX = 2.4818              # (c1_ls_jacobian_fd/results/n1_freeze_fix_sweep.csv)
GB15_WALL = 515.3               # its committed wall clock (different session)
A1_PRECOND_PCT = 42.6           # A1: LS Newton M6 medium M0.84 precond share

checks = CheckList("B14 -- Schur+AMG structural preconditioner")
ab_rows = []
AB_HEADER = ("part,case,m_inf,method,wall_s,seed_s,assembly_s,precond_s,"
             "linsolve_s,precond_pct,n_steps,n_lin_iters,n_stalled,"
             "n_fallback,n_refactor,gamma,dgamma_vs_ref,m_max,converged,"
             "target_reached")


# ---------------------------------------------------------------------------
# Case construction (the TE polyline MUST come from the authoritative
# geometry -- a hand-rolled x_te off by ~2e-4 matches 0 TE nodes -> NaN).
# ---------------------------------------------------------------------------

def mvop_naca(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def mvop_m6(mesh):
    a = np.radians(ALPHA_M6)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def fused_free_system(mesh, mvop, alpha_deg):
    """Freestream-state wake_ls matrix, far-field-Dirichlet-reduced."""
    phi_ext = np.zeros(mvop.n_total)
    phi_ext[:mvop.n_main] = freestream_phi(mesh.nodes, alpha_deg, 1.0)
    cut_nodes = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    phi_ext[mvop.cm.ext_dof_of_node[cut_nodes]] = phi_ext[cut_nodes]
    A = mvop.assemble_matrix(closure="wake_ls", te_kutta="pressure",
                             phi_ext=phi_ext).tocsr()
    ff = np.unique(mesh.boundary_faces["farfield"])
    is_dir = np.zeros(mvop.n_total, dtype=bool)
    is_dir[ff] = True
    free = np.flatnonzero(~is_dir)
    return A[free][:, free], free


def add_row(part, case, m_inf, method, r, wall, gamma_ref=None):
    """One schur_ab.csv row from a solver result (lifting or ramp)."""
    tt = r.get("timings_total", r.get("timings", {}))
    if "levels" in r:      # ramp: totals over levels
        lv = r["levels"]
        n_steps = sum(l["n_newton"] for l in lv)
        n_iters = sum(l["n_lin_iters"] for l in lv)
        n_refac = sum(l["n_refactor"] for l in lv)
        n_fall = sum(l.get("n_schur_fallback", 0) for l in lv)
        n_stall = n_fall  # per-level stall counts are not recorded; fallbacks
        conv = all(l["converged"] for l in lv)
        target = r.get("target_reached", "")
    else:
        n_steps = r.get("n_outer", r.get("n_newton", 0))
        n_iters = r["n_gmres_total"]
        n_refac = r["n_refactor"]
        n_fall = r.get("n_schur_fallback", 0)
        n_stall = r["n_gmres_stalled"]
        conv = r["converged"]
        target = ""
    pre = tt.get("precond", 0.0)
    dg = "" if gamma_ref is None else f"{abs(r['gamma'] - gamma_ref):.3e}"
    ab_rows.append((part, case, m_inf, method, f"{wall:.1f}",
                    f"{tt.get('seed', 0.0):.1f}",
                    f"{tt.get('assembly', 0.0):.1f}", f"{pre:.1f}",
                    f"{tt.get('linsolve', 0.0):.1f}",
                    f"{100.0 * pre / wall:.1f}", n_steps, n_iters, n_stall,
                    n_fall, n_refac, f"{r['gamma']:.8f}", dg,
                    f"{float(np.sqrt(r['mach2_max'])):.4f}", int(conv),
                    target))
    return ab_rows[-1]


# ---------------------------------------------------------------------------
# Part 1 -- GB14.1: the pre-registered J_aa diagnostic.
# ---------------------------------------------------------------------------
print("=== Part 1 (GB14.1): J_aa invertibility/conditioning -- measured ===")
diag_rows = []
diag_cases = [("naca 2.5D", "coarse", NACA_DIR / "coarse.msh", mvop_naca,
               ALPHA_NACA),
              ("naca 2.5D", "medium", NACA_DIR / "medium.msh", mvop_naca,
               ALPHA_NACA),
              ("M6 wake-free", "coarse", M6_DIR / "coarse.msh", mvop_m6,
               ALPHA_M6)]
if GATED:
    diag_cases.append(("M6 wake-free", "medium", M6_DIR / "medium.msh",
                       mvop_m6, ALPHA_M6))
for case, level, path, factory, alpha in diag_cases:
    if not path.exists():
        print(f"  [{case} {level}] SKIP (mesh not generated)")
        continue
    mesh = read_mesh(path)
    mv = factory(mesh)
    A_free, free = fused_free_system(mesh, mv, alpha)
    t0 = time.perf_counter()
    schur = SchurReducedSystem(A_free, free, mv.n_main,
                               n_aux_expected=mv.n_ext)
    t_ms = 1e3 * (time.perf_counter() - t0)
    d = jaa_diagnostic(schur)
    diag_rows.append((case, level, schur.n_mf, schur.n_aux,
                      f"{d['cond1']:.3e}", 1, f"{t_ms:.1f}"))
    print(f"  [{case} {level}] n_aux={schur.n_aux} cond1={d['cond1']:.3e} "
          f"split+splu {t_ms:.0f} ms")
write_csv(OUT, "jaa_diag.csv",
          "case,level,n_main_free,n_aux,cond1_est,factor_ok,split_factor_ms",
          diag_rows)
checks.add("GB14.1", "J_aa factors on every measured case, cond1 finite",
           "; ".join(f"{r[0]} {r[1]}: {r[4]}" for r in diag_rows),
           "splu succeeds, cond1 < 1e14",
           len(diag_rows) >= 3 and all(float(r[4]) < 1e14 for r in diag_rows))

# ---------------------------------------------------------------------------
# Part 2 -- GB14.2: 2.5D coarse correctness A/B (the surrogate's stall case).
# ---------------------------------------------------------------------------
print("\n=== Part 2 (GB14.2): 2.5D coarse -- schur vs spsolve ===")
naca_coarse = read_mesh(NACA_DIR / "coarse.msh")

t0 = time.perf_counter()
ref = solve_multivalued_lifting(mvop_naca(naca_coarse), naca_coarse, 0.5,
                                alpha_deg=ALPHA_NACA)
add_row(2, "naca coarse lifting", 0.5, "spsolve", ref,
        time.perf_counter() - t0)
t0 = time.perf_counter()
s = solve_multivalued_lifting(mvop_naca(naca_coarse), naca_coarse, 0.5,
                              alpha_deg=ALPHA_NACA, precond="schur")
add_row(2, "naca coarse lifting", 0.5, "schur", s, time.perf_counter() - t0,
        gamma_ref=ref["gamma"])
dg_lift = abs(s["gamma"] - ref["gamma"])
ok_lift = (ref["converged"] and s["converged"] and dg_lift < 1e-8
           and s["n_gmres_stalled"] == 0 and s["n_schur_fallback"] == 0)
print(f"  lifting: gamma {ref['gamma']:.8f} -> {s['gamma']:.8f} "
      f"(|dgamma|={dg_lift:.2e})")

NEWTON_COMMON = dict(m_inf=0.7, alpha_deg=ALPHA_NACA, farfield="neumann",
                     n_seed=20, n_newton_max=25)
t0 = time.perf_counter()
ref_n = solve_multivalued_newton(mesh=naca_coarse,
                                 mvop=mvop_naca(naca_coarse), **NEWTON_COMMON)
add_row(2, "naca coarse newton M0.7", 0.7, "spsolve", ref_n,
        time.perf_counter() - t0)
t0 = time.perf_counter()
s_n = solve_multivalued_newton(mesh=naca_coarse, mvop=mvop_naca(naca_coarse),
                               precond="schur", **NEWTON_COMMON)
add_row(2, "naca coarse newton M0.7", 0.7, "schur", s_n,
        time.perf_counter() - t0, gamma_ref=ref_n["gamma"])
dg_newton = abs(s_n["gamma"] - ref_n["gamma"])
ok_newton = (ref_n["converged"] and s_n["converged"] and dg_newton < 1e-8
             and s_n["n_gmres_stalled"] == 0 and s_n["n_schur_fallback"] == 0)
print(f"  newton : gamma {ref_n['gamma']:.8f} -> {s_n['gamma']:.8f} "
      f"(|dgamma|={dg_newton:.2e})")
checks.add("GB14.2", "coarse lifting: schur == spsolve gamma, no fallback",
           f"|dgamma|={dg_lift:.2e}", "< 1e-8, 0 stalls/fallbacks", ok_lift,
           note="the B11 surrogate stalled HERE to gamma 0.0033 vs 0.139")
checks.add("GB14.2", "coarse Newton M0.7: schur == spsolve gamma",
           f"|dgamma|={dg_newton:.2e}", "< 1e-8, 0 stalls/fallbacks",
           ok_newton,
           note="supersonic pocket: Terms 2/3 invisible to the Term-1 AMG")

# ---------------------------------------------------------------------------
# Part 3 -- GB14.3: the discriminating tier (2.5D medium lifting, ILU's
# grave). This is the pre-registered "what a real escape means" case.
# ---------------------------------------------------------------------------
print("\n=== Part 3 (GB14.3): 2.5D MEDIUM lifting -- the ILU-divergence "
      "tier ===")
naca_medium = read_mesh(NACA_DIR / "medium.msh")
t0 = time.perf_counter()
ref_m = solve_multivalued_lifting(mvop_naca(naca_medium), naca_medium, 0.5,
                                  alpha_deg=ALPHA_NACA)
add_row(3, "naca medium lifting", 0.5, "spsolve", ref_m,
        time.perf_counter() - t0)
t0 = time.perf_counter()
s_m = solve_multivalued_lifting(mvop_naca(naca_medium), naca_medium, 0.5,
                                alpha_deg=ALPHA_NACA, precond="schur")
t_sm = time.perf_counter() - t0
add_row(3, "naca medium lifting", 0.5, "schur", s_m, t_sm,
        gamma_ref=ref_m["gamma"])
dg_med = abs(s_m["gamma"] - ref_m["gamma"])
ok_med = (ref_m["converged"] and s_m["converged"] and dg_med < 1e-8
          and s_m["n_gmres_stalled"] == 0 and s_m["n_schur_fallback"] == 0)
print(f"  gamma {ref_m['gamma']:.8f} -> {s_m['gamma']:.8f} "
      f"(|dgamma|={dg_med:.2e}); ILU here: {ILU_MEDIUM_GAMMA} (diverged)")
checks.add("GB14.3", "2.5D medium lifting: schur CONVERGES where ILU diverged",
           f"gamma={s_m['gamma']:.8f}; |dgamma|={dg_med:.2e}",
           "converged, < 1e-8 vs spsolve, 0 stalls/fallbacks", ok_med,
           note=f"ILU committed: gamma={ILU_MEDIUM_GAMMA} 77 stalls "
                "(b11_ls_infra/solver_ab.csv)")

# ---------------------------------------------------------------------------
# Parts 4+5 -- GB14.4 (gated): ONERA M6 3D capability + the timing A/B.
# ---------------------------------------------------------------------------
SUB_KW = dict(m_inf=M6_SUB, alpha_deg=ALPHA_M6, farfield="neumann",
              n_outer_max=60, tol_residual=1e-7)
TRANS_KW = dict(m_target=M6_TRANS, alpha_deg=ALPHA_M6, farfield="neumann",
                n_seed=40, n_newton_max=80, tol_residual=1e-10)


def m6_pair(level, part):
    """Fresh same-session A/B (lagged-LU committed recipe vs schur) for one
    M6 level: subsonic lifting + transonic Newton ramp. The medium arms are
    cached to gitignored npz (the P5 policy: committed evidence = CSV/PNG)."""
    mesh = read_mesh(M6_DIR / f"{level}.msh")
    out = {}
    for regime, tag in (("sub", f"M6 {level} lifting"),
                        ("trans", f"M6 {level} ramp M0.84")):
        for method in ("lagged", "schur"):
            key = f"{level}_{regime}_{method}"
            cache = OUT / f"{key}.npz" if level == "medium" else None
            if cache is not None and cache.exists():
                d = np.load(cache, allow_pickle=True)
                r, wall = d["result"].item(), float(d["wall_s"])
                print(f"  [{key}] cached ({wall:.1f}s recorded)")
            else:
                t0 = time.perf_counter()
                if regime == "sub":
                    kw = dict(SUB_KW)
                    if method == "lagged":
                        kw["direct_refactor_every"] = 1000
                    else:
                        kw["precond"] = "schur"
                    r = solve_multivalued_lifting(mvop_m6(mesh), mesh, **kw)
                else:
                    kw = dict(TRANS_KW, **B_NEWTON_M6_DEFAULTS)
                    if method == "schur":
                        kw["precond"] = "schur"
                    r = solve_multivalued_newton_transonic(
                        mvop=mvop_m6(mesh), mesh=mesh, **kw)
                wall = time.perf_counter() - t0
                print(f"  [{key}] solved {wall:.1f}s gamma={r['gamma']:.6f}",
                      flush=True)
                if cache is not None:
                    slim = {k: v for k, v in r.items()
                            if k not in ("phi", "rho_tilde")}
                    np.savez(cache, result=np.array(slim, dtype=object),
                             wall_s=wall)
            out[f"{regime}_{method}"] = (r, wall)
        ref, _ = out[f"{regime}_lagged"]
        r, wall = out[f"{regime}_schur"]
        add_row(part, tag, M6_SUB if regime == "sub" else M6_TRANS, "lagged",
                ref, out[f"{regime}_lagged"][1])
        add_row(part, tag, M6_SUB if regime == "sub" else M6_TRANS, "schur",
                r, wall, gamma_ref=ref["gamma"])
    out["mesh"] = mesh
    return out


if GATED and (M6_DIR / "coarse.msh").exists():
    print("\n=== Part 4 (GB14.4): ONERA M6 COARSE A/B (gated) ===")
    c = m6_pair("coarse", 4)
    dg = abs(c["sub_schur"][0]["gamma"] - c["sub_lagged"][0]["gamma"])
    dgt = abs(c["trans_schur"][0]["gamma"] - c["trans_lagged"][0]["gamma"])
    checks.add("GB14.4", "M6 coarse: schur solves M0.5 AND the M0.84 ramp",
               f"sub |dgamma|={dg:.2e}; ramp |dgamma|={dgt:.2e}",
               "converged + target_reached, |dgamma| <= 1e-4 vs lagged",
               c["sub_schur"][0]["converged"]
               and c["trans_schur"][0]["target_reached"]
               and dg <= 1e-4 and dgt <= 1e-4)
else:
    print("\n=== Parts 4-5 SKIPPED -- set PYFP3D_TRANSONIC_GATES=1 (and "
          "generate cases/meshes/onera_m6_wakefree) ===")

if GATED and (M6_DIR / "medium.msh").exists():
    print("\n=== Part 5 (GB14.4 headline): ONERA M6 MEDIUM A/B (gated) ===")
    m = m6_pair("medium", 5)
    r_sub_l, w_sub_l = m["sub_lagged"]
    r_sub_s, w_sub_s = m["sub_schur"]
    r_trn_l, w_trn_l = m["trans_lagged"]
    r_trn_s, w_trn_s = m["trans_schur"]
    dg_sub = abs(r_sub_s["gamma"] - r_sub_l["gamma"])
    dg_trn = abs(r_trn_s["gamma"] - r_trn_l["gamma"])
    mmax_s = float(np.sqrt(r_trn_s["mach2_max"]))
    checks.add("GB14.4", "M6 medium: schur solves M0.5 AND the M0.84 ramp",
               f"sub |dgamma|={dg_sub:.2e}; ramp |dgamma|={dg_trn:.2e}",
               "converged + target_reached, |dgamma| <= 1e-4 vs lagged",
               r_sub_s["converged"] and r_trn_s["target_reached"]
               and dg_sub <= 1e-4 and dg_trn <= 1e-4)
    checks.add("GB14.4", "M6 medium M0.84 schur state == the committed GB15.4"
               " physics",
               f"gamma={r_trn_s['gamma']:.6f}, M_max={mmax_s:.4f}",
               f"within 1% of gamma {GB15_GAMMA}, M_max {GB15_MMAX}",
               abs(r_trn_s["gamma"] - GB15_GAMMA) / GB15_GAMMA <= 0.01
               and abs(mmax_s - GB15_MMAX) / GB15_MMAX <= 0.01)
    n_fall = sum(l.get("n_schur_fallback", 0) for l in r_trn_s["levels"])
    tt = r_trn_s["timings_total"]
    pre_pct = 100.0 * tt["precond"] / w_trn_s
    print(f"  RECORDED: ramp wall {w_trn_l:.0f}s (lagged) vs {w_trn_s:.0f}s "
          f"(schur) = {w_trn_l / w_trn_s:.2f}x; precond share "
          f"{A1_PRECOND_PCT}% (A1 lagged) -> {pre_pct:.1f}% (schur); "
          f"{n_fall} fallbacks")

# ---------------------------------------------------------------------------
# Figures.
# ---------------------------------------------------------------------------
print("\n=== Figures ===")


def _rows_for(part, case, method):
    for r in ab_rows:
        if r[0] == part and r[1] == case and r[3] == method:
            return r
    return None


# (1) wall-clock A/B bars over every case that ran.
fig, ax = plt.subplots(figsize=(9.5, 4.6))
cases_plot = []
for part, case in dict.fromkeys((r[0], r[1]) for r in ab_rows):
    base = (_rows_for(part, case, "spsolve") or
            _rows_for(part, case, "lagged"))
    schur_r = _rows_for(part, case, "schur")
    if base and schur_r:
        cases_plot.append((case, base, schur_r))
xs = np.arange(len(cases_plot))
b_base = ax.bar(xs - 0.19, [float(c[1][4]) for c in cases_plot], 0.38,
                color=BASELINE, label="baseline (spsolve / lagged-LU)")
b_sch = ax.bar(xs + 0.19, [float(c[2][4]) for c in cases_plot], 0.38,
               color=S1_BLUE, label='precond="schur" (B14)')
for x, (case, base, sch) in zip(xs, cases_plot):
    wb, ws = float(base[4]), float(sch[4])
    ax.text(x, max(wb, ws) * 1.03, f"{wb / ws:.2f}x", ha="center",
            fontsize=9, color=INK_2)
ax.set_xticks(xs)
ax.set_xticklabels([c[0].replace(" lifting", "\nlifting")
                    .replace(" ramp", "\nramp") for c in cases_plot],
                   fontsize=8.5)
ax.set_ylabel("wall clock (s)")
ax.set_yscale("log")
ax.grid(axis="y")
ax.set_title("B14: wall-clock A/B -- schur vs the committed direct/lagged-LU "
             "recipes\n(annotation = baseline/schur speedup; >1 = schur "
             "faster; TIMING IS RECORDED, NOT GATED)")
ax.legend()
finish(fig, OUT, "b14_wall_ab.png")

# (2) phase fractions (the A1 "precond-bound" claim, re-measured per arm).
PHASE_KEYS = ["seed", "assembly", "precond", "linsolve", "residual", "kutta"]
PHASE_COLORS = [S5_VIOLET, S3_YELLOW, CRITICAL, S1_BLUE, S2_AQUA, S4_ROSE,
                BASELINE]
fig, ax = plt.subplots(figsize=(9.0, 4.6))
frac_cases = [(c, b, s) for c, b, s in cases_plot
              if "M6" in c or "medium" in c]
labels, bottoms = [], None
xs = np.arange(len(frac_cases) * 2)
mat = []
for case, base, sch in frac_cases:
    for row in (base, sch):
        wall = float(row[4])
        secs = [float(row[5]), float(row[6]), float(row[7]), float(row[8])]
        seed, asm, pre, lin = secs
        other = max(wall - sum(secs), 0.0)
        mat.append(np.array([seed, asm, pre, lin, 0, 0, other]) / wall)
        labels.append(f"{case}\n{row[3]}")
mat = np.array(mat)
bottom = np.zeros(len(mat))
for j, (key, color) in enumerate(zip(PHASE_KEYS[:4] + ["", "", "other"],
                                     PHASE_COLORS)):
    if not key:
        continue
    ax.bar(np.arange(len(mat)), mat[:, j], 0.6, bottom=bottom, color=color,
           label=key)
    bottom += mat[:, j]
for i, row in enumerate(mat):
    ax.text(i, 1.02, f"precond\n{100 * row[2]:.0f}%", ha="center", fontsize=8,
            color=CRITICAL)
ax.set_xticks(np.arange(len(mat)))
ax.set_xticklabels(labels, fontsize=7, rotation=20, ha="right")
ax.grid(axis="y")
ax.set_ylabel("fraction of wall clock")
ax.set_ylim(0, 1.18)
ax.set_title("B14: per-phase wall-clock fractions -- the A1 precond "
             "bottleneck (42.6% at M6 medium M0.84)\nvs the schur arm "
             "(factorizations replaced by thin-strip LU + AMG V-cycles)")
ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.18))
finish(fig, OUT, "b14_phase_fractions.png")

if GATED and (M6_DIR / "medium.msh").exists():
    # (3) GMRES iterations per Newton step across the schur medium ramp.
    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    x0 = 0
    lv_colors = [S1_BLUE, S2_AQUA, S3_YELLOW, S4_ROSE, S5_VIOLET, INK_2]
    for li, lv in enumerate(r_trn_s["levels"]):
        iters = [sr.get("n_lin_iters", 0) for sr in lv["step_records"]]
        xs = np.arange(x0, x0 + len(iters))
        ax.plot(xs, iters, "o-", ms=3.5, lw=1.2,
                color=lv_colors[li % len(lv_colors)],
                label=f"M{lv['m_inf']:.2f} ({lv['n_newton']} steps)")
        x0 += len(iters)
    ax.set_xlabel("Newton step (cumulative across ramp levels)")
    ax.set_ylabel("GMRES iterations (reduced system)")
    ax.set_title("B14: reduced-GMRES iterations per Newton step, M6 medium "
                 "M0.84 schur ramp\n(no full-size factorization anywhere; "
                 "0 fallbacks = no step ever fell back to spsolve)")
    ax.legend(fontsize=8, ncol=3)
    finish(fig, OUT, "b14_gmres_iters.png")

    # (4) final-level residual histories, lagged vs schur.
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for r, label, color, ls in ((r_trn_l, "lagged-LU (committed recipe)",
                                 BASELINE, "-"),
                                (r_trn_s, 'precond="schur"', S1_BLUE, "--")):
        hist = r["levels"][-1]["residual_history"]
        ax.semilogy(hist, ls, color=color, lw=2.0, label=label)
    ax.set_xlabel("Newton step (final ramp level, M0.84)")
    ax.set_ylabel("|R|_inf (free DOFs)")
    ax.set_title("B14: the two arms walk the same Newton path at M0.84\n"
                 "(same freeze recipe; the linear solver is the only change)")
    ax.legend()
    finish(fig, OUT, "b14_residual_ab.png")

    # (5) the computed field: spanwise circulation A/B + tip Mach (schur).
    mesh = m["mesh"]
    mv_fig = mvop_m6(mesh)
    te_z = mesh.nodes[np.asarray(mv_fig.cm.te_nodes)][:, 2]
    o = np.argsort(te_z)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax = axes[0]
    for r, label, color, ls, m_val in (
            (r_sub_l, "M0.5 lagged", BASELINE, "-", M6_SUB),
            (r_sub_s, "M0.5 schur", S2_AQUA, "--", M6_SUB),
            (r_trn_l, "M0.84 lagged", "#888888", "-", M6_TRANS),
            (r_trn_s, "M0.84 schur", S1_BLUE, "--", M6_TRANS)):
        gz = mv_fig.te_jump(r["phi_ext"])
        ax.plot(te_z[o] / B_SEMI, gz[o], ls, color=color, lw=1.8,
                label=f"{label} (gamma={r['gamma']:.4f})")
    ax.axhline(0, color=BASELINE, lw=0.8)
    ax.set_xlabel("z / b_semi")
    ax.set_ylabel("gamma(z) = TE jump")
    ax.set_title("spanwise circulation: the arms overlay")
    ax.legend(fontsize=8)
    ax = axes[1]
    cen = mesh.nodes[mesh.elements].mean(axis=1)
    mach = np.sqrt(mv_fig.element_mach2(r_trn_s["phi_ext"], M6_TRANS,
                                        u_inf=1.0))
    msk = ((np.abs(cen[:, 1]) < 0.025) & (cen[:, 2] > -0.05)
           & (cen[:, 2] < 1.1 * B_SEMI) & (cen[:, 0] > -0.4)
           & (cen[:, 0] < 2.0))
    norm = TwoSlopeNorm(vmin=0.3, vcenter=1.0,
                        vmax=max(1.05, min(mach[msk].max(), 3.0)))
    sc = ax.scatter(cen[msk, 0], cen[msk, 2], c=mach[msk], s=7,
                    cmap="RdBu_r", norm=norm)
    zz = np.linspace(0.0, B_SEMI, 40)
    ax.plot([x_le(z) for z in zz], zz, color="0.35", lw=1.0)
    ax.plot([x_te(z) for z in zz], zz, color="0.35", lw=1.0)
    ax.set_xlabel("x")
    ax.set_ylabel("z (span)")
    mmax_s = float(np.sqrt(r_trn_s["mach2_max"]))
    ax.set_title(f"local Mach, schur M0.84 solution (M_max={mmax_s:.3f})")
    fig.colorbar(sc, ax=ax, shrink=0.85, label="local Mach (red = M>1)")
    fig.suptitle("B14: the M6 medium fields computed BY the schur path", y=1.03)
    finish(fig, OUT, "b14_field.png")

write_csv(OUT, "schur_ab.csv", AB_HEADER, ab_rows)
code = checks.report(OUT)
sys.exit(code)
