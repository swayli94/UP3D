"""
Track A / A2 -- GATED intervention leg: the H1 discriminator (GA2.2).

*** NOT RUN in the 2026-07-16 scaffold session (user-arbitrated). ***

Question. run_a2.py's zero-solve legs cannot decide H1 on their own: the
conforming field was SOLVED with the jittery per-station Gamma prescribed
across the wake, so re-reading any estimator on that self-consistent field
inherits the jitter. The decisive test is the T3/E fixed-Gamma methodology
(cases/demo/p5_onera_m6/INVESTIGATION_kutta_closure.md): re-solve the density
field warm-started from the cache with Gamma FIXED to a smooth diagnostic
profile, then re-measure the Kutta targets on the new field.

  D = roughness_d2( kutta_targets(phi') )  /  roughness_d2( Gamma_smooth )

Pre-registered verdict (roadmap/track_a.md A2, H1):
  D > 3    => H1 CONFIRMED -- the probe estimator REGENERATES the jitter
              from a spanwise-smooth field; the noise lives in the target
              definition (probe placement), not in the flow.
  D < 1.5  => H1 REFUTED  -- smooth Gamma yields smooth targets; the cached
              jitter was self-consistent field content (enforcement/solution).
  else     => mixed; report the split honestly.

Guard rails:
  * Gamma_smooth (a local-quadratic fit of the cached Gamma) is a DIAGNOSTIC
    INPUT in a fixed-Gamma solve -- the exact operation the P5 T-tests and
    test E performed. It is NOT the refuted smoothing-as-fix route
    (INVESTIGATION_gamma_smoothing.md): nothing here feeds a smoother into
    the live secant, and A2 proposes no fix.
  * Each state runs a CONTROL solve first (gamma_fixed = the cached Gamma):
    the warm start must reproduce the cached state's Kutta targets before
    the smooth-Gamma leg means anything (the T2/T3 "path check").
  * Warm-start fixed-Gamma solves only -- no from-scratch continuation
    (the "warm-start evidence does not validate a from-scratch recipe"
    lesson does not bite here: A2 makes no recipe claim).

Cost (T-series precedent): coarse ~1-3 min/solve, medium ~2-5 min/solve
warm-started; 2 solves per state (control + smooth). Run at the 16-thread
cap (agent-rules rule 1):

  NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
      PYFP3D_TRANSONIC_GATES=1 \
      python cases/analysis/a2_te_kutta_fidelity/run_a2_interventions.py

The coarse leg always runs; add PYFP3D_A2_MEDIUM=1 for the medium leg
(needs the P5 medium_solution.npz local cache, ~10 min total).
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(HERE))

import matplotlib.pyplot as plt                                    # noqa: E402

import _metrics as M                                               # noqa: E402
from cases.demo._common import (                                   # noqa: E402
    CheckList, S1_BLUE, S2_AQUA, CRITICAL, apply_style, finish, write_csv,
)
from pyfp3d.constraints.wake import kutta_targets                  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                           # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                          # noqa: E402
from pyfp3d.solve.continuation import TRANSONIC_DEFAULTS           # noqa: E402
from pyfp3d.solve.picard import solve_subsonic_lifting             # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
RUN_MEDIUM = os.environ.get("PYFP3D_A2_MEDIUM", "0") == "1"

M_INF, ALPHA = 0.84, 3.06
P5_RES = REPO_ROOT / "cases/demo/p5_onera_m6/results"
MESHES = REPO_ROOT / "cases/meshes"

STATES = [("coarse", P5_RES / "coarse_solution.npz",
           MESHES / "onera_m6/coarse.msh")]
if RUN_MEDIUM:
    STATES.append(("medium", P5_RES / "medium_solution.npz",
                   MESHES / "onera_m6/medium.msh"))

# Verdict thresholds -- PRE-REGISTERED in roadmap/track_a.md A2 (H1). Do not
# tune them after seeing the numbers.
D_CONFIRM, D_REFUTE = 3.0, 1.5


def fixed_gamma_solve(mc, wc, phi_init, gamma_fixed, n_max=400):
    """Warm-started fixed-Gamma density re-solve -- a verbatim mirror of
    continuation.solve_transonic_lifting's `_density_solve` (the operation
    every P5 T-test ran), at the P5 recipe's rtol=1e-7 and the landed
    polish's omega_rho=0.5 (T2: under-relaxed density converges the eval
    ~4x deeper at identical fixed points)."""
    return solve_subsonic_lifting(
        mc, wc, m_inf=M_INF, alpha_deg=ALPHA, u_inf=1.0,
        omega=1.0, upwind_c=TRANSONIC_DEFAULTS["upwind_c"],
        m_crit=TRANSONIC_DEFAULTS["m_crit"],
        damping_theta=TRANSONIC_DEFAULTS["damping_theta"],
        tol_rho=1e-8, n_picard_max=n_max, forcing=0.0,
        phi_init=phi_init, gamma_fixed=gamma_fixed,
        rtol=1e-7, maxiter=3000,
        farfield_spanwise_gamma=True, omega_rho=0.5,
    )


def run_state(cl, level, npz, mesh_path, rows, curves):
    if not npz.exists() or not mesh_path.exists():
        print(f"[skip] {level}: missing "
              f"{npz.name if not npz.exists() else mesh_path}")
        return
    d = np.load(npz, allow_pickle=True)
    gamma = np.asarray(d["gamma"], dtype=np.float64)
    z = np.asarray(d["station_z"], dtype=np.float64)
    phi = np.asarray(d["phi"], dtype=np.float64)
    mc, wc = cut_wake(read_mesh(str(mesh_path)))
    assert np.allclose(wc.station_z, z), "cut_wake no longer matches cache"
    o = np.argsort(z)

    g_smooth = M.local_fit(z, gamma)
    r_cached = M.roughness_d2(gamma[o])
    r_smooth = M.roughness_d2(g_smooth[o])

    # -- control: the warm start must reproduce the cached state's targets
    print(f"  [{level}] control solve (gamma_fixed = cached) ...", flush=True)
    t0 = time.perf_counter()
    rc = fixed_gamma_solve(mc, wc, phi, gamma)
    tgt_ctrl = kutta_targets(rc["phi"], wc)
    ctrl_dev = float(np.max(np.abs(tgt_ctrl - kutta_targets(phi, wc))))
    cl.add("GA2.2", f"{level}: control reproduces cached Kutta targets",
           f"max dev {ctrl_dev:.2e} ({time.perf_counter() - t0:.0f}s)",
           "small vs max|Gamma| (path check, T2/T3b)",
           ctrl_dev < 0.05 * float(np.max(np.abs(gamma))))

    # -- intervention: smooth Gamma in, what roughness comes back?
    print(f"  [{level}] intervention solve (gamma_fixed = smooth) ...",
          flush=True)
    t0 = time.perf_counter()
    ri = fixed_gamma_solve(mc, wc, phi, g_smooth)
    tgt_back = kutta_targets(ri["phi"], wc)
    r_back = M.roughness_d2(tgt_back[o])
    D = r_back / r_smooth if r_smooth > 0 else np.inf
    verdict = ("H1 CONFIRMED" if D > D_CONFIRM
               else "H1 REFUTED" if D < D_REFUTE else "MIXED")
    cl.add("GA2.2", f"{level}: discriminator D = r(targets back)/r(smooth in)",
           f"D={D:.2f} (r {r_smooth:.4f} -> {r_back:.4f}; cached {r_cached:.4f};"
           f" {time.perf_counter() - t0:.0f}s)",
           f">{D_CONFIRM} confirm / <{D_REFUTE} refute (pre-registered)",
           True, note=verdict)

    for i in range(len(z)):
        rows.append((level, f"{z[i]:.5f}", f"{gamma[i]:.6f}",
                     f"{g_smooth[i]:.6f}", f"{tgt_back[i]:.6f}",
                     f"{tgt_ctrl[i]:.6f}"))
    curves[level] = (z[o], gamma[o], g_smooth[o], tgt_back[o])


def main():
    if not GATED:
        print("A2 intervention leg is gated (2 transonic fixed-Gamma solves "
              "per level). Set PYFP3D_TRANSONIC_GATES=1 to run. Nothing done.")
        return 0
    apply_style()
    cl = CheckList("Track A / A2 -- GA2.2 fixed-Gamma intervention")
    rows, curves = [], {}
    for level, npz, mesh_path in STATES:
        run_state(cl, level, npz, mesh_path, rows, curves)
    if rows:
        write_csv(OUT, "a2_intervention.csv",
                  "level,z,gamma_cached,gamma_smooth_in,target_back_probe,"
                  "target_back_control", rows)
        fig, axes = plt.subplots(1, len(curves), figsize=(6.4 * len(curves), 4.8),
                                 squeeze=False)
        for ax, (level, (z, g, gs, tb)) in zip(axes[0], curves.items()):
            ax.plot(z, g, color=S1_BLUE, marker=".", ms=3, lw=1.2,
                    label="cached Γ (jittery)")
            ax.plot(z, gs, color=S2_AQUA, lw=2.0,
                    label="smooth Γ (diagnostic input)")
            ax.plot(z, tb, color=CRITICAL, marker=".", ms=3, lw=1.2,
                    label="probe targets read back on φ(Γ_smooth)")
            ax.set(xlabel="z", title=level)
            ax.legend(fontsize=8.5)
        axes[0][0].set_ylabel("Γ / Kutta target")
        fig.suptitle("A2/GA2.2 — does the probe estimator regenerate the "
                     "jitter from a smooth field?", fontweight="bold")
        finish(fig, OUT, "a2_intervention.png")
    return cl.report(OUT, fname="checks_interventions.csv")


if __name__ == "__main__":
    sys.exit(main())
