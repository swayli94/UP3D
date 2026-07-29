"""Development-time bit-identity A/B (phase two GS0.3, decision D1).

Bit identity is a strong and useful check -- but only where it is meaningful:
SAME machine, SAME environment, before vs after an edit. As a permanent test it
is a load- and version-sensitive false alarm (audit 2026-07-28 §6.2: 5 of 10
suite failures were 1-ULP differences). So it lives here as a tool instead.

Workflow:

    python bench/bitcheck.py --save bench/results/bit_before.npz   # before editing
    # ... make the change ...
    python bench/bitcheck.py --save bench/results/bit_after.npz
    python bench/bitcheck.py --diff bench/results/bit_before.npz \
                                    bench/results/bit_after.npz

The diff reports, per probe: bitwise equal or not, the number of differing
entries, the max absolute/relative difference, and the max difference in ULPs
of the largest entry -- so a round's report can say "subsonic path bit-identical,
supersonic path moved by design" with numbers instead of adjectives.

Probes (cheap, ~30 s total, coarse meshes only):

    res_laplace     P1 residual assembly, rho == 1
    mat_laplace     stiffness matrix data, rho == 1
    rho_tilde_sub   artificial density at a SUBCRITICAL state (must stay
                    bit-identical through any S1 change: nu == 0 there)
    rho_tilde_sup   artificial density at a SUPERCRITICAL state (this is what
                    S1 changes on purpose)
    phi_laplace     lifting Laplace solve (phi, Gamma)
    phi_transonic   coarse NACA0012 M0.80 transonic Newton (phi, Gamma)
    phi_wing        coarse ONERA M6 M0.8395 transonic Newton (phi, Gamma)
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from pyfp3d.kernels.jacobian import PicardOperator            # noqa: E402
from pyfp3d.kernels.upwind import UpwindOperator              # noqa: E402
from pyfp3d.mesh.reader import read_mesh                      # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                     # noqa: E402
from pyfp3d.physics.isentropic import density_field           # noqa: E402
from pyfp3d.solve.newton import solve_newton_transonic        # noqa: E402
from pyfp3d.solve.picard import solve_laplace_lifting         # noqa: E402

NACA = REPO / "cases/meshes/naca0012_2.5d/coarse.msh"
M6 = REPO / "cases/meshes/onera_m6/coarse.msh"


def collect():
    """Run every probe, return {name: float array}."""
    out = {}
    mc, wc = cut_wake(read_mesh(NACA))
    op = PicardOperator(mc.nodes, mc.elements)
    upw = UpwindOperator(mc.nodes, mc.elements, weighted=False)

    # a smooth, non-degenerate potential; scaled so q^2 is subcritical
    x = mc.nodes[:, 0]
    phi = x + 0.05 * x ** 2
    rho1 = np.ones(len(mc.elements))
    out["res_laplace"] = op.assemble_residual(phi, rho1).copy()
    out["mat_laplace"] = op.assemble_matrix(rho1).tocsr().data.copy()

    grad, q2 = op.velocities(phi)
    grad, q2 = grad.copy(), q2.copy()
    for tag, m_inf in (("sub", 0.50), ("sup", 0.85)):
        rho = density_field(q2, m_inf)
        out[f"rho_tilde_{tag}"] = upw.rho_tilde(
            grad, q2, rho, m_inf, 1.5, 0.95).copy()

    r = solve_laplace_lifting(mc, wc, alpha_deg=2.0)
    out["phi_laplace"] = np.asarray(r["phi"]).copy()
    out["gamma_laplace"] = np.atleast_1d(np.asarray(r["gamma"])).copy()

    r = solve_newton_transonic(
        mc, wc, m_inf=0.80, alpha_deg=1.25, m_start=0.70, dm=0.025,
        dm_min=0.003, freeze_tol=1e-6,
        newton_kw=dict(freeze_refresh_max=8, precond="direct",
                       n_newton_max=60))
    out["phi_transonic"] = np.asarray(r["phi"]).copy()
    out["gamma_transonic"] = np.atleast_1d(np.asarray(r["gamma"])).copy()

    if M6.exists():
        mcw, wcw = cut_wake(read_mesh(M6))
        r = solve_newton_transonic(
            mcw, wcw, m_inf=0.8395, alpha_deg=3.06, dm=0.05, dm_min=0.01,
            freeze_tol=1e-6, intermediate_tol=1e-5,
            newton_kw=dict(freeze_refresh_max=8, precond="amg",
                           n_newton_max=60, farfield_spanwise_gamma=True))
        out["phi_wing"] = np.asarray(r["phi"]).copy()
        out["gamma_wing"] = np.atleast_1d(np.asarray(r["gamma"])).copy()
    return out


def diff(before_npz, after_npz):
    a = np.load(before_npz)
    b = np.load(after_npz)
    # ★ GS1b.3: the THREAD COUNT is part of the environment, not a detail.
    # Measured 2026-07-29: comparing an 8-thread run against a 16-thread
    # reference reports `gamma_wing` and `phi_wing` as moved (79 / 6470 values,
    # 16 / 2 ULPs) with NO code change at all -- those probes run parallel
    # kernels whose summation order follows the thread count. Attributing that
    # to a code change is exactly the false attribution this tool exists to
    # prevent, so a mismatch is now reported before the table.
    ta = str(a["_n_threads"]) if "_n_threads" in a.files else "unrecorded"
    tb = str(b["_n_threads"]) if "_n_threads" in b.files else "unrecorded"
    if ta != tb or "unrecorded" in (ta, tb):
        print(f"!! thread counts: before = {ta}, after = {tb}. Bit-identity is "
              f"only meaningful at the SAME count -- 2-16 ULP differences in the "
              f"parallel-kernel probes (phi_wing / gamma_wing) are expected "
              f"across counts and are NOT a code change.\n")
    keys = sorted(k for k in set(a.files) | set(b.files)
                  if not k.startswith("_"))
    print(f"{'probe':18s} {'bitwise':>8s} {'n_diff':>9s} "
          f"{'max|d|':>11s} {'max rel':>10s} {'ULPs':>7s}")
    n_moved = 0
    for k in keys:
        if k not in a.files or k not in b.files:
            print(f"{k:18s} {'MISSING':>8s}")
            n_moved += 1
            continue
        x, y = np.asarray(a[k], float).ravel(), np.asarray(b[k], float).ravel()
        if x.shape != y.shape:
            print(f"{k:18s} {'SHAPE':>8s}  {x.shape} vs {y.shape}")
            n_moved += 1
            continue
        same = bool(np.array_equal(x, y))
        d = np.abs(x - y)
        nd = int(np.count_nonzero(d))
        scale = np.maximum(np.abs(x), np.abs(y))
        with np.errstate(invalid="ignore", divide="ignore"):
            rel = np.nanmax(d / np.where(scale > 0, scale, 1.0)) if d.size else 0.0
        big = float(np.max(np.abs(x))) if x.size else 1.0
        ulp = d.max() / np.spacing(big) if d.size and big > 0 else 0.0
        print(f"{k:18s} {'yes' if same else 'NO':>8s} {nd:9d} "
              f"{d.max() if d.size else 0.0:11.3e} {rel:10.3e} {ulp:7.1f}")
        if not same:
            n_moved += 1
    print(f"\n{n_moved}/{len(keys)} probe(s) moved")
    return n_moved


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--save", metavar="NPZ")
    ap.add_argument("--diff", nargs=2, metavar=("BEFORE", "AFTER"))
    args = ap.parse_args()
    if args.diff:
        sys.exit(1 if diff(*args.diff) else 0)
    if not args.save:
        ap.error("pass --save NPZ or --diff BEFORE AFTER")
    out = Path(args.save)
    out.parent.mkdir(parents=True, exist_ok=True)
    data = collect()
    data["_n_threads"] = np.array(
        [int(os.environ.get("NUMBA_NUM_THREADS", 0)),
         int(os.environ.get("OMP_NUM_THREADS", 0))])
    np.savez(out, **data)
    print(f"wrote {out} ({len(data)} probes)")
    for k, v in sorted(data.items()):
        print(f"  {k:18s} shape {np.shape(v)}")


if __name__ == "__main__":
    main()
