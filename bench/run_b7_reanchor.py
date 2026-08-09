"""Measure everything b7's gates assert, on ROUND vs FLAT, before re-specifying them.

The 2026-08-06 gated run failed three b7 tests, and the round-tip attribution is already
settled for the [M1] subsonic leg by a decisive A/B (round: not converged, residual
7.81e-03 at the 60-outer cap; flat: converged 5.86e-08 in 20 outers).

But the transonic gate asserts FIVE things and pytest stops at the first -- the clamp
counts, (0,3) and (0,2). So it is not yet known whether the OTHER four still hold on the
round mesh, and "3 floored cells" versus "M_max left its band" or "Gamma(z) broke" call
for completely different re-specs. This measures all of them before anything is changed.

★ b7's transonic gate deliberately does NOT assert `converged` -- the transonic Picard
tail parks near 1e-6 -- so its physicality requirement IS the 0/0 clamp assertion. That
is why the clamp counts cannot simply be relaxed: they are the whole guarantee.

Caches phi + the diagnostics BEFORE reporting.

Outputs (TRACKED): bench/gate_results/b7_reanchor.csv
"""
import csv, os, sys, time
os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np                                                  # noqa: E402
from pathlib import Path                                            # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
from tests.test_b7_onera_m6 import (ALPHA, M_INF, _gamma_of_z,      # noqa: E402
                                    _setup)
from pyfp3d.solve.picard_ls import (solve_multivalued_lifting,       # noqa: E402
                                    solve_multivalued_transonic)

CSV = os.path.join(HERE, "gate_results", "b7_reanchor.csv")
SC = os.environ.get("PYFP3D_SCRATCH", "/tmp/claude-1000/-home-lrz-codes-UP3D/"
                    "3c5b43c4-b62c-4a09-b4da-9b9c7128d43e/scratchpad")
#: (tag, family, mesh file). M4's family has no _flat variant, hence no control there.
LEGS = (("M1_round", "onera_m6", "coarse.msh"),
        ("M1_flat", "onera_m6", "coarse_flat.msh"),
        ("M4_round", "onera_m6_wakefree", "coarse.msh"))


def report(tag, kind, mesh, cm, mvop, r, wall, rows):
    phi = r["phi_ext"]
    np.savez_compressed(f"{SC}/b7_{tag}_{kind}.npz", phi=np.asarray(phi),
                        residual_norm=float(r["residual_norm"]), wall_s=wall)
    z, g = _gamma_of_z(mesh, cm, phi, mvop)
    m_max = float(np.sqrt(r["mach2_max"]))
    half = len(g) // 2
    row = dict(tag=tag, kind=kind,
               converged=bool(r.get("converged")),
               residual_norm=float(r["residual_norm"]),
               n_limited=int(r["n_limited"]), n_floored=int(r["n_floored"]),
               m_max=round(m_max, 5),
               # the four other things the transonic gate checks
               a_res_lt_1e4=bool(r["residual_norm"] < 1e-4),
               b_mmax_in_band=bool(1.0 < m_max < 2.5),
               c_gamma_min=round(float(g.min()), 6),
               c_gamma_min_ok=bool(g.min() > -1e-3),
               d_gamma_tip=round(float(abs(g[-1])), 6),
               d_gamma_tip_ok=bool(abs(g[-1]) < 0.02),
               e_outboard_rises=int((np.diff(g[half:]) > 1e-3).sum()),
               e_outboard_ok=bool(int((np.diff(g[half:]) > 1e-3).sum()) == 0),
               wall_s=round(wall, 1))
    rows.append(row)
    ok = [k for k in ("a_res_lt_1e4", "b_mmax_in_band", "c_gamma_min_ok",
                      "d_gamma_tip_ok", "e_outboard_ok") if row[k]]
    print(f"  {tag:10} {kind:9} conv={str(row['converged']):5} "
          f"|R|={row['residual_norm']:.3e} lim/flr={row['n_limited']}/{row['n_floored']} "
          f"M_max={row['m_max']:.4f} gmin={row['c_gamma_min']:+.5f} "
          f"gtip={row['d_gamma_tip']:.5f} rises={row['e_outboard_rises']} "
          f"| 其余四条通过 {len(ok)}/5 -> {ok} ({wall:.0f}s)", flush=True)


def main():
    rows = []
    print("b7 re-anchor: everything the gates assert, ROUND vs FLAT\n")
    for tag, fam, name in LEGS:
        p = Path(REPO) / "cases" / "meshes" / fam / name
        if not p.exists():
            print(f"  {tag}: {fam}/{name} missing"); continue
        mesh, cm, mvop = _setup(p)
        t0 = time.perf_counter()
        r = solve_multivalued_lifting(mvop, mesh, 0.5, alpha_deg=ALPHA,
                                      farfield="neumann", n_outer_max=60,
                                      tol_residual=1e-7)
        report(tag, "subsonic", mesh, cm, mvop, r, time.perf_counter() - t0, rows)
        mesh, cm, mvop = _setup(p)
        t0 = time.perf_counter()
        r = solve_multivalued_transonic(
            mvop, mesh, M_INF, alpha_deg=ALPHA, farfield="neumann",
            m_start=0.60, dm=0.04, n_outer_seed=120, n_outer_level=600,
            tol_residual=1e-7)
        report(tag, "transonic", mesh, cm, mvop, r, time.perf_counter() - t0, rows)
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
