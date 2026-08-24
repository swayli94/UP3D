"""
COST PROBE ONLY -- xcoarse / coarse wing-body legs at M0.88.

★★★ THIS ROUND HAS NO CRITERIA AND PRODUCES NO VERDICT. It exists because I
asserted that "raise the sample independence" was expensive WITHOUT PRICING IT,
which is the exact error this project logged in phase 3 round 21 ("I deferred an
almost-free obligation four rounds on a cost I had made up") and made a standing
rule in the GS4.1 document §7.6: before restating a debt a second time, spend the
zero compute needed to pin down what it costs.

So: measure the wall clock, the convergence flags and the clamp counts of the two
cheaper mesh levels, and CACHE phi so the real round pays nothing. Nothing here
may be read as a result about the recovery operator -- that is what the next
pre-registered round is for.

★ The recipe is COPIED FROM `failure_modes.py` verbatim (same seed call,
same CONF_SEED_KW / CONF_RAMP_NK, same taper, same dm / dm_min / freeze_tol /
intermediate_tol), because the whole point is that the new levels sit on the SAME
recipe as the medium states already on disk. A drifted recipe would make the
levels incomparable -- the 5th question, cross-recipe.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "bench"))

import capability_matrix as cap                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)

MDIR = REPO / "cases/meshes/onera_m6_wingbody_conforming"
OUTD = REPO / "bench/gate_results/gs40d_levels"
CSV = REPO / "bench/gate_results/gs40d_level_cost.csv"
M_TARGET, ALPHA, R_C = 0.88, 3.06, 0.05        # le14's production leg
LEVELS = ("xcoarse", "coarse")


def leg(level):
    msh = MDIR / f"{level}.msh"
    if not msh.exists():
        return dict(level=level, note="mesh absent -- regenerate", converged=None)
    man = read_manifest(msh)
    sha = mesh_fingerprint(msh)["sha256"]
    assert man is None or sha == man["sha256"], f"{level}: mesh moved from manifest"

    t0 = time.perf_counter()
    mc, wc = cut_wake(read_mesh(msh))
    t_mesh = time.perf_counter() - t0

    t = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", R_C * B_SEMI)
    t1 = time.perf_counter()
    try:
        seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART, alpha_deg=ALPHA,
                                    **cap.CONF_SEED_KW)
        nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
                  phi_init=seed["phi"], gamma_init=seed["gamma"],
                  n_picard_seed=0, tip_taper=t)
        r = solve_newton_transonic(mc, wc, m_inf=M_TARGET, alpha_deg=ALPHA,
                                   m_start=cap.WB_MSTART, dm=cap.DM, dm_min=0.01,
                                   freeze_tol=1e-5, intermediate_tol=1e-4,
                                   newton_kw=nk)
    except Exception as exc:                                       # noqa: BLE001
        return dict(level=level, n_nodes=len(mc.nodes), n_tets=len(mc.elements),
                    wall_s=round(time.perf_counter() - t1, 1),
                    converged=False, note=f"{type(exc).__name__}: {exc}")
    wall = time.perf_counter() - t1

    #: ★ cache FIRST, report after -- "cache before you report" is now code order,
    #: not advice: two reporting-layer defects this season each destroyed a
    #: finished solve. Everything the next round could want goes in.
    OUTD.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUTD / f"{level}.npz", phi=r["phi"], gamma=r["gamma"],
                        conv=bool(r["converged"]),
                        nlim=int(r["n_limited"]), nflr=int(r["n_floored"]),
                        res=float(r["residual_history"][-1]),
                        hist=np.asarray(r["residual_history"], dtype=float),
                        m_final=float(r["m_final"] if r["m_final"] is not None else np.nan),
                        mesh_sha=sha)

    #: ★ the GS4.0 honesty fields, used here for the first time by a new consumer
    return dict(level=level, n_nodes=len(mc.nodes), n_tets=len(mc.elements),
                mesh_sha=sha[:12], t_mesh_s=round(t_mesh, 1),
                wall_s=round(wall, 1), converged=bool(r["converged"]),
                n_limited=int(r["n_limited"]), n_floored=int(r["n_floored"]),
                res_final=float(r["residual_history"][-1]),
                m_final=r["m_final"], m_last_converged=r["m_last_converged"],
                target_reached=r["target_reached"],
                n_levels=len(r["level_results"]), note="")


def main():
    print("★ COST PROBE -- no criteria, no verdict. Wall clock + convergence only.")
    print(f"  recipe copied verbatim from failure_modes (M{M_TARGET} / "
          f"alpha {ALPHA} / r_c {R_C}); medium is already on disk.\n")
    rows = []
    for lv in LEVELS:
        print(f"[{lv}] running ...", flush=True)
        row = leg(lv)
        rows.append(row)
        print("   " + "  ".join(f"{k}={v}" for k, v in row.items()
                                if k not in ("note",)) + f"  {row.get('note','')}",
              flush=True)
    CSV.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV.name}; phi cached under {OUTD.name}/ so the next round "
          "pays nothing.")
    print("★ Nothing here is a result about the recovery operator.")


if __name__ == "__main__":
    main()
