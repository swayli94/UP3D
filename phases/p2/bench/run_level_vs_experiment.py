"""Which mesh level is actually CLOSEST to the committed M6 experiment?

The convergence reading says refinement does not reach first order on any 3-D combination.
That does NOT license shipping a coarse mesh: "not converging" means refinement does not
approach a limit, not that every level is equally right. To rank the levels you need a
REFERENCE, and the M6 wing has one committed -- TEST 2308, the only condition in the file
(M 0.8395, alpha 3.06, Re_MAC 11.72e6).

Two design points that decide whether this measurement means anything:

★ The M_max < 1.4 screen is DROPPED here, deliberately. It was introduced to screen
  CONVERGENCE safety, not accuracy. Under it, medium cannot reach the experimental Mach at
  all (conforming stops at 0.65-0.683, level-set at 0.60) while xcoarse and coarse both
  reach 0.84 -- and they reach it precisely BECAUSE a coarse mesh under-resolves the peak
  and so reports a lower M_max. Keeping the screen would therefore rig the accuracy
  comparison in favour of the coarse meshes. Any converged solution is admitted; M_max is
  recorded alongside, not used as a filter.

★ The reference is VISCOUS (Re 11.72e6) and this solver is inviscid, so there is a model
  floor no level can beat -- GV5.3 measured pooled RMS 0.1288 at medium k=0, and S2
  measured an inviscid model floor of 0.052-0.071. So the absolute RMS is NOT the result.
  The result is the RANKING and its DIRECTION: does the error against experiment fall as the
  mesh refines? That question is answerable even though the floor is nonzero, and it needs
  no convergence model at all -- which is the point, since p proved unreliable and R only
  ever speaks about a limit, never about correctness.

Reuses GV5.3's committed machinery by import (parse_experiment, ETAS, the station stack)
rather than re-deriving it, so these numbers sit on the same experiment parsing and the same
station set as the committed evidence.

Outputs (TRACKED): bench/gate_results/level_vs_experiment.csv
"""

import csv
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "bench", "studies", "v5_3_m6_cp"))

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.post.unified import section_cp                          # noqa: E402

#: GV5.3's committed experiment parsing and station set -- imported, not re-derived.
from run import ETAS, N_UNMASKED, parse_experiment                  # noqa: E402

CSV = os.path.join(_GATE, "level_vs_experiment.csv")
M_INF, ALPHA = 0.8395, 3.06
B_SEMI = cap.B_SEMI

#: every M6 wing level available on each path. `coarse` is the SHIPPED (h_far-clamped)
#: level and `coarse_ss` its self-similar twin -- both kept, because they are different
#: meshes and the whole question is which mesh is closest.
CASES = [
    ("conforming", "onera_m6", "xcoarse_ss", cap.conf_wing),
    ("conforming", "onera_m6", "coarse_ss", cap.conf_wing),
    ("conforming", "onera_m6", "coarse", cap.conf_wing),
    ("conforming", "onera_m6", "medium", cap.conf_wing),
    ("level-set", "onera_m6_wakefree", "xcoarse_ss", cap.ls_wing),
    ("level-set", "onera_m6_wakefree", "coarse_ss", cap.ls_wing),
    ("level-set", "onera_m6_wakefree", "coarse", cap.ls_wing),
    ("level-set", "onera_m6_wakefree", "medium", cap.ls_wing),
]


def station_rms(mesh, kw, m_eff, eta, exp):
    """RMS of (computed - experimental) Cp at one station, sampled at the
    experiment's own x/c points. Same construction as GV5.3's station_rms."""
    sec = section_cp(mesh, eta=eta, b_semi=B_SEMI, m_inf=m_eff, **kw)
    e = exp[eta]
    tot, n = 0.0, 0
    for want_upper in (True, False):
        side = "upper" if want_upper else "lower"
        m = e["upper"] == want_upper
        if not np.any(m):
            continue
        cp_i = np.interp(e["x"][m], sec[f"x_{side}"], sec[f"cp_{side}"])
        tot += float(np.sum((cp_i - e["cp"][m]) ** 2))
        n += int(m.sum())
    return (tot / max(n, 1)) ** 0.5, n


def main():
    exp = parse_experiment()
    rows = []
    for path, mdir, level, fn in CASES:
        mesh_path = os.path.join(REPO, "cases", "meshes", mdir, f"{level}.msh")
        if not os.path.exists(mesh_path):
            print(f"  {path}/{level}: mesh missing", flush=True)
            continue
        print(f"\n=== {path} / {level} @ M{M_INF} alpha {ALPHA} ===", flush=True)
        import time
        t0 = time.perf_counter()
        try:
            mesh, op, r, phi, mvop = fn(mesh_path, M_INF, ALPHA)
        except Exception as exc:                                   # noqa: BLE001
            print(f"  ERROR {type(exc).__name__}: {exc}", flush=True)
            rows.append(dict(path=path, level=level, status="ERROR",
                             note=f"{type(exc).__name__}: {exc}"))
            continue
        wall = time.perf_counter() - t0
        if mvop is None:
            kw, m_eff = dict(phi=phi), M_INF
            conv = bool(r["converged"])
            res = float(r["residual_history"][-1])
            m_att = float(r.get("m_last_converged", r.get("m_final", M_INF)))
        else:
            mf = r.get("m_final", M_INF)
            kw, m_eff = dict(mvop=mvop, phi_ext=phi), mf
            conv = bool(r.get("target_reached", False))
            res = float(r["levels"][-1]["residual_norm"])
            m_att = float(r.get("m_last_converged", mf))
        per, tot, ntot = {}, 0.0, 0
        for eta in ETAS[:N_UNMASKED]:          # tip-masked stations excluded, as in GV5.3
            try:
                rms, n = station_rms(mesh, kw, m_eff, eta, exp)
            except Exception as exc:                               # noqa: BLE001
                print(f"    eta {eta}: {type(exc).__name__}: {exc}", flush=True)
                continue
            per[eta] = rms
            tot += rms ** 2 * n
            ntot += n
        pooled = (tot / ntot) ** 0.5 if ntot else float("nan")
        print(f"  conv={conv} m_attained={m_att} |R|={res:.3e} "
              f"pooled_RMS={pooled:.4f} ({wall:.0f}s)", flush=True)
        print("    " + "  ".join(f"eta{k:.2f}:{v:.4f}" for k, v in per.items()),
              flush=True)
        rows.append(dict(path=path, level=level, n_nodes=len(mesh.nodes),
                         status="OK" if conv else "NOT_CONVERGED",
                         converged=conv, m_attained=m_att, res_final=res,
                         pooled_rms=round(pooled, 6), wall_s=round(wall, 1),
                         **{f"rms_eta{k:.2f}": round(v, 6) for k, v in per.items()},
                         note=""))
    if rows:
        keys = sorted({k for r in rows for k in r})
        with open(CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        print(f"\nwrote {CSV}")
    print("\n=== DOES THE ERROR AGAINST EXPERIMENT FALL WITH REFINEMENT? ===")
    for path in ("conforming", "level-set"):
        sel = [r for r in rows if r.get("path") == path and r.get("pooled_rms")]
        if not sel:
            continue
        print(f"  {path}:")
        for r in sorted(sel, key=lambda r: r["n_nodes"]):
            print(f"    {r['level']:11s} {r['n_nodes']:>7d} nodes  "
                  f"pooled RMS {r['pooled_rms']:.4f}  "
                  f"{'' if r['converged'] else '(NOT CONVERGED)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
