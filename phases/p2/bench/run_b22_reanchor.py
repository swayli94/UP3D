"""Measure b22's M6 LS 3-D ramp anchors on the ROUND-tip meshes.

Why: test_b22_ls_3d_anchors.py holds ABSOLUTE anchors (gamma at rtol 1e-4, deliberately
tight -- "20x the measured run-to-run spread, four orders below a B20-sized move") that
were measured on the FLAT cap. `cases/meshes/onera_m6_wakefree/` was regenerated ROUND
on 2026-08-04, and the 2026-08-06 gated run failed both tests:

    coarse   target_reached False, message printed m_final=0.84
    medium   target_reached False, m_final=0.675

★ The medium 0.675 is not a surprise: LE-15 independently measured the round-tip LS M6
medium envelope failing at M0.6750, mode CLAMPING, with the clamp count rising
monotonically with Mach (3 -> 4 -> 5 -> 9). Same number from a different route.

★ The coarse message is self-contradicting -- "no longer reaches M0.84 (m_final=0.84)"
-- because the assertion is on `target_reached` while the message prints `m_final`, the
Mach ATTEMPTED rather than the Mach CONVERGED. Same key confusion that cost 40 minutes
on the G8.2 ceiling measurement. This script reports both, plus the per-level converged
flags, so the coarse case can actually be read.

Caches phi/gamma/histories BEFORE reporting (the discipline, violated four times).

Outputs (TRACKED): bench/gate_results/b22_reanchor.csv
"""
import csv, os, sys, time
os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np                                                  # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te                      # noqa: E402
from pyfp3d.solve.newton_ls import (B_NEWTON_M6_DEFAULTS,           # noqa: E402
                                    solve_multivalued_newton_transonic)
from pyfp3d.wake import (CutElementMap, MultivaluedOperator,        # noqa: E402
                         WakeLevelSet)
from tests.test_b22_ls_3d_anchors import ANCHORS, RAMP, M6_DIR      # noqa: E402

CSV = os.path.join(_GATE, "b22_reanchor.csv")
SC = os.environ.get("PYFP3D_SCRATCH", "/tmp/claude-1000/-home-lrz-codes-UP3D/"
                    "3c5b43c4-b62c-4a09-b4da-9b9c7128d43e/scratchpad")
ALPHA = 3.06


def main():
    rows = []
    print("b22 re-anchor: the M6 LS ramps on the ROUND-tip wakefree meshes\n")
    for level in ("coarse", "medium"):
        p = M6_DIR / f"{level}.msh"
        if not p.exists():
            print(f"  {level}: mesh missing"); continue
        mesh = read_mesh(p)
        a = np.radians(ALPHA)
        wls = WakeLevelSet(
            np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
            direction=(np.cos(a), np.sin(a), 0.0))
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]))
        mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
        t0 = time.perf_counter()
        r = solve_multivalued_newton_transonic(mvop=mvop, mesh=mesh, **RAMP,
                                               **B_NEWTON_M6_DEFAULTS)
        wall = time.perf_counter() - t0
        #: ★ FIFTH caching lesson, 2026-08-06: caching is not binary. The first run
        #: saved phi/gamma/residual/M_max and NOT `level_results` -- and the re-anchor
        #: needs the PER-LEVEL clamp counts, because the terminal states here are
        #: clamped and GS1.4 says a clamped state is not a solution, so the anchor
        #: belongs on the highest level that converged with ZERO clamps. Saving "the
        #: obvious fields" cost a 37-minute re-run.
        import json
        with open(f"{SC}/b22_levels_{level}.json", "w") as _fh:
            json.dump(r.get("levels") or [], _fh, indent=1, default=float)
        np.savez_compressed(f"{SC}/b22_reanchor_{level}.npz",
            phi=np.asarray(r["phi_ext"] if "phi_ext" in r else r.get("phi", [])),
            gamma=np.atleast_1d(np.asarray(r["gamma"])),
            residual_history=np.asarray(r.get("residual_history", []), float),
            mach2_max=float(r["mach2_max"]), wall_s=wall)
        lv = r.get("levels") or []
        ok = [l for l in lv if l.get("converged")]
        #: ★ the level key is `m_inf` on the LEVEL-SET driver and `m` on the
        #: CONFORMING one (newton_ls.py:1107 vs newton.py's level_results) -- the same
        #: quantity under two names across the two paths, which is the cross-path
        #: inconsistency discipline #9 exists for. Reading `m` here raised KeyError
        #: AFTER the coarse solve; the cache-before-report rule meant it survived.
        m_ok = max((float(l["m_inf"]) for l in ok), default=None)
        row = dict(level=level, target_reached=bool(r.get("target_reached")),
                   m_final=r.get("m_final"), m_highest_converged=m_ok,
                   gamma=float(np.atleast_1d(r["gamma"])[0]) if np.size(r["gamma"]) else None,
                   m_max=round(float(np.sqrt(r["mach2_max"])), 6),
                   res_final=float(np.asarray(r["residual_history"], float)[-1]),
                   n_limited=int(r["n_limited"]), n_floored=int(r["n_floored"]),
                   n_levels=len(lv), n_levels_converged=len(ok),
                   m_highest_clean=max((float(l["m_inf"]) for l in ok
                                        if not l.get("n_limited")
                                        and not l.get("n_floored")), default=None),
                   anchor_gamma=ANCHORS[level]["gamma"],
                   anchor_m_max=ANCHORS[level]["m_max"], wall_s=round(wall, 1))
        rows.append(row)
        print(f"  {level:7} target_reached={row['target_reached']} "
              f"m_final={row['m_final']} highest_converged={m_ok} "
              f"gamma={row['gamma']} (anchor {row['anchor_gamma']}) "
              f"M_max={row['m_max']} (anchor {row['anchor_m_max']}) "
              f"|R|={row['res_final']:.2e} lim/flr={row['n_limited']}/{row['n_floored']} "
              f"levels {len(ok)}/{len(lv)} ({wall:.0f}s)  [cached]", flush=True)
        print(f"     highest CLEAN (converged, 0 clamps) = {row['m_highest_clean']}")
        #: built outside the f-string: nesting quotes inside an f-string does not
        #: parse on 3.11, and this is the second time this phase.
        per = [(round(float(l["m_inf"]), 4), bool(l["converged"]),
                str(l.get("n_limited")) + "/" + str(l.get("n_floored")))
               for l in lv]
        print(f"     per-level (m, conv, lim/flr): {per}",
              flush=True)
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
