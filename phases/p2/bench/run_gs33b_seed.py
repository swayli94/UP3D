"""GS3.3b: can a cheap linearised seed replace the 14.44 s Picard warm start?

Pre-registered in docs/dev_phase_two/20260802-2000-gs33b-prereg.md, committed before this
file. Criteria P1 (the ramp must converge), P2 (the RAMP's wall must beat the baseline --
not the seed line, because a worse seed can hand its saving back as Newton steps), P3
(invariance RECORDED rather than vetoed, since a seed is a legitimate algorithmic choice
and this project has measured the solution to be non-unique at 1e-4), P4 (record level-0
steps and GMRES counts, to tell "saved time" from "moved the cost").

Outputs (TRACKED): bench/gate_results/gs33b_seed.csv
"""
import csv, os, sys, time

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
sys.path.insert(0, REPO); sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,           # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve import newton as N                                # noqa: E402
from pyfp3d.solve.picard import solve_laplace_lifting               # noqa: E402
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

OUT = os.path.join(_GATE); os.makedirs(OUT, exist_ok=True)
ALPHA, M_TGT = 3.06, 0.84

mc, wc = cut_wake(read_mesh(os.path.join(REPO, "cases/meshes/onera_m6/medium.msh")))
s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
lvl = []
_orig = N.solve_newton_lifting
def _timed(*a, **k):
    t = time.perf_counter(); r = _orig(*a, **k)
    lvl.append((k.get("m_inf"), time.perf_counter() - t, r["n_newton"],
                r.get("n_gmres_total"))); return r
N.solve_newton_lifting = _timed


def leg(tag, over, pre=None):
    lvl.clear()
    t_pre = 0.0
    if pre is not None:
        t0 = time.perf_counter(); over = dict(over, **pre()); t_pre = time.perf_counter() - t0
    kw = dict(NEWTON_M6_RECIPE); kw["newton_kw"] = dict(kw["newton_kw"], **over)
    t0 = time.perf_counter()
    try:
        r = N.solve_newton_transonic(mc, wc, m_inf=M_TGT, alpha_deg=ALPHA, **kw)
    except Exception as exc:                                        # noqa: BLE001
        print(f"  {tag:12s} FAILED: {type(exc).__name__}: {exc}", flush=True)
        return dict(leg=tag, note=f"{type(exc).__name__}: {exc}")
    w = time.perf_counter() - t0 + t_pre
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"],
                                np.asarray(r["phi"]), alpha_deg=ALPHA,
                                s_ref=s_ref, m_inf=M_TGT)
    o = np.argsort(wc.station_z)
    kj = float(cl_kj_3d(np.atleast_1d(r["gamma"])[o], wc.station_z[o], s_ref, B_SEMI))
    tt = r.get("timings_total", {})
    row = dict(leg=tag, wall_s=round(w, 2), pre_s=round(t_pre, 2),
               converged=bool(r["converged"]),
               seed_s=round(tt.get("seed", 0.0), 2),
               linsolve_s=round(tt.get("linsolve", 0.0), 2),
               lvl0_steps=lvl[0][2] if lvl else None,
               lvl0_wall=round(lvl[0][1], 2) if lvl else None,
               gmres_total=sum(g or 0 for _, _, _, g in lvl),
               cl_p=f["cl"], cl_kj=kj)
    print(f"  {tag:12s} {w:6.2f}s conv={row['converged']} seed={row['seed_s']:5.2f} "
          f"pre={t_pre:5.2f} lvl0={row['lvl0_wall']}s/{row['lvl0_steps']}st "
          f"gmres={row['gmres_total']} cl_p={f['cl']:.9f}", flush=True)
    return row


def laplace_seed():
    r0 = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA)
    return dict(phi_init=np.asarray(r0["phi"]), n_picard_seed=0,
                gamma_init=np.atleast_1d(np.asarray(r0["gamma"])))


print(f"GS3.3b -- M6 medium M{M_TGT}, seed study")
rows = [leg("S-base", dict(n_picard_seed=5)),
        leg("S-none", dict(n_picard_seed=0)),
        leg("S-p2", dict(n_picard_seed=2)),
        leg("S-lap", dict(), pre=laplace_seed)]
with open(os.path.join(OUT, "gs33b_seed.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
    w.writeheader(); w.writerows(rows)
print("\nwrote", os.path.join(OUT, "gs33b_seed.csv"))

base = rows[0]
print("\n=== the registered reading (P1 -> P2 -> P3) ===")
for r in rows[1:]:
    if r.get("note"):
        print(f"  {r['leg']:8s} P1 FAIL ({r['note']})"); continue
    d = base["wall_s"] - r["wall_s"]
    p3 = max(abs(r[k] - base[k]) / abs(base[k]) for k in ("cl_p", "cl_kj"))
    verdict = ("P1 FAIL (not converged)" if not r["converged"] else
               f"saves {d:+.1f}s ({100*d/base['wall_s']:+.0f} %)" +
               ("  P2 FAIL" if d <= 0 else
                "  under the 3 s worth-it bar" if d < 3 else "  P2 OK"))
    band = ("same root (<=1e-8)" if p3 <= 1e-8 else
            "inside the known 1e-8..1e-4 non-uniqueness band -> USER" if p3 <= 1e-4
            else "OUTSIDE the known band -> leg FAILS")
    print(f"  {r['leg']:8s} {verdict};  P3 {p3:.2e} = {band}")
