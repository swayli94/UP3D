"""GS1b.6: the phase-one ON/OFF table -- batch 1 (gates and benchmarks).

Protocol (pre-registered in docs/dev_phase_two/20260730-0000-s1b-onoff-recompute
.md): every in-scope quantity is computed BOTH ways on the same mesh, the same
recipe and the same thread count, and written with the phase-one COMMITTED value
beside it. Criterion R1 is the load-bearing one: the OFF leg must reproduce the
committed number, otherwise something else drifted and the whole table is
untrustworthy until that is found.

Scope is provably bounded (round file sec 2): no shock => ON and OFF are
bit-identical (locked by a test), and the level-set path is unwired (grep = 0), so
only the CONFORMING path with a supersonic zone can move.

Outputs: bench/gate_results/entropy_ab.csv  (TRACKED -- bench/results/ is
gitignored, which is how the GS1.5 close-out ended up claiming a committed artifact
that was never committed)
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = HERE.parent
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                          # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                         # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve                 # noqa: E402
from pyfp3d.post.shock import shock_report                        # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients           # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting              # noqa: E402

OUT = _GATE
OUT.mkdir(exist_ok=True)

#: phase-one / phase-two committed values, with where each is written down. These
#: are what criterion R1 checks the OFF leg against.
COMMITTED = {
    # M1a lock dict in tests/test_s1_m1a_envelope.py (GS1.5a)
    ("m1a", "coarse"): 0.242797,
    ("m1a", "medium"): 0.253351,
    ("m1a", "fine"): 0.255662,
    # the M1 gate FAIL recorded in the GS1.5 close-out (coarse, C = 1.5)
    ("m1", "coarse"): 0.458974,
    # bench/baseline_2026-07-28.csv, wing/coarse
    ("wing", "coarse"): 0.25579327,
}

CASES = (
    # (group, level, m_inf, alpha, upwind_c, quantity note)
    ("m1a", "coarse", 0.72, 1.25, 1.5),
    ("m1a", "medium", 0.72, 1.25, 1.5),
    ("m1a", "fine", 0.72, 1.25, 1.5),
    ("m1", "coarse", 0.80, 1.25, 1.5),
    ("m1", "medium", 0.80, 1.25, 1.5),
)


def solve(mc, wc, m, a, C, ent, phi=None, gam=None):
    kw = dict(m_inf=m, alpha_deg=a, upwind_c=C, m_crit=0.95, freeze_tol=1e-6,
              freeze_refresh_max=8, precond="direct", direct_refactor_every=4,
              n_newton_max=80, entropy_correction=ent)
    if phi is not None:
        kw.update(phi_init=phi, gamma_init=gam, n_picard_seed=0)
    return solve_newton_lifting(mc, wc, **kw)


def usable(r):
    return bool(r["converged"]) and not r.get("clamped", False)


def continue_to(mc, wc, m_target, a, C, ent, m_start=0.60, dm0=0.04,
                max_halvings=5):
    """Adaptive continuation (the GS1b.2 Q5 machinery). Returns (result,
    m_reached, halvings, solves) -- the halvings are part of the RECORD: criterion
    R4 wants the cases that needed smaller Mach steps named."""
    m, phi, gam, halv, n, last = m_start, None, None, 0, 0, None
    while halv <= max_halvings:
        m_next = min(m + dm0 / (2 ** halv), m_target)
        r = solve(mc, wc, m_next, a, C, ent, phi, gam)
        n += 1
        if usable(r):
            phi, gam, m, last = r["phi"], r["gamma"], m_next, r
            if abs(m - m_target) < 1e-12:
                return last, m, halv, n
        else:
            halv += 1
    return last, m, halv - 1, n


def main():
    rows = []
    for group, level, m_inf, alpha, C in CASES:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {group}/{level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        dz = float(np.ptp(mc.nodes[:, 2]))
        for ent in (False, True):
            t0 = time.perf_counter()
            r, m_reached, halv, nsolve = continue_to(mc, wc, m_inf, alpha, C,
                                                     ent)
            wall = time.perf_counter() - t0
            if r is None:
                rows.append(dict(group=group, level=level, m_inf=m_inf,
                                 entropy=ent, reached=None, usable=False,
                                 note="no usable state", wall_s=round(wall, 1)))
                print(f"  {group}/{level} entropy={ent}: no usable state "
                      f"({wall:.0f}s)", flush=True)
                continue
            rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                            m_inf=m_reached), m_reached)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], r["phi"],
                                        alpha_deg=alpha, s_ref=dz,
                                        m_inf=m_reached)
            com = COMMITTED.get((group, level))
            rel = (None if (com is None or ent)
                   else round(f["cl"] / com - 1.0, 6))
            rows.append(dict(
                group=group, level=level, m_inf=m_inf, entropy=ent,
                reached=round(m_reached, 6),
                at_target=abs(m_reached - m_inf) < 1e-12, usable=usable(r),
                cl_p=round(f["cl"], 6), gamma=round(float(r["gamma"][0]), 8),
                x_shock=rep["upper"].get("x_shock"),
                m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                res_final=r["residual_history"][-1], n_newton=r["n_newton"],
                sigma_min=r.get("sigma_min"), m1_detected=r.get("m1_max"),
                n_shock_cells=r.get("n_shock_cells"),
                committed_cl=com, off_vs_committed_rel=rel,
                n_halvings=halv, n_solves=nsolve, note="",
                wall_s=round(wall, 1)))
            q = rows[-1]
            tag = "" if com is None or ent else (
                f"  vs committed {com}: {100 * rel:+.3f} %")
            print(f"  {group}/{level:7s} entropy={str(ent):5s} "
                  f"reached={q['reached']}"
                  f"{'' if q['at_target'] else ' (SHORT)'} "
                  f"cl={q['cl_p']:.6f} x_sh={q['x_shock']} "
                  f"M_max={q['m_max']} halv={halv} ({wall:.0f}s){tag}",
                  flush=True)

    with open(OUT / "entropy_ab.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "entropy_ab.csv")

    print("\n=== R1: does the OFF leg reproduce the committed numbers? ===")
    for (g, lv), com in COMMITTED.items():
        off = next((r for r in rows if r.get("group") == g
                    and r["level"] == lv and r["entropy"] is False), None)
        if off is None or off.get("cl_p") is None:
            print(f"  {g}/{lv}: no OFF value (not run or not usable)")
            continue
        rel = off["cl_p"] / com - 1.0
        print(f"  {g}/{lv:7s} OFF cl {off['cl_p']:.6f} vs committed {com}: "
              f"{100*rel:+.3f} %  {'OK' if abs(rel) < 5e-3 else 'DRIFT'}")

    print("\n=== the ON/OFF table (cl_p, x_shock) ===")
    for g, lv, m, a, C in CASES:
        pair = {r["entropy"]: r for r in rows
                if r.get("group") == g and r["level"] == lv}
        if len(pair) < 2 or any(p.get("cl_p") is None for p in pair.values()):
            print(f"  {g}/{lv}: incomplete pair -- recorded, not hidden")
            continue
        o, n = pair[False], pair[True]
        print(f"  {g}/{lv:7s} M{m}: cl {o['cl_p']:.6f} -> {n['cl_p']:.6f} "
              f"({100*(n['cl_p']/o['cl_p']-1):+.2f} %)   "
              f"x_shock {o['x_shock']:.4f} -> {n['x_shock']:.4f}   "
              f"1-sigma {100*(1-n['sigma_min']):.3f} %")


if __name__ == "__main__":
    main()
