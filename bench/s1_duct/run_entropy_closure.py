"""GS1b.3 criterion I: what does the entropy correction do to the closure's
0/0 structure?

GS1.7 measured that the probe closure is

    F = T(phi(Gamma)) - Gamma = T_0 + (b - 1) Gamma,   root Gamma* = T_0/(1 - b)

with BOTH T_0 and (1 - b) vanishing like ~h^0.6 -- the root is a ratio of two
quantities that both go to zero. That is a structural risk independent of the
shock physics, so criterion I records (no threshold) whether the entropy
correction changes it. If cl becomes h-convergent while (1 - b) keeps collapsing,
the closure is a SEPARATE disease still waiting.

b is taken from the matrices, not from a Gamma sweep (GS1.3b's route, one linear
solve instead of a sweep of field solves):

    F = T(phi) - Gamma,  K = dT/dphi_free,  dphi/dGamma = -J_ff^-1 B
    => b = dT/dGamma = -K J_ff^-1 B

and T_0 = (1 - b) * Gamma* follows from the root relation at the converged state,
so nothing else needs solving.

Isentropic references (GS1.7, M0.72): b = 0.9477 / 0.9657 / 0.9771,
1 - b = 0.0523 / 0.0343 / 0.0229 (ratio 0.66 per level = h^0.60); at
M0.7875 coarse the sweep fit gave b = 0.9742 against a matrix value of 0.933.
This run measures M0.7875 on BOTH legs with one protocol, so the ON/OFF
comparison is like for like.

Outputs: results/gs1b_3_closure.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.solve.newton import NewtonWorkspace                # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting           # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
M_TARGET = 0.7875
M_START = 0.72
DM0 = {"coarse": 0.01, "medium": 0.01, "fine": 0.02}
MAX_HALVINGS = 4
LEVELS = ("coarse", "medium", "fine")
C = 1.5


def solve_at(mc, wc, m, ent, phi=None, gam=None):
    kw = dict(m_inf=m, alpha_deg=ALPHA, upwind_c=C, m_crit=0.95,
              freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
              direct_refactor_every=4, n_newton_max=80,
              entropy_correction=ent)
    if phi is not None:
        kw.update(phi_init=phi, gamma_init=gam, n_picard_seed=0)
    return solve_newton_lifting(mc, wc, **kw)


def usable(r):
    return bool(r["converged"]) and not r.get("clamped", False)


def continue_to(mc, wc, level, ent):
    """The GS1b.2 Q5 adaptive continuation (repeated halving), verbatim."""
    dm0 = DM0[level]
    m, phi, gam, last, halv, n = M_START, None, None, None, 0, 0
    while True:
        m_next = min(m + dm0 / (2 ** halv), M_TARGET)
        r = solve_at(mc, wc, m_next, ent, phi, gam)
        n += 1
        if usable(r):
            phi, gam, m, last = r["phi"], r["gamma"], m_next, r
            if abs(m - M_TARGET) < 1e-12:
                return last, m, n
        else:
            halv += 1
            if halv > MAX_HALVINGS:
                return last, m, n


def closure_b(mc, wc, r, ent):
    """b = -K J_ff^-1 B at the converged state (probe rendering)."""
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws.set_mach(M_TARGET)
    phi_free = np.asarray(r["phi"], dtype=np.float64)[:ws.n_red][ws.free].copy()
    gamma = np.atleast_1d(np.asarray(r["gamma"], dtype=np.float64)).copy()
    if ent:
        _, _, st = ws.eval_residual(phi_free, gamma, C, 0.95, 3.0, 0.05)
        ws.refresh_sigma(st, frozen=None)
    _, F, st = ws.eval_residual(phi_free, gamma, C, 0.95, 3.0, 0.05)
    J_ff, B = ws.assemble_coupled(st, C, 0.95, 0.05)
    lu = spla.splu(J_ff.tocsc())
    dphi_dgamma = -lu.solve(np.asarray(B.todense()).ravel())
    K = ws.K                                   # dF/dphi_free, constant (probe)
    b = float(np.asarray(K @ dphi_dgamma).ravel()[0])
    return b, float(F[0]), ws.sigma_frozen


def main():
    rows = []
    for level in LEVELS:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        for ent in (False, True):
            t0 = time.perf_counter()
            r, m_reached, n_solve = continue_to(mc, wc, level, ent)
            if r is None:
                print(f"  {level} entropy={ent}: no usable state")
                continue
            b, f0, sig = closure_b(mc, wc, r, ent)
            g = float(np.atleast_1d(r["gamma"])[0])
            rows.append(dict(
                level=level, n_dof=len(mc.nodes), entropy=ent, upwind_c=C,
                m_reached=round(m_reached, 6),
                at_target=abs(m_reached - M_TARGET) < 1e-12,
                gamma_star=round(g, 8), b=round(b, 6),
                one_minus_b=round(1.0 - b, 6),
                T0=round((1.0 - b) * g, 8),
                F_at_state=f0,
                sigma_min=(round(float(sig.min()), 8) if sig is not None
                           else None),
                wall_s=round(time.perf_counter() - t0, 1), n_solves=n_solve))
            q = rows[-1]
            print(f"  {level:7s} entropy={str(ent):5s} "
                  f"M={m_reached:.5f}{'' if q['at_target'] else ' (SHORT)'} "
                  f"Gamma*={q['gamma_star']:.6f} b={q['b']:.6f} "
                  f"1-b={q['one_minus_b']:.6f} T0={q['T0']:.6f} "
                  f"({q['wall_s']}s)", flush=True)

    with open(OUT / "gs1b_3_closure.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1b_3_closure.csv")

    print("\n=== criterion I: the 0/0 closure structure, ON vs OFF ===")
    for ent in (False, True):
        sub = [r for r in rows if r["entropy"] == ent and r["at_target"]]
        if len(sub) < 2:
            print(f"  entropy={ent}: fewer than two levels at the target")
            continue
        print(f"  entropy={str(ent):5s}  "
              + "  ".join(f"{r['level']}: 1-b={r['one_minus_b']:.5f} "
                          f"T0={r['T0']:.6f}" for r in sub))
        for a, c in zip(sub[:-1], sub[1:]):
            print(f"      {a['level']} -> {c['level']}: "
                  f"(1-b) ratio {c['one_minus_b'] / a['one_minus_b']:.3f}, "
                  f"T0 ratio {c['T0'] / a['T0']:.3f}")


if __name__ == "__main__":
    main()
