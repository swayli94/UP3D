"""GV2.1 transpiration channel exactness (Track V2 gate).

Binding text: docs/roadmap/track_v.md GV2.1(a)-(c) (2026-07-22 re-spec);
pre-registered bands: cases/analysis/v2_transpiration_channel/
PRE_REGISTRATION.md (committed before the first execution). Regenerates
every CSV/PNG in results/ and exits 0 iff all pre-registered assertions
hold (honest FAIL otherwise).

  (a) manufactured Fourier blowing on the M0 cylinder vs the analytic
      exterior Laplace solution: relmax strict decrease + measured order
      >= 1.0 on coarse/medium/fine;
  (b) m_dot = 0 bit-identity on ALL five driver legs (Picard x3,
      conforming Newton, LS b_base);
  (c) conforming-Newton Jacobian EXACT under a lagged m_dot (structural
      bit-invariance + FD check + residual identity).

Run:  python cases/analysis/v2_transpiration_channel/run.py
"""

import os
import sys

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard import (
    solve_laplace,
    solve_subsonic,
    solve_subsonic_lifting,
)
from pyfp3d.viscous.transpiration import assemble_transpiration_rhs
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

from tests.mesh_utils import (
    cylinder_blowing_m_dot,
    cylinder_blowing_phi_exact,
    cylinder_phi_exact,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CYL_DIR = os.path.join(REPO, "cases", "meshes", "cylinder_2.5d")
NACA_DIR = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d")

V0, N_MODE = 0.1, 2
M_SUB, ALPHA = 0.3, 2.0
M_NEWTON = 0.5
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR = 1.5, 0.95, 3.0, 0.05
CASE_ARGS = dict(upwind_c=UPWIND_C, m_crit=M_CRIT, m_cap=M_CAP,
                 rho_floor=RHO_FLOOR)
M_DOT_NEWTON = 0.01   # uniform blowing for the (c) load
M_DOT_LS = 0.005

SUMMARY = []  # (gate, metric, band, measured, verdict)


def _record(gate, metric, band, measured, ok):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append((gate, metric, band, measured, verdict))
    print(f"  [{verdict:8s}] {gate} {metric}: measured={measured} "
          f"(band: {band})")


def _write_csv(name, header, rows):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"  wrote {path}")


# ---------------------------------------------------------------------------
# Gate (a) -- Fourier-mode blowing on the M0 cylinder
# ---------------------------------------------------------------------------

def gate_a():
    print("GV2.1(a) manufactured Fourier blowing on the M0 cylinder")
    rows = []
    for level in ("coarse", "medium", "fine"):
        mesh = read_mesh(os.path.join(CYL_DIR, f"{level}.msh"))
        nodes, elements = mesh.nodes, mesh.elements
        wall_faces = mesh.boundary_faces["wall"]
        wall_nodes = np.unique(wall_faces)
        ff_nodes = np.unique(mesh.boundary_faces["farfield"])
        phi_ex = cylinder_blowing_phi_exact(nodes, v0=V0, n_mode=N_MODE)
        rhs = assemble_transpiration_rhs(
            nodes, wall_faces,
            cylinder_blowing_m_dot(nodes, v0=V0, n_mode=N_MODE))
        res = solve_laplace(nodes, elements, ff_nodes, phi_ex[ff_nodes],
                            body_source_rhs=rhs, rtol=1e-11, maxiter=3000)
        err = np.abs(res["phi"] - phi_ex)
        scale = np.abs(phi_ex).max()
        relmax = err.max() / scale
        relmax_w = np.abs(res["phi"][wall_nodes] - phi_ex[wall_nodes]).max() / scale
        rel_l2 = float(np.sqrt(np.mean(err ** 2))) / scale
        rows.append((level, len(nodes), len(wall_nodes), relmax, relmax_w,
                     rel_l2, res["n_cg_iterations"], res["residual_norm"]))
        _record("GV2.1(a)", f"relmax {level}", "recorded",
                f"{relmax:.6e}", None)
        _record("GV2.1(a)", f"relmax_wall {level}", "recorded",
                f"{relmax_w:.6e}", None)
        _record("GV2.1(a)", f"relL2 {level}", "recorded",
                f"{rel_l2:.6e}", None)
        _record("GV2.1(a)", f"n_cg {level}", "recorded",
                res["n_cg_iterations"], None)
        _record("GV2.1(a)", f"residual {level}", "recorded",
                f"{res['residual_norm']:.3e}", None)

    e = [r[3] for r in rows]
    strict = e[0] > e[1] > e[2]
    _record("GV2.1(a)", "relmax strict decrease c>m>f", "binding",
            f"{e[0]:.4e} > {e[1]:.4e} > {e[2]:.4e}", strict)
    p1 = float(np.log2(e[0] / e[1]))
    p2 = float(np.log2(e[1] / e[2]))
    _record("GV2.1(a)", "order coarse->medium", ">= 1.0", f"{p1:.3f}",
            p1 >= 1.0)
    _record("GV2.1(a)", "order medium->fine", ">= 1.0", f"{p2:.3f}",
            p2 >= 1.0)
    _write_csv(
        "gv2_1a_refinement.csv",
        "level,n_nodes,n_wall,relmax,relmax_wall,relL2,n_cg,residual",
        rows)
    return rows, (p1, p2)


# ---------------------------------------------------------------------------
# Gate (b) -- m_dot = 0 bit-identity on all drivers
# ---------------------------------------------------------------------------

def _leg(gate, name, run_a, run_b, keys):
    ra, rb = run_a(), run_b()
    ok = True
    for k in keys:
        va, vb = ra[k], rb[k]
        same = (np.array_equal(va, vb) if isinstance(va, np.ndarray)
                else va == vb)
        ok &= bool(same)
        _record(gate, f"{name} {k} identical", "np.array_equal",
                same, same)
    return ok


def gate_b():
    print("GV2.1(b) m_dot = 0 bit-identity on all drivers")

    # leg 1: solve_laplace (cylinder coarse, base flow)
    mesh = read_mesh(os.path.join(CYL_DIR, "coarse.msh"))
    nodes, elements = mesh.nodes, mesh.elements
    ff = np.unique(mesh.boundary_faces["farfield"])
    phi_ff = cylinder_phi_exact(nodes)[ff]
    _leg("GV2.1(b)", "solve_laplace",
         lambda: solve_laplace(nodes, elements, ff, phi_ff,
                               rtol=1e-11, maxiter=3000),
         lambda: solve_laplace(nodes, elements, ff, phi_ff,
                               body_source_rhs=np.zeros(len(nodes)),
                               rtol=1e-11, maxiter=3000),
         ["phi", "residual_norm"])

    # leg 2: solve_subsonic (cylinder coarse, M0.3)
    phi_ff2 = nodes[ff, 0].copy()
    _leg("GV2.1(b)", "solve_subsonic",
         lambda: solve_subsonic(nodes, elements, ff, phi_ff2, m_inf=M_SUB),
         lambda: solve_subsonic(nodes, elements, ff, phi_ff2, m_inf=M_SUB,
                                body_source_rhs=np.zeros(len(nodes))),
         ["phi", "residual_history"])

    # leg 3: solve_subsonic_lifting (NACA coarse cut, M0.3)
    naca = read_mesh(os.path.join(NACA_DIR, "coarse.msh"))
    mc, wc = cut_wake(naca)
    n_cut = len(mc.nodes)
    _leg("GV2.1(b)", "solve_subsonic_lifting",
         lambda: solve_subsonic_lifting(mc, wc, m_inf=M_SUB, alpha_deg=ALPHA),
         lambda: solve_subsonic_lifting(mc, wc, m_inf=M_SUB, alpha_deg=ALPHA,
                                        body_source_rhs=np.zeros(n_cut)),
         ["phi", "gamma", "residual_history"])

    # leg 4: conforming Newton (NACA coarse, M0.5, external_rhs=zeros)
    _leg("GV2.1(b)", "newton_lifting",
         lambda: solve_newton_lifting(mc, wc, m_inf=M_NEWTON,
                                      alpha_deg=ALPHA, **CASE_ARGS),
         lambda: solve_newton_lifting(mc, wc, m_inf=M_NEWTON,
                                      alpha_deg=ALPHA,
                                      external_rhs=np.zeros(n_cut),
                                      **CASE_ARGS),
         ["phi", "gamma", "residual_history"])

    # leg 5: LS Newton (NACA coarse, M0.3, wall_rhs=zeros)
    z = naca.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(naca.nodes, naca.elements, wls,
                       wall_nodes=np.unique(naca.boundary_faces["wall"]))
    mvop = MultivaluedOperator(naca.nodes, naca.elements, cm, levelset=wls)
    ls_args = dict(m_inf=M_SUB, alpha_deg=ALPHA, farfield="neumann",
                   n_seed=10, n_newton_max=15)
    _leg("GV2.1(b)", "ls_newton",
         lambda: solve_multivalued_newton(mvop=mvop, mesh=naca, **ls_args),
         lambda: solve_multivalued_newton(mvop=mvop, mesh=naca,
                                          wall_rhs=np.zeros(len(naca.nodes)),
                                          **ls_args),
         ["phi_ext", "residual_history"])


# ---------------------------------------------------------------------------
# Gate (c) -- conforming Newton Jacobian EXACT under lagged m_dot
# ---------------------------------------------------------------------------

def gate_c():
    print("GV2.1(c) conforming-Newton Jacobian exact under lagged m_dot")
    naca = read_mesh(os.path.join(NACA_DIR, "coarse.msh"))
    mc, wc = cut_wake(naca)
    base = solve_newton_lifting(mc, wc, m_inf=M_NEWTON, alpha_deg=ALPHA,
                                **CASE_ARGS)
    assert base["converged"]

    mdot = np.zeros(len(mc.nodes))
    mdot[np.unique(mc.boundary_faces["wall"])] = M_DOT_NEWTON
    b_ext = assemble_transpiration_rhs(mc.nodes, mc.boundary_faces["wall"],
                                       mdot)

    ws0 = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws0.set_mach(M_NEWTON)
    wsb = NewtonWorkspace(mc, wc, alpha_deg=ALPHA, external_rhs=b_ext)
    wsb.set_mach(M_NEWTON)
    phi_free = np.asarray(base["phi"])[: ws0.n_red][ws0.free].copy()
    gamma = np.asarray(base["gamma"], dtype=np.float64).copy()

    # (1) structural bit-invariance of the assembled Jacobian
    _, _, state0 = ws0.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                     M_CAP, RHO_FLOOR)
    _, _, stateb = wsb.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                     M_CAP, RHO_FLOOR)
    J0, B0 = ws0.assemble_coupled(state0, UPWIND_C, M_CRIT, RHO_FLOOR)
    Jb, Bb = wsb.assemble_coupled(stateb, UPWIND_C, M_CRIT, RHO_FLOOR)
    same_J = bool(np.array_equal(J0.data, Jb.data))
    same_B = bool(np.array_equal(B0.data, Bb.data))
    _record("GV2.1(c)", "J_ff bit-invariant under b_ext", "np.array_equal",
            same_J, same_J)
    _record("GV2.1(c)", "B bit-invariant under b_ext", "np.array_equal",
            same_B, same_B)

    # (2) residual identity R'_free = R_free - (T^T b_ext)[free]
    R0, F0, _ = ws0.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                  M_CAP, RHO_FLOOR)
    Rb, Fb, _ = wsb.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                  M_CAP, RHO_FLOOR)
    ident = float(np.max(np.abs(Rb - (R0 - (ws0.con.T.T @ b_ext)[ws0.free]))))
    _record("GV2.1(c)", "residual identity max abs err", "< 1e-12",
            f"{ident:.3e}", ident < 1e-12)
    _record("GV2.1(c)", "F invariant (np.array_equal)", "np.array_equal",
            bool(np.array_equal(F0, Fb)), bool(np.array_equal(F0, Fb)))

    # (3) FD exactness with b_ext ACTIVE (test_p8 protocol)
    rng = np.random.default_rng(11)
    eps = 1e-5
    fd_rows = []
    worst = 0.0
    for i in range(3):
        delta = rng.standard_normal(wsb.n_free)
        delta /= np.abs(delta).max()
        R_p, _, _ = wsb.eval_residual(phi_free + eps * delta, gamma,
                                      UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
        R_m, _, _ = wsb.eval_residual(phi_free - eps * delta, gamma,
                                      UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
        fd = (R_p - R_m) / (2.0 * eps)
        exact = Jb @ delta
        scale = max(np.abs(fd).max(), 1e-30)
        rel = float(np.abs(exact - fd).max() / scale)
        worst = max(worst, rel)
        fd_rows.append((f"J_ff dir {i}", rel))
        _record("GV2.1(c)", f"FD J_ff dir {i}", "< 1e-5", f"{rel:.3e}",
                rel < 1e-5)
    for j in range(wsb.n_st):
        dg = np.zeros(wsb.n_st)
        dg[j] = eps
        R_p, _, _ = wsb.eval_residual(phi_free, gamma + dg, UPWIND_C,
                                      M_CRIT, M_CAP, RHO_FLOOR)
        R_m, _, _ = wsb.eval_residual(phi_free, gamma - dg, UPWIND_C,
                                      M_CRIT, M_CAP, RHO_FLOOR)
        fd = (R_p - R_m) / (2.0 * eps)
        col = np.asarray(Bb[:, j].todense()).ravel()
        scale = max(np.abs(fd).max(), 1e-30)
        rel = float(np.abs(col - fd).max() / scale)
        worst = max(worst, rel)
        fd_rows.append((f"B gamma col {j}", rel))
        _record("GV2.1(c)", f"FD B gamma col {j}", "< 1e-5", f"{rel:.3e}",
                rel < 1e-5)
    _write_csv("gv2_1c_fd.csv", "check,rel_err", fd_rows)

    # RECORDED: nonzero-load consistent solve
    r = solve_newton_lifting(mc, wc, m_inf=M_NEWTON, alpha_deg=ALPHA,
                             external_rhs=b_ext, **CASE_ARGS)
    _record("GV2.1(c)", "nonzero load converged", "recorded",
            f"{r['converged']} (|R| {r['residual_history'][-1]:.2e}, "
            f"|dgamma| {np.max(np.abs(base['gamma'] - r['gamma'])):.3e})",
            None)
    return worst


# ---------------------------------------------------------------------------
# panel + summary
# ---------------------------------------------------------------------------

def _panel(a_rows, a_orders):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    lv = [r[0] for r in a_rows]
    h = [0.10, 0.05, 0.025]
    e = [r[3] for r in a_rows]
    ew = [r[4] for r in a_rows]
    el2 = [r[5] for r in a_rows]
    ax = axes[0]
    ax.loglog(h, e, "o-", label="relmax (all nodes)")
    ax.loglog(h, ew, "s--", label="relmax (wall)")
    ax.loglog(h, el2, "^:", label="relL2")
    ax.loglog(h, [e[0] * (hh / h[0]) for hh in h], "k-.", alpha=0.5,
              label="O(h) reference")
    ax.loglog(h, [e[0] * (hh / h[0]) ** 2 for hh in h], "k:", alpha=0.5,
              label="O(h^2) reference")
    ax.set_xlabel("h_wall")
    ax.set_ylabel("phi error vs analytic")
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    ax.set_title(f"GV2.1(a) cylinder Fourier blowing "
                 f"(orders {a_orders[0]:.2f} / {a_orders[1]:.2f})")
    ax.grid(True, which="both", alpha=0.3)
    for x, y, s in zip(h, e, lv):
        ax.annotate(s, (x, y), textcoords="offset points", xytext=(6, 6),
                    fontsize=8)

    ax = axes[1]
    ax.axis("off")
    lines = ["GV2.1(b) bit-identity legs:",
             "  solve_laplace / solve_subsonic",
             "  solve_subsonic_lifting",
             "  newton_lifting (external_rhs=0)",
             "  ls_newton (wall_rhs=0)",
             "",
             "GV2.1(c) Jacobian under lagged mdot:",
             "  J_ff / B bit-invariant",
             "  FD rel err < 1e-5 (see summary.csv)",
             "  residual identity < 1e-12"]
    ax.text(0.02, 0.98, "\n".join(lines), va="top", family="monospace",
            fontsize=9)
    ax.set_title("GV2.1(b)/(c) -- see summary.csv")
    fig.tight_layout()
    path = os.path.join(RESULTS, "gv2_1_panels.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}")


def main():
    a_rows, a_orders = gate_a()
    gate_b()
    gate_c()
    _panel(a_rows, a_orders)
    _write_csv("summary.csv", "gate,metric,band,measured,verdict", SUMMARY)
    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV2.1: {n_pass} PASS / {n_fail} FAIL / {n_rec} RECORDED")
    if n_fail:
        print("HONEST FAIL -- see summary.csv")
        sys.exit(1)


if __name__ == "__main__":
    main()
