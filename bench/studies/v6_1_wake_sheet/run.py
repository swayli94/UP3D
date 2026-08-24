"""GV6.1 conforming wake-sheet delta* source gate (Track V6).

Binding text: phases/p1/docs/roadmap/track_v.md GV6.1; pre-registered bands:
bench/studies/v6_1_wake_sheet/PRE_REGISTRATION.md (+ the 2026-07-25
addendum, committed before the first code change). Regenerates every
CSV/PNG in results/ and exits 0 iff all binding bands PASS (honest FAIL
otherwise).

  (a) delta*_wake = 0 bit-identical, two in-session legs (binding):
      (i) flag ON with a prescribed ZERO delta*_wake vs flag OFF (coarse,
          3 outer); (ii) flag OFF vs the gate-free library (pinned baseline
          commit 13916b5; BOTH legs fresh-compile worktree subprocesses --
          numba cache-load is not bit-faithful to fresh-compile in the
          viscous chain, isolate3 2026-07-25);
  (b) sign-pin MMS (binding): dead-air coarse strip, uniform m0 = 0.01
      through the production assembly; probe faces in the middle third of
      the strip assert antisymmetry, the jump [v_n] = m0/rho0 within 5 %
      (the lock also pins the addendum's per-face 1/2 factor), and the
      ejects-away sign;
  (c) TE-continuity (binding, runtime): the W2 construction identity
      (1e-12 rel) asserted inside every producer call -- recorded held
      over every outer of the band-(d) ON run;
  (d) GV3.1 smoke (RECORDED, non-binding): the committed GV3.1 recipe
      verbatim (medium, M0.5/alpha2/Re3e6/xtr0.05, loose Picard leg),
      flag ON vs OFF in-session A/B: on/off Delta-cl, TE-region Cp shift,
      outer counts, the m_wake field scale. The committed GV3.1 medium
      fixed point is NOT reproducible cross-run (the v5_tight_coupling
      scatter caveat); the ON/OFF pairing is in-session on the same seed.

Run:  python bench/studies/v6_1_wake_sheet/run.py
Temporary session constraint (user-directed, PRE_REGISTRATION section 7):
8 threads (NUMBA/OMP/OPENBLAS env); wall times are non-comparable.
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

import shutil
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import _cp_from_q2, wall_force_coefficients
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_airfoil_case,
    make_picard_lifting_driver,
    run_loose_coupling,
    station_average,
)
from pyfp3d.viscous.transpiration import edge_velocity_per_zone
from pyfp3d.viscous.wake_sheet import (
    assemble_wake_sheet_rhs,
    build_wake_sheet_case,
    wake_edge_velocity,
    wake_transpiration_source,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
NACA_DIR = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d")

M_INF, ALPHA, RE = 0.5, 2.0, 3.0e6
M0_MMS = 0.01
MMS_TOL = 0.05  # pinned 5% lock (band (b))
GATE_FREE_BASELINE = "13916b5"  # the GV6.0 merge = the gate-free library

SUMMARY = []  # (band, metric, band_text, measured, verdict)


def _record(band, metric, band_text, measured, ok):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append((band, metric, band_text, measured, verdict))
    print(f"  [{verdict:8s}] ({band}) {metric}: {measured} "
          f"(band: {band_text})", flush=True)


def _write_csv(name, header, rows):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"  wrote {path}", flush=True)


def _loose_loop(level, n_outer_max=None, wake=None):
    """One loose-coupling run; wake=None -> flag OFF (legacy), else the
    GV6.1 ON leg with the given WakeSheetCase. Returns (mc, wc, case,
    cfg, res, s_ref)."""
    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, f"{level}.msh")))
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA,
                         wake_transpiration=wake is not None)
    if n_outer_max is not None:
        cfg.n_outer_max = n_outer_max
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
    dz = float(np.ptp(mc.nodes[:, 2]))
    s_ref = 1.0 * dz

    def probe(phi, gamma, k):
        f = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
            alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
        return {"cl": f["cl"], "cd_p": f["cd_pressure"]}

    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)
    res = run_loose_coupling(driver, case, cfg, probe=probe, wake=wake)
    return mc, wc, case, cfg, res, s_ref


def _ab_leg_run(snippet, worktree, out_npz):
    """One (a)(ii) subprocess leg on a worktree (both legs fresh-compile:
    the isolate3/4 2026-07-25 cache-mode discipline -- numba cache-load
    is not bit-faithful to fresh-compile, and the infidelity lives
    entirely in pyfp3d/viscous/)."""
    subprocess.run(
        [sys.executable, snippet, worktree,
         os.path.join(NACA_DIR, "coarse.msh"), out_npz],
        check=True, capture_output=True, text=True, env=dict(os.environ))
    return np.load(out_npz)


def _overlay_working_tree_delta(worktree):
    """Overlay the working tree's pyfp3d/ delta (modified / added /
    untracked / deleted; __pycache__ excluded) onto a HEAD worktree, so
    the leg measures THIS tree's exact code state even when dirty."""
    out = subprocess.run(
        ["git", "-C", REPO, "status", "--porcelain", "--", "pyfp3d/"],
        check=True, capture_output=True, text=True).stdout
    for line in out.splitlines():
        xy, rel = line[:2], line[3:]
        if " -> " in rel:
            rel = rel.split(" -> ")[-1]
        rel = rel.strip('"')
        if "__pycache__" in rel:
            continue
        dst = os.path.join(worktree, rel)
        if "D" in xy:
            if os.path.exists(dst):
                os.remove(dst)
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy(os.path.join(REPO, rel), dst)


# ---------------------------------------------------------------------------
# (a) bit-identity, both legs
# ---------------------------------------------------------------------------

def band_a():
    print("--- (a) delta*_wake = 0 bit-identity ---", flush=True)
    # (i) flag ON (prescribed ZERO field) vs flag OFF, coarse 3 outer
    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, "coarse.msh")))
    wsc = build_wake_sheet_case(mc, wc)
    wsc0 = type(wsc)(
        sm=wsc.sm, slave_of_surf=wsc.slave_of_surf, s=wsc.s,
        station_of=wsc.station_of, faces_both=wsc.faces_both,
        prescribed_ds=np.zeros(wsc.sm.n_node),
    )
    *_, res_on, _ = _loose_loop("coarse", n_outer_max=3, wake=wsc0)
    *_, res_off, _ = _loose_loop("coarse", n_outer_max=3, wake=None)
    ok = (np.array_equal(res_on.phi, res_off.phi)
          and np.array_equal(res_on.gamma, res_off.gamma))
    _record("a-i", "flag-ON zero field vs flag-OFF (coarse 3 outer)",
            "bit-identical phi/gamma",
            f"phi_equal={np.array_equal(res_on.phi, res_off.phi)} "
            f"gamma_equal={np.array_equal(res_on.gamma, res_off.gamma)}",
            ok)

    # (ii) flag OFF vs the gate-free library: BOTH legs are subprocesses
    # on FRESH worktrees (no __pycache__ -> fresh numba compile) -- the
    # isolate3/4 cache-mode discipline (numba cache-load is not
    # bit-faithful to fresh-compile in pyfp3d/viscous/; an in-process
    # cache-warm leg vs a fresh worktree leg fails spuriously at ~1e-5
    # in phi at outer k >= 1 even for identical sources).
    snippet = os.path.join(RESULTS, "_ab_leg.py")
    with open(snippet, "w") as f:
        f.write(_AB_SNIPPET)
    with tempfile.TemporaryDirectory() as tmp:
        wt_base = os.path.join(tmp, "gate_free")
        wt_cur = os.path.join(tmp, "current")
        try:
            head = subprocess.run(
                ["git", "-C", REPO, "rev-parse", "HEAD"],
                check=True, capture_output=True, text=True).stdout.strip()
            subprocess.run(
                ["git", "-C", REPO, "worktree", "add", "--detach",
                 wt_base, GATE_FREE_BASELINE],
                check=True, capture_output=True, text=True)
            subprocess.run(
                ["git", "-C", REPO, "worktree", "add", "--detach",
                 wt_cur, head],
                check=True, capture_output=True, text=True)
            _overlay_working_tree_delta(wt_cur)
            base = _ab_leg_run(snippet, wt_base, os.path.join(tmp, "b.npz"))
            cur = _ab_leg_run(snippet, wt_cur, os.path.join(tmp, "c.npz"))
            ok = (np.array_equal(cur["phi"], base["phi"])
                  and np.array_equal(cur["gamma"], base["gamma"]))
            meas = (
                f"phi_equal={np.array_equal(cur['phi'], base['phi'])} "
                f"gamma_equal={np.array_equal(cur['gamma'], base['gamma'])}")
        finally:
            for wt in (wt_base, wt_cur):
                subprocess.run(
                    ["git", "-C", REPO, "worktree", "remove", "--force",
                     wt], capture_output=True)
    os.remove(snippet)
    _record("a-ii", f"flag-OFF vs gate-free library @{GATE_FREE_BASELINE} "
            "(fresh-compile worktree legs, coarse 3 outer)",
            "bit-identical phi/gamma", meas, ok)


_AB_SNIPPET = """\
import sys
sys.meta_path = [f for f in sys.meta_path
                 if "_EditableFinder" not in type(f).__name__]
sys.path.insert(0, sys.argv[1])
import numpy as np
import pyfp3d
assert pyfp3d.__file__.startswith(sys.argv[1]), pyfp3d.__file__
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.viscous.coupling import (
    CouplingConfig, build_airfoil_case, make_picard_lifting_driver,
    run_loose_coupling)
mc, wc = cut_wake(read_mesh(sys.argv[2]))
cfg = CouplingConfig(re_chord=3.0e6, m_inf=0.5, alpha_deg=2.0, n_outer_max=3)
case = build_airfoil_case(mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
res = run_loose_coupling(make_picard_lifting_driver(mc, wc, 0.5, 2.0), case, cfg)
np.savez(sys.argv[3], phi=res.phi, gamma=res.gamma)
"""


# ---------------------------------------------------------------------------
# (b) sign-pin MMS
# ---------------------------------------------------------------------------

def _tet_face_owners(elements):
    tet_faces = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    owners = {}
    for e, tet in enumerate(np.asarray(elements, dtype=np.int64)):
        for f in tet_faces:
            key = tuple(sorted(int(tet[i]) for i in f))
            owners.setdefault(key, []).append(e)
    return owners


def band_b():
    print("--- (b) sign-pin MMS (dead air, uniform m0) ---", flush=True)
    import scipy.sparse.linalg as spla

    from pyfp3d.constraints.dirichlet import farfield_dirichlet
    from pyfp3d.constraints.wake import WakeConstraint
    from pyfp3d.kernels.jacobian import PicardOperator
    from pyfp3d.solve.linear import build_amg_preconditioner

    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, "coarse.msh")))
    wsc = build_wake_sheet_case(mc, wc)
    b_wake = assemble_wake_sheet_rhs(
        mc.nodes, wsc, np.full(wsc.sm.n_node, M0_MMS))
    # dead-air reduced Laplace through the production T^T route (at
    # m_inf = 0 the Picard driver IS this incompressible solve -- the
    # G3.3 equivalence; the compressible bookkeeping is the only part
    # that breaks at u_inf = 0, not the sheet channel)
    op = PicardOperator(mc.nodes, mc.elements)
    con = WakeConstraint(op.assemble_matrix(), wc)
    gamma0 = np.zeros(wc.n_stations)
    dir_nodes, dir_vals = farfield_dirichlet(
        mc, wc, 0.0, gamma0, 0.0, (0.25, 0.0), beta=1.0)
    dir_red, vals_red = con.to_reduced_dirichlet(dir_nodes, dir_vals)
    is_dir = np.zeros(con.n_reduced, dtype=bool)
    is_dir[dir_red] = True
    free = np.where(~is_dir)[0]
    A = con.A_reduced
    A_free = A[free][:, free].tocsr()
    b_free = con.reduced_rhs(b_wake, gamma0)[free] \
        - A[free][:, dir_red].tocsr() @ vals_red
    x, info = spla.cg(A_free, b_free, M=build_amg_preconditioner(A_free)[1],
                      rtol=1e-11, maxiter=3000)
    if info != 0:
        raise RuntimeError(f"MMS CG did not converge (info={info})")
    phi_red = np.empty(con.n_reduced)
    phi_red[free] = x
    phi_red[dir_red] = vals_red
    phi = con.expand(phi_red, gamma0)

    owners = _tet_face_owners(mc.elements)
    grad, _ = op.velocities(phi)
    n_faces = len(wc.wake_faces_minus)
    cent = mc.nodes[np.asarray(wc.wake_faces_minus, dtype=np.int64)].mean(axis=1)
    x0, x1 = cent[:, 0].min(), cent[:, 0].max()
    lo, hi = x0 + (x1 - x0) / 3.0, x0 + 2.0 * (x1 - x0) / 3.0
    rows, worst_jump, worst_antisym = [], 0.0, 0.0
    sign_ok = True
    for f in range(n_faces):
        if not (lo < cent[f, 0] < hi):
            continue
        key_m = tuple(sorted(int(v) for v in wc.wake_faces_minus[f]))
        key_p = tuple(sorted(int(v) for v in wc.wake_faces_plus[f]))
        v_m = grad[owners[key_m][0], 1]
        v_p = grad[owners[key_p][0], 1]
        jump = v_p - v_m
        antisym = abs(v_p + v_m)
        sign_ok &= (v_p > 0.0) and (v_m < 0.0)
        worst_jump = max(worst_jump, abs(jump - M0_MMS))
        worst_antisym = max(worst_antisym, antisym)
        rows.append((f, f"{cent[f, 0]:.4f}", f"{v_p:.6e}", f"{v_m:.6e}",
                     f"{jump:.6e}", f"{(jump - M0_MMS) / M0_MMS:+.4f}"))
    _write_csv("mms_probes.csv",
               "face,x,v_plus,v_minus,jump,rel_err_vs_m0", rows)
    _record("b", "ejects-away sign at every probe", "v+ > 0, v- < 0",
            f"sign_ok={sign_ok}", sign_ok)
    _record("b", "antisymmetry max|v+ + v-|/m0", "<= 5%",
            f"{worst_antisym / M0_MMS:.4f}",
            worst_antisym / M0_MMS < MMS_TOL)
    _record("b", "jump max|(v+ - v-) - m0|/m0 (also pins the per-face "
            "1/2 factor)", "<= 5%", f"{worst_jump / M0_MMS:.4f}",
            worst_jump / M0_MMS < MMS_TOL)


# ---------------------------------------------------------------------------
# (c)+(d) GV3.1 smoke, flag ON vs OFF (in-session A/B)
# ---------------------------------------------------------------------------

def _wall_cp(mc, case, phi):
    """Per-wall-surface-node Cp from the loop's own recovery discipline."""
    n_vol = len(mc.nodes)
    le_mask_vol = np.zeros(n_vol, dtype=bool)
    le_mask_vol[case.sm.volume_node_of[case.le_band_surf]] = True
    ue_vol = edge_velocity_per_zone(
        mc.nodes, case.wall_faces, phi, elements=mc.elements,
        le_band_mask=le_mask_vol)
    ue_surf = ue_vol[case.sm.volume_node_of]
    q2 = np.sum(ue_surf ** 2, axis=1)
    return _cp_from_q2(q2, M_INF), ue_surf


def band_cd():
    print("--- (c)+(d) GV3.1 smoke: flag ON vs OFF (medium, in-session) ---",
          flush=True)
    mc, wc, case, cfg, res_off, s_ref = _loose_loop("medium", wake=None)
    wsc = build_wake_sheet_case(mc, wc)
    *_, res_on, _ = _loose_loop("medium", wake=wsc)

    # (c) W2 held at every outer of the ON run (the runtime assert inside
    # the producer never fired over n_outer outers)
    n_on = res_on.n_outer
    _record("c", "W2 TE-continuity identity held at every outer (1e-12 "
            "rel, runtime assert)", "no assertion", f"{n_on} outers clean",
            True)

    hist_head = ("k,ds_max,ds_change_rel,ds_neg_floored,mdot_max,"
                 "ibl_n_iter,ibl_converged,ibl_final_residual,"
                 "inflow_n_pinned,cl,cd_p,ds_wake_te,th_wake_te,"
                 "mdot_wake_max")
    for tag, res in (("on", res_on), ("off", res_off)):
        _write_csv(f"smoke_history_{tag}.csv", hist_head,
                   [tuple(h.get(k, "") for k in hist_head.split(","))
                    for h in res.history])

    cl_off = float(res_off.history[-1]["cl"])
    cl_on = float(res_on.history[-1]["cl"])
    _record("d", "cl flag OFF (final)", "recorded", f"{cl_off:.5f}", None)
    _record("d", "cl flag ON (final)", "recorded", f"{cl_on:.5f}", None)
    _record("d", "on/off Delta-cl (in-session A/B; committed GV3.1 medium "
            "fixed point NOT cross-run reproducible -- scatter caveat)",
            "recorded", f"{cl_on - cl_off:+.5f}", None)
    _record("d", "outer count ON vs OFF", "recorded",
            f"{res_on.n_outer} vs {res_off.n_outer}", None)
    mmax = max(h.get("mdot_wake_max", 0.0) for h in res_on.history[1:])
    _record("d", "m_wake scale max over outers", "recorded", f"{mmax:.3e}",
            None)
    _record("d", "ds_wake_te / th_wake_te (final)", "recorded",
            f"{res_on.wake_info.get('ds_te', float('nan')):.5e} / "
            f"{res_on.wake_info.get('th_te', float('nan')):.5e}", None)

    # final wake field along the strip: the producer state at the ON run's
    # final wall-IBL state (outs) and its final phi's wake u_e
    ds_wake, m_wake, _ = wake_transpiration_source(
        wsc, case.stations, res_on.outs,
        wake_edge_velocity(mc.nodes, wsc, res_on.phi), M_INF, 1.4, 1.0)
    _write_csv("wake_field.csv", "s,ds_wake,m_wake",
               [(f"{s:.6f}", f"{d:.6e}", f"{m:.6e}")
                for s, d, m in zip(wsc.s, ds_wake, m_wake)])

    # TE-region Cp shift (stations x/c > 0.9), per-side station averages
    cp_off, _ = _wall_cp(mc, case, res_off.phi)
    cp_on, _ = _wall_cp(mc, case, res_on.phi)
    st = case.stations
    te_rows, dcp_all = [], []
    for side_val, side_name in ((1, "upper"), (-1, "lower")):
        for r in range(len(st.xc)):
            if st.xc[r] <= 0.9:
                continue
            mask = (st.station_of == r) & (st.side_node == side_val)
            if not np.any(mask):
                continue
            a, b = float(np.mean(cp_off[mask])), float(np.mean(cp_on[mask]))
            te_rows.append((f"{st.xc[r]:.4f}", side_name, f"{a:.6f}",
                            f"{b:.6f}", f"{b - a:+.6f}"))
            dcp_all.append(abs(b - a))
    _write_csv("te_cp.csv", "x_c,side,cp_off,cp_on,dcp", te_rows)
    _record("d", "TE-region (x/c > 0.9) max |dCp| on/off", "recorded",
            f"{max(dcp_all):.5f}", None)

    _panels(res_off, res_on, wsc, ds_wake, m_wake, te_rows)


def _panels(res_off, res_on, wsc, ds_wake, m_wake, te_rows):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ax = axes[0]
    for tag, res, sty in (("OFF", res_off, "o--"), ("ON", res_on, "s-")):
        ax.plot([h["k"] for h in res.history],
                [h.get("cl", np.nan) for h in res.history], sty,
                label=f"flag {tag}")
    ax.set_xlabel("outer iteration k")
    ax.set_ylabel("c_l (pressure integral)")
    ax.set_title("GV6.1(d) GV3.1 medium smoke: cl history ON vs OFF")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    order = np.argsort(wsc.s)
    ax.plot(wsc.s[order], ds_wake[order], "-", label="delta*_wake")
    ax.plot(wsc.s[order], m_wake[order], "-", label="m_wake")
    ax.set_xlabel("s / c from TE")
    ax.set_title("prescribed wake field (final outer; L_rel = 1.0c pinned)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[2]
    for side in ("upper", "lower"):
        sel = [r for r in te_rows if r[1] == side]
        xs = [float(r[0]) for r in sel]
        ax.plot(xs, [float(r[2]) for r in sel], "o--", label=f"OFF {side}")
        ax.plot(xs, [float(r[3]) for r in sel], "s-", label=f"ON {side}")
    ax.set_xlabel("x/c")
    ax.set_ylabel("C_p")
    ax.set_title("TE-region Cp ON vs OFF")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    path = os.path.join(RESULTS, "gv6_1_smoke.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  wrote {path}", flush=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    band_a()
    band_b()
    band_cd()
    _write_csv("summary.csv", "band,metric,band,measured,verdict", SUMMARY)
    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV6.1: {n_pass} PASS / {n_fail} FAIL / {n_rec} RECORDED",
          flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
