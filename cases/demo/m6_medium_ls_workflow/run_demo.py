"""
M6 MEDIUM LEVEL-SET WORKFLOW VALIDATION — methods × meshes × regimes.

The capability this demo evidences (unlocked by B11/B12/B13's linear-algebra
work): the ONERA M6 MEDIUM level-set solve, both SUBSONIC (M∞0.5) and TRANSONIC
(M∞0.84), on the WAKE-FREE workflow mesh — at a wall-clock now comparable to the
conforming path (B13: lifting 447→68 s, seed+Newton ~330→112 s).

It cross-checks three axes and visualizes physical reasonableness:

  AXIS 1 — MESH (level-set path): wake-free (M4, the workflow target; the wake
    is built analytically from the TE polyline, nothing embedded) vs
    wake-embedded (M1, the dual-mesh A/B partner). B7 established they agree
    ~2% at coarse; this checks it at medium.
  AXIS 2 — METHOD: level-set vs conforming (body-fitted wake). Conforming
    subsonic is solved here; conforming transonic is the committed P5 Picard
    field (`cases/demo/p5_onera_m6/results/medium_solution.npz`).
  AXIS 3 — REGIME: subsonic M0.5 vs transonic M0.84.

Figures (committed evidence):
  spanwise_loading.png — Γ(z) root→tip (all methods/meshes, both regimes); the
    tip must → 0 (the spanwise clip / free-edge rule, discretely).
  section_cp.png       — Cp(x/c) at η = 0.44 / 0.65 / 0.90, LS vs conforming;
    transonic panels carry the Cp* critical line and the shock.
  wake_potential.png    — perturbation velocity potential φ' on a spanwise slice
    (η≈0.65), showing the circulation pattern and the wake trailing from the TE.
  tip_mach.png          — local Mach on a chord-plane slab over the outer span,
    subsonic vs transonic, showing the tip acceleration (and, transonically,
    the supersonic pocket).

COST / CACHING (P5 policy): every heavy solve caches its field to a GITIGNORED
`results/*.npz`; the committed evidence is the PNGs + `summary.csv` +
`checks.csv`. The solve stage is gated — set PYFP3D_TRANSONIC_GATES=1 (or pass
--solve) to (re)compute missing caches; without it the demo plots from existing
caches and SKIPs anything missing. Cap threads at 16 INCLUDING BLAS/OMP:
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \\
    PYFP3D_TRANSONIC_GATES=1 python cases/demo/m6_medium_ls_workflow/run_demo.py

Runtime (cold, all caches missing): ~30–60 min (2 LS subsonic ~70 s, 2 LS
transonic ramps, 1 conforming subsonic; conforming transonic is loaded, not
solved). Requires the M1 + M4 medium meshes (gitignored; generate via
cases/meshes/onera_m6{,_wakefree}/generate_*.py).
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, INK, MUTED, S1_BLUE, S2_AQUA, S3_YELLOW, S4_ROSE,
    CheckList, MESH_DIR, REPO_ROOT, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_le, x_te  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve  # noqa: E402
from pyfp3d.post.shock import shock_report  # noqa: E402
from pyfp3d.post.surface import cl_kj_3d, planform_area  # noqa: E402
from pyfp3d.post.surface_ls import section_cp_curve_levelset  # noqa: E402
from pyfp3d.solve.picard import solve_subsonic_lifting  # noqa: E402
from pyfp3d.solve.picard_ls import (  # noqa: E402
    solve_multivalued_lifting, solve_multivalued_transonic,
)
from pyfp3d.wake import (  # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

OUT = Path(__file__).resolve().parent / "results"
WF_MESH = MESH_DIR / "onera_m6_wakefree" / "medium.msh"
EMB_MESH = MESH_DIR / "onera_m6" / "medium.msh"
P5_CONF = REPO_ROOT / "cases" / "demo" / "p5_onera_m6" / "results" / \
    "medium_solution.npz"

ALPHA = 3.06
M_SUB, M_TRANS = 0.5, 0.84
ETAS = (0.44, 0.65, 0.90)
LAGGED_K = 1000

GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1" or \
    "--solve" in sys.argv

# ---------------------------------------------------------------------------
# Level-set solves (cached).
# ---------------------------------------------------------------------------

def _build_ls(mesh):
    a = np.radians(ALPHA)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def ls_solution(mesh_path, regime, key):
    """Return dict {mesh, mvop, phi_ext, gamma, cl_kj, m_max, converged,
    n_limited, n_floored, wall_s, levels}. Caches phi_ext to results/<key>.npz;
    metrics are re-derived from phi_ext (cheap post) so the cache is small."""
    if not mesh_path.exists():
        return None
    mesh = read_mesh(mesh_path)
    mvop = _build_ls(mesh)
    cache = OUT / f"{key}.npz"
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        phi_ext = d["phi_ext"]
        wall_s = float(d["wall_s"]) if "wall_s" in d else 0.0
        levels = (d["levels"].tolist()
                  if "levels" in d and d["levels"].size else None)
    elif GATED:
        OUT.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        if regime == "sub":
            r = solve_multivalued_lifting(
                mvop, mesh, M_SUB, alpha_deg=ALPHA, farfield="neumann",
                n_outer_max=60, tol_residual=1e-7,
                direct_refactor_every=LAGGED_K)
            levels = None
        else:
            # tol_residual=1e-5 is ABOVE the transonic Picard residual plateau
            # (the P4/B6/N5 slow shock-position tail; diagnosed on this mesh:
            # levels park at 1e-5..4e-4 with M_max bounded and ≤3 clamped cells
            # of 329k). 1e-7 would burn every level's full budget on the
            # plateau (~1 h). The accepted state is the B7-semantics bounded /
            # engineering-converged solution, not a strict fixed point.
            r = solve_multivalued_transonic(
                mvop, mesh, M_TRANS, alpha_deg=ALPHA, farfield="neumann",
                m_start=0.60, dm=0.04, n_outer_seed=120, n_outer_level=200,
                tol_residual=1e-5, direct_refactor_every=LAGGED_K)
            levels = [(L["m_inf"], L["mach_max"], int(L["converged"]),
                       L["gamma"], L["n_limited"], L["n_floored"])
                      for L in r["levels"]]
        wall_s = time.time() - t0
        phi_ext = r["phi_ext"]
        np.savez(cache, phi_ext=phi_ext, wall_s=wall_s,
                 levels=np.array(levels, dtype=object) if levels else
                 np.array([]))
        print(f"  [{key}] solved {wall_s:.1f}s", flush=True)
    else:
        print(f"  [{key}] SKIP (no cache; set PYFP3D_TRANSONIC_GATES=1)")
        return None

    m_inf = M_SUB if regime == "sub" else M_TRANS
    gamma_st = mvop.te_jump(phi_ext)
    z = mesh.nodes[np.asarray(mvop.cm.te_nodes)][:, 2]
    o = np.argsort(z)
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    return {
        "path": "levelset", "mesh": mesh, "mvop": mvop, "phi_ext": phi_ext,
        "z": z[o], "gamma_z": gamma_st[o],
        "cl_kj": cl_kj_3d(gamma_st, z, s_ref=s_ref, b_semi=B_SEMI),
        "m_max": float(np.sqrt(np.max(mvop.element_mach2(
            phi_ext, m_inf, u_inf=1.0)))),
        "wall_s": wall_s, "levels": levels, "m_inf": m_inf, "s_ref": s_ref,
    }


# ---------------------------------------------------------------------------
# Conforming solves (subsonic here; transonic = committed P5 field).
# ---------------------------------------------------------------------------

def conforming_solution(regime, key):
    if regime == "trans":
        if not P5_CONF.exists() or not EMB_MESH.exists():
            print(f"  [{key}] SKIP (P5 conforming cache or mesh missing)")
            return None
        d = np.load(P5_CONF)
        mc, wc = cut_wake(read_mesh(EMB_MESH))
        phi, gamma, station_z = d["phi"], d["gamma"], d["station_z"]
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        return {
            "path": "conforming", "mc": mc, "phi": phi,
            "z": station_z, "gamma_z": gamma,
            "cl_kj": cl_kj_3d(gamma, station_z, s_ref=s_ref, b_semi=B_SEMI),
            "m_max": float(np.sqrt(d["mach2_max"])),
            "converged": bool(d["converged"]), "m_inf": M_TRANS,
        }
    # subsonic conforming
    if not EMB_MESH.exists():
        return None
    mc, wc = cut_wake(read_mesh(EMB_MESH))
    cache = OUT / f"{key}.npz"
    if cache.exists():
        d = np.load(cache)
        phi = d["phi"]
    elif GATED:
        OUT.mkdir(parents=True, exist_ok=True)
        r = solve_subsonic_lifting(mc, wc, M_SUB, alpha_deg=ALPHA)
        phi = r["phi"]
        np.savez(cache, phi=phi, gamma=r["gamma"], station_z=wc.station_z)
        print(f"  [{key}] solved (conforming subsonic)", flush=True)
    else:
        print(f"  [{key}] SKIP (no cache)")
        return None
    d = np.load(cache) if cache.exists() else None
    gamma = d["gamma"] if d is not None else r["gamma"]
    station_z = d["station_z"] if d is not None else wc.station_z
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    return {
        "path": "conforming", "mc": mc, "phi": phi,
        "z": station_z, "gamma_z": gamma,
        "cl_kj": cl_kj_3d(gamma, station_z, s_ref=s_ref, b_semi=B_SEMI),
        "m_inf": M_SUB,
    }


def section(sol, eta):
    """Unified section Cp(x/c) at eta for either path."""
    if sol["path"] == "levelset":
        return section_cp_curve_levelset(sol["mesh"], sol["mvop"],
                                         sol["phi_ext"], eta=eta,
                                         b_semi=B_SEMI, m_inf=sol["m_inf"])
    return section_cp_curve(sol["mc"], sol["phi"], eta=eta, b_semi=B_SEMI,
                            m_inf=sol["m_inf"])


# ---------------------------------------------------------------------------
# Figures.
# ---------------------------------------------------------------------------

STYLE = {"LS wake-free": (S1_BLUE, "-"), "LS embedded": (S2_AQUA, "--"),
         "conforming": (S4_ROSE, ":")}


def fig_spanwise_loading(sols):
    """Γ(z) root→tip for every available (method,mesh), subsonic + transonic."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, regime, title in ((axes[0], "sub", f"subsonic  M∞={M_SUB}"),
                              (axes[1], "trans", f"transonic  M∞={M_TRANS}")):
        for label, sol in sols[regime].items():
            if sol is None:
                continue
            c, ls = STYLE[label]
            ax.plot(sol["z"] / B_SEMI, sol["gamma_z"], ls, color=c, lw=1.8,
                    label=f"{label}  (cl={sol['cl_kj']:.3f})")
        ax.axhline(0, color=BASELINE, lw=0.8)
        ax.set_xlabel("z / b_semi"), ax.set_ylabel("Γ(z)")
        ax.set_title(title), ax.legend(fontsize=8)
    fig.suptitle("ONERA M6 medium: spanwise circulation Γ(z), tip → 0", y=1.02)
    finish(fig, OUT, "spanwise_loading.png")


def fig_section_cp(sols):
    """Cp(x/c) at ETAS: LS wake-free vs conforming, subsonic + transonic."""
    fig, axes = plt.subplots(2, len(ETAS), figsize=(4 * len(ETAS), 7.4),
                             sharex=True)
    for row, regime, mtag in ((0, "sub", M_SUB), (1, "trans", M_TRANS)):
        for col, eta in enumerate(ETAS):
            ax = axes[row, col]
            for label in ("LS wake-free", "LS embedded", "conforming"):
                sol = sols[regime].get(label)
                if sol is None:
                    continue
                c, ls = STYLE[label]
                cv = section(sol, eta)
                ax.plot(cv["x_upper"], cv["cp_upper"], ls, color=c, lw=1.4,
                        label=label if col == 0 else None)
                ax.plot(cv["x_lower"], cv["cp_lower"], ls, color=c, lw=1.4)
            if regime == "trans":
                rep = shock_report(section(sols[regime]["LS wake-free"], eta),
                                   mtag)
                ax.axhline(-rep["cp_critical"], color=CRITICAL, lw=0.8,
                           ls="--", label="−Cp*" if col == 0 else None)
            ax.invert_yaxis()
            ax.set_title(f"η={eta}  ({'sub' if row == 0 else 'trans'})")
            if col == 0:
                ax.set_ylabel("Cp"), ax.legend(fontsize=7)
            if row == 1:
                ax.set_xlabel("x / c")
    fig.suptitle("ONERA M6 medium: section Cp — level-set vs conforming", y=1.0)
    finish(fig, OUT, "section_cp.png")


def _slab(coords, axis, plane, half):
    return np.abs(coords[:, axis] - plane) < half


def fig_wake_potential(sols):
    """Perturbation velocity potential φ' on a spanwise slice (η≈0.65),
    level-set wake-free, subsonic + transonic. Shows the circulation pattern
    and the wake trailing downstream from the TE."""
    z0 = 0.65 * B_SEMI
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, regime, title in ((axes[0], "sub", f"subsonic M∞={M_SUB}"),
                              (axes[1], "trans", f"transonic M∞={M_TRANS}")):
        sol = sols[regime].get("LS wake-free")
        if sol is None:
            ax.set_title(f"{title} (no cache)")
            continue
        mesh, mvop = sol["mesh"], sol["mvop"]
        phi = mvop.main_potential(sol["phi_ext"])
        n = mesh.nodes
        a = np.radians(ALPHA)
        phi_inf = n[:, 0] * np.cos(a) + n[:, 1] * np.sin(a)   # u_inf=1
        # near-field window at the η=0.65 station (root->downstream a few chords)
        m = (_slab(n, 2, z0, 0.04 * B_SEMI) & (n[:, 0] > -0.6) &
             (n[:, 0] < 3.2) & (np.abs(n[:, 1]) < 1.4))
        sc = ax.tricontourf(n[m, 0], n[m, 1], (phi - phi_inf)[m], levels=28,
                            cmap="RdBu_r")
        ax.plot([x_le(z0), x_te(z0)], [0, 0], color="k", lw=2.2)  # section
        ax.plot([x_te(z0), 3.2], [0, 0], color=INK, lw=1.0, ls="--")  # wake
        ax.set_title(f"{title}  (η=0.65)")
        ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_aspect("equal")
        ax.set_xlim(-0.6, 3.2), ax.set_ylim(-1.4, 1.4)
        fig.colorbar(sc, ax=ax, shrink=0.8, label="φ'")
    fig.suptitle("ONERA M6 medium: perturbation potential near the wake "
                 "(level-set wake-free)", y=1.02)
    finish(fig, OUT, "wake_potential.png")


def fig_tip_mach(sols):
    """Local Mach on a chord-plane slab (|y| small) over the outer span,
    level-set wake-free, subsonic vs transonic — the tip acceleration."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, regime, title in ((axes[0], "sub", f"subsonic M∞={M_SUB}"),
                              (axes[1], "trans", f"transonic M∞={M_TRANS}")):
        sol = sols[regime].get("LS wake-free")
        if sol is None:
            ax.set_title(f"{title} (no cache)")
            continue
        mesh, mvop = sol["mesh"], sol["mvop"]
        cen = mesh.nodes[mesh.elements].mean(axis=1)
        mach = np.sqrt(mvop.element_mach2(sol["phi_ext"], sol["m_inf"],
                                          u_inf=1.0))
        # tip window: outer span → just past the tip, over the local chord
        m = (_slab(cen, 1, 0.0, 0.025) & (cen[:, 2] > 0.72 * B_SEMI) &
             (cen[:, 2] < 1.10 * B_SEMI) & (cen[:, 0] > 0.3) &
             (cen[:, 0] < 1.7))
        mtop = min(mach[m].max(), 3.0)
        norm = TwoSlopeNorm(vmin=0.3, vcenter=1.0, vmax=max(1.05, mtop))
        sc = ax.scatter(cen[m, 0], cen[m, 2], c=mach[m], s=10, cmap="RdBu_r",
                        norm=norm)
        ax.axhline(B_SEMI, color="k", lw=1.0, ls="--", label="tip (z=b)")
        # local LE / TE across this window (swept)
        zz = np.linspace(0.72 * B_SEMI, B_SEMI, 30)
        ax.plot([x_le(z) for z in zz], zz, color="0.35", lw=1.0)
        ax.plot([x_te(z) for z in zz], zz, color="0.35", lw=1.0)
        ax.set_title(f"{title}  (M_max={sol['m_max']:.3f})")
        ax.set_xlabel("x"), ax.set_ylabel("z (span)")
        ax.legend(fontsize=7, loc="upper right")
        fig.colorbar(sc, ax=ax, shrink=0.8, label="local Mach (red = M>1)")
    fig.suptitle("ONERA M6 medium: local Mach near the tip "
                 "(level-set wake-free, chord-plane slab)", y=1.02)
    finish(fig, OUT, "tip_mach.png")


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------

def main():
    apply_style()
    checks = CheckList("M6 medium level-set workflow — methods × meshes × "
                       "regimes")
    print("solving / loading (level-set wake-free, embedded; conforming)...")
    sols = {"sub": {}, "trans": {}}
    for regime in ("sub", "trans"):
        sols[regime]["LS wake-free"] = ls_solution(WF_MESH, regime,
                                                   f"ls_wf_{regime}")
        sols[regime]["LS embedded"] = ls_solution(EMB_MESH, regime,
                                                  f"ls_emb_{regime}")
        sols[regime]["conforming"] = conforming_solution(regime,
                                                         f"conf_{regime}")

    have_any = any(v is not None for r in sols.values() for v in r.values())
    if not have_any:
        print("SKIP: no solutions available (set PYFP3D_TRANSONIC_GATES=1 to "
              "solve). Nothing to plot.")
        return 0

    # ---- figures -------------------------------------------------------
    fig_spanwise_loading(sols)
    fig_section_cp(sols)
    fig_wake_potential(sols)
    fig_tip_mach(sols)

    # ---- summary CSV ---------------------------------------------------
    rows = []
    for regime in ("sub", "trans"):
        for label, sol in sols[regime].items():
            if sol is None:
                continue
            mm = sol.get("m_max")
            rows.append((regime, label, sol.get("path", ""),
                         f"{sol['cl_kj']:.5f}",
                         f"{mm:.4f}" if mm is not None else "n/a",
                         f"{sol['gamma_z'][0]:.5f}",
                         f"{sol['gamma_z'][-1]:.5f}",
                         f"{sol.get('wall_s', 0.0):.1f}"))
    write_csv(OUT, "summary.csv",
              "regime,label,path,cl_kj,m_max,gamma_root,gamma_tip,wall_s", rows)

    # ---- self-checks ---------------------------------------------------
    _run_checks(checks, sols)
    return checks.report(OUT, "checks.csv")


def _run_checks(checks, sols):
    # (1) tip circulation -> 0 (spanwise clip / free-edge rule), all LS sols
    for regime in ("sub", "trans"):
        for label in ("LS wake-free", "LS embedded"):
            sol = sols[regime].get(label)
            if sol is None:
                continue
            gt = abs(sol["gamma_z"][-1]) / (abs(sol["gamma_z"]).max() + 1e-30)
            checks.add("WF", f"{regime}/{label}: Γ(tip)/Γ_max → 0",
                       f"{gt:.4f}", "< 0.10", gt < 0.10)
    # (2) mesh A/B (LS wake-free vs embedded) cl agreement
    for regime in ("sub", "trans"):
        wf, emb = sols[regime].get("LS wake-free"), sols[regime].get(
            "LS embedded")
        if wf and emb:
            rel = abs(wf["cl_kj"] - emb["cl_kj"]) / abs(emb["cl_kj"])
            checks.add("MESH", f"{regime}: LS wake-free vs embedded cl",
                       f"{rel * 100:.2f}%", "< 5%", rel < 0.05)
    # (3) method A/B (LS vs conforming) cl agreement
    for regime in ("sub", "trans"):
        wf, cf = sols[regime].get("LS wake-free"), sols[regime].get(
            "conforming")
        if wf and cf:
            rel = abs(wf["cl_kj"] - cf["cl_kj"]) / abs(cf["cl_kj"])
            checks.add("METHOD", f"{regime}: LS vs conforming cl",
                       f"{rel * 100:.2f}%", "< 15%", rel < 0.15)
    # (4) transonic M_max physical/bounded (< M_cap = 3)
    for label in ("LS wake-free", "LS embedded"):
        sol = sols["trans"].get(label)
        if sol is None:
            continue
        checks.add("TRANS", f"{label}: M_max bounded (< 3)",
                   f"{sol['m_max']:.3f}", "< 3.0", sol["m_max"] < 3.0)


if __name__ == "__main__":
    sys.exit(main())
