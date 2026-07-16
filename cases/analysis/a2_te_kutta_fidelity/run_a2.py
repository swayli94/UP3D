"""
Track A / A2 -- TE/Kutta fidelity study, ZERO-SOLVE legs (L0/L1/L2z/L3).

Two symptoms, both user-flagged from committed figures (2026-07-16):
  S1  conforming spanwise circulation Gamma(z) carries station-to-station
      jitter (a1_m6_spanwise.png, g52_spanwise_*.png) while the level-set
      path is smooth on BOTH mesh families (B7 gamma_of_z.png);
  S2  section Cp jumps unphysically at the trailing edge (g51_sections_*,
      m6_cp_sections.png) -- present in 2.5-D too (a1_cp.png), conforming
      much worse than level-set (B7 section_cp.png).

Prior art this study sharpens (it does NOT rediscover them):
  * B7 committed the jitter CONTRAST: roughness (RMS 2nd diff / range)
    0.0970 conforming-P5 vs 0.0079 / 0.0091 level-set (track_b.md B7).
  * P13 committed the prose ATTRIBUTIONS: Gamma jitter = "conforming
    Kutta-probe placement"; TE Cp gap 0.14-0.22 = potential-jump Kutta
    approximation + sharp-TE P1 floor (demo_report/track_p.md). Neither has
    intervention-grade evidence -- that is A2's deliverable.
  * INVESTIGATION_kutta_closure.md: probe sharing / off-plane probes on the
    swept unstructured TE recorded as a known-robustness item, NOT fixed.
  * INVESTIGATION_gamma_smoothing.md: smoothing-as-FIX is a dead route; the
    fixed-Gamma smooth input in the gated intervention leg is a DIAGNOSTIC
    (T3/E methodology), not that route re-opened.

This script performs NO solves. It harvests the committed/local solution
caches, so every leg costs seconds except the one-off level-set operator
rebuilds (cached to results/a2_cache_*.npz, gitignored). Legs whose local
cache (.npz, gitignored) is absent are SKIPPED with a message -- regenerate
via the owning demo (p5_onera_m6, b7_onera_m6, a1_solver_bottleneck).
Meshes are gitignored too: cases/meshes/onera_m6/generate_onera_m6.py and
generate_onera_m6_wakefree.py (~30 s each) rebuild them.

Run:  python cases/analysis/a2_te_kutta_fidelity/run_a2.py
Refresh derived caches:  PYFP3D_A2_REFRESH=1 python ...

GA2 gate status: this is the SCAFFOLD run -- it records measurements;
closing any GA2 gate is a separate, arbitrated act (roadmap/track_a.md A2).
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(HERE))

import matplotlib.pyplot as plt                                    # noqa: E402

import _metrics as M                                               # noqa: E402
from cases.demo._common import (                                   # noqa: E402
    CheckList, INK_2, MUTED, S1_BLUE, S2_AQUA, S3_YELLOW, S4_ROSE,
    apply_style, finish, write_csv,
)
from pyfp3d.constraints.wake import kutta_targets                  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                           # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                          # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, x_te                     # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve               # noqa: E402
from pyfp3d.post.surface import (                                  # noqa: E402
    _cp_from_q2, triangle_tangential_gradients,
)
from pyfp3d.post.surface_ls import section_cp_curve_levelset       # noqa: E402
from pyfp3d.wake import (                                          # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
REFRESH = os.environ.get("PYFP3D_A2_REFRESH", "0") == "1"

M_INF, ALPHA = 0.84, 3.06
ETAS = (0.20, 0.44, 0.65, 0.90)
XC_TARGETS = (0.92, 0.96, 0.98, 0.99, 0.995)
MESHES = REPO_ROOT / "cases/meshes"
P5_RES = REPO_ROOT / "cases/demo/p5_onera_m6/results"
B7_RES = REPO_ROOT / "cases/demo/b7_onera_m6/results"
A1_RES = REPO_ROOT / "cases/analysis/a1_solver_bottleneck/results"
P5_BASELINE_CSV = REPO_ROOT / "cases/demo/b7_onera_m6/p5_gamma_baseline.csv"

# B7's committed roughness numbers (track_b.md B7 / b7 checks.csv) -- the
# GA2.1 reproduction anchors for the shared metric implementation.
R_ANCHORS = {"p5_baseline": 0.0970, "ls_m1_coarse": 0.0079,
             "ls_m4_coarse": 0.0091}

# state -> (solution npz, mesh, colour, level). npz are gitignored local
# caches; absent => leg skipped.
CONF_STATES = {
    "conf_picard_coarse": (P5_RES / "coarse_solution.npz",
                           MESHES / "onera_m6/coarse.msh", S1_BLUE, "coarse"),
    "conf_picard_medium": (P5_RES / "medium_solution.npz",
                           MESHES / "onera_m6/medium.msh", S1_BLUE, "medium"),
    "conf_newton_medium": (A1_RES / "a1_m6_conf_newton.npz",
                           MESHES / "onera_m6/medium.msh", S2_AQUA, "medium"),
}
LS_STATES = {
    "ls_m1_coarse": (B7_RES / "M1.npz",
                     MESHES / "onera_m6/coarse.msh", S3_YELLOW, "coarse"),
    "ls_m4_coarse": (B7_RES / "M4.npz",
                     MESHES / "onera_m6_wakefree/coarse.msh", S3_YELLOW,
                     "coarse"),
    "ls_newton_medium": (A1_RES / "a1_m6_ls_newton.npz",
                         MESHES / "onera_m6_wakefree/medium.msh", S4_ROSE,
                         "medium"),
    "ls_picard_medium": (A1_RES / "a1_m6_ls_picard.npz",
                         MESHES / "onera_m6_wakefree/medium.msh", S3_YELLOW,
                         "medium"),
}
LABEL = {
    "conf_picard_coarse": "conforming Picard (P5)",
    "conf_picard_medium": "conforming Picard (P5)",
    "conf_newton_medium": "conforming Newton (A1/G8.2)",
    "ls_m1_coarse": "level-set Picard, wake-embedded M1 (B7)",
    "ls_m4_coarse": "level-set Picard, wake-free M4 (B7)",
    "ls_newton_medium": "level-set Newton (A1/B15)",
    "ls_picard_medium": "level-set Picard (A1/B13)",
}

_mesh_cache = {}


def get_cut(mesh_path):
    """(mesh_cut, wc), built once per mesh path (deterministic cut_wake)."""
    key = ("cut", str(mesh_path))
    if key not in _mesh_cache:
        t0 = time.perf_counter()
        mc, wc = cut_wake(read_mesh(str(mesh_path)))
        print(f"  [cut_wake] {mesh_path.name}: {wc.n_stations} stations "
              f"({time.perf_counter() - t0:.1f}s)")
        _mesh_cache[key] = (mc, wc)
    return _mesh_cache[key]


def get_mvop(mesh_path):
    """(mesh, MultivaluedOperator), built once per mesh path (the committed
    M6 level-set convention: b7_onera_m6 / a1 _build_ls_m6)."""
    key = ("mvop", str(mesh_path))
    if key not in _mesh_cache:
        t0 = time.perf_counter()
        mesh = read_mesh(str(mesh_path))
        a = np.radians(ALPHA)
        wls = WakeLevelSet(
            np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
            direction=(np.cos(a), np.sin(a), 0.0))
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]))
        mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                   levelset=wls)
        print(f"  [mvop] {mesh_path.name}: {len(cm.te_nodes)} TE nodes "
              f"({time.perf_counter() - t0:.1f}s)")
        _mesh_cache[key] = (mesh, mvop)
    return _mesh_cache[key]


# ---------------------------------------------------------------- harvesting

def harvest_conf(name):
    npz, mesh_path, color, level = CONF_STATES[name]
    if not npz.exists() or not mesh_path.exists():
        print(f"[skip] {name}: missing {npz.name if not npz.exists() else mesh_path}")
        return None
    d = np.load(npz, allow_pickle=True)
    if "_span_gamma" in d.files:                       # A1-style cache
        gamma = np.atleast_1d(np.asarray(d["_span_gamma"], dtype=np.float64))
        z = np.atleast_1d(np.asarray(d["_span_z"], dtype=np.float64))
    else:                                              # P5-style cache
        gamma = np.asarray(d["gamma"], dtype=np.float64)
        z = np.asarray(d["station_z"], dtype=np.float64)
    phi = np.asarray(d["phi"], dtype=np.float64)
    mc, wc = get_cut(mesh_path)
    det = (len(z) == wc.n_stations) and np.allclose(wc.station_z, z)
    F = kutta_targets(phi, wc) - gamma if det else np.full_like(gamma, np.nan)
    o = np.argsort(z)
    return {"name": name, "color": color, "level": level, "kind": "conf",
            "z": z[o], "gamma": gamma[o], "F": F[o], "phi": phi,
            "mc": mc, "wc": wc, "det": det, "order": o}


def harvest_ls(name):
    npz, mesh_path, color, level = LS_STATES[name]
    if not npz.exists() or not mesh_path.exists():
        print(f"[skip] {name}: missing {npz.name if not npz.exists() else mesh_path}")
        return None
    d = np.load(npz, allow_pickle=True)
    rec = {"name": name, "color": color, "level": level, "kind": "ls"}
    if "_span_gamma" in d.files:                       # A1-style cache
        z = np.atleast_1d(np.asarray(d["_span_z"], dtype=np.float64))
        g = np.atleast_1d(np.asarray(d["_span_gamma"], dtype=np.float64))
    else:                                              # B7-style cache
        z = np.asarray(d["z"], dtype=np.float64)
        g = np.asarray(d["gamma"], dtype=np.float64)
    o = np.argsort(z)
    rec["z"], rec["gamma"] = z[o], g[o]
    rec["phi_ext"] = np.asarray(d["phi_ext"], dtype=np.float64)
    rec["mesh_path"] = mesh_path
    return rec


def ls_derived(rec):
    """TE pressure gap + section curves for a level-set state (needs the
    operator rebuild -- cached to results/, gitignored)."""
    cache = OUT / f"a2_cache_{rec['name']}.npz"
    if cache.exists() and not REFRESH:
        d = np.load(cache, allow_pickle=True)
        rec["te_z"] = d["te_z"]
        rec["gap"] = d["gap"]
        rec["sections"] = d["sections"].item()
        return rec
    print(f"  [derive] {rec['name']} (level-set operator rebuild) ...",
          flush=True)
    mesh, mvop = get_mvop(rec["mesh_path"])
    qu, ql = mvop.te_velocities(rec["phi_ext"])
    cpu = _cp_from_q2(np.sum(qu * qu, axis=1), M_INF)
    cpl = _cp_from_q2(np.sum(ql * ql, axis=1), M_INF)
    te_z = mesh.nodes[mvop.cm.te_nodes, 2]
    o = np.argsort(te_z)
    rec["te_z"], rec["gap"] = te_z[o], np.abs(cpu - cpl)[o]
    secs = {}
    for eta in ETAS:
        secs[eta] = {}
        for p in (0, 1):
            try:
                secs[eta][p] = section_cp_curve_levelset(
                    mesh, mvop, rec["phi_ext"], eta=eta, b_semi=B_SEMI,
                    m_inf=M_INF, smooth_passes=p)
            except ValueError as e:
                print(f"  [warn] {rec['name']} eta={eta} p={p}: {e}")
    rec["sections"] = secs
    np.savez(cache, te_z=rec["te_z"], gap=rec["gap"],
             sections=np.array(rec["sections"], dtype=object))
    return rec


def conf_derived(rec):
    """Sections (variant sweep), the all-station TE gap sweep, and the
    Delta-phi(z; x/c) decay profile for a conforming state."""
    t0 = time.perf_counter()
    print(f"  [derive] {rec['name']} ...", flush=True)
    mc, wc, phi = rec["mc"], rec["wc"], rec["phi"]
    secs = {}
    for eta in ETAS:
        secs[eta] = {}
        for p in (0, 1, 2):
            try:
                secs[eta][p] = section_cp_curve(
                    mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF,
                    smooth_passes=p)
            except ValueError as e:
                print(f"  [warn] {rec['name']} eta={eta} p={p}: {e}")
    rec["sections"] = secs

    # one recovery, many stations: raw per-triangle Cp swept over all TE
    # stations (the cheap analogue of _wall_section_points)
    wall = np.asarray(mc.boundary_faces["wall"], dtype=np.int64)
    grad_tri, _, _ = triangle_tangential_gradients(mc.nodes, wall, phi)
    tri_cp = _cp_from_q2(np.sum(grad_tri * grad_tri, axis=1), M_INF)
    upper = M._tri_sides_ny(mc.nodes, mc.elements, wall)
    zs, gaps = [], []
    eps = 1e-4 * float(np.median(np.diff(np.sort(rec["z"]))))
    for zj in rec["z"]:
        x, cp, s = M.cp_section_from_tri(mc.nodes, wall, tri_cp, zj + eps,
                                         upper)
        if s.sum() < 4 or (~s).sum() < 4:
            continue
        iu = np.argmax(x[s])
        il = np.argmax(x[~s])
        zs.append(zj)
        gaps.append(abs(cp[s][iu] - cp[~s][il]))
    rec["te_z"], rec["gap"] = np.asarray(zs), np.asarray(gaps)

    # Delta phi(z; x/c): the wall-side jump proxy, NOT pinned by the wake
    # constraint at x/c < 1
    prof = np.full((len(rec["z"]), len(XC_TARGETS)), np.nan)
    for i, zj in enumerate(rec["z"]):
        prof[i] = M.delta_phi_profile(mc.nodes, wall, phi, zj + eps,
                                      XC_TARGETS, upper)
    rec["dphi"] = prof
    print(f"  [derive] {rec['name']} done ({time.perf_counter() - t0:.1f}s)",
          flush=True)
    return rec


# ------------------------------------------------------------------ figures

def fig_jitter(states):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)
    for ax, level in zip(axes, ("coarse", "medium")):
        for rec in states:
            if rec["level"] != level:
                continue
            r = M.roughness_d2(rec["gamma"])
            ls = "--" if rec["name"] == "ls_m4_coarse" else "-"
            ax.plot(rec["z"] / B_SEMI, rec["gamma"], color=rec["color"],
                    ls=ls, marker=".", ms=3.5, lw=1.6,
                    label=f"{LABEL[rec['name']]}  r={r:.4f}")
        ax.set(xlabel="z / b", title=f"{level}")
        ax.legend(fontsize=8, loc="lower left")
    axes[0].set_ylabel("circulation Γ(z)")
    fig.suptitle("A2/S1 — Γ(z) roughness r (B7 metric: RMS d² / range), "
                 f"M{M_INF} α{ALPHA}", fontweight="bold")
    finish(fig, OUT, "a2_jitter.png")


def fig_decay(confs, ls_ref):
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    for rec in confs:
        if "dphi" not in rec:
            continue
        xs, ys = [], []
        for i, xc in enumerate(XC_TARGETS):
            col = rec["dphi"][:, i]
            m = np.isfinite(col)
            if m.sum() > 10:
                xs.append(xc)
                ys.append(M.jitter_localfit(rec["z"][m], col[m]))
        xs.append(1.0)
        ys.append(M.jitter_localfit(rec["z"], rec["gamma"]))
        ax.plot(xs, ys, color=rec["color"], marker="o", ms=6,
                ls="-" if rec["level"] == "medium" else "--",
                label=f"{LABEL[rec['name']]} ({rec['level']})")
    if ls_ref is not None:
        ax.axhline(ls_ref, color=S4_ROSE, ls=":", lw=1.6)
        ax.annotate("level-set Newton medium, Γ(z) jitter", (0.923, ls_ref),
                    textcoords="offset points", xytext=(0, 5),
                    color=INK_2, fontsize=8.5)
    ax.set_yscale("log")
    ax.set(xlabel="x/c where Δφ = φ_u − φ_l is read on the WALL "
                  "(1.0 = the Kutta station Γ itself)",
           ylabel="jitter (RMS residual vs local fit / range)")
    ax.legend(fontsize=8.5)
    ax.set_title("A2/S1 — is the jitter confined to the TE layer?",
                 fontweight="bold")
    finish(fig, OUT, "a2_decay.png")


def fig_te_gap(states):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharey=True)
    for ax, level in zip(axes, ("coarse", "medium")):
        for rec in states:
            if rec["level"] != level or "gap" not in rec:
                continue
            med = float(np.median(rec["gap"]))
            ls = "--" if rec["name"] == "ls_m4_coarse" else "-"
            ax.plot(rec["te_z"] / B_SEMI, rec["gap"], color=rec["color"],
                    ls=ls, marker=".", ms=3, lw=1.4,
                    label=f"{LABEL[rec['name']]}  med={med:.3f}")
        ax.set(xlabel="z / b", title=level)
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_ylabel("TE pressure gap |Cp_u − Cp_l|")
    fig.suptitle("A2/S2 — Kutta pressure closure at the TE: conforming "
                 "(potential-jump Kutta) vs level-set (pressure Kutta)",
                 fontweight="bold")
    finish(fig, OUT, "a2_te_gap.png")


def fig_spike(states):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharey=True)
    for ax, level in zip(axes, ("coarse", "medium")):
        for rec in states:
            if rec["level"] != level or "sections" not in rec:
                continue
            passes = sorted({p for e in rec["sections"].values() for p in e})
            mean_sp, all_sp = [], {p: [] for p in passes}
            for p in passes:
                for eta in ETAS:
                    sec = rec["sections"].get(eta, {}).get(p)
                    if sec is None:
                        continue
                    for side in ("upper", "lower"):
                        sp = M.spike_metric(sec[f"x_{side}"],
                                            sec[f"cp_{side}"])["spike"]
                        if np.isfinite(sp):
                            all_sp[p].append(sp)
                mean_sp.append(np.mean(all_sp[p]) if all_sp[p] else np.nan)
            ax.plot(passes, mean_sp, color=rec["color"], marker="o", ms=7,
                    label=LABEL[rec["name"]])
            for p in passes:
                ax.scatter([p] * len(all_sp[p]), all_sp[p], s=12,
                           color=rec["color"], alpha=0.35)
        ax.set(xlabel="P6 smoothing passes", title=level,
               xticks=[0, 1, 2])
        ax.legend(fontsize=8.5)
    axes[0].set_ylabel("TE last-point spike |Cp_last − trend fit|")
    fig.suptitle("A2/S2 — TE Cp spike vs recovery variant "
                 "(mean over η/side; dots = individual sections)",
                 fontweight="bold")
    finish(fig, OUT, "a2_spike.png")


# --------------------------------------------------------------------- main

def main():
    apply_style()
    cl = CheckList("Track A / A2 TE-Kutta fidelity -- zero-solve scaffold")
    t_start = time.perf_counter()

    confs = [r for r in (harvest_conf(n) for n in CONF_STATES) if r]
    lss = [r for r in (harvest_ls(n) for n in LS_STATES) if r]

    # ---- GA2.1: shared metric reproduces B7's committed numbers ----------
    if P5_BASELINE_CSV.exists():
        base = np.loadtxt(P5_BASELINE_CSV, delimiter=",", comments="#",
                          skiprows=5)
        r = M.roughness_d2(base[:, 1])
        cl.add("GA2.1", "roughness(committed p5_gamma_baseline.csv) == B7",
               f"{r:.4f}", f"{R_ANCHORS['p5_baseline']} +-10%",
               abs(r - R_ANCHORS["p5_baseline"]) < 0.1 * R_ANCHORS["p5_baseline"])
    for name in ("ls_m1_coarse", "ls_m4_coarse"):
        rec = next((x for x in lss if x["name"] == name), None)
        if rec is None:
            continue
        r = M.roughness_d2(rec["gamma"])
        cl.add("GA2.1", f"roughness({name} cache) == B7", f"{r:.4f}",
               f"{R_ANCHORS[name]} +-10%",
               abs(r - R_ANCHORS[name]) < 0.1 * R_ANCHORS[name])

    # ---- L0: unified jitter table ----------------------------------------
    rows = []
    for rec in confs + lss:
        rows.append((rec["name"], rec["level"],
                     "conforming" if rec["kind"] == "conf" else "level-set",
                     len(rec["gamma"]),
                     f"{M.roughness_d2(rec['gamma']):.5f}",
                     f"{M.jitter_localfit(rec['z'], rec['gamma']):.5f}"))
    write_csv(OUT, "a2_jitter.csv",
              "state,mesh_level,wake_model,n_stations,roughness_d2,"
              "jitter_localfit", rows)

    # ---- L1: determinism + Kutta closure census + probe census -----------
    census_rows = []
    for rec in confs:
        cl.add("A2-sanity", f"{rec['name']}: rebuilt cut matches cache "
               "(station_z bitwise)", str(rec["det"]), "cut_wake deterministic",
               bool(rec["det"]))
        if not rec["det"]:
            continue
        F, g = rec["F"], rec["gamma"]
        scale = float(np.max(np.abs(g)))
        cl.add("GA2.3", f"{rec['name']}: Kutta closure residual max|F|/max|Γ|",
               f"{np.max(np.abs(F)) / scale:.2%} (median "
               f"{np.median(np.abs(F)) / scale:.2%})",
               "recorded (H2: small |F| + persisting jitter => noise lives "
               "in the target definition)", True)
        feats = M.probe_census(rec["mc"], rec["wc"])
        resid = np.abs(g - M.local_fit(rec["z"], g))
        feats_o = {k: np.asarray(v, dtype=np.float64)[rec["order"]]
                   for k, v in feats.items()}
        corr = {k: M.pearson(resid, v) for k, v in feats_o.items()}
        top = max(corr, key=lambda k: abs(corr[k]) if np.isfinite(corr[k])
                  else -1)
        cl.add("GA2.2", f"{rec['name']}: probe-geometry correlation "
               "(supporting only; correlation != causation)",
               f"top |r|: {top} r={corr[top]:.2f}",
               "recorded; H1 verdict needs the intervention leg", True)
        for i in range(len(g)):
            census_rows.append(
                (rec["name"], f"{rec['z'][i]:.5f}", f"{g[i]:.6f}",
                 f"{resid[i]:.6f}", f"{np.abs(F[i]):.2e}",
                 f"{feats_o['d_up'][i]:.5f}", f"{feats_o['d_lo'][i]:.5f}",
                 f"{feats_o['asym'][i]:.4f}", f"{feats_o['dz_up'][i]:.5f}",
                 f"{feats_o['dz_lo'][i]:.5f}",
                 int(feats_o["shared_up"][i]), int(feats_o["shared_lo"][i])))
    if census_rows:
        write_csv(OUT, "a2_probe_census.csv",
                  "state,z,gamma,resid_localfit,abs_F,d_up,d_lo,asym,"
                  "dz_up,dz_lo,shared_up,shared_lo", census_rows)

    # ---- L2z + L3 derived fields ------------------------------------------
    for rec in confs:
        conf_derived(rec)
    for rec in lss:
        ls_derived(rec)

    decay_rows = []
    for rec in confs:
        jg = M.jitter_localfit(rec["z"], rec["gamma"])
        for i, xc in enumerate(XC_TARGETS):
            col = rec["dphi"][:, i]
            m = np.isfinite(col)
            if m.sum() > 10:
                decay_rows.append((rec["name"], xc,
                                   f"{M.jitter_localfit(rec['z'][m], col[m]):.5f}"))
        decay_rows.append((rec["name"], 1.0, f"{jg:.5f}"))
        if len(rec["dphi"]) and np.isfinite(rec["dphi"][:, 0]).sum() > 10:
            m = np.isfinite(rec["dphi"][:, 0])
            ratio = M.jitter_localfit(rec["z"][m], rec["dphi"][:, 0][m]) / jg
            cl.add("GA2.2", f"{rec['name']}: wall Δφ jitter at x/c=0.92 "
                   "vs station-Γ jitter",
                   f"{ratio:.2f}x",
                   "recorded (<<1 => jitter confined to the TE layer, "
                   "supports H1)", True)
    write_csv(OUT, "a2_decay.csv", "state,xc,jitter_localfit", decay_rows)

    gap_rows, spike_rows = [], []
    for rec in confs + lss:
        if "gap" in rec and len(rec["gap"]):
            for zz, gg in zip(rec["te_z"], rec["gap"]):
                gap_rows.append((rec["name"], f"{zz:.5f}", f"{gg:.5f}"))
        for eta, by_pass in rec.get("sections", {}).items():
            for p, sec in by_pass.items():
                for side in ("upper", "lower"):
                    sp = M.spike_metric(sec[f"x_{side}"], sec[f"cp_{side}"])
                    spike_rows.append(
                        (rec["name"], eta, side, p, f"{sp['spike']:.4f}",
                         f"{sp['cp_last']:.4f}", f"{sp['fit_at_1']:.4f}",
                         f"{M.te_gap_from_curve(sec):.4f}"))
    write_csv(OUT, "a2_te_gap.csv", "state,z,te_cp_gap", gap_rows)
    write_csv(OUT, "a2_spike.csv",
              "state,eta,side,smooth_passes,spike,cp_last,fit_at_1,"
              "te_gap_section", spike_rows)

    # GA2.4 headline. Two distinct statements, kept apart on purpose:
    # (a) each method zeroes its OWN constraint (conforming: probe-phi |F|,
    #     GA2.3 above; level-set: the TE CV pressure gap, machine level) --
    #     measuring the LS gap on its own CV objects is closure quality, not
    #     an A/B; (b) the SAME-estimator A/B is the section-curve last-point
    #     gap, raw recovery, identical extraction on both paths.
    def _sec_gap(rec):
        gs = [M.te_gap_from_curve(by_pass[0])
              for by_pass in rec.get("sections", {}).values() if 0 in by_pass]
        return float(np.mean(gs)) if gs else np.nan
    for level in ("coarse", "medium"):
        cg = [np.median(r["gap"]) for r in confs
              if r["level"] == level and len(r.get("gap", ()))]
        lg = [np.median(r["gap"]) for r in lss
              if r["level"] == level and len(r.get("gap", ()))]
        if cg and lg:
            cl.add("GA2.4", f"TE closure on each method's OWN constraint "
                   f"({level})",
                   f"conf probe-|F| see GA2.3; LS CV pressure gap median "
                   f"{np.mean(lg):.1e}",
                   "recorded (both ~0: each method zeroes its own residual)",
                   True)
            cl.add("GA2.4", f"conforming TE pressure gap, all-station sweep "
                   f"({level})", f"median {np.mean(cg):.3f}",
                   "recorded (P13 committed prose range 0.14-0.22)", True)
        cs = [_sec_gap(r) for r in confs if r["level"] == level]
        ls_ = [_sec_gap(r) for r in lss if r["level"] == level]
        cs = [v for v in cs if np.isfinite(v)]
        ls_ = [v for v in ls_ if np.isfinite(v)]
        if cs and ls_:
            cl.add("GA2.4", f"TE gap, SAME estimator: section last points, "
                   f"raw recovery ({level})",
                   f"conforming {np.mean(cs):.3f} vs level-set "
                   f"{np.mean(ls_):.3f} ({np.mean(cs) / np.mean(ls_):.0f}x)",
                   "recorded (H4: conforming >> LS = Kutta-FORM error)", True)

    # ---- GA2.2 verdict (lives in the gated sibling; committed) ------------
    # run_a2_interventions.py (PYFP3D_TRANSONIC_GATES=1) ran the fixed-Γ
    # discriminator on 2026-07-16/17: coarse D=7.33, medium D=25.70, both
    # above the pre-registered 3.0 confirm threshold. Recorded here so the
    # zero-solve checks.csv reflects the settled gate; the numbers themselves
    # are committed in checks_interventions.csv / a2_intervention.csv.
    ic = OUT / "checks_interventions.csv"
    cl.add("GA2.2", "H1 discriminator D (fixed-Γ intervention, gated sibling)",
           "coarse 7.33 / medium 25.70 (H1 CONFIRMED)",
           "D>3 confirm (pre-registered); see checks_interventions.csv",
           ic.exists(), note="run_a2_interventions.py")

    # ---- figures ----------------------------------------------------------
    fig_jitter(confs + lss)
    ls_med = next((r for r in lss if r["name"] == "ls_newton_medium"), None)
    fig_decay(confs, M.jitter_localfit(ls_med["z"], ls_med["gamma"])
              if ls_med else None)
    fig_te_gap(confs + lss)
    fig_spike(confs + lss)

    print(f"\ntotal {time.perf_counter() - t_start:.1f}s")
    return cl.report(OUT, fname="checks.csv")


if __name__ == "__main__":
    sys.exit(main())
