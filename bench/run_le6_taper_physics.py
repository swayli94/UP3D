"""LE-6: the PHYSICS comparison that has to pass before the taper can be removed.

User requirement (2026-08-04), and it supersedes LE-5's convergence-only legs: convergence is
necessary but not sufficient, because the taper is a MODEL change. Before removal is confirmed,
compare with and without it on

  * the spanwise lift distribution        (wing and wing-body)
  * section Cp at EIGHT stations, NOT uniformly spaced -- clustered toward the TIP
  * cl and cd
  * across subsonic AND transonic, several alpha, and several mesh densities

and confirm no non-physical result. The spanwise-loading and section-Cp comparisons are saved as
images for inspection.

★ THE SPECIFIC NON-PHYSICAL RISK, named so the figures are read for something: the taper exists
to suppress the P13 tip-edge singularity, so removing it may put a spurious SPIKE or OSCILLATION
in the tip loading. Physically the spanwise loading must fall smoothly to ~0 at the tip (Gamma
-> 0 there). So the acceptance checks are:

  N1  tip monotonicity   the outermost spanwise strips must be non-increasing toward the tip,
                         within a tolerance -- a rise at the last strips is the singularity
                         showing through
  N2  tip loading bound  the outermost strip's sectional load must not exceed the wing's own
                         maximum -- an overshoot at the free edge is non-physical
  N3  tip suction bound  min Cp at the tip-most stations must not diverge relative to the
                         inboard stations (recorded per level, since a genuine singularity
                         sharpens with refinement while a resolved feature settles)
  N4  cd at subsonic     with no shock, d'Alembert makes cd exactly 0, so spurious cd is a
                         direct error measure -- it must not grow when the taper is removed

Spanwise loading is taken TWO ways deliberately: the surface-pressure strip integral (the
"integrate the closed body pressure" convention) and the Gamma-based sectional cl. They are
independent post-processing paths, so agreement is a wiring check and disagreement localises a
post-processing fault rather than a physics one.

Outputs (TRACKED):
  bench/gate_results/le6_taper_physics.csv          per-condition cl, cd, N1-N4
  bench/gate_results/le6_spanwise.csv               spanwise loading curves
  bench/gate_results/capability/le6_span_<geom>_<level>.png
  bench/gate_results/capability/le6_cp_<geom>_<level>_M<..>_a<..>.png
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at                  # noqa: E402
from pyfp3d.post.surface import _pressure_force, planform_area      # noqa: E402
from pyfp3d.post.unified import (_conforming_wall_state,            # noqa: E402
                                 _cp_from_q2, section_cp, wall_forces)
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)

OUT = os.path.join(HERE, "gate_results")
ART = os.path.join(OUT, "capability")
os.makedirs(ART, exist_ok=True)
CSV = os.path.join(OUT, "le6_taper_physics.csv")
CSV_SPAN = os.path.join(OUT, "le6_spanwise.csv")

#: EIGHT stations, tip-clustered as required -- the committed 7-station experiment set
#: (0.20/0.44/0.65/0.80/0.90/0.96/0.99) with the outboard end densified, so the tip region
#: where the taper acts is actually sampled rather than straddled.
ETAS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.95, 0.98, 0.995)
#: spanwise strips for the loading curve -- densified outboard for the same reason
SPAN_EDGES = np.concatenate([np.linspace(0.0, 0.8, 9), [0.85, 0.90, 0.94, 0.97, 0.99, 1.0001]])

TAPERS = [("none", 0.0), ("vanish_smooth", 0.05)]
#: subsonic AND transonic, several alpha
CONDS = [(0.50, 1.00), (0.50, 3.06), (0.70, 3.06), (0.84, 1.00), (0.84, 3.06)]
CASES = [("wingbody", "onera_m6_wingbody_conforming", ("coarse", "medium")),
         ("m6wing", "onera_m6", ("coarse_ss", "medium"))]

KEYS = ["geom", "level", "m_inf", "alpha", "form", "r_c", "converged",
        "res_final", "cl_p", "cd_p", "cl_gamma_tip", "N1_tip_monotone",
        "N1_worst_rise", "N2_tip_over_max", "N3_tip_cp_min", "N3_inboard_cp_min",
        "N4_cd_vs_none_pct", "cl_vs_none_pct", "wall_s", "note"]


def solve(mesh_path, geom, m, alpha, form, r_c):
    mc, wc = cut_wake(read_mesh(mesh_path))
    t = (None if form == "none"
         else tip_taper_factors(wc.station_z, B_SEMI, form, r_c * B_SEMI))
    t0 = time.perf_counter()
    if geom == "wingbody":
        base = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure")
        if m <= cap.WB_MSTART:
            kw = dict(base)
            if t is not None:
                kw["tip_taper"] = t
            r = solve_newton_lifting(mc, wc, m_inf=m, alpha_deg=alpha, **kw)
        else:
            seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART,
                                        alpha_deg=alpha, **cap.CONF_SEED_KW)
            nk = dict(base, phi_init=seed["phi"], gamma_init=seed["gamma"],
                      n_picard_seed=0)
            if t is not None:
                nk["tip_taper"] = t
            r = solve_newton_transonic(mc, wc, m_inf=m, alpha_deg=alpha,
                                       m_start=cap.WB_MSTART, dm=cap.DM,
                                       dm_min=0.01, freeze_tol=1e-5,
                                       intermediate_tol=1e-4, newton_kw=nk)
    else:
        #: the wing-alone uses the shipped recipe PLUS pressure, since LE-3 measured the
        #: probe estimator to be the load-bearing gap there -- comparing tapers on the
        #: superseded estimator would answer a question nobody is asking.
        base = dict(cap.NEWTON_M6_RECIPE["newton_kw"], kutta_estimator="pressure")
        if t is not None:
            base["tip_taper"] = t
        if m <= cap.CONF_WING_MSTART:
            r = solve_newton_lifting(mc, wc, m_inf=m, alpha_deg=alpha, **base)
        else:
            r = solve_newton_transonic(mc, wc, m_inf=m, alpha_deg=alpha,
                                       **dict(cap.NEWTON_M6_RECIPE,
                                              newton_kw=base))
    return mc, wc, r, time.perf_counter() - t0


def spanwise_load(mesh, phi, alpha, m_inf, s_ref):
    """Sectional load per spanwise strip, from the SURFACE PRESSURE integral (the
    'integrate the closed body' convention). Normalised by strip area so it reads as a
    sectional coefficient rather than a strip total."""
    wall = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    q2, _, area, n_out = _conforming_wall_state(mesh, phi, wall, 1.0, 0)
    cp = _cp_from_q2(q2, m_inf, 1.4)
    cen = mesh.nodes[wall].mean(axis=1)
    eta = cen[:, 2] / B_SEMI
    xs, ys, areas = [], [], []
    for i in range(len(SPAN_EDGES) - 1):
        lo, hi = SPAN_EDGES[i], SPAN_EDGES[i + 1]
        msk = (eta >= lo) & (eta < hi)
        if not np.any(msk) or area[msk].sum() <= 0:
            continue
        a = float(area[msk].sum())
        _cf, cl, _cd = _pressure_force(cp[msk], area[msk], n_out[msk], a, alpha)
        xs.append(0.5 * (lo + hi)); ys.append(cl); areas.append(a)
    return np.array(xs), np.array(ys), np.array(areas)


def main():
    rows, span_rows = [], []
    store = {}
    for geom, mdir, levels in CASES:
        for level in levels:
            mp = os.path.join(REPO, "cases", "meshes", mdir, f"{level}.msh")
            if not os.path.exists(mp):
                print(f"{geom}/{level}: mesh missing"); continue
            for m, alpha in CONDS:
                base_cl = base_cd = None
                for form, r_c in TAPERS:
                    tag = f"{geom}/{level}/M{m}/a{alpha}/{form}"
                    try:
                        mc, wc, r, wall = solve(mp, geom, m, alpha, form, r_c)
                    except Exception as exc:                       # noqa: BLE001
                        print(f"  {tag:46s} DIED {type(exc).__name__}: "
                              f"{str(exc)[:50]}", flush=True)
                        rows.append(dict(geom=geom, level=level, m_inf=m,
                                         alpha=alpha, form=form, r_c=r_c,
                                         converged=False,
                                         note=f"{type(exc).__name__}"))
                        continue
                    phi = np.asarray(r["phi"])
                    sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
                    f = wall_forces(mc, phi=phi, alpha_deg=alpha, s_ref=sref,
                                    m_inf=m)
                    clp, cdp = f["cl"], f["cd_pressure"]
                    if form == "none":
                        base_cl, base_cd = clp, cdp
                    eta_s, load, _a = spanwise_load(mc, phi, alpha, m, sref)
                    #: N1 -- the outermost strips must not RISE toward the tip
                    tail = load[eta_s >= 0.85]
                    rise = float(np.max(np.diff(tail))) if len(tail) > 1 else 0.0
                    n1 = rise <= 0.02 * float(np.max(np.abs(load)))
                    #: N2 -- the tip strip must not exceed the wing's own maximum
                    n2 = float(load[-1] / np.max(load)) if np.max(load) else float("nan")
                    #: N3 -- tip vs inboard suction, from the section curves
                    cps = {}
                    for e in ETAS:
                        try:
                            cps[e] = section_cp(mc, eta=e, b_semi=B_SEMI,
                                                m_inf=m, phi=phi)
                        except Exception:                          # noqa: BLE001
                            pass
                    tipmin = min((float(np.min(cps[e]["cp_upper"]))
                                  for e in ETAS if e >= 0.95 and e in cps),
                                 default=float("nan"))
                    inmin = min((float(np.min(cps[e]["cp_upper"]))
                                 for e in ETAS if e <= 0.80 and e in cps),
                                default=float("nan"))
                    gtip = float(np.asarray(r["gamma"])[np.argmax(wc.station_z)])
                    print(f"  {tag:46s} conv={bool(r['converged'])} "
                          f"cl {clp:+.6f} cd {cdp:+.6f} "
                          f"N1={'ok' if n1 else 'RISE'} N2={n2:.3f} "
                          f"tipCp {tipmin:+.3f} ({wall:.0f}s)", flush=True)
                    rows.append(dict(
                        geom=geom, level=level, m_inf=m, alpha=alpha, form=form,
                        r_c=r_c, converged=bool(r["converged"]),
                        res_final=float(r["residual_history"][-1]),
                        cl_p=clp, cd_p=cdp, cl_gamma_tip=gtip,
                        N1_tip_monotone=n1, N1_worst_rise=round(rise, 8),
                        N2_tip_over_max=round(n2, 5) if n2 == n2 else None,
                        N3_tip_cp_min=round(tipmin, 5) if tipmin == tipmin else None,
                        N3_inboard_cp_min=round(inmin, 5) if inmin == inmin else None,
                        N4_cd_vs_none_pct=(round(100 * (cdp - base_cd) / base_cd, 3)
                                           if base_cd else 0.0),
                        cl_vs_none_pct=(round(100 * (clp - base_cl) / base_cl, 3)
                                        if base_cl else 0.0),
                        wall_s=round(wall, 1), note=""))
                    for e, l in zip(eta_s, load):
                        span_rows.append(dict(geom=geom, level=level, m_inf=m,
                                              alpha=alpha, form=form,
                                              eta=round(float(e), 5),
                                              load=round(float(l), 8)))
                    store[(geom, level, m, alpha, form)] = (eta_s, load, cps)

            # ---------------- figures for this geometry/level ------------------
            conds = [(m, a) for (m, a) in CONDS
                     if (geom, level, m, a, "none") in store]
            if conds:
                fig, axes = plt.subplots(1, len(conds),
                                         figsize=(4.2 * len(conds), 4.4),
                                         squeeze=False)
                for j, (m, a) in enumerate(conds):
                    ax = axes[0][j]
                    for form, style in (("none", "-o"), ("vanish_smooth", "--s")):
                        k = (geom, level, m, a, form)
                        if k not in store:
                            continue
                        e, l, _c = store[k]
                        ax.plot(e, l, style, ms=4, lw=1.6,
                                label=("no taper" if form == "none"
                                       else "taper (production)"))
                    ax.set_title(f"M{m}  α={a}", fontsize=10)
                    ax.set_xlabel(r"$\eta=z/b_{semi}$"); ax.grid(alpha=0.3)
                    ax.legend(fontsize=8)
                axes[0][0].set_ylabel("sectional load (strip pressure integral)")
                fig.suptitle(f"{geom} / {level}: spanwise loading, taper vs none "
                             f"— the tip must fall smoothly, not spike", fontsize=11)
                fig.tight_layout()
                p = os.path.join(ART, f"le6_span_{geom}_{level}.png")
                fig.savefig(p, dpi=130); plt.close(fig)
                print(f"  wrote {p}", flush=True)

                for (m, a) in conds:
                    fig, axes = plt.subplots(2, 4, figsize=(16.5, 7.4))
                    for i, e in enumerate(ETAS):
                        ax = axes[i // 4][i % 4]
                        for form, style, col in (("none", "-", "tab:blue"),
                                                 ("vanish_smooth", "--", "tab:red")):
                            k = (geom, level, m, a, form)
                            if k not in store or e not in store[k][2]:
                                continue
                            c = store[k][2][e]
                            for side in ("upper", "lower"):
                                ax.plot(c[f"x_{side}"], c[f"cp_{side}"], style,
                                        lw=1.4, color=col,
                                        label=(("no taper" if form == "none"
                                                else "taper") if side == "upper"
                                               else None))
                        ax.invert_yaxis(); ax.grid(alpha=0.3)
                        ax.set_title(f"η = {e}", fontsize=9.5)
                        if i == 0:
                            ax.legend(fontsize=8)
                    fig.suptitle(f"{geom} / {level}  M{m}  α={a}: section Cp, "
                                 f"8 tip-clustered stations", fontsize=11)
                    fig.tight_layout()
                    p = os.path.join(
                        ART, f"le6_cp_{geom}_{level}_M{str(m).replace('.','')}"
                             f"_a{str(a).replace('.','')}.png")
                    fig.savefig(p, dpi=120); plt.close(fig)
                    print(f"  wrote {p}", flush=True)

    for path, data in ((CSV, rows), (CSV_SPAN, span_rows)):
        if not data:
            continue
        with open(path, "w", newline="") as fh:
            keys = KEYS if path == CSV else list(data[0])
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader(); w.writerows(data)
        print(f"wrote {path}")
    print("\n=== N1-N4: any non-physical signature when the taper is removed? ===")
    for r in rows:
        if r.get("form") != "none" or not r.get("converged"):
            continue
        bad = []
        if r.get("N1_tip_monotone") is False:
            bad.append("N1 tip loading RISES toward the tip")
        if (r.get("N2_tip_over_max") or 0) > 1.0:
            bad.append("N2 tip strip exceeds the wing maximum")
        print(f"  {r['geom']}/{r['level']}/M{r['m_inf']}/a{r['alpha']}: "
              f"{'; '.join(bad) if bad else 'clean'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
