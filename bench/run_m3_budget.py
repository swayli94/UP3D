"""M3 error budget: is wall accuracy the binding term, or is it shock position?

Pre-registered in docs/dev_phase_two/20260731-1600-m3-budget-prereg.md and committed
before this file was written. Read that for the protocol and, more importantly, for the
interpretation rules R1-R4 -- in particular R2, the validity gate that decides whether
the entropy-ON leg is a physics reading at all or the G8.2 3-D donor-cycle defect
reproducing.

Extractor reused verbatim from cases/analysis/v5_3_m6_cp/run.py (same-extractor
discipline): the same committed experiment file, the same station_rms, the same seven
stations with the first five outside the tip mask.

Outputs (TRACKED): bench/gate_results/m3_budget.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_report                          # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,           # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.newton import (solve_newton_lifting,               # noqa: E402
                                 solve_newton_transonic)
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

OUT = os.path.join(HERE, "gate_results")
os.makedirs(OUT, exist_ok=True)

M_INF, ALPHA = 0.8395, 3.06            # TEST 2308 dataset label verbatim
MAC = 0.64607
ETAS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)
N_UNMASKED = 5
EXP_FILE = os.path.join(REPO, "cases", "reference_data",
                        "onera_m6_experiment", "experiment-Cp.dat")
#: GV5.3's committed same-extractor k = 0 pooled RMS -- R1's reference at medium
GV53_K0_POOLED = {"medium": 0.1288}
#: P14 pressure-Kutta anchors (cl_p, cl_KJ), GV5.3's P14_ANCHOR verbatim
P14_ANCHOR = {"coarse": (0.262778, 0.268813),
              "medium": (0.277628, 0.282263)}
M6_NEWTON_KW = dict(farfield_spanwise_gamma=True, precond="direct",
                    direct_refactor_every=1000, n_newton_max=60)
#: (entropy, kutta) legs. The probe legs are the first round's; the
#: pressure legs are what the committed P14 anchor actually used.
LEGS = ((False, "probe"), (True, "probe"),
        (False, "pressure"), (True, "pressure"))


def parse_experiment(path=EXP_FILE):
    """GV5.3's parser verbatim."""
    zones = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("ZONE"):
                zones.append([])
            elif s and not s.startswith(("TITLE", "VARIABLES")):
                p = s.split()
                if len(p) == 5:
                    zones[-1].append((float(p[1]), float(p[2]), float(p[3]),
                                      float(p[4])))
    exp = {}
    for z in zones:
        a = np.asarray(z)
        exp[round(float(a[0, 1]), 2)] = {"x": a[:, 0], "zl": a[:, 2],
                                         "cp": a[:, 3], "upper": a[:, 2] > 0.0}
    return exp


def station_rms(curves, exp, eta):
    """GV5.3's station_rms verbatim (flip=False, the W2-validated mapping)."""
    e = exp[eta]
    tot, n = 0.0, 0
    for want_upper in (True, False):
        side = "upper" if want_upper else "lower"
        m = e["upper"] == want_upper
        if not np.any(m):
            continue
        cp_i = np.interp(e["x"][m], curves[eta][f"x_{side}"],
                         curves[eta][f"cp_{side}"])
        tot += float(np.sum((cp_i - e["cp"][m]) ** 2))
        n += int(m.sum())
    return (tot / max(n, 1)) ** 0.5, n


#: the x/c decomposition bands. Registered here rather than in the pre-registration
#: because the pre-registration's R3 has already been ANSWERED (shock position is not a
#: leading term at coarse) -- this is a follow-up decomposition of what IS left, it
#: carries no pass/fail, and it is free post-processing on states already solved.
BANDS = (("LE", 0.0, 0.15), ("MID", 0.15, 0.85), ("TE", 0.85, 1.01))


def band_rms(curves, exp, eta):
    """Sum-of-squares of the same Cp residual, split by x/c band and by side.

    Same interpolation and same side mapping as station_rms -- so the bands add up
    to it by construction (asserted at the call site), which is the only way a
    decomposition can be trusted not to be measuring something else.
    """
    e = exp[eta]
    out = {}
    for name, lo, hi in BANDS:
        for want_upper in (True, False):
            side = "upper" if want_upper else "lower"
            m = (e["upper"] == want_upper) & (e["x"] >= lo) & (e["x"] < hi)
            if not np.any(m):
                out[f"{name}_{side}"] = (0.0, 0)
                continue
            cp_i = np.interp(e["x"][m], curves[eta][f"x_{side}"],
                             curves[eta][f"cp_{side}"])
            out[f"{name}_{side}"] = (float(np.sum((cp_i - e["cp"][m]) ** 2)),
                                     int(m.sum()))
    return out


def solve(mc, wc, entropy, kutta="probe", n_newton_max=None):
    """The P14 transonic recipe verbatim, entropy (and now the Kutta form) variable.

    ★ 2026-07-31: `kutta` was added after the first budget round measured its own cl
    4.8 % below the committed P14 anchor and traced that to the estimator, not to the
    seed. NEWTON_M6_RECIPE inherits the library default "probe", while the committed
    P14 M0.84 anchor was produced with "pressure" -- and P14's own headline was that
    the probe path sits 4.5 % / 4.3 % BELOW the pressure and level-set answers at
    exactly this condition. The pressure leg follows P14's recipe: level 0 is seeded
    from a PROBE Newton solve at M0.70, because the quadratic Kutta row has a smaller
    basin and the transonic driver only cold-seeds level 0.

    NEWTON_M6_RECIPE already carries its own newton_kw, so the flag is MERGED into
    it rather than passed alongside -- passing a second newton_kw is a TypeError,
    and silently replacing the recipe's would have changed precond /
    direct_refactor_every / freeze_refresh_max underneath the comparison. That is
    also why M6_NEWTON_KW below is only an assertion target: the recipe is the
    single source of truth for the solver settings.
    """
    kw = dict(NEWTON_M6_RECIPE)
    # ★ the drift guard runs on the RECIPE, before any intentional override --
    # otherwise a deliberate, recorded deviation (n_newton_max below) trips the
    # very check that exists to catch UNintended drift. It did exactly that on
    # 2026-08-01 and cost a four-leg sweep, which is the argument for the order.
    for k, v in M6_NEWTON_KW.items():
        assert kw["newton_kw"][k] == v, (
            f"the P14 recipe's newton_kw[{k}] = {kw['newton_kw'][k]} no longer "
            f"matches GV5.3's recorded {v} -- the comparison basis moved")
    kw["newton_kw"] = dict(kw["newton_kw"], entropy_correction=entropy)
    if n_newton_max is not None:
        # recorded deviation from "the recipe verbatim" -- see the caller
        kw["newton_kw"]["n_newton_max"] = int(n_newton_max)
    if kutta == "pressure":
        r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA,
                                  entropy_correction=entropy, **M6_NEWTON_KW)
        kw["newton_kw"].update(kutta_estimator="pressure", phi_init=r0["phi"],
                               gamma_init=r0["gamma"], n_picard_seed=0)
    return solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)


def main(levels=("coarse",)):
    exp = parse_experiment()
    # GV5.3's W2 experiment-side guard, kept: a station whose max Cp is not at the
    # LE would mean the side mapping is broken and every RMS below is meaningless.
    for eta, e in exp.items():
        i_le = int(np.argmax(e["cp"]))
        if e["x"][i_le] >= 0.05:
            raise RuntimeError(f"W2 guard: eta={eta} max-Cp at x/c="
                               f"{e['x'][i_le]:.4f}, not the LE")
    rows, pooled = [], {}
    for level in levels:
        path = os.path.join(REPO, "cases", "meshes", "onera_m6", f"{level}.msh")
        if not os.path.exists(path):
            print(f"skip {level}: mesh missing (regenerate via "
                  f"cases/meshes/onera_m6/generate_onera_m6.py)")
            continue
        mc, wc = cut_wake(read_mesh(path))
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        for entropy, kutta in LEGS:
            tag = ("ON" if entropy else "OFF") + f"/{kutta}"
            t0 = time.perf_counter()
            r = solve(mc, wc, entropy, kutta)
            wall = time.perf_counter() - t0
            phi = np.asarray(r["phi"])
            gamma = np.atleast_1d(np.asarray(r["gamma"]))
            res = float(r.get("residual_history", [np.nan])[-1])
            print(f"\n{level} entropy {tag}: conv={r.get('converged')} "
                  f"|R|={res:.2e} m_final={r.get('m_final', M_INF)} "
                  f"({wall:.0f}s)", flush=True)

            # ---- B3: the R2 validity gate, read BEFORE any RMS is quoted -----
            smin = r.get("sigma_min")
            sconv = r.get("sigma_converged", r.get("entropy_converged"))
            nsh, m1 = r.get("n_shock_cells"), r.get("m1_max")
            nflr = r.get("n_floored")
            print(f"  B3 entropy health: sigma_min={smin} n_shock={nsh} "
                  f"m1_max={m1} transport_converged={sconv} n_floored={nflr}")
            r2_ok = True
            if entropy:
                r2_ok = (smin is not None and float(smin) > 0.0
                         and (sconv is None or bool(sconv)))
                print(f"  R2 validity gate: "
                      f"{'PASS -- this IS an entropy reading' if r2_ok else 'FAIL -- the G8.2 3-D donor-cycle defect reproducing; its RMS is NOT an entropy effect'}")

            # ---- B1/B2 ------------------------------------------------------
            curves = {}
            for eta in ETAS:
                curves[eta] = section_cp_curve(mc, phi, eta=eta,
                                               b_semi=B_SEMI, m_inf=M_INF)
            per = {}
            for eta in ETAS:
                rms, npts = station_rms(curves, exp, eta)
                sh = shock_report(curves[eta], M_INF)
                xs = sh["upper"].get("x_shock")
                per[eta] = (rms, npts, xs)
                mask = "" if eta < 0.95 else "  (tip-masked, RECORDED-only)"
                print(f"    eta {eta:.2f}: RMS {rms:.4f} ({npts} pts), "
                      f"x_shock {'missing' if xs is None else f'{xs:.4f}'}"
                      f"{mask}")
            # -- the band decomposition over the unmasked stations --------
            acc = {}
            for eta in ETAS[:N_UNMASKED]:
                for kk, (ss, nn) in band_rms(curves, exp, eta).items():
                    a0, n0 = acc.get(kk, (0.0, 0))
                    acc[kk] = (a0 + ss, n0 + nn)
                # the decomposition must add back up to the station RMS
                tot_ss = sum(v[0] for v in band_rms(curves, exp, eta).values())
                tot_n = sum(v[1] for v in band_rms(curves, exp, eta).values())
                assert abs((tot_ss / max(tot_n, 1)) ** 0.5
                           - per[eta][0]) < 1e-12, (
                    f"band decomposition does not reconstruct the station RMS "
                    f"at eta={eta} -- it is measuring something else")
            band_out = {}
            print("  band decomposition (RMS over the 5 unmasked stations):")
            for name, _, _ in BANDS:
                for side in ("upper", "lower"):
                    ss, nn = acc[f"{name}_{side}"]
                    v = (ss / max(nn, 1)) ** 0.5 if nn else float("nan")
                    band_out[f"band_{name}_{side}"] = round(v, 6)
                    print(f"    {name:3s} {side:5s}: RMS {v:.4f} ({nn} pts)")
            all_ss = sum(v[0] for v in acc.values())
            all_n = sum(v[1] for v in acc.values())
            print(f"    all-point pooled RMS = "
                  f"{(all_ss/max(all_n,1))**0.5:.4f} ({all_n} pts)")
            pool = float(np.mean([per[e][0] for e in ETAS[:N_UNMASKED]]))
            pooled[(level, tag)] = pool
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], phi,
                                        alpha_deg=ALPHA, s_ref=s_ref,
                                        m_inf=M_INF)
            o = np.argsort(wc.station_z)
            clkj = cl_kj_3d(gamma[o], wc.station_z[o], s_ref, B_SEMI)
            print(f"  B1 pooled 5-station RMS = {pool:.4f}")
            print(f"  B4 cl_p {f['cl']:.6f}  cl_KJ {float(clkj):.6f}  "
                  f"(P14 anchors {P14_ANCHOR.get(level)})")
            rows.append(dict(
                level=level, entropy=tag, kutta=kutta,
                converged=bool(r.get("converged")),
                res_final=res, wall_s=round(wall, 1),
                pooled_rms_5=round(pool, 6),
                allpoint_rms_5=round((all_ss / max(all_n, 1)) ** 0.5, 6),
                **band_out,
                **{f"rms_eta{e:.2f}": round(per[e][0], 6) for e in ETAS},
                **{f"xshock_eta{e:.2f}": per[e][2] for e in ETAS},
                cl_p=round(f["cl"], 6), cl_kj=round(float(clkj), 6),
                p14_cl_p=P14_ANCHOR.get(level, (None, None))[0],
                p14_cl_kj=P14_ANCHOR.get(level, (None, None))[1],
                sigma_min=smin, n_shock_cells=nsh, m1_max=m1,
                transport_converged=sconv, n_floored=nflr,
                r2_validity=("n/a" if not entropy
                             else ("PASS" if r2_ok else "FAIL"))))

    with open(os.path.join(OUT, "m3_budget.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", os.path.join(OUT, "m3_budget.csv"))

    print("\n=== the registered reading (R1 -> R2 -> R3) ===")
    for level in levels:
        off = pooled.get((level, "OFF/probe"))
        on = pooled.get((level, "ON/probe"))
        if off is None or on is None:
            print(f"  {level}: incomplete pair -- recorded, not hidden")
            continue
        ref = GV53_K0_POOLED.get(level)
        if ref is None:
            print(f"  R1 {level}: no committed same-extractor RMS -- OFF leg "
                  f"{off:.4f} RECORDED")
        else:
            rel = off / ref - 1.0
            print(f"  R1 {level}: OFF {off:.4f} vs committed {ref}: "
                  f"{100*rel:+.2f} %  {'OK' if abs(rel) < 0.05 else 'DRIFT'}")
        r2 = next(r["r2_validity"] for r in rows
                  if r["level"] == level and r["entropy"] == "ON/probe")
        d = on - off
        print(f"  R2 {level}: {r2}")
        print(f"  {level} pooled RMS  OFF {off:.4f} -> ON {on:.4f}  "
              f"(delta {d:+.4f})")
        for k in ("OFF/pressure", "ON/pressure"):
            v = pooled.get((level, k))
            if v is not None:
                print(f"  ★ {level} {k}: pooled RMS {v:.4f}  "
                      f"(vs OFF/probe {off:.4f}: {v - off:+.4f})")
        if r2 != "PASS":
            print("  R3: NOT applied -- the ON leg failed the validity gate, so "
                  "this delta is the 3-D donor-cycle defect, not an entropy "
                  "effect. The round's deliverable is that diagnosis.")
        elif d <= -0.02:
            print("  R3: shock position IS a leading term => the 3-D "
                  "donor-cycle repair outranks prism layers")
        elif abs(d) < 0.01:
            print("  R3: shock position is NOT a leading term => S2 proceeds "
                  "with prisms as planned")
        else:
            print("  R3: in between => RECORDED, order unchanged, medium leg "
                  "decides")
    return 0


if __name__ == "__main__":
    sys.exit(main(tuple(sys.argv[1:]) or ("coarse",)))
