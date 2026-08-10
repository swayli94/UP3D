"""GB30.3 -- the López dissipation lever (B30 PRE_REGISTRATION.md).

Zero library change: `upwind_c` and `upwind_c_post` are existing kwargs in
both solvers (never used by a production recipe -- capability review,
2026-07-14). The staging is driven MANUALLY (re-solve at the same Mach
with decreasing c), which is exactly what `upwind_c_post` does internally.

Protocol (pre-registered):
  1. L1 climb: re-run GB30.2's first failing level at upwind_c=2.0 (the
     López M6 climb value; UP3D default 1.5), same seed, rest of the recipe
     verbatim. If it converges (production semantics), keep climbing at
     dm=0.0125 (<= 0.84).
  2. L1 staging: at the new top converged level re-solve strictly at
     c = 1.8, then 1.6 (the López Table 4.13 2.0 -> 1.6 sharpening).
  3. Cost: cl_p / gamma / Mmax drift vs the committed c=1.5 anchors
     (RECORDED); corridor corrM guardrail (<= 1.3).

Verdicts (pre-registered):
  B30-L1-PASS  ceiling +>= 1 rung AND the 1.6 staging survives
  B30-L1-PART  climbs only at elevated dissipation (staging fails) ->
               cost table to the user
  B30-L1-FAIL  ceiling unmoved -> dissipation is not the constraint; the
               C-class tip cure becomes the named next candidate

Legs follow the GB30.2 tree: SAME -> both legs; SPLIT -> LS only.
Run:  python cases/analysis/b30_transonic_ceiling/run_g3.py [LS CONF]
Artifacts: results/g3_levers.csv, results/g3_levers.png
"""

import csv
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb30 import (ALPHA, DEMO_OUT, LS_MESH_DIR, OUT, build_conf,
                  build_ls_flat, corridor_mmax, load_mesh, mach2_ls,
                  cl_p_ls, cl_p_conf, run_conf_level, run_ls_level)
from pyfp3d.solve.newton import NewtonWorkspace
from wb30 import mach2_conf
from pyfp3d.constraints.wake import kutta_targets

C_CLIMB = 2.0          # the López M6 climb value (Table 4.13)
C_STAGES = [1.8, 1.6]  # the López post-ramp sharpening
DM_CLIMB = 0.0125
M_TOP = 0.84
SEED = {"LS": ("ls_flat_medium_084.npz", "phi_ext"),
        "CONF": ("conf_medium_079.npz", "phi")}
# ^ demo FALLBACK seeds, used only when the GB30.2 chain has no converged
# state below the first failing level (LS case: 0.7875 dies at once).
ANCHOR = {"LS": dict(m_last=0.775, cl_p=0.2475),
          "CONF": dict(m_last=0.79, cl_p=0.2579)}


def g2_first_fail(leg):
    """GB30.2's census Mach for this leg (fallback = the B29 anchor)."""
    try:
        with open(OUT / "g2_census.csv") as fh:
            for row in csv.DictReader(fh):
                if row["leg"] == leg:
                    return float(row["m_census"]), row["converged"] == "True"
    except FileNotFoundError:
        pass
    return {"LS": 0.7875, "CONF": 0.80}[leg], False


def g2_best_converged(leg, m_fail):
    """Best GB30.2 seed for this leg: the HIGHEST converged chain state
    below m_fail. The G2 strict warm-start chain can converge past the
    demo anchor (T2, measured: CONF 0.80 -> 0.82) -- the lever climb must
    continue the chain, not restart from the older demo seed.
    Returns (m_inf, phi, gamma) or None (-> demo fallback)."""
    best = None
    for f in sorted(OUT.glob(f"g2_{leg}_m*.npz")):
        d = np.load(f, allow_pickle=True)
        s = json.loads(str(d["summary"]))
        if not s.get("converged") or s["m_inf"] >= m_fail - 1e-12:
            continue
        if best is None or s["m_inf"] > best[0]:
            best = (s["m_inf"], d["phi"], d["gamma"])
    return best


def leg_ls(mesh):
    _wls, _cm, mvop = build_ls_flat(mesh)
    m_fail, reached = g2_first_fail("LS")
    best = g2_best_converged("ls", m_fail)
    if best is not None:
        phi0 = best[1]
        print(f"  [LS] seed = GB30.2 chain top m={best[0]}", flush=True)
    else:
        d = np.load(DEMO_OUT / SEED["LS"][0])
        phi0 = d[SEED["LS"][1]]
        print(f"  [LS] seed = demo {SEED['LS'][0]} (no converged G2 "
              f"state below {m_fail})", flush=True)
    if reached:
        print("  [LS] GB30.2 reached 0.84 at c=1.5 -- lever moot", flush=True)
        return [], []
    # ---- climb at c=2.0
    climb, m, phi = [], m_fail, phi0
    while True:
        rec = run_ls_level(mesh, mvop, m, phi, upwind_c=C_CLIMB,
                           tag=f"g3c{C_CLIMB}")
        climb.append(rec)
        if not rec["converged"] or m >= M_TOP - 1e-12:
            break
        phi = rec["phi"]
        m = round(m + DM_CLIMB, 4)
    top = [r for r in climb if r["converged"]]
    if not top:
        return climb, []        # the lever did not move the ceiling
    m_top = top[-1]["m_inf"]
    phi_top = top[-1]["phi"]
    # ---- staging 2.0 -> 1.8 -> 1.6 at m_top
    stages, phi = [], phi_top
    for c in C_STAGES:
        rec = run_ls_level(mesh, mvop, m_top, phi, upwind_c=c,
                           tag=f"g3c{c}")
        stages.append(rec)
        if not rec["converged"]:
            break
        phi = rec["phi"]
    return climb, stages


def leg_conf():
    mc, wc = build_conf("medium")
    m_fail, reached = g2_first_fail("CONF")
    best = g2_best_converged("conf", m_fail)
    if best is not None:
        phi0, gamma0 = best[1], best[2]
        print(f"  [CONF] seed = GB30.2 chain top m={best[0]}", flush=True)
    else:
        d = np.load(DEMO_OUT / SEED["CONF"][0])
        phi0 = d[SEED["CONF"][1]]
        gamma0 = kutta_targets(phi0, wc)  # T6: zeros-Gamma diverges (G2)
        print(f"  [CONF] seed = demo {SEED['CONF'][0]} + reconstructed "
              f"Gamma (no converged G2 state below {m_fail})", flush=True)
    if reached:
        print("  [CONF] GB30.2 reached 0.84 at c=1.5 -- lever moot",
              flush=True)
        return [], []
    climb, m, phi, gamma = [], m_fail, phi0, gamma0
    while True:
        rec = run_conf_level(mc, wc, m, phi, gamma_init=gamma,
                             upwind_c=C_CLIMB, tag=f"g3c{C_CLIMB}")
        climb.append(rec)
        if not rec["converged"] or m >= M_TOP - 1e-12:
            break
        phi, gamma = rec["phi"], rec["gamma"]
        m = round(m + DM_CLIMB, 4)
    top = [r for r in climb if r["converged"]]
    if not top:
        return climb, []
    m_top = top[-1]["m_inf"]
    phi_top, gamma_top = top[-1]["phi"], top[-1]["gamma"]
    stages, phi, gamma = [], phi_top, gamma_top
    for c in C_STAGES:
        rec = run_conf_level(mc, wc, m_top, phi, gamma_init=gamma,
                             upwind_c=c, tag=f"g3c{c}")
        stages.append(rec)
        if not rec["converged"]:
            break
        phi, gamma = rec["phi"], rec["gamma"]
    return climb, stages


def measure(mesh, mvop, leg, rec, mc=None, wc=None):
    """cl_p / mmax / corridor on a returned state (cost + guardrail)."""
    phi = rec["phi"]
    m_inf = rec["m_inf"]
    if leg == "LS":
        m2 = mach2_ls(mvop, phi, m_inf)
        cents = mesh.nodes[mesh.elements].mean(axis=1)
        cl = cl_p_ls(mesh, mvop, phi, m_inf)
    else:
        ws = NewtonWorkspace(mc, wc, ALPHA, 1.0, 1.4, (0.25, 0.0), True,
                             kutta_estimator="pressure")
        ws.set_mach(m_inf)
        m2 = mach2_conf(ws, phi, rec["gamma"], m_inf)
        cents = mc.nodes[mc.elements].mean(axis=1)
        cl = cl_p_conf(mc, phi, m_inf)
    corr_m, _xyz = corridor_mmax(m2, cents)
    return dict(cl_p=cl, mmax=float(np.sqrt(m2.max())), corr_mmax=corr_m)


def main():
    only = sys.argv[1:] or ["LS", "CONF"]
    rows = []
    if "LS" in only:
        mesh = load_mesh(LS_MESH_DIR / "medium.msh")
        print(f"=== GB30.3 LS leg: climb at c={C_CLIMB} ===", flush=True)
        climb, stages = leg_ls(mesh)
        _wls, _cm, mvop = build_ls_flat(mesh)
        rows += summarize("LS", climb, stages,
                          lambda r: measure(mesh, mvop, "LS", r))
    if "CONF" in only:
        print(f"=== GB30.3 CONF leg: climb at c={C_CLIMB} ===", flush=True)
        climb, stages = leg_conf()
        mc, wc = build_conf("medium")
        rows += summarize("CONF", climb, stages,
                          lambda r: measure(None, None, "CONF", r,
                                            mc=mc, wc=wc))
    _write(rows)
    _fig(rows)
    _verdict_hint(rows)


def summarize(leg, climb, stages, meas):
    rows = []
    for rec, phase in ([(r, "climb") for r in climb]
                       + [(r, "stage") for r in stages]):
        row = dict(leg=leg, phase=phase, m=rec["m_inf"],
                   upwind_c=rec["upwind_c"], converged=rec["converged"],
                   accept_reason=rec["accept_reason"],
                   res=rec["residual_norm"],
                   n_limited=rec["n_limited"], n_floored=rec["n_floored"],
                   n_newton=rec["n_newton"], wall_s=round(rec["wall_s"], 1))
        if rec["converged"]:
            row.update(meas(rec))
        rows.append(row)
        print(f"  [{leg} {phase} c={rec['upwind_c']} m={rec['m_inf']}] "
              f"conv={rec['converged']} res={rec['residual_norm']:.1e} "
              f"clamps={rec['n_limited']}+{rec['n_floored']}", flush=True)
    return rows


def _write(rows):
    keys = ["leg", "phase", "m", "upwind_c", "converged", "accept_reason",
            "res", "n_limited", "n_floored", "n_newton", "wall_s",
            "cl_p", "mmax", "corr_mmax"]
    with open(OUT / "g3_levers.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g3_levers.csv'}")


def _fig(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for leg, ax in zip(("LS", "CONF"), axes):
        rr = [r for r in rows if r["leg"] == leg and r["converged"]]
        if not rr:
            ax.set_title(f"{leg}: no converged lever level")
            continue
        labels = [f"c={r['upwind_c']}\nM{r['m']}" for r in rr]
        cl = [r.get("cl_p", np.nan) for r in rr]
        mm = [r.get("mmax", np.nan) for r in rr]
        x = np.arange(len(rr))
        ax.bar(x - 0.18, cl, width=0.36, label="cl_p")
        ax.bar(x + 0.18, mm, width=0.36, label="M_max")
        anch = ANCHOR[leg]
        ax.axhline(anch["cl_p"], color="tab:blue", ls=":", lw=1.0,
                   label=f"c=1.5 anchor cl_p {anch['cl_p']}")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(f"{leg}: converged lever levels (anchor m_last "
                     f"{anch['m_last']} @ c=1.5)")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8)
    fig.suptitle("GB30.3: the López dissipation lever -- converged-state "
                 "cost readout")
    fig.tight_layout()
    fig.savefig(OUT / "g3_levers.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT / 'g3_levers.png'}")


def _verdict_hint(rows):
    print("\n--- verdict hint (pre-reg GB30.3) ---")
    for leg in ("LS", "CONF"):
        rr = [r for r in rows if r["leg"] == leg]
        if not rr:
            continue
        anch = ANCHOR[leg]
        conv = [r for r in rr if r["converged"]]
        climbed = [r for r in conv if r["phase"] == "climb"
                   and r["m"] > anch["m_last"] + 1e-9]
        staged16 = [r for r in conv if r["phase"] == "stage"
                    and abs(r["upwind_c"] - 1.6) < 1e-9]
        if climbed and staged16:
            print(f"  {leg}: ceiling {anch['m_last']} -> "
                  f"{climbed[-1]['m']} at c=2.0 AND 1.6 staging survives "
                  f"-> B30-L1-PASS candidate (refresh the demo anchors; "
                  f"cl_p drift {climbed[-1].get('cl_p')} vs "
                  f"{anch['cl_p']} RECORDED)")
        elif climbed:
            print(f"  {leg}: climbs at c=2.0 but staging incomplete -> "
                  f"B30-L1-PART candidate (elevated dissipation only; "
                  f"cost table above, user's call)")
        else:
            print(f"  {leg}: ceiling unmoved at c=2.0 -> B30-L1-FAIL "
                  f"candidate (dissipation is not the constraint; C-class "
                  f"tip cure becomes the named next candidate)")


if __name__ == "__main__":
    main()
