"""The capability matrix's cross-cell deliverable: envelope tables + force-coefficient plots.

Reads only committed CSVs (bench/gate_results/capability_matrix.csv and
capability_alpha.csv) -- it never re-solves, so it is cheap to re-run after any row is
re-measured.

Three things the per-cell plots cannot show:
  1. cl_p(M) for every cell on shared axes, grouped by geometry, with the two wake paths
     overlaid -- the cross-path agreement (or disagreement) is the point.
  2. the mesh-refinement axis: the envelope FALLS with refinement on several cells, which
     is a property of the criteria binding earlier, not of the solver degrading, and it is
     only visible with the levels side by side.
  3. the stop-reason split. "Envelope M0.75" means two completely different things
     depending on whether convergence or the M_max < 1.4 criterion bound first, and a
     single number hides that.

Outputs (TRACKED): bench/gate_results/capability/summary_cl.png
                   bench/gate_results/capability/summary_envelope.png
                   bench/gate_results/capability_summary.md
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "gate_results")
ART = os.path.join(OUT, "capability")
MMAX_LIMIT = 1.4

GEOMS = [("naca2.5d", "NACA0012 2.5-D"), ("m6wing", "ONERA M6 wing"),
         ("wingbody", "ONERA M6 wing-body")]
LEVELS = ["coarse", "medium", "fine"]
PATHS = [("conforming", "-", "o"), ("level-set", "--", "s")]
LEVEL_COLOR = {"coarse": "tab:blue", "medium": "tab:red", "fine": "tab:green"}


def load(name):
    p = os.path.join(OUT, name)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def fnum(r, k):
    v = r.get(k)
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def main():
    rows = load("capability_matrix.csv")
    arows = load("capability_alpha.csv")
    if not rows:
        raise SystemExit("no matrix CSV")

    # ---------------------------------------------------------- cl_p(M) plots
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0))
    for ax, (geom, title) in zip(axes, GEOMS):
        for path, ls, mk in PATHS:
            for level in LEVELS:
                sel = [r for r in rows if r["geom"] == geom
                       and r["path"] == path and r["level"] == level
                       and r["status"].startswith("CLEAN")]
                if not sel:
                    continue
                sel.sort(key=lambda r: fnum(r, "m_inf"))
                m = [fnum(r, "m_inf") for r in sel]
                cl = [fnum(r, "cl_p") for r in sel]
                ax.plot(m, cl, ls, marker=mk, ms=4.5,
                        color=LEVEL_COLOR[level], lw=1.4,
                        label=f"{path} {level}")
                #: ★ mark the points that converged cleanly but exceed M_max < 1.4.
                #: They are real solutions -- excluded by the criterion, not by the
                #: solver -- so hiding them would misreport the capability.
                over = [(fnum(r, "m_inf"), fnum(r, "cl_p")) for r in sel
                        if r["status"] == "CLEAN_OVER_MMAX"]
                if over:
                    ax.plot([o[0] for o in over], [o[1] for o in over], "x",
                            color=LEVEL_COLOR[level], ms=10, mew=2)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"$M_\infty$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.5, loc="upper left")
    axes[0].set_ylabel(r"$c_{l,p}$  (surface-pressure integral)")
    fig.suptitle(r"Capability matrix: $c_{l,p}(M_\infty)$ — "
                 r"$\times$ = converged but $M_{max}\geq1.4$", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(ART, "summary_cl.png"), dpi=130)
    plt.close(fig)

    # -------------------------------------------------- envelope bar chart
    env = {}
    for r in rows:
        c = r["cell"]
        if r["status"] == "CLEAN":
            env.setdefault(c, {"m": None, "stop": None})
            env[c]["m"] = fnum(r, "m_inf")
        elif c in env or r["status"] != "CLEAN":
            env.setdefault(c, {"m": None, "stop": None})
            if env[c]["stop"] is None:
                env[c]["stop"] = r["status"]
    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    cells = sorted(env, key=lambda c: (c.split("_")[0], c))
    xs = range(len(cells))
    colors = {"NOT_CONVERGED": "tab:red", "CLEAN_OVER_MMAX": "tab:orange",
              None: "tab:green"}
    for x, c in zip(xs, cells):
        m, stop = env[c]["m"], env[c]["stop"]
        ax.bar(x, m if m else 0.0, color=colors.get(stop, "tab:gray"),
               edgecolor="k", lw=0.5)
        if m:
            ax.text(x, m + 0.005, f"{m:.2f}", ha="center", fontsize=8)
        else:
            ax.text(x, 0.01, "none", ha="center", fontsize=8, rotation=90)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(cells, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(r"envelope $M_\infty$ (clean, $M_{max}<1.4$)")
    ax.set_title("green = ladder exhausted · orange = stopped by $M_{max}$ "
                 "criterion · red = stopped by convergence")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(ART, "summary_envelope.png"), dpi=130)
    plt.close(fig)

    # ------------------------------------------------------------ markdown
    L = ["# Capability matrix — summary (generated by bench/make_capability_summary.py)",
         "",
         "Generated from the committed CSVs; no solve. `envelope` = highest ladder Mach "
         "that converged cleanly with `M_max < 1.4`.", "",
         "## Envelope by cell", "",
         "| cell | nodes | envelope M | stopped by | M_max there | cl_p | cl_KJ |",
         "|---|---|---|---|---|---|---|"]
    for c in cells:
        sel = [r for r in rows if r["cell"] == c]
        cl = [r for r in sel if r["status"] == "CLEAN"]
        last = cl[-1] if cl else None
        stop = env[c]["stop"] or "ladder exhausted"
        L.append(f"| {c} | {sel[0]['n_nodes']} | "
                 f"{('M%.2f' % fnum(last,'m_inf')) if last else '**none**'} | {stop} | "
                 f"{last['m_max'] if last else '-'} | "
                 f"{last['cl_p'] if last else '-'} | "
                 f"{last['cl_kj'] if last else '-'} |")
    if arows:
        L += ["", "## Alpha-relaxed points (M_max < 1.4 kept; only alpha lowered)", "",
              "| cell | M | alpha | status | M_max | cl_p |", "|---|---|---|---|---|---|"]
        for r in arows:
            L.append(f"| {r['cell']} | {r['m_inf']} | {r['alpha']} | {r['status']} | "
                     f"{r['m_max']} | {r['cl_p']} |")
    L += ["", "## How to read this", "",
          "- **Two different stop reasons.** `CLEAN_OVER_MMAX` rows are genuine "
          "converged solutions excluded by the `M_max < 1.4` criterion, not solver "
          "failures. `NOT_CONVERGED` rows are solver failures. An envelope number "
          "alone conflates them.",
          "- **The envelope falls with refinement** on several cells. Refinement "
          "resolves a sharper local peak, so whichever criterion binds first binds "
          "earlier -- this is not the solver degrading.",
          "- **Ladder resolution is 0.02-0.03 near the boundary** and there is no rung "
          "below M0.50, so each envelope is located only to that precision.",
          "- `descent10` separates a genuine stall from an exhausted iteration budget; "
          "`m_attained` guards against a row claiming a Mach the solver never reached.",
          "- **`cl_KJ` is blank on the level-set rows and on every NACA row, by "
          "construction, not by omission.** The level-set path carries no per-station "
          "circulation DOF, so the spanwise Kutta-Joukowski integral has no input "
          "there; and the 2.5-D NACA cases have no span to integrate over. "
          "Cross-path lift comparison therefore uses `cl_p`, which is how B9/B27 did "
          "their cross-model comparisons too.",
          "- **`conf_wing_coarse` exhausted the ladder at M0.84 with M_max 1.3976** -- "
          "just under the limit. Read as 'the ladder ran out', not as headroom: one "
          "more rung would very likely have tripped the criterion.",
          ""]
    with open(os.path.join(OUT, "capability_summary.md"), "w") as fh:
        fh.write("\n".join(L))
    print("wrote summary_cl.png, summary_envelope.png, capability_summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
