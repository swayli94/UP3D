"""R13 -- join the three committed m1_gate CSVs. Zero solves.

Binding text: phases/p5/docs/dev_phase_five/20260823-0500-r13-prereg.md (committed first).

Both axes were already measured and committed. What had never been done is the
comparison: putting the 16-thread and 8-thread runs side by side per (level, C, seed).

Run:  PYTHONNOUSERSITE=1 python bench/studies/r13_thread_seed_join/run.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
GR = os.path.join(ROOT, "bench/gate_results")

CANON = ("xcoarse", "coarse", "medium")      # G-SCOPE
TOL_X, TOL_CL = 0.0055, 0.03                 # M1's requirement, chord / relative
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}")


def load(name, threads):
    rows = list(csv.DictReader(open(os.path.join(GR, name + ".csv"))))
    dropped = [r for r in rows if r["level"] not in CANON]
    rows = [r for r in rows if r["level"] in CANON]
    conv = [r for r in rows
            if r["converged"] == "True"
            and int(float(r.get("n_limited") or 0)) == 0
            and int(float(r.get("n_floored") or 0)) == 0]
    print(f"  G-PROV  {name}.csv  threads={threads}  rows={len(rows)} "
          f"(dropped {len(dropped)} non-canonical) usable={len(conv)}")
    if dropped:
        print(f"          ★ G-SCOPE excluded levels: "
              f"{sorted({r['level'] for r in dropped})}")
    return {(r["level"], r["C"], r["n_picard_seed"]): r for r in rows}, \
           {(r["level"], r["C"], r["n_picard_seed"]) for r in conv}


def main():
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    assert "pyfp3d.solve.newton" not in sys.modules, "G-NOSOLVE"
    print("  G-NOSOLVE  no solver module imported; three committed CSVs only")
    a16, u16 = load("m1_gate_default", 16)
    a8, u8 = load("m1_gate_default_8t", 8)
    a3, u3 = load("m1_gate_default_3seed", 16)

    # ---- R-THREAD: the comparison nobody had made -------------------------
    keys = sorted(set(a16) & set(a8))
    if not keys:
        _record("R-THREAD", "alignment", "must align on (level, C, seed)",
                "no shared keys", "★ STOP -- kill criterion 5"); return 2
    both, split, rows = [], [], []
    for k in keys:
        r16, r8 = a16[k], a8[k]
        rec = {"level": k[0], "C": k[1], "seed": k[2],
               "usable_16t": k in u16, "usable_8t": k in u8,
               "x16": r16["x_shock"] or "", "x8": r8["x_shock"] or "",
               "cl16": r16["cl_p"] or "", "cl8": r8["cl_p"] or ""}
        if k in u16 and k in u8:
            dx = abs(float(r16["x_shock"]) - float(r8["x_shock"]))
            dcl = abs(float(r16["cl_p"]) - float(r8["cl_p"])) / abs(float(r8["cl_p"]))
            rec.update(dx=dx, dcl_rel=dcl,
                       over=bool(dx > TOL_X or dcl > TOL_CL))
            both.append(rec)
        elif (k in u16) != (k in u8):
            split.append(rec)
        rows.append(rec)
    with open(os.path.join(HERE, "results", "thread_join.csv"), "w", newline="") as f:
        ks = sorted({x for d in rows for x in d})
        w = csv.DictWriter(f, fieldnames=ks); w.writeheader(); w.writerows(rows)
    print(f"\n  {'level':8}{'C':>5}{'seed':>6}  {'16t usable':>11}{'8t usable':>11}"
          f"  {'|dx_shock|':>11}{'|dcl| rel':>11}")
    for r in rows:
        d = (f"{r['dx']:11.6f}{100*r['dcl_rel']:10.2f}%"
             if "dx" in r else f"{'--':>11}{'--':>11}")
        print(f"  {r['level']:8}{r['C']:>5}{r['seed']:>6}  {str(r['usable_16t']):>11}"
              f"{str(r['usable_8t']):>11}  {d}"
              + ("   ★ OVER" if r.get("over") else ""))
    over = [r for r in both if r["over"]]
    _record("R-THREAD", "per-cell 16t vs 8t, both-usable cells only",
            f"any |dx_shock| > {TOL_X} c or |dcl| > {TOL_CL:.0%} => the ANSWER is "
            "thread-dependent, not just convergence",
            f"{len(both)} both-usable cells, {len(over)} over tolerance"
            + (f" ({', '.join(r['level']+'/C'+r['C']+'/s'+r['seed'] for r in over)})"
               if over else "")
            + f";  {len(split)} cells where CONVERGENCE ITSELF is thread-dependent"
            + (f" ({', '.join(r['level']+'/C'+r['C']+'/s'+r['seed'] for r in split)})"
               if split else ""),
            "★ R-THREAD: the ANSWER is thread-dependent" if over else
            "R-THREAD: thread dependence shows up as convergence, NOT as which root")

    # ---- R-SEED: cross-seed spread, the script's own UNDEFINED rule -------
    cells, srows = {}, []
    for (lv, C, sd), r in a3.items():
        if (lv, C, sd) in u3:
            cells.setdefault((lv, C), {})[sd] = float(r["cl_p"])
    worst, worstk = -1.0, None
    for (lv, C), d in sorted(cells.items()):
        if len(d) < 2:
            srows.append({"level": lv, "C": C, "n_conv": len(d), "spread_rel": None})
            print(f"  cross-seed {lv:7} C={C:>4}: {len(d)} converged seed(s) "
                  f"-> UNDEFINED (the script's own rule)")
            continue
        v = list(d.values()); sp = (max(v) - min(v)) / abs(min(v))
        srows.append({"level": lv, "C": C, "n_conv": len(d),
                      "spread_rel": round(sp, 6)})
        if sp > worst:
            worst, worstk = sp, (lv, C)
        print(f"  cross-seed {lv:7} C={C:>4}: {len(d)} seeds {sorted(d)} "
              f"-> spread {100*sp:.2f} %" + ("  ★ > 3 %" if sp > 0.03 else ""))
    with open(os.path.join(HERE, "results", "seed_spread.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["level", "C", "n_conv", "spread_rel"])
        w.writeheader(); w.writerows(srows)
    n_und = sum(1 for r in srows if r["spread_rel"] is None)
    _record("R-SEED", "cross-seed cl_p spread from the committed 3-seed run",
            "coarse C=1.0 should reproduce ~2.76 %; any cell > 3 % means the seed "
            "alone spends the whole (b)(c) budget",
            f"worst {100*worst:.2f} % at {worstk};  {n_und}/{len(srows)} cells "
            f"UNDEFINED (<2 converged seeds)",
            f"R-SEED: worst {100*worst:.2f} %"
            + (" -- OVER the 3 % budget" if worst > 0.03 else " -- inside 3 %"))

    # ---- R-BLIND: what the gate PRINTS -----------------------------------
    src = open(os.path.join(ROOT, "bench/run_m1_gate.py")).read()
    prints_seed = "cross-seed cl_p spread" in src
    prints_thread = ("8t" in src and "16" in src and
                     any(s in src for s in ("thread_join", "vs 8", "8 vs 16")))
    _record("R-BLIND", "does the gate PRINT a cross-seed / cross-thread comparison",
            "cross-seed printed + cross-thread not => the capability boundary's "
            "'the gate never compares across seeds' is STALE and the real gap is threads",
            f"cross-seed printed: {prints_seed};  cross-thread compared: {prints_thread}",
            "★ R-BLIND: cross-seed IS printed, cross-thread is NOT -- my previous "
            "statement to the user was wrong and the boundary sentence is stale"
            if prints_seed and not prints_thread else
            "R-BLIND: neither is printed" if not prints_seed else "R-BLIND: both")

    with open(os.path.join(HERE, "results", "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
