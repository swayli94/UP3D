"""Offline comparison of the four baseline-provenance arms. Reads only; safe to re-run.

usage: python compare.py [OUT_DIR]   (default: this file's directory)

Writes OUT_DIR/compare_values.csv (one row per locked value, every arm side by side)
and OUT_DIR/compare_summary.md. Pairings are the ones registered BEFORE the run:
  ISA, cold:  native_cold  vs generic_cold
  ISA, warm:  native_warm  vs generic_warm     <- same state as the 12 reds
  cache, native:  native_cold vs native_warm
  cache, generic: generic_cold vs generic_warm
and every arm vs the committed summary.csv values.
Scope: NUMBA_CPU_NAME controls only numba/LLVM JIT codegen, not OpenBLAS/SciPy/PyAMG.
"""
import csv
import json
import math
import os
import sys
import xml.etree.ElementTree as ET

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
ARMS = ["native_cold", "native_warm", "generic_cold", "generic_warm"]
PAIRS = [("ISA cold", "native_cold", "generic_cold"),
         ("ISA warm", "native_warm", "generic_warm"),
         ("cache native", "native_cold", "native_warm"),
         ("cache generic", "generic_cold", "generic_warm")]


def fnum(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def load_arm(arm):
    path = os.path.join(OUT, arm, "fresh.jsonl")
    vals, meta = {}, {}
    if not os.path.exists(path):
        return None, None
    count = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            rec = json.loads(line)
            base = (rec["nodeid"], rec["fn"], rec["gate_dir"])
            count[base] = count.get(base, 0) + 1
            call = count[base]
            for k, row in rec["fresh"].items():
                for col in rec["cols"]:
                    if col not in row:
                        continue
                    key = (*base, call, k, col)
                    vals[key] = (row[col], (rec["committed"].get(k) or {}).get(col))
                    meta[key] = {"spread_tol": (rec.get("spread_tol") or {}).get(col),
                                 "rel_tol": rec.get("rel_tol")}
    return vals, meta


def junit(arm):
    path = os.path.join(OUT, arm, "junit.xml")
    if not os.path.exists(path):
        return None
    root = ET.parse(path).getroot()
    failed, counts = [], {"tests": 0, "failures": 0, "errors": 0, "skipped": 0}
    suites = [root] if root.tag == "testsuite" else root.findall("testsuite")
    for s in suites:
        for c in counts:
            counts[c] += int(s.get(c, 0))
        for tc in s.iter("testcase"):
            if tc.find("failure") is not None or tc.find("error") is not None:
                failed.append(f"{tc.get('classname')}::{tc.get('name')}")
    return counts, sorted(failed)


def rel(a, b):
    if a is None or b is None:
        return None
    return abs(a - b) / max(abs(b), 1e-30)


def main():
    data = {a: load_arm(a) for a in ARMS}
    present = [a for a in ARMS if data[a][0] is not None]
    keys = sorted(set().union(*[set(data[a][0]) for a in present])) if present else []
    rows = []
    for key in keys:
        r = {"nodeid": key[0], "fn": key[1], "gate_dir": key[2], "call": key[3],
             "row_key": key[4], "col": key[5]}
        committed = None
        for a in present:
            v = data[a][0].get(key)
            r[a] = v[0] if v else ""
            if v and v[1] not in (None, ""):
                committed = v[1]
            m = data[a][1].get(key) if v else None
            if m:
                r["spread_tol"] = m["spread_tol"]
                r["rel_tol"] = m["rel_tol"]
        r["committed"] = committed or ""
        for a in present:
            r[f"rel_{a}_vs_committed"] = rel(fnum(r[a]), fnum(committed))
        for name, x, y in PAIRS:
            if x in present and y in present:
                r[f"bitsame[{name}]"] = (r[x] == r[y]) if r[x] and r[y] else ""
                r[f"rel[{name}]"] = rel(fnum(r[y]), fnum(r[x]))
        rows.append(r)

    if rows:
        cols = list(dict.fromkeys(c for r in rows for c in r))
        with open(os.path.join(OUT, "compare_values.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: ("" if r.get(c) is None else
                                (f"{r[c]:.6e}" if isinstance(r.get(c), float) else r.get(c)))
                            for c in cols})

    L = ["# baseline-provenance arms: comparison", "",
         f"source: `{OUT}`; values compared at full precision (`repr`).",
         "Scope: `NUMBA_CPU_NAME=generic` changes numba/LLVM JIT codegen only.", ""]
    L += ["## Arms", "", "| arm | tests | failures | errors | skipped | numba hits | numba misses | dispatchers | wall |",
          "|---|---|---|---|---|---|---|---|---|"]
    progress = {}
    pp = os.path.join(OUT, "progress.txt")
    if os.path.exists(pp):
        for line in open(pp, encoding="utf-8"):
            parts = dict(p.split("=", 1) for p in line.split()[1:] if "=" in p)
            progress[line.split()[0]] = parts
    fails = {}
    for a in ARMS:
        j = junit(a)
        st = {}
        sp = os.path.join(OUT, a, "numba_stats.json")
        if os.path.exists(sp):
            st = json.load(open(sp, encoding="utf-8"))["total"]
        if j is None:
            L.append(f"| {a} | (no junit) | | | | | | | |")
            continue
        c, f = j
        fails[a] = f
        L.append(f"| {a} | {c['tests']} | {c['failures']} | {c['errors']} | {c['skipped']} | "
                 f"{st.get('hits', '')} | {st.get('misses', '')} | {st.get('dispatchers', '')} | "
                 f"{progress.get(a, {}).get('wall', '')} s |")
    L += ["", "## Failed nodes per arm", ""]
    allf = sorted(set().union(*fails.values())) if fails else []
    if allf:
        L += ["| node | " + " | ".join(a for a in ARMS if a in fails) + " |",
              "|---|" + "---|" * len(fails)]
        for n in allf:
            L.append(f"| `{n}` | " + " | ".join("FAIL" if n in fails[a] else "pass"
                                                 for a in ARMS if a in fails) + " |")
    L += ["", "## Pairwise (all locked values)", "",
          "| pairing | values | bit-identical | max rel diff |", "|---|---|---|---|"]
    for name, x, y in PAIRS:
        b = [r.get(f"bitsame[{name}]") for r in rows if r.get(f"bitsame[{name}]") not in ("", None)]
        d = [r.get(f"rel[{name}]") for r in rows if r.get(f"rel[{name}]") is not None]
        if b:
            L.append(f"| {name} ({x} vs {y}) | {len(b)} | {sum(1 for v in b if v)} | "
                     f"{max(d) if d else float('nan'):.3e} |")
    L += ["", "## Per gate: max rel diff vs committed, per arm", "",
          "| gate | values | " + " | ".join(present) + " |", "|---|---|" + "---|" * len(present)]
    gates = sorted({r["gate_dir"] for r in rows})
    for g in gates:
        gr = [r for r in rows if r["gate_dir"] == g]
        cells = []
        for a in present:
            d = [r[f"rel_{a}_vs_committed"] for r in gr if r.get(f"rel_{a}_vs_committed") is not None]
            cells.append(f"{max(d):.3e}" if d else "")
        L.append(f"| {g} | {len(gr)} | " + " | ".join(cells) + " |")
    L += ["", "## Per gate: pairwise bit-identical / max rel diff", "",
          "| gate | " + " | ".join(p[0] for p in PAIRS) + " |", "|---|" + "---|" * len(PAIRS)]
    for g in gates:
        gr = [r for r in rows if r["gate_dir"] == g]
        cells = []
        for name, _, _ in PAIRS:
            b = [r.get(f"bitsame[{name}]") for r in gr if r.get(f"bitsame[{name}]") not in ("", None)]
            d = [r.get(f"rel[{name}]") for r in gr if r.get(f"rel[{name}]") is not None]
            cells.append(f"{sum(1 for v in b if v)}/{len(b)}, {max(d):.2e}" if b and d else "")
        L.append(f"| {g} | " + " | ".join(cells) + " |")

    # D13 runs fixture (own comparison code in the gate)
    try:
        import numpy as np
        arrs = {a: dict(np.load(os.path.join(OUT, a, "d13_runs.npz")))
                for a in ARMS if os.path.exists(os.path.join(OUT, a, "d13_runs.npz"))}
        if arrs:
            L += ["", "## D13 `runs` fixture (npz), pairwise", "",
                  "| pairing | arrays | bit-identical | max abs diff |", "|---|---|---|---|"]
            for name, x, y in PAIRS:
                if x in arrs and y in arrs:
                    ks = sorted(set(arrs[x]) & set(arrs[y]))
                    same = sum(1 for k in ks if np.array_equal(arrs[x][k], arrs[y][k]))
                    mx = max((float(np.max(np.abs(arrs[x][k].astype(float) - arrs[y][k].astype(float))))
                              for k in ks if arrs[x][k].shape == arrs[y][k].shape and arrs[x][k].size), default=float("nan"))
                    L.append(f"| {name} | {len(ks)} | {same} | {mx:.3e} |")
            for a, d in arrs.items():
                asym = []
                for lv in ("coarse", "medium"):
                    u, lo = d.get(f"{lv}/upper/ds"), d.get(f"{lv}/lower/ds")
                    if u is not None and lo is not None:
                        asym.append(f"{lv} max asym {float(np.max(u / lo))!r}")
                L.append(f"- {a}: " + "; ".join(asym))
    except Exception as e:  # noqa: BLE001
        L.append(f"\n(D13 npz comparison failed: {e!r})")

    with open(os.path.join(OUT, "compare_summary.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
