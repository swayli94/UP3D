"""Per-gate classification of the migrated-machine reds from the four committed arms.

usage: python classify.py            (reads ./arms/, ./compare_values.csv, the round-10 console log)
writes: classification.csv, premise_b.csv, numba_cache.csv  (all numbers quoted in the round-11 verdict)

Registered BEFORE the run (round-10 audit §3.4 / §12): the four arms and the four pairings
  ISA cold = native_cold vs generic_cold   ISA warm = native_warm vs generic_warm
  cache native = native_cold vs native_warm   cache generic = generic_cold vs generic_warm
NOT registered before the run (post hoc, descriptive only): the per-gate summary statistics and
the labels below. They describe; they do not decide a disposition (that is a user ruling).
Scope: NUMBA_CPU_NAME=generic changes numba/LLVM JIT codegen only, not OpenBLAS/SciPy/PyAMG.
"""
import csv
import json
import os
import re
import statistics as st
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = ["native_cold", "native_warm", "generic_cold", "generic_warm"]
LOG = os.path.join(HERE, "..", "devflow_adoption_audit", "targeted_12red_20261104_console.txt")


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def rel(a, b):
    a, b = f(a), f(b)
    if a is None or b is None:
        return None
    return abs(a - b) / max(abs(b), 1e-30)


def failed(arm):
    root = ET.parse(os.path.join(HERE, "arms", arm, "junit.xml")).getroot()
    out = set()
    for tc in root.iter("testcase"):
        if tc.find("failure") is not None or tc.find("error") is not None:
            out.add(f"{tc.get('classname')}::{tc.get('name')}")
    return out


def label(gap_max, isa_max, cache_max):
    if gap_max <= 1e-6:
        return "no gap above 1e-6"
    if cache_max > 1e-6:
        return "cache-state sensitive (cold vs warm moves it) + codegen"
    if isa_max >= 0.1 * gap_max:
        return "codegen-sensitive at the gap's order (ISA max >= 0.1 x gap max)"
    return "NOT explained by numba JIT codegen (ISA max < 0.1 x gap max)"


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "compare_values.csv"), encoding="utf-8")))
    fails = {a: failed(a) for a in ARMS}
    out = []
    for g in sorted({r["gate_dir"] for r in rows}):
        gr = [r for r in rows if r["gate_dir"] == g and f(r["committed"]) is not None]
        gap = [x for x in (rel(r["native_warm"], r["committed"]) for r in gr) if x is not None]
        isa = [x for x in (rel(r["generic_warm"], r["native_warm"]) for r in gr) if x is not None]
        cac = [x for x in (rel(r["native_cold"], r["native_warm"]) for r in gr) if x is not None]
        tag = g.split("/")[-1].split("_")[0]
        out.append({
            "gate": g.split("/")[-1], "n_values": len(gap),
            "n_gap_gt_1e-6": sum(x > 1e-6 for x in gap),
            "gap_median": f"{st.median(gap):.3e}", "gap_max": f"{max(gap):.3e}",
            "isa_warm_median": f"{st.median(isa):.3e}", "isa_warm_max": f"{max(isa):.3e}",
            "cache_native_median": f"{st.median(cac):.3e}", "cache_native_max": f"{max(cac):.3e}",
            "failed_in": " ".join(a for a in ARMS if any(f".test_{tag}_" in n for n in fails[a])) or "-",
            "label_post_hoc": label(max(gap), max(isa), max(cac)),
        })
    # D13 has its own comparison code: summarise from the npz dumps
    import numpy as np
    asym = {}
    for a in ARMS:
        d = dict(np.load(os.path.join(HERE, "arms", a, "d13_runs.npz")))
        asym[a] = float(np.max(d["coarse/upper/ds"] / d["coarse/lower/ds"]))
    out.append({
        "gate": "D13_ibl_vs_xfoil (npz)", "n_values": "", "n_gap_gt_1e-6": "",
        "gap_median": "", "gap_max": "",
        "isa_warm_median": "", "isa_warm_max": f"{abs(asym['generic_warm'] - asym['native_warm']):.3e} (abs, coarse max asym)",
        "cache_native_median": "", "cache_native_max": f"{abs(asym['native_cold'] - asym['native_warm']):.3e} (abs, coarse max asym)",
        "failed_in": " ".join(a for a in ARMS if any("test_D13_" in n for n in fails[a])),
        "label_post_hoc": "coarse max asym per arm " + ", ".join(f"{a}={v!r}" for a, v in asym.items())
                          + "; limit 1.15, historical 0.90-1.12",
    })
    with open(os.path.join(HERE, "classification.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    # premise (b) of ruling A: two warm runs, same codegen, different cache dirs
    log = open(LOG, encoding="utf-8").read()
    pat = re.compile(r"E\s+AssertionError: (.+?) 的 (\w+)[：:]本次算出 ([-0-9.e+]+)")
    pb = []
    for key, col, val in pat.findall(log):
        k = key.strip().replace("'", "")
        cand = [r for r in rows if r["col"] == col and r["row_key"].replace("'", "") == k]
        nw = cand[0]["native_warm"] if len(cand) == 1 else ""
        pb.append({"row_key": key.strip(), "col": col, "targeted_run_9e": val,
                   "native_warm_repr": nw,
                   "equal_at_9e": (f"{float(nw):.9e}" == f"{float(val):.9e}") if nw else "unmatched"})
    m = re.search(r"assert ([0-9.]+) <= 1\.15", log)
    pb.append({"row_key": "D13 coarse max asym", "col": "asym", "targeted_run_9e": m.group(1),
               "native_warm_repr": repr(asym["native_warm"]),
               "equal_at_9e": f"repr-equal={m.group(1) == repr(asym['native_warm'])}"})
    with open(os.path.join(HERE, "premise_b.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pb[0]))
        w.writeheader()
        w.writerows(pb)

    nc = []
    for a in ARMS:
        s = json.load(open(os.path.join(HERE, "arms", a, "numba_stats.json"), encoding="utf-8"))
        per = s["per_dispatcher"]
        nc.append({"arm": a, "dispatchers": len(per),
                   "dispatchers_used": sum(1 for v in per.values() if v["n_signatures"]),
                   "hits": s["total"]["hits"], "misses": s["total"]["misses"],
                   "cpu": json.load(open(os.path.join(HERE, "arms", a, "codegen.json")))["effective_cpu_name"]})
    with open(os.path.join(HERE, "numba_cache.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(nc[0]))
        w.writeheader()
        w.writerows(nc)
    for name in ("classification.csv", "premise_b.csv", "numba_cache.csv"):
        print(f"== {name}")
        print(open(os.path.join(HERE, name), encoding="utf-8").read())


if __name__ == "__main__":
    main()
