"""G27 -- B27 consistency diff: the refreshed B18 demo vs the committed anchors
(B27 PRE_REGISTRATION.md, 2026-07-20; gates GB27.1/GB27.2).

Three comparisons, all item-by-item into results/g27_consistency.csv:

  1. conforming legs (demo L1/L2) vs the OLD committed B18 demo CSVs (read
     from git HEAD: cl_vs_mach.csv 4-decimal strings, cross_model.csv, and
     the checks.csv M_max anchors). B21/B22 touched only the LS path, so bit
     agreement is expected; a drift is RECORDED as an independent B21/B22
     finding (pre-reg T1), not a B27 failure.
  2. LS A/C ceiling legs (demo L4/L5 npz caches) vs the committed B26
     g1_summary.csv (m_last / m_final / cls / die res / nlim / nflr / die
     Mmax / end-state cl_p / reached). Same code, same 16 threads -> bit
     agreement expected; drift is disclosed (pre-reg T5).
  3. the same legs' per-level records vs the committed B26 g1_levels.csv
     (per level: res / mach_max / nlim / nflr / cl_kj / n_newton / tag /
     converged).

verdict column: "bit" = exactly equal (4-decimal string equality for the
conforming cl_p rows), "drift" = any difference (abs_diff recorded).

Run:  python cases/analysis/b27_b18_demo_refresh/run_g27.py
Artifacts: results/g27_consistency.csv
"""

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEMO = REPO_ROOT / "cases/demo/b18_wingbody_transonic/results"
DEMO_REL = "cases/demo/b18_wingbody_transonic/results"
B26 = REPO_ROOT / "cases/analysis/b26_ls_transonic_ceiling/results"
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ROWS = []


def add(family, leg, item, committed, rerun):
    """Record one compared quantity; verdict bit iff exactly equal."""
    c, r = str(committed), str(rerun)
    try:
        diff = abs(float(c) - float(r))
        bit = diff == 0.0
    except ValueError:
        diff = ""
        bit = c == r
    ROWS.append(dict(family=family, leg=leg, item=item, committed=c,
                     rerun=r, abs_diff=diff, verdict="bit" if bit else "drift"))
    return bit


def git_show(rel):
    return subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=REPO_ROOT,
                          capture_output=True, text=True, check=True).stdout


def read_csv_text(text):
    return list(csv.DictReader(text.splitlines()))


def read_csv_path(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- conforming
def diff_conforming():
    old_cl = {(r["mach"], r["path"], r["resolution"]): r["cl_p"]
              for r in read_csv_text(git_show(f"{DEMO_REL}/cl_vs_mach.csv"))}
    new_cl = {(r["mach"], r["path"], r["resolution"]): r["cl_p"]
              for r in read_csv_path(DEMO / "cl_vs_mach.csv")}
    for key, old_v in sorted(old_cl.items()):
        leg = f"conf {key[1]} {key[2]} M{key[0]}"
        if key in new_cl:
            add("conforming_cl_p", leg, "cl_p(4dec)", old_v, new_cl[key])
        else:
            add("conforming_cl_p", leg, "cl_p(4dec)", old_v, "<missing>")

    old_x = read_csv_text(git_show(f"{DEMO_REL}/cross_model.csv"))
    new_x = read_csv_path(DEMO / "cross_model.csv")
    for ro in old_x:
        if ro["resolution"] != "coarse":
            continue  # the medium 0.5 row is a hardcoded anchor pair
        rn = [r for r in new_x if r["mach"] == ro["mach"]
              and r["resolution"] == "coarse"]
        add("conforming_cl_p", f"conf coarse cross M{ro['mach']}",
            "cross conf cl_p(4dec)", ro["conf"],
            rn[0]["conf"] if rn else "<missing>")

    # the checks.csv scalar anchors: conf coarse M_max and the reached flags
    old_checks = git_show(f"{DEMO_REL}/checks.csv")
    m = re.search(r"conforming_coarse_M084,reached=(\w+) cl_p=[\d.]+ \(M_max ([\d.]+)\)",
                  old_checks)
    if m:
        d = np.load(DEMO / "conf_coarse_084.npz")
        add("conforming_anchor", "conf coarse M0.84", "reached",
            m.group(1), str(bool(d["reached"])))
        add("conforming_anchor", "conf coarse M0.84", "M_max(2dec)",
            m.group(2), f"{float(d['mmax']):.2f}")
    m = re.search(r"conforming_medium_M079,reached=(\w+)", old_checks)
    if m:
        d = np.load(DEMO / "conf_medium_079.npz")
        add("conforming_anchor", "conf medium M0.79", "reached",
            m.group(1), str(bool(d["reached"])))


# ---------------------------------------------------------------- LS A/C
SUMMARY_ITEMS = [  # (demo npz field, B26 g1_summary column)
    ("reached", "target_reached"), ("m_last", "m_last_converged"),
    ("m_final", "die_m"), ("cls", "cls"), ("res", "die_res"),
    ("nlim", "die_nlim"), ("nflr", "die_nflr"), ("mmax", "die_mmax"),
    ("clp", "cl_p"),
]
LEVEL_ITEMS = ["residual_norm", "mach_max", "n_limited", "n_floored",
               "cl_kj", "n_newton", "tag", "converged"]


def _norm(v):
    """Normalize a value for exact comparison across npz/json/csv."""
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, float):
        return repr(v)
    return str(v)


def diff_ls():
    summary = {(r["side"], r["level"]): r for r in read_csv_path(B26 / "g1_summary.csv")}
    levels = {}
    for r in read_csv_path(B26 / "g1_levels.csv"):
        levels.setdefault((r["side"], r["level"]), []).append(r)
    for side in "AC":
        for level in ("coarse", "medium"):
            leg = f"LS {side} {level} ceiling"
            cache = DEMO / f"ls_{side}_{level}_084.npz"
            if not cache.exists():
                add("ls_summary", leg, "<cache>", "present", "<missing>")
                continue
            d = np.load(cache, allow_pickle=True)
            demo = dict(reached=bool(d["reached"]), m_last=d["m_last"].item(),
                        m_final=float(d["m_final"]), cls=str(d["cls"]),
                        res=float(d["res"]), nlim=int(d["nlim"]),
                        nflr=int(d["nflr"]), mmax=float(d["mmax"]),
                        clp=float(d["clp"]))
            s = summary[(side, level)]
            for dk, sk in SUMMARY_ITEMS:
                add("ls_summary", leg, sk, s[sk], _norm(demo[dk]))
            # per-level diff vs g1_levels
            dl = json.loads(str(d["levels_json"]))
            bl = levels[(side, level)]
            add("ls_levels", leg, "n_levels", len(bl), len(dl))
            for i, (bv, dv) in enumerate(zip(bl, dl)):
                add("ls_levels", leg, f"level[{i}].m_inf", bv["m_inf"],
                    repr(float(dv["m_inf"])))
                for k in LEVEL_ITEMS:
                    add("ls_levels", leg, f"level[{i}@{float(dv['m_inf']):.4g}].{k}",
                        bv[k], _norm(dv[k]))


# ---------------------------------------------------------------- main
def main():
    diff_conforming()
    diff_ls()
    keys = ["family", "leg", "item", "committed", "rerun", "abs_diff", "verdict"]
    with open(OUT / "g27_consistency.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(ROWS)
    n_bit = sum(r["verdict"] == "bit" for r in ROWS)
    drifts = [r for r in ROWS if r["verdict"] == "drift"]
    print(f"wrote {OUT/'g27_consistency.csv'}: {n_bit}/{len(ROWS)} bit-identical, "
          f"{len(drifts)} drift(s)")
    fam = {}
    for r in ROWS:
        f = fam.setdefault(r["family"], [0, 0])
        f[r["verdict"] != "bit"] += 1
    print("\n--- gate summary (pre-reg 2.4) ---")
    for k, (b, d) in sorted(fam.items()):
        print(f"  {k}: {b} bit, {d} drift")
    gb271 = all(r["verdict"] == "bit" for r in ROWS
                if r["family"].startswith("conforming"))
    gb272 = all(r["verdict"] == "bit" for r in ROWS
                if r["family"] in ("ls_summary", "ls_levels"))
    print(f"  GB27.1 (conforming bit-reproduces the committed B18 anchors): "
          f"{'PASS' if gb271 else 'RECORDED -- drift, independent B21/B22 finding'}")
    print(f"  GB27.2 (LS A/C bit-reproduces the committed B26 anchors): "
          f"{'PASS' if gb272 else 'RECORDED -- drift disclosed above'}")
    for r in drifts[:40]:
        print(f"    DRIFT {r['family']} {r['leg']} {r['item']}: "
              f"{r['committed']} -> {r['rerun']}")


if __name__ == "__main__":
    main()
