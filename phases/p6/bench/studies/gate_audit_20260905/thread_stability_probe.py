"""pytest plugin: replace the evidence lock with a COLLECTOR.

The real `assert_matches_committed` raises on the FIRST mismatch, so a red run
names one column and hides the rest.  This patch reuses each gate's own
fixtures (no reimplementation of the solve) and reports every column's worst
relative difference against the committed CSV.
"""
import csv, os, collections
import tests._gate_evidence as GE

WORST = collections.defaultdict(float)      # (gate, col) -> worst rel
COUNT = collections.defaultdict(int)


def _patched(gate_dir, fresh, measured, rel_tol=GE.DEFAULT_REL_TOL,
             key_of=None, refresh_hint=None, filename="summary.csv"):
    path = os.path.join(str(gate_dir), filename)
    gate = os.path.basename(str(gate_dir).rstrip("/"))
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if key_of is None:
        def key_of(r):
            return tuple(r[k] for k in ("level", "surface", "case", "x_c", "nx", "leg")
                         if k in r)
    n = 0
    for r in rows:
        k = key_of(r)
        if k not in fresh:
            continue
        for col in measured:
            if col not in r or r[col] in ("", "-", "nan"):
                continue
            got = fresh[k].get(col)
            if got is None:
                continue
            want = float(r[col]); n += 1
            rel = abs(float(got) - want) / max(abs(want), 1e-30)
            COUNT[(gate, col)] += 1
            if rel > WORST[(gate, col)]:
                WORST[(gate, col)] = rel
    return n


GE.assert_matches_committed = _patched
for _m in ("test_D05_euler_naca0012", "test_D06_euler_rae2822",
           "test_D07_euler_onera_m6", "test_D08_rans_naca0012"):
    try:
        mod = __import__(f"tests.D.{_m}", fromlist=["x"])
        mod.assert_matches_committed = _patched
    except Exception as e:                                       # noqa: BLE001
        print(f"  [plugin] could not patch {_m}: {e}")


def pytest_sessionfinish(session, exitstatus):
    print("\n\n=== per-column worst relative difference vs the committed CSV ===")
    print(f"{'gate':28} {'column':26} {'n':>4} {'worst rel':>11}  verdict")
    for (gate, col), rel in sorted(WORST.items(), key=lambda kv: -kv[1]):
        v = "STABLE" if rel <= 1e-6 else ("**UNSTABLE**" if rel > 1e-6 else "")
        print(f"{gate:28} {col:26} {COUNT[(gate,col)]:4d} {rel:11.2e}  {v}")
