"""pytest plugin for the baseline-provenance arms: RECORDS ONLY, never changes a verdict.

Loaded with ``-p prov_plugin`` (module on PYTHONPATH). It:

1. wraps ``tests._gate_evidence.assert_matches_committed`` / ``assert_within_spread``
   BEFORE any test module imports them (``pytest_configure`` runs before collection), and
   appends to ``$PROV_OUT_DIR/fresh.jsonl`` the fresh values at full precision (``repr``)
   together with the committed CSV values for the same keys; then calls the original;
2. dumps D13's module-scoped ``runs`` fixture (D13 has its own comparison code) to
   ``$PROV_OUT_DIR/d13_runs.npz``;
3. at session end, writes per-dispatcher numba cache hit/miss counts for every pyfp3d
   module-level dispatcher to ``numba_stats.json`` (pilot for ruling A, premise 1 --
   recorded, NOT asserted), and where ``pyfp3d`` was imported from.

Every recording step is wrapped so that a recording error can only lose a record, never
alter an outcome. Recording errors go to ``plugin_errors.txt``.
"""
import csv
import inspect
import json
import os
import sys
import traceback

import pytest

OUT = os.environ["PROV_OUT_DIR"]
WT = os.path.realpath(os.environ["PROV_WORKTREE"])
_state = {"nodeid": None, "d13_done": False}


def _err(where):
    with open(os.path.join(OUT, "plugin_errors.txt"), "a", encoding="utf-8") as fh:
        fh.write(f"--- {where} @ {_state['nodeid']}\n{traceback.format_exc()}\n")


def _ser(v):
    try:
        return repr(float(v))
    except Exception:
        return repr(v)


def _default_key(r):
    return tuple(r[k] for k in ("level", "surface", "case", "x_c", "nx", "leg") if k in r)


def _record(fn_name, sig, args, kwargs):
    ba = sig.bind(*args, **kwargs)
    ba.apply_defaults()
    a = ba.arguments
    gate_dir, fresh = str(a["gate_dir"]), a["fresh"]
    cols = list(a["measured"]) if "measured" in a else list(a["spread_tol"])
    key_of = a.get("key_of") or _default_key
    committed = {}
    path = os.path.join(gate_dir, a.get("filename") or "summary.csv")
    if os.path.exists(path):
        with open(path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                k = key_of(r)
                if k in fresh:
                    committed[repr(k)] = {c: r.get(c) for c in cols}
    rec = {
        "nodeid": _state["nodeid"], "fn": fn_name,
        "gate_dir": os.path.relpath(os.path.realpath(gate_dir), WT),
        "cols": cols,
        "fresh": {repr(k): {c: _ser(v) for c, v in row.items()} for k, row in fresh.items()},
        "committed": committed,
    }
    if "spread_tol" in a:
        rec["spread_tol"] = {c: float(t) for c, t in a["spread_tol"].items()}
    if "rel_tol" in a:
        rec["rel_tol"] = float(a["rel_tol"])
    with open(os.path.join(OUT, "fresh.jsonl"), "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")


def _wrap(mod, name):
    orig = getattr(mod, name)
    sig = inspect.signature(orig)

    def wrapper(*args, **kwargs):
        try:
            _record(name, sig, args, kwargs)
        except Exception:
            _err(f"record {name}")
        return orig(*args, **kwargs)

    wrapper.__wrapped__ = orig
    wrapper.__doc__ = orig.__doc__
    setattr(mod, name, wrapper)


def pytest_configure(config):
    import tests._gate_evidence as ge
    assert os.path.realpath(os.path.dirname(os.path.dirname(ge.__file__))) == WT, (
        f"tests._gate_evidence imported from {ge.__file__}, not the worktree {WT}")
    _wrap(ge, "assert_matches_committed")
    _wrap(ge, "assert_within_spread")


def pytest_collection_finish(session):
    mod = sys.modules.get("pyfp3d")
    where = os.path.realpath(mod.__file__) if mod is not None else None
    with open(os.path.join(OUT, "pyfp3d_origin.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"{where}\n")
    if where is None or not where.startswith(WT + os.sep):
        pytest.exit(f"pyfp3d imported from {where}, not the worktree {WT}", returncode=3)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    _state["nodeid"] = item.nodeid
    yield
    _state["nodeid"] = None


def _flatten(prefix, obj, out):
    import numpy as np
    if isinstance(obj, dict):
        for k, v in obj.items():
            _flatten(f"{prefix}/{k}" if prefix else str(k), v, out)
    elif isinstance(obj, (list, tuple)) and obj and all(isinstance(x, dict) for x in obj):
        for i, v in enumerate(obj):
            _flatten(f"{prefix}/{i}", v, out)
    else:
        try:
            arr = np.asarray(obj)
            if arr.dtype.kind in "fiub":
                out[prefix] = arr
        except Exception:
            pass


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    yield
    try:
        if (not _state["d13_done"] and "test_D13_ibl_vs_xfoil" in item.nodeid
                and "runs" in getattr(item, "funcargs", {})):
            import numpy as np
            flat = {}
            _flatten("", item.funcargs["runs"], flat)
            np.savez(os.path.join(OUT, "d13_runs.npz"), **flat)
            _state["d13_done"] = True
    except Exception:
        _err("d13 dump")


def pytest_sessionfinish(session, exitstatus):
    try:
        from numba.core.dispatcher import Dispatcher
        seen, per = set(), {}
        for name, mod in list(sys.modules.items()):
            if not (name == "pyfp3d" or name.startswith("pyfp3d.")) or mod is None:
                continue
            for attr, obj in list(vars(mod).items()):
                if isinstance(obj, Dispatcher) and id(obj) not in seen:
                    seen.add(id(obj))
                    st = obj.stats
                    per[f"{name}.{attr}"] = {
                        "hits": int(sum(st.cache_hits.values())),
                        "misses": int(sum(st.cache_misses.values())),
                        "n_signatures": len(obj.signatures),
                    }
        tot = {"dispatchers": len(per),
               "hits": sum(v["hits"] for v in per.values()),
               "misses": sum(v["misses"] for v in per.values()),
               "exitstatus": int(exitstatus)}
        with open(os.path.join(OUT, "numba_stats.json"), "w", encoding="utf-8") as fh:
            json.dump({"total": tot, "per_dispatcher": per}, fh, indent=1, sort_keys=True)
    except Exception:
        _err("numba stats")
