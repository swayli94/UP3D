"""Stamp the GEOMETRY AXES into every mesh manifest's `extra`, so a level name
never has to be decoded.

WHY. The M6 wing family encodes TWO INDEPENDENT axes in one name, and only one of
them is written down:

    `_ss`    -> the far-field clamp is OFF (self-similar).  NOT about the tip.
    `_flat`  -> flat tip cap.
    (absent) -> ROUND tip cap ... but ONLY SINCE 2026-08-04.

The second axis is carried by the ABSENCE of a suffix, and that absence changed
meaning once. That change is the 2.0 % that stood for five days and was restated
as "unexplained" in eight documents before phase 3 round 36 pinned it to the
flat->round mesh-family switch. So the name is not a safe thing to reason from.

WHAT THIS DOES INSTEAD. It reads each generator's OWN `LEVELS` table -- the
authoritative source, not my re-derivation of the naming rule -- and writes the
axes next to the mesh, where they are read. Renaming the levels would be the
other fix, but `coarse_ss` alone appears in 28 py/md files and a committed CSV,
so a rename either rewrites committed evidence (the G3 hazard) or leaves stale
CSVs behind. This is the same information at a fraction of the risk.

★ `extra` is quarantined from `sha256` by construction (locked by
tests/test_gs40_provenance.py), so stamping cannot change a mesh's identity.
★ It records h_far and h_far/h_wall rather than a boolean "clamped": the ratio
IS the question ("is this level on the refinement ray?"), and 120 versus less
answers it without trusting my reading of the clamp logic.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest, write_manifest  # noqa: E402


def _load(gen_path, names):
    """Import a generator module by path and pull the named attributes."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_gen", gen_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return {n: getattr(mod, n, None) for n in names}


def main():
    m6 = _load(REPO / "cases/meshes/onera_m6/generate_onera_m6.py",
               ("LEVELS", "RICHARDSON_LADDER"))
    rt = _load(REPO / "cases/meshes/onera_m6_roundtip/generate_onera_m6_roundtip.py",
               ("LEVELS", "RICHARDSON_LADDER"))
    print(f"onera_m6 LEVELS: {len(m6['LEVELS'])}   RICHARDSON_LADDER "
          f"{m6['RICHARDSON_LADDER']}")
    print(f"roundtip LEVELS: {len(rt['LEVELS'])}   RICHARDSON_LADDER "
          f"{rt['RICHARDSON_LADDER']}\n")

    FAM = (("cases/meshes/onera_m6", m6, False),
           ("cases/meshes/onera_m6_roundtip", rt, True))
    n = 0
    for rel, gen, always_round in FAM:
        d = REPO / rel
        for msh in sorted(d.glob("*.msh")):
            level = msh.stem
            p = gen["LEVELS"].get(level)
            if p is None:
                print(f"  {rel}/{msh.name}: level not in the generator's LEVELS "
                      f"-- SKIPPED (recording nothing beats guessing)")
                continue
            h_wall, h_far = float(p["h_wall"]), float(p["h_far"])
            ratio = h_far / h_wall
            #: ★ derived from the generator's OWN rule, quoted inline
            tip = "round" if always_round else (
                "flat" if level.endswith("_flat") else "round")
            on_ray = abs(ratio - 120.0) < 1e-9
            axes = dict(
                tip_cap=tip,
                h_wall=h_wall, h_far=h_far,
                h_far_over_h_wall=round(ratio, 4),
                on_refinement_ray=on_ray,
                in_richardson_ladder=level in (gen["RICHARDSON_LADDER"] or ()),
                axes_note=(
                    "TWO INDEPENDENT AXES, only one of which the level name states. "
                    "`_ss` = the far-field clamp h_far=min(2.5,120*h_wall) is OFF "
                    "(self-similar) and says NOTHING about the tip. `_flat` = flat "
                    "tip cap; ABSENCE of `_flat` = ROUND, but only since 2026-08-04 "
                    "-- before that date the same name meant FLAT, and that silent "
                    "flip is the 2.0 % that stood five days (phase 3 round 36, "
                    "F-MESH). h_far_over_h_wall == 120 means the clamp did not bite, "
                    "i.e. this level IS on the refinement ray. Source: each "
                    "generator's own LEVELS table, read at stamping time."),
            )
            prev = read_manifest(msh) or {}
            keep = {k: v for k, v in (prev.get("extra") or {}).items()
                    if k not in axes}
            write_manifest(msh, **{**keep, **axes})
            n += 1
            print(f"  {rel}/{msh.name:22} tip={tip:5} h_wall={h_wall:<6} "
                  f"h_far={h_far:<5} ratio={ratio:6.1f} "
                  f"{'ON-RAY' if on_ray else 'CLAMPED (off the ray)'}"
                  f"{'  [RICHARDSON]' if axes['in_richardson_ladder'] else ''}")

    #: G-IDENTITY: stamping must not move any sha256 -- extra is quarantined
    bad = [str(m) for rel, _g, _r in FAM
           for m in (REPO / rel).glob("*.msh")
           if (read_manifest(m) or {}).get("sha256") != mesh_fingerprint(m)["sha256"]]
    assert not bad, f"G-IDENTITY: stamping changed a fingerprint: {bad}"
    print(f"\n{n} manifests stamped; G-IDENTITY PASS (no sha256 moved).")


if __name__ == "__main__":
    main()
