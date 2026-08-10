"""Is the loose viscous loop run-to-run reproducible? (GS1.4 by-catch)

`tests/test_v6_wake_sheet.py::test_ab_bit_identity_gate_free_library` compares
the loose-loop phi of the current tree against a pinned pre-V6 commit and
asserted bit identity. It began failing with an O(1) difference (max relative
1.02 over 6104/6106 nodes). Attribution, in order:

  1. NOT caused by the GS1.4 clamp change -- reverting it changes nothing.
  2. NOT caused by this branch -- the same test fails on `main`.
  3. The premise itself is false: running the SAME COMMIT twice, on the same
     machine with the same thread count, gives max relative 1.024 over
     6104/6106 nodes -- exactly the magnitude the test reports.

So the loose loop at coarse is not run-to-run reproducible, and a bit-identity
(or any tight-tolerance) A/B on it cannot be a regression guard. This script is
the standing evidence for that; it runs the same leg twice and prints the
difference. The inviscid path IS deterministic (bench/bitcheck.py: 10/10), which
localises the non-determinism to pyfp3d/viscous/ -- consistent with phase one's
own discipline #12.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import tests.test_v6_wake_sheet as V     # noqa: E402


def main():
    tmp = Path(tempfile.mkdtemp())
    snip = tmp / "ab_leg.py"
    snip.write_text(V._AB_SNIPPET)
    outs = []
    for i in (1, 2):
        d = V._ab_leg(snip, tmp / f"wt{i}", V.GATE_FREE_BASELINE,
                      tmp / f"o{i}.npz")
        outs.append(np.asarray(d["phi"], dtype=np.float64))
        print(f"leg {i} (ref {V.GATE_FREE_BASELINE}): phi[:3] = "
              f"{outs[-1][:3]}", flush=True)
    a, b = outs
    d = np.abs(a - b)
    sc = np.maximum(np.abs(a), np.abs(b))
    rel = float(np.nanmax(d / np.where(sc > 0, sc, 1.0)))
    print(f"\nsame commit, twice: bitwise equal = {bool(np.array_equal(a, b))}")
    print(f"  differing entries {int((d > 0).sum())}/{a.size}   "
          f"max|d| {d.max():.3e}   max relative {rel:.3e}")
    print("\nreading: the loose loop is NOT reproducible run to run, so no "
          "A/B on it can serve as a regression guard.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
