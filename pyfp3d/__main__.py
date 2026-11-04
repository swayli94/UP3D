from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence

from .freestream import run_freestream


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Apply a uniform inviscid free stream to a tetrahedron JSON mesh.",
    )
    parser.add_argument("input", help="path to a format_version 1 tetrahedron JSON file")
    args = parser.parse_args(argv)

    try:
        result = run_freestream(args.input)
    except (OSError, ValueError) as error:
        print(f"pyfp3d: {error}", file=sys.stderr)
        return 1

    print(json.dumps(result, ensure_ascii=False, allow_nan=False, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
