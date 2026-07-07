#!/usr/bin/env python3
"""Regenerate the committed golden outputs used by tests/test_golden.py.

Runs each case in ``golden_cases.CASES`` against the pinned pharokka database
and copies the curated ``KEY_FILES`` into ``tests/test_data/golden/<case>/``,
overwriting whatever is there.  Run this whenever pharokka's output is meant to
change (a real behaviour change, or a deliberate database version bump) — then
review the git diff of the golden files before committing.

Usage (from the repo root)::

    PHAROKKA_DB=/path/to/pharokka_db python tests/generate_golden.py
    # optionally limit to specific cases:
    PHAROKKA_DB=/path/to/pharokka_db python tests/generate_golden.py standard meta
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from golden_cases import CASES, GOLDEN_DIR, KEY_FILES, run_pharokka

THREADS = 8


def main(argv: list[str]) -> int:
    db_dir = Path(os.environ.get("PHAROKKA_DB", "tests/test_data/database"))
    if not db_dir.is_dir():
        print(
            f"ERROR: pharokka database not found at {db_dir}. "
            f"Set PHAROKKA_DB to your database directory."
        )
        return 2

    wanted = argv or sorted(CASES)
    unknown = [c for c in wanted if c not in CASES]
    if unknown:
        print(
            f"ERROR: unknown case(s): {', '.join(unknown)}. "
            f"Known: {', '.join(sorted(CASES))}"
        )
        return 2

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    for case in wanted:
        print(f"\n══ generating golden: {case} ══")
        with tempfile.TemporaryDirectory(prefix=f"golden_{case}_") as tmp:
            tmp_out = Path(tmp) / "out"
            run_pharokka(CASES[case], tmp_out, db_dir, THREADS)

            dest = GOLDEN_DIR / case
            if dest.exists():
                shutil.rmtree(dest)
            dest.mkdir(parents=True)

            copied = []
            for name in KEY_FILES:
                src = tmp_out / name
                if src.is_file():
                    shutil.copy2(src, dest / name)
                    copied.append(name)
            print(f"   copied {len(copied)} file(s): {', '.join(copied)}")

    print(
        f"\nDone. Golden outputs under {GOLDEN_DIR}. Review the diff before committing."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
