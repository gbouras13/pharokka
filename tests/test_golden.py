"""Golden-output regression tests.

For each case in ``golden_cases.CASES`` this runs pharokka end-to-end and diffs
the curated key output files against the committed reference copies under
``test_data/golden/<case>/``.  Unlike test_overall.py (which only asserts a zero
exit code), these assert that the *content* of pharokka's output is unchanged.

The comparison (compare_outputs.compare_against_golden) tolerates timestamps,
float-formatting differences, and the documented non-deterministic mash
tie-break, so a failure here means a genuine change in pharokka's output.

These are marked ``golden`` so they can be deselected with ``-m "not golden"``.
Regenerate the reference outputs with ``tests/generate_golden.py`` after an
intended output change.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_outputs import compare_against_golden
from golden_cases import CASES, GOLDEN_DIR, run_pharokka

DATABASE_DIR = Path(os.environ.get("PHAROKKA_DB", "tests/test_data/database"))
THREADS = 8


@pytest.fixture(scope="session")
def ensure_database():
    """Ensure the pinned pharokka database is present before running cases.

    ``pharokka install`` is idempotent — databases.instantiate_install only
    downloads when files are missing — so this is a fast no-op when the DB is
    already in place.  It exists to guard test ordering: pytest may collect
    test_golden before test_overall::test_download, so we cannot rely on that
    test having installed the DB first.
    """
    subprocess.run(["pharokka", "install", "-o", str(DATABASE_DIR)], check=True)
    return DATABASE_DIR


@pytest.mark.golden
@pytest.mark.parametrize("case", sorted(CASES))
def test_golden(case, tmp_path, ensure_database):
    """Run a pharokka case and diff key outputs against the committed golden copy."""
    golden = GOLDEN_DIR / case
    if not golden.is_dir() or not any(golden.iterdir()):
        pytest.skip(
            f"no golden outputs committed for case {case!r} "
            f"(generate with tests/generate_golden.py)"
        )

    run_pharokka(CASES[case], tmp_path, DATABASE_DIR, THREADS)

    diffs, notices = compare_against_golden(tmp_path, golden)
    for n in notices:
        print(n)
    assert not diffs, (
        f"Golden output mismatch for case {case!r} "
        f"(regenerate with tests/generate_golden.py if intended):\n" + "\n".join(diffs)
    )
