"""Shared definitions for the golden-output regression tests.

``CASES`` maps a case name to the pharokka CLI arguments (everything that
follows ``pharokka`` on the command line, minus ``-d``/``-o``/``-t``/``-f``,
which the runner appends).  ``KEY_FILES`` is the curated set of user-facing,
deterministic output files that are committed under ``test_data/golden/<case>/``
and diffed against a fresh run by ``tests/test_golden.py``.

Reproducibility notes
---------------------
* The pharokka database is version-pinned (``databases.VERSION = "1.8.0"`` with
  an md5 check), so the database-dependent outputs (mash / mmseqs / phrog) are
  reproducible until someone deliberately bumps the database — at which point a
  golden failure is the *intended* signal to regenerate.
* Every case passes ``-l PHARTEST``.  Without an explicit locus tag pharokka
  generates a random 8-character prefix (see ``post_processing.create_gff``),
  which would make the locus tags — and therefore every output file — differ on
  each run.  (In ``-m``/meta mode the locus tag is derived from the contig name
  regardless, but passing it is harmless and keeps the cases uniform.)
* The mash top-hits file has a documented non-deterministic tie-break that the
  comparator (``compare_outputs._mash_tiebreak``) tolerates.
* The GenBank ``LOCUS`` line carries today's date; the comparator skips it.

Regenerate the golden outputs with::

    PHAROKKA_DB=/path/to/pharokka_db python tests/generate_golden.py
"""

import shlex
import subprocess
from pathlib import Path

TEST_DATA = Path(__file__).resolve().parent / "test_data"
OVERALL = TEST_DATA / "overall"
GOLDEN_DIR = TEST_DATA / "golden"

CASES = {
    # Standard single-contig phage, full pipeline (phanotate + mmseqs + mash +
    # tRNA/CRISPR/tmRNA scans).
    "standard": f"run -i {OVERALL}/Standard_examples/SAOMS1.fasta -l PHARTEST",
    # --fast path: pyhmmer instead of mmseqs2 (mash still runs).
    "fast": f"run -i {OVERALL}/Standard_examples/SAOMS1.fasta --fast -l PHARTEST",
    # Meta mode over multiple contigs (per-contig locus tags).
    "meta": f"run -i {OVERALL}/Meta_example/combined_meta.fasta -m -l PHARTEST",
    # CRISPR detection (MinCED).
    "crispr": f"run -i {OVERALL}/CRISPR_example/Biggiephage_A_fullcontig_CasΦ1.fasta -l PHARTEST",
    # tmRNA detection (Aragorn).
    "tmrna": f"run -i {OVERALL}/tmRNA_example/NC_051700.fasta -l PHARTEST",
    # GenBank input path (--genbank).
    "genbank": f"run -i {OVERALL}/genbank_examples/SAOMS1.gbk --genbank -l PHARTEST",
}

# Curated, user-facing, deterministic output files.  Logs, binaries (.msh),
# and volatile intermediates (*.faa/*.ffn, tmp dirs) are intentionally excluded.
# A file that a given case does not produce is simply not committed for that
# case, and golden comparison only checks files that exist under the golden dir.
KEY_FILES = [
    "pharokka.gff",
    "pharokka.gbk",
    "pharokka.tbl",
    "pharokka_cds_final_merged_output.tsv",
    "pharokka_cds_functions.tsv",
    "pharokka_length_gc_cds_density.tsv",
    "pharokka_top_hits_mash_inphared.tsv",
    "pharokka_minced.gff",
    "pharokka_aragorn.gff",
    "trnascan_out.gff",
]


def run_pharokka(case_args: str, out_dir: Path, db_dir: Path, threads: int = 8) -> None:
    """Run a single pharokka case into ``out_dir`` (``-f`` overwrites).

    Uses an argv list (no shell) so the non-ASCII CRISPR input filename and any
    quoted sub-arguments are passed through verbatim.  Raises RuntimeError with
    captured output on a non-zero exit so test failures are diagnosable.
    """
    cmd = [
        "pharokka",
        *shlex.split(case_args),
        "-d",
        str(db_dir),
        "-o",
        str(out_dir),
        "-t",
        str(threads),
        "-f",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"pharokka failed (exit {proc.returncode}):\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  stdout:\n{proc.stdout}\n  stderr:\n{proc.stderr}"
        )
