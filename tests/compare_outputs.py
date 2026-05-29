#!/usr/bin/env python3
"""
Compare pharokka outputs between two runs (e.g. refactored dev vs bioconda 1.9.1).
Ignores timestamp-dependent content (log files, GenBank LOCUS date, GFF pragma dates).
For TSV/TXT files, floating-point fields are compared with numeric tolerance so that
formatting differences (e.g. "82.0" vs "82", "3.473e-22" vs "3.4729999999999998e-22",
"1.741e-7" vs "1.741e-07") do not count as mismatches.

Usage:
    python tests/compare_outputs.py <dir_new> <dir_ref>

Exit code 0 = identical (modulo timestamps / float formatting), non-zero = differences found.
"""
import math
import os
import re
import sys
from pathlib import Path

# ── patterns that are timestamp/run-specific and should be ignored ──────────
SKIP_LINE_PATTERNS = [
    re.compile(r"^\d{4}-\d{2}-\d{2}"),          # log lines: 2026-05-26 ...
    re.compile(r"^#.*pharokka.*run"),             # GFF pragma with run info
    re.compile(r"^LOCUS\s+\S+\s+\d+ bp"),        # GenBank LOCUS date field
    re.compile(r"^##date"),                        # any ##date pragma
    re.compile(r"Creation Date"),
    re.compile(r"Time to find repeats:"),         # MinCED timing in minced_spacers.txt
]

# file extensions to skip entirely
SKIP_EXTENSIONS = {".log"}

# files to skip by name
SKIP_FILENAMES = {"pharokka.log"}

# directory components to skip entirely (any file under these dirs is ignored)
SKIP_DIRS = {"logs"}

# ── float-tolerant TSV comparison ────────────────────────────────────────────

# Matches a field that is entirely a number (int or float, with optional sign/exponent).
# Anchored so "1000/1000" or "No_PHROGs_HMM" never match.
_NUMERIC_FIELD_RE = re.compile(
    r'^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$'
)


def _numeric_equal(a: str, b: str) -> bool:
    """Return True if a and b represent the same number within tolerance.

    Handles:
    - identical strings            → True (fast path)
    - int vs float:  "82" == "82.0"
    - exponent padding: "1.741e-7" == "1.741e-07"
    - float precision: "3.473e-22" == "3.4729999999999998e-22"
      (same IEEE 754 double, different string repr from different parsers)
    """
    if a == b:
        return True
    if _NUMERIC_FIELD_RE.match(a) and _NUMERIC_FIELD_RE.match(b):
        try:
            fa, fb = float(a), float(b)
            if math.isnan(fa) and math.isnan(fb):
                return True
            if math.isinf(fa) or math.isinf(fb):
                return fa == fb
            # Use a tight relative tolerance: these should be the *same* number
            # just serialised differently, not genuinely different values.
            return math.isclose(fa, fb, rel_tol=1e-9, abs_tol=1e-300)
        except (ValueError, OverflowError):
            pass
    return False


def _tsv_lines_match(a: str, b: str) -> bool:
    """Compare two tab-separated lines with per-field float tolerance."""
    if a == b:
        return True
    fields_a = a.split('\t')
    fields_b = b.split('\t')
    if len(fields_a) != len(fields_b):
        return False
    return all(_numeric_equal(fa, fb) for fa, fb in zip(fields_a, fields_b))


# Matches a contig name that is purely an integer (unicycler-style numeric headers).
_INTEGER_CONTIG_RE = re.compile(r'^\d+$')


def _cds_functions_v191_intcontig_bugs(new_lines: list, ref_lines: list) -> tuple:
    """Detect v1.9.1 pandas int/str contig-name bug in pharokka_cds_functions.tsv.

    When contig names are purely numeric (e.g. unicycler integer headers like '1','2','4'),
    pharokka v1.9.1 reads minced.gff / aragorn.gff with pandas which infers the 'contig'
    column as int64.  The loop variable is a str (from length_df["contig"].astype("string")).
    pandas comparison int(4) == str("4") silently returns False, so CRISPRs and tmRNAs are
    counted as 0 even though the raw GFF files contain the annotations.

    The NEW (polars) code correctly counts them.  We treat these as known v1.9.1 bugs and
    demote them to notices rather than hard failures.

    Returns (bug_count: int, non_bug_diff_count: int).
    """
    def parse_cds_functions(lines):
        """Return {(Description, contig): count_str} skipping the header."""
        d = {}
        for line in lines:
            if not line or line.startswith("Description"):
                continue
            parts = line.split('\t')
            if len(parts) != 3:
                continue
            desc, count_str, contig = parts
            d[(desc, contig)] = count_str
        return d

    new_d = parse_cds_functions(new_lines)
    ref_d = parse_cds_functions(ref_lines)

    bug_count = 0
    non_bug_diff_count = 0

    for key in set(new_d) | set(ref_d):
        desc, contig = key
        new_val = new_d.get(key)
        ref_val = ref_d.get(key)
        if new_val == ref_val:
            continue
        # v1.9.1 bug: CRISPRs or tmRNAs on purely-integer contig names counted as 0 in ref
        if (desc in {"CRISPRs", "tmRNAs"}
                and _INTEGER_CONTIG_RE.match(contig)
                and new_val is not None
                and new_val != "0"
                and (ref_val is None or ref_val == "0")):
            bug_count += 1
        else:
            non_bug_diff_count += 1

    return bug_count, non_bug_diff_count


def _mash_tiebreak(a: str, b: str) -> bool:
    """Return True if a and b differ ONLY because of a non-deterministic mash tie-break.

    mash dist output order is non-deterministic when multiple database entries share
    the same minimum Hamming distance (distance=0.0, same 1000/1000 hashes).  When
    this happens, the Accession chosen as the 'top hit' may differ between runs even
    though both answers are equally valid.  We treat such differences as harmless.

    Columns of pharokka_top_hits_mash_inphared.tsv:
      0: contig   1: Accession   2: mash_distance   3: mash_pval
      4: mash_matching_hashes   5+: INPHARED metadata

    A pair of rows (a, b) is a tie-break if columns 0, 2, 3, 4 are numerically
    equal (same contig, distance, pval, hash-ratio) but column 1 differs.
    """
    if a == b:
        return False  # identical → not a tie-break, handled elsewhere
    fa = a.split('\t')
    fb = b.split('\t')
    if len(fa) < 5 or len(fb) < 5:
        return False
    # contig (col 0) must be identical
    if fa[0] != fb[0]:
        return False
    # Accession (col 1) must differ
    if fa[1] == fb[1]:
        return False
    # mash_distance (col 2), mash_pval (col 3), mash_matching_hashes (col 4) must match
    for col in (2, 3, 4):
        if not _numeric_equal(fa[col], fb[col]):
            return False
    return True


# ── file helpers ─────────────────────────────────────────────────────────────

def filter_lines(path: Path) -> list[str]:
    """Read a file and return lines with timestamp-like content removed."""
    try:
        lines = path.read_text(errors="replace").splitlines()
    except Exception as e:
        return [f"<ERROR reading {path}: {e}>"]
    filtered = []
    for line in lines:
        if any(p.search(line) for p in SKIP_LINE_PATTERNS):
            continue
        filtered.append(line)
    return filtered


def compare_dirs(dir_new: Path, dir_ref: Path, prefix: str = "") -> tuple[list[str], list[str]]:
    """Recursively compare two directories.

    Returns (diffs, notices) where diffs are real mismatches (cause non-zero exit)
    and notices are informational skips (e.g. non-deterministic mash tie-breaks).
    """
    diffs = []
    notices = []

    new_files = {f.relative_to(dir_new) for f in dir_new.rglob("*") if f.is_file()}
    ref_files = {f.relative_to(dir_ref) for f in dir_ref.rglob("*") if f.is_file()}

    only_new = new_files - ref_files
    only_ref = ref_files - new_files
    common   = new_files & ref_files

    def should_skip(rel: Path) -> bool:
        return (
            rel.suffix in SKIP_EXTENSIONS
            or rel.name in SKIP_FILENAMES
            or bool(SKIP_DIRS.intersection(rel.parts))
        )

    for f in sorted(only_new):
        if should_skip(f):
            continue
        diffs.append(f"  ONLY IN NEW : {f}")

    for f in sorted(only_ref):
        if should_skip(f):
            continue
        diffs.append(f"  ONLY IN REF : {f}")


    for rel in sorted(common):
        if should_skip(rel):
            continue

        fn = dir_new / rel
        fr = dir_ref / rel

        ln = filter_lines(fn)
        lr = filter_lines(fr)

        # For TSV/TXT files where row order may differ, compare sorted with
        # per-field float tolerance so formatting differences (e.g. "82" vs "82.0",
        # "3.473e-22" vs "3.4729999999999998e-22") do not count as mismatches.
        if rel.suffix in {".tsv", ".txt"}:
            sn = sorted(ln)
            sr = sorted(lr)
            if len(sn) != len(sr):
                diffs.append(f"  DIFFER (sorted) : {rel}")
                diffs.append(f"    line count: new={len(sn)} ref={len(sr)}")
                continue
            mismatches = [
                (i, a, b)
                for i, (a, b) in enumerate(zip(sn, sr))
                if not _tsv_lines_match(a, b)
            ]
            if mismatches:
                # For the mash top-hits file, filter out non-deterministic tie-breaks.
                # mash dist output order is non-deterministic when multiple database entries
                # share the same minimum distance (e.g. distance=0.0, 1000/1000 hashes).
                # The chosen Accession may differ between runs even though both are valid.
                is_mash_file = rel.name == "pharokka_top_hits_mash_inphared.tsv"
                if is_mash_file:
                    real_mismatches = [
                        (i, a, b) for (i, a, b) in mismatches
                        if not _mash_tiebreak(a, b)
                    ]
                    tiebreak_count = len(mismatches) - len(real_mismatches)
                    if tiebreak_count > 0:
                        notices.append(
                            f"  MASH TIE-BREAK SKIP : {rel} "
                            f"({tiebreak_count} row(s) differ only in tied-distance Accession — "
                            f"mash output order is non-deterministic for equal distances)"
                        )
                    if not real_mismatches:
                        # All differences were non-deterministic tie-breaks → no real diff
                        continue
                    # Mix: only report real differences
                    mismatches = real_mismatches

                # For pharokka_cds_functions.tsv, check for the v1.9.1 pandas int/str
                # contig-name bug: when contig names are purely numeric (unicycler integer
                # headers), pandas infers the GFF contig column as int64 but the loop
                # variable is str, so CRISPRs/tmRNAs are counted as 0 in the reference
                # even though the raw GFF contains the annotations.  Demote to notices.
                is_cds_functions_file = rel.name == "pharokka_cds_functions.tsv"
                if is_cds_functions_file:
                    bug_count, non_bug_count = _cds_functions_v191_intcontig_bugs(ln, lr)
                    if bug_count > 0:
                        notices.append(
                            f"  V1.9.1 INT-CONTIG BUG SKIP : {rel} "
                            f"({bug_count} CRISPRs/tmRNAs row(s) differ because pharokka v1.9.1 "
                            f"fails to count features on purely-numeric contig names — pandas "
                            f"reads the GFF contig column as int64 vs str in the loop)"
                        )
                    if non_bug_count == 0:
                        # All differences are the known v1.9.1 bug → no real diff
                        continue
                    # Some real mismatches remain; fall through to report them

                diffs.append(f"  DIFFER (sorted) : {rel}")
                for i, a, b in mismatches[:11]:
                    diffs.append(f"    new[{i}]: {a[:120]}")
                    diffs.append(f"    ref[{i}]: {b[:120]}")
                if len(mismatches) > 11:
                    diffs.append("    ... (truncated)")
        else:
            if ln != lr:
                diffs.append(f"  DIFFER : {rel}")
                for i, (a, b) in enumerate(zip(ln, lr)):
                    if a != b:
                        diffs.append(f"    new[{i}]: {a[:120]}")
                        diffs.append(f"    ref[{i}]: {b[:120]}")
                        if i > 10:
                            diffs.append("    ... (truncated)")
                            break
                if len(ln) != len(lr):
                    diffs.append(f"    line count: new={len(ln)} ref={len(lr)}")

    return diffs, notices


def main():
    if len(sys.argv) != 3:
        print("Usage: compare_outputs.py <dir_new> <dir_ref>")
        sys.exit(2)

    dir_new = Path(sys.argv[1])
    dir_ref = Path(sys.argv[2])

    if not dir_new.is_dir():
        print(f"ERROR: {dir_new} is not a directory")
        sys.exit(2)
    if not dir_ref.is_dir():
        print(f"ERROR: {dir_ref} is not a directory")
        sys.exit(2)

    print(f"Comparing:\n  NEW: {dir_new}\n  REF: {dir_ref}\n")

    diffs, notices = compare_dirs(dir_new, dir_ref)

    if notices:
        for n in notices:
            print(n)

    if diffs:
        print(f"DIFFERENCES FOUND ({len(diffs)} issues):")
        for d in diffs:
            print(d)
        sys.exit(1)
    else:
        if notices:
            print("ALL OUTPUTS MATCH (modulo timestamps / float formatting / known v1.9.1 quirks).")
        else:
            print("ALL OUTPUTS MATCH (modulo timestamps / float formatting).")
        sys.exit(0)


if __name__ == "__main__":
    main()
