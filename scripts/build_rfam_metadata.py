#!/usr/bin/env python3
"""Builds Rfam_metadata.tsv for the pharokka database directory.

pharokka needs a small, fast lookup from Rfam accession to family name, type
and description so that cmscan hits can be annotated without parsing the
329 MB Rfam.cm flatfile at runtime.

Inputs (all from the Rfam FTP for the pinned release):
  family.txt.gz           - family table (accession, id, description, type)
  clan_membership.txt.gz  - clan accession -> family accession
  Rfam.cm                 - used only to cross-check that every model in the
                            database has a metadata row

Usage:
    python build_rfam_metadata.py \
        --family family.txt.gz \
        --clan-membership clan_membership.txt.gz \
        --cm Rfam.cm \
        --out Rfam_metadata.tsv
"""

import argparse
import gzip
import re
import sys

# 1-indexed columns of the Rfam family.txt MySQL dump that we care about
_F_ACC = 0
_F_ID = 1
_F_DESCRIPTION = 3
_F_GA = 6
_F_TYPE = 18
_F_MIN_FIELDS = 19

_ACC_RE = re.compile(r"^RF\d{5}\t")


def _open(path):
    return (
        gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path)
    )


def parse_family(path):
    """Parses family.txt.

    Some fields (comment, description) contain embedded newlines, so a family
    record is only started by a line beginning with an RF accession followed by
    a tab.  Continuation lines are ignored - every field we want precedes the
    free-text comment.
    """
    families = {}
    with _open(path) as fh:
        for line in fh:
            if not _ACC_RE.match(line):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < _F_MIN_FIELDS:
                continue
            families[fields[_F_ACC]] = {
                "rfam_id": fields[_F_ID],
                "description": fields[_F_DESCRIPTION],
                "ga_threshold": fields[_F_GA],
                "type": fields[_F_TYPE],
            }
    return families


def parse_clan_membership(path):
    """Parses clan_membership.txt -> {rfam_acc: clan_acc}."""
    clans = {}
    with _open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0].startswith("CL"):
                clans[fields[1]] = fields[0]
    return clans


def accessions_in_cm(path):
    """Every ACC present in the CM flatfile, to verify metadata completeness."""
    accessions = set()
    with _open(path) as fh:
        for line in fh:
            if line.startswith("ACC "):
                accessions.add(line.split(None, 1)[1].strip())
    return accessions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", required=True)
    parser.add_argument("--clan-membership", required=True)
    parser.add_argument("--cm", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    families = parse_family(args.family)
    clans = parse_clan_membership(args.clan_membership)
    cm_accessions = accessions_in_cm(args.cm)

    print(f"parsed {len(families)} families, {len(clans)} clan memberships")
    print(f"{len(cm_accessions)} models in {args.cm}")

    missing = cm_accessions - set(families)
    if missing:
        # not fatal - pharokka falls back to the cmscan description - but it
        # means family.txt and Rfam.cm are from different releases
        print(
            f"WARNING: {len(missing)} model(s) have no family.txt row, "
            f"e.g. {sorted(missing)[:5]}",
            file=sys.stderr,
        )

    with open(args.out, "w") as out:
        out.write("rfam_acc\trfam_id\ttype\tdescription\tga_threshold\tclan_acc\n")
        for acc in sorted(cm_accessions):
            info = families.get(acc, {})
            row = [
                acc,
                info.get("rfam_id", ""),
                info.get("type", ""),
                info.get("description", ""),
                info.get("ga_threshold", ""),
                clans.get(acc, ""),
            ]
            # guard against stray tabs/newlines in the free-text fields
            out.write("\t".join(f.replace("\t", " ").replace("\n", " ") for f in row))
            out.write("\n")

    print(f"wrote {len(cm_accessions)} rows to {args.out}")


if __name__ == "__main__":
    main()
