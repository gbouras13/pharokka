#!/bin/bash
# Builds the Rfam component of the pharokka database.
#
# The Rfam release is PINNED - pharokka ships a fixed Rfam version inside the
# database tarball so that annotations are reproducible across runs and
# machines.  To move to a new Rfam release you must also bump RFAM_VERSION and
# the database VERSION in src/pharokka/databases.py, and upload a new tarball.
#
# Requires: Infernal (cmpress) on $PATH.
#
# Output files, which belong at the top level of the pharokka database dir:
#   Rfam.cm.i1f  Rfam.cm.i1i  Rfam.cm.i1m  Rfam.cm.i1p
#   Rfam.clanin
#   Rfam_metadata.tsv
#
# The Rfam.cm flatfile is deleted after pressing: cmscan reads the .i1* files
# and only uses the flatfile path as a base name, so shipping it would add
# 329 MB to every user's database for nothing.

set -euo pipefail

RFAM_VERSION="15.1"
RFAM_FTP="https://ftp.ebi.ac.uk/pub/databases/Rfam/${RFAM_VERSION}"
OUTDIR="${1:-rfam_db}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "==> Downloading Rfam ${RFAM_VERSION}"
curl -sSL -O "${RFAM_FTP}/Rfam.cm.gz"
curl -sSL -O "${RFAM_FTP}/Rfam.clanin"
curl -sSL -O "${RFAM_FTP}/database_files/family.txt.gz"
curl -sSL -O "${RFAM_FTP}/database_files/clan_membership.txt.gz"

echo "==> Decompressing covariance models"
gunzip -f Rfam.cm.gz

echo "==> Pressing covariance models (this takes a few minutes)"
cmpress -F Rfam.cm

echo "==> Building metadata table"
python "${SCRIPT_DIR}/build_rfam_metadata.py" \
    --family family.txt.gz \
    --clan-membership clan_membership.txt.gz \
    --cm Rfam.cm \
    --out Rfam_metadata.tsv

echo "==> Cleaning up intermediates"
# Rfam.cm is only needed to build the pressed files and the metadata table.
# cmscan reads the .i1* files, so the 329 MB flatfile is not shipped.
rm -f family.txt.gz clan_membership.txt.gz Rfam.cm

echo "==> Done. Files in ${OUTDIR}:"
ls -la Rfam.cm.i1? Rfam.clanin Rfam_metadata.tsv
echo
echo "Copy these into the pharokka database directory before creating the"
echo "database tarball, then update db_url and md5 in src/pharokka/databases.py."
