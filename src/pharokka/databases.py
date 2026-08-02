"""
to tar DBs
export GZIP=-9

To remake the MMseqs2 PHROG profile DB with MMseqs2 v18
wget https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz
tar -xzf MSA_phrogs.tar.gz
# to rename for the headers
for f in MSA_Phrogs_M50_FASTA/*.fma; do     mv "$f" "${f%.fma}"; done
tar -czf MSA_phrogs_renamed.tar.gz MSA_Phrogs_M50_FASTA/
mmseqs tar2db  MSA_phrogs_renamed.tar.gz  msa --output-dbtype 11
mmseqs msa2profile msa phrogs_profile_db


tar cvzf pharokka_v1.8.0_databases.tar.gz pharokka_v1.8.0_databases

to create the mmseqs2 databases (CARD)
with the latest download v 3.2.7 on August 21 2023
mmseqs createdb protein_fasta_protein_homolog_model.fasta CARD

# also need aro.tsv for this one

# for VFDB, only need the FASTA

# VFDB update as of August 18 2023 (not versioned)
# clustered

mmseqs easy-cluster VFDB_setB_pro_form.fas VFDBclusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1
mmseqs createdb VFDBclusterRes_rep_seq.fasta vfdb

"""

import hashlib
import os
import shutil
import tarfile
from pathlib import Path

import requests
from alive_progress import alive_bar
from loguru import logger

from .util import remove_directory

VERSION = "1.11.0"

# Rfam release pinned into the pharokka database tarball.  Bump this (and the
# database VERSION above) only when a new Rfam is packaged and uploaded.
RFAM_VERSION = "15.1"

# to hold information about the different DBs
VERSION_DICTIONARY = {
    "1.2.0": {
        "md5": "0014b7a982dbf071f8856a5a29a95e30",
        "major": 1,
        "minor": 2,
        "minorest": 0,
        "db_url": "https://zenodo.org/record/7563578/files/pharokka_v1.2.0_database.tar.gz",
        "dir_name": "pharokka_v1.2.0_databases",
        "inphared_mash": "5Jan2023_genomes.fa.msh",
        "inphared_annot": "5Jan2023_data.tsv",
    },
    "1.4.0": {
        "md5": "c21144209b993c06fae2dac906d73b96",
        "major": 1,
        "minor": 4,
        "minorest": 0,
        "db_url": "https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz",
        "dir_name": "pharokka_v1.4.0_databases",
        "inphared_mash": "1Aug2023_genomes.fa.msh",
        "inphared_annot": "1Aug2023_data.tsv",
    },
    "1.8.0": {
        "md5": "a63c485241b900a11989bd1821bfbb09",
        "major": 1,
        "minor": 8,
        "minorest": 0,
        "db_url": "https://zenodo.org/record/17110353/files/pharokka_v1.8.0_databases.tar.gz",
        "dir_name": "pharokka_v1.8.0_databases",
        "inphared_mash": "9Aug2025_genomes.fa.msh",
        "inphared_annot": "9Aug2025_data.tsv",
    },
    # v1.11.0 adds the Rfam 15.1 covariance models (cmpress'd) for --rfam
    # ncRNA annotation.  See scripts/build_rfam_db.sh for how they are built.
    # TODO(gbouras13): fill in db_url + md5 once the tarball is uploaded to Zenodo.
    "1.11.0": {
        "md5": "TODO_MD5_AFTER_UPLOAD",
        "major": 1,
        "minor": 11,
        "minorest": 0,
        "db_url": "https://zenodo.org/record/TODO/files/pharokka_v1.11.0_databases.tar.gz",
        "dir_name": "pharokka_v1.11.0_databases",
        "inphared_mash": "9Aug2025_genomes.fa.msh",
        "inphared_annot": "9Aug2025_data.tsv",
    },
}


PHROG_DB_NAMES = [
    # version marker file shipped inside the database tarball; derived so it
    # tracks VERSION rather than silently checking for a stale marker
    f"VERSION_{VERSION.replace('.', '_')}",
    "phrogs_profile_db",
    "phrogs_profile_db.dbtype",
    "phrogs_profile_db.index",
    "phrogs_profile_db.lookup",
    "phrogs_profile_db.source",
    "phrogs_profile_db_h",
    "phrogs_profile_db_h.dbtype",
    "phrogs_profile_db_h.index",
]

PHROG_HMM_NAMES = ["all_phrogs.h3m"]

VFDB_DB_NAMES = [
    "vfdb",
    "vfdb.dbtype",
    "vfdb.index",
    "vfdb.lookup",
    "vfdb.source",
    "vfdb_h",
    "vfdb_h.dbtype",
    "vfdb_h.index",
    "VFDBclusterRes_cluster.tsv",
    "VFDBclusterRes_rep_seq.fasta",
]

CARD_DB_NAMES = [
    "CARD",
    "CARD.dbtype",
    "CARD.index",
    "CARD.lookup",
    "CARD.source",
    "CARD_h",
    "CARD_h.dbtype",
    "CARD_h.index",
]

# Rfam covariance models, cmpress'd.  Only required when --rfam is used, so
# these are checked separately (check_rfam_installation) rather than in
# check_db_installation, which gates the whole-database download.
RFAM_DB_NAMES = [
    "Rfam.cm",
    "Rfam.cm.i1f",
    "Rfam.cm.i1i",
    "Rfam.cm.i1m",
    "Rfam.cm.i1p",
    "Rfam.clanin",
    "Rfam_metadata.tsv",
]


def instantiate_install(db_dir):
    instantiate_dir(db_dir)
    # check the database is installed
    logger.info(f"Checking Pharokka database installation in {db_dir}.")
    downloaded_flag = check_db_installation(db_dir)
    if downloaded_flag:
        logger.info("All Pharokka Databases have already been Downloaded and Checked.")
    else:
        logger.info("Some Databases are missing.")

        db_url = VERSION_DICTIONARY[VERSION]["db_url"]
        requiredmd5 = VERSION_DICTIONARY[VERSION]["md5"]

        logger.info(f"Downloading Pharokka Databases from {db_url}.")

        # derived from VERSION, not hardcoded - otherwise every database
        # version bump silently keeps writing the old tarball name
        tarball = f"{VERSION_DICTIONARY[VERSION]['dir_name']}.tar.gz"
        tarball_path = Path(f"{db_dir}/{tarball}")

        download(db_url, tarball_path)

        md5_sum = calc_md5_sum(tarball_path)

        if md5_sum == requiredmd5:
            logger.info(f"Database file download OK: {md5_sum}")
        else:
            logger.error(
                f"Error: corrupt database file! MD5 should be '{requiredmd5}' but is '{md5_sum}'"
            )

        logger.info(f"Extracting DB tarball: file={tarball_path}, output={db_dir}")
        untar(tarball_path, db_dir)
        tarball_path.unlink()


"""
lots of this code from the marvellous bakta https://github.com/oschwengers/bakta, db.py specifically
"""


def download(db_url: str, tarball_path: Path):
    headers = {
        "User-Agent": f"pharokka/{VERSION} (contact: george.bouras@adelaide.edu.au)"
    }

    try:
        with tarball_path.open("wb") as fh_out, requests.get(
            db_url, stream=True, headers=headers
        ) as resp:
            total_length = resp.headers.get("content-length")
            if total_length is not None:
                total_length = int(total_length)
            with alive_bar(total=total_length, scale="SI") as bar:
                for data in resp.iter_content(chunk_size=1024 * 1024):
                    fh_out.write(data)
                    bar(count=len(data))
    except IOError:
        logger.error(
            f"ERROR: Could not download file from Zenodo! url={db_url}, path={tarball_path}"
        )


def calc_md5_sum(tarball_path: Path, buffer_size: int = 1024 * 1024) -> str:
    md5 = hashlib.md5()
    with tarball_path.open("rb") as fh:
        data = fh.read(buffer_size)
        while data:
            md5.update(data)
            data = fh.read(buffer_size)
    return md5.hexdigest()


def untar(tarball_path: Path, output_path: Path):
    try:
        with tarball_path.open("rb") as fh_in, tarfile.open(
            fileobj=fh_in, mode="r:gz"
        ) as tar_file:
            try:
                tar_file.extractall(path=str(output_path), filter="data")
            except TypeError:
                tar_file.extractall(path=str(output_path))

        tarpath = os.path.join(output_path, VERSION_DICTIONARY[VERSION]["dir_name"])

        files_to_move = [
            f for f in os.listdir(tarpath) if os.path.isfile(os.path.join(tarpath, f))
        ]

        for file_name in files_to_move:
            source_path = os.path.join(tarpath, file_name)
            destination_path = os.path.join(output_path, file_name)
            shutil.move(source_path, destination_path)
        remove_directory(tarpath)

    except OSError:
        logger.error(f"Could not extract {tarball_path} to {output_path}")


def instantiate_dir(db_dir):
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)


def check_db_installation(db_dir):
    downloaded_flag = True
    # PHROGS files
    for file_name in PHROG_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if not os.path.isfile(path):
            logger.info(f"PHROGs Database file {file_name} is missing.")
            downloaded_flag = False
            break
    # VFDB
    for file_name in VFDB_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if not os.path.isfile(path):
            logger.info("VFDB Databases are missing.")
            downloaded_flag = False
            break
    # CARD
    for file_name in CARD_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if not os.path.isfile(path):
            logger.info("CARD Databases are missing.")
            downloaded_flag = False
            break
    # annot.tsv
    path = os.path.join(db_dir, "phrog_annot_v4.tsv")
    if not os.path.isfile(path):
        logger.info("PHROGs Annotation File is missing.")
        downloaded_flag = False

    # mash files
    path = os.path.join(db_dir, VERSION_DICTIONARY[VERSION]["inphared_annot"])
    if not os.path.isfile(path):
        logger.info("INPHARED Mash Annotation File is missing.")
        downloaded_flag = False

    path = os.path.join(db_dir, VERSION_DICTIONARY[VERSION]["inphared_mash"])
    if not os.path.isfile(path):
        logger.info("INPHARED Mash Sketch File is missing.")
        downloaded_flag = False

    return downloaded_flag


def check_rfam_installation(db_dir):
    """Checks that the cmpress'd Rfam database is present.

    Only called when --rfam is specified.  Rfam ships inside the pharokka
    database tarball from v1.11.0 onwards, so a user with an older database
    directory will have everything else but not these files - hence the
    explicit, actionable error rather than a bare missing-file traceback.

    :param db_dir: pharokka database directory
    :return: True if every Rfam file is present
    """
    missing = [
        file_name
        for file_name in RFAM_DB_NAMES
        if not os.path.isfile(os.path.join(db_dir, file_name))
    ]

    if missing:
        logger.error(
            f"--rfam was specified but the Rfam database is missing from {db_dir}."
        )
        logger.error(f"Missing file(s): {', '.join(missing)}")
        logger.error(
            "Rfam was added to the pharokka database in v1.11.0. Please re-run "
            "'pharokka install' to download the updated database."
        )
        return False

    return True
