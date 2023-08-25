"""
to tar DBs
export GZIP=-9
tar cvzf pharokka_v1.4.0_databases.tar.gz pharokka_v1.4.0_databases

to create the mmseqs2 databases (CARD)
with the latest download v 3.2.7 on August 21
mmseqs createdb protein_fasta_protein_homolog_model.fasta CARD

# also need aro.tsv for this one

# for VFDB, only need the FASTA

# VFDB update as of August 18 2023 (not versioned)
# clustered 

mmseqs easy-cluster VFDB_setB_pro_form.fas VFDBclusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1
mmseqs createdb VFDBclusterRes_rep_seq.fasta vfdb

"""

#!/usr/bin/env python3
import hashlib
import os
import shutil
import tarfile
from pathlib import Path

import requests
from alive_progress import alive_bar
from loguru import logger
from post_processing import remove_directory

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
    }
}

VERSION_DICTIONARY = {
    "1.4.0": {
        "md5": "c21144209b993c06fae2dac906d73b96",
        "major": 1,
        "minor": 4,
        "minorest": 0,
        "db_url": "https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz",
        "dir_name": "pharokka_v1.4.0_databases",
        "inphared_mash": "1Aug2023_genomes.fa.msh",
        "inphared_annot": "1Aug2023_data.tsv",
    }
}


PHROG_DB_NAMES = [
    "phrogs_db",
    "phrogs_db.dbtype",
    "phrogs_db.index",
    "phrogs_profile_db",
    "phrogs_profile_db.dbtype",
    "phrogs_profile_db.index",
    "phrogs_profile_db_consensus",
    "phrogs_profile_db_consensus.dbtype",
    "phrogs_profile_db_consensus.index",
    "phrogs_profile_db_h",
    "phrogs_profile_db_h.index",
    "phrogs_profile_db_seq",
    "phrogs_profile_db_seq.dbtype",
    "phrogs_profile_db_seq.index",
    "phrogs_profile_db_seq_h",
    "phrogs_profile_db_seq_h.index",
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


def instantiate_install(db_dir):
    instantiate_dir(db_dir)
    # check the database is installed
    logger.info(f"Checking Pharokka database installation in {db_dir}.")
    downloaded_flag = check_db_installation(db_dir)
    if downloaded_flag == True:
        logger.info("All Pharokka Databases have already been Downloaded and Checked.")
    else:
        logger.info("Some Databases are missing.")

        db_url = VERSION_DICTIONARY["1.4.0"]["db_url"]
        requiredmd5 = VERSION_DICTIONARY["1.4.0"]["md5"]

        logger.info(f"Downloading Pharokka Databases from {db_url}.")

        tarball_path = Path(f"{db_dir}/pharokka_v1.4.0_databases.tar.gz")

        download(db_url, tarball_path)

        # get_database_zenodo(db_dir)

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
    try:
        with tarball_path.open("wb") as fh_out, requests.get(
            db_url, stream=True
        ) as resp:
            total_length = resp.headers.get("content-length")
            if total_length is not None:  # content length header is set
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
            tar_file.extractall(path=str(output_path))

        tarpath = os.path.join(output_path, VERSION_DICTIONARY["1.4.0"]["dir_name"])

        # Get a list of all files in the source directory
        files_to_move = [
            f for f in os.listdir(tarpath) if os.path.isfile(os.path.join(tarpath, f))
        ]

        # Move each file to the destination directory
        for file_name in files_to_move:
            source_path = os.path.join(tarpath, file_name)
            destination_path = os.path.join(output_path, file_name)
            shutil.move(source_path, destination_path)
        # remove the directort
        remove_directory(tarpath)

    except OSError:
        logger.error(f"Could not extract {tarball_path} to {output_path}")


def instantiate_dir(db_dir):
    if os.path.isdir(db_dir) == False:
        os.mkdir(db_dir)


def check_db_installation(db_dir):
    downloaded_flag = True
    # PHROGS files
    for file_name in PHROG_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if os.path.isfile(path) == False:
            logger.info("PHROGs Databases are missing.")
            downloaded_flag = False
            break
    # VFDB
    for file_name in VFDB_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if os.path.isfile(path) == False:
            logger.info("VFDB Databases are missing.")
            downloaded_flag = False
            break
    # CARD
    for file_name in CARD_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if os.path.isfile(path) == False:
            logger.info("CARD Databases are missing.")
            downloaded_flag = False
            break
    # annot.tsv
    path = os.path.join(db_dir, "phrog_annot_v4.tsv")
    if os.path.isfile(path) == False:
        logger.info("PHROGs Annotation File is missing.")
        downloaded_flag = False

    # mash files
    path = os.path.join(db_dir, VERSION_DICTIONARY["1.4.0"]["inphared_annot"])
    if os.path.isfile(path) == False:
        logger.info("INPHARED Mash Annotation File is missing.")
        downloaded_flag = False

    # mash files
    path = os.path.join(db_dir, VERSION_DICTIONARY["1.4.0"]["inphared_mash"])
    if os.path.isfile(path) == False:
        logger.info("INPHARED Mash Sketch File is missing.")
        downloaded_flag = False

    return downloaded_flag
