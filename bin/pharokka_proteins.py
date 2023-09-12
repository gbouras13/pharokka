#!/usr/bin/env python3

import os
import sys
import time
from pathlib import Path

from databases import check_db_installation
from input_commands import (check_dependencies, instantiate_dirs,
                            validate_fasta, validate_threads)
from loguru import logger
from post_processing import remove_directory, remove_file
from proteins import (Pharok_Prot, get_input_proteins, run_mmseqs_proteins,
                      run_pyhmmer_proteins)
from util import get_version


def main():
    # get the args
    args = get_input_proteins()

    logger.add(lambda _: sys.exit(1), level="ERROR")

    if args.citation == True:
        logger.info("If you use Pharokka in your research, please cite:")
        logger.info(
            "George Bouras, Roshan Nepal, Ghais Houtak, Alkis James Psaltis, Peter-John Wormald, Sarah Vreugde."
        )
        logger.info("Pharokka: a fast scalable bacteriophage annotation tool.")
        logger.info("Bioinformatics, Volume 39, Issue 1, January 2023, btac776.")
        logger.info("https://doi.org/10.1093/bioinformatics/btac776.")
        logger.info(
            "You should also cite the full list of tools Pharokka uses, which can be found at https://github.com/gbouras13/pharokka#citation."
        )
        sys.exit()

    # set the prefix
    if args.prefix == "Default":
        prefix = "pharokka_proteins"
    else:
        prefix = args.prefix

    # instantiate outdir
    out_dir = instantiate_dirs(args.outdir, False, args.force)  # args.meta always false

    # set log dir
    logdir = Path(f"{out_dir}/logs")

    # set the db dir
    if args.database == "Default":
        db_dir = os.path.join(os.path.realpath(__file__), "../", "databases/")
    else:
        db_dir = args.database

    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(args.outdir, f"pharokka_proteins_{start_time}.log")
    # adds log file
    logger.add(log_file)

    # preamble
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Starting Pharokka v{get_version()}")
    logger.info("Running pharokka_proteins.py to annotate proteins.")
    logger.info("Command executed: {}", args)
    logger.info("Repository homepage is https://github.com/gbouras13/pharokka")
    logger.info("Written by George Bouras: george.bouras@adelaide.edu.au")

    # check the database is installed
    logger.info("Checking database installation.")

    database_installed = check_db_installation(db_dir)
    if database_installed == True:
        logger.info("All databases have been successfully checked.")
    else:
        logger.error(
            "The database directory was unsuccessfully checked. Please run install_databases.py"
        )

    # dependencies
    logger.info("Checking dependencies.")
    check_dependencies()

    # instantiation/checking fasta and gene_predictor
    validate_fasta(args.infile)
    validate_threads(args.threads)

    # define input - overwrite if terminase reorienter is true
    input_fasta = args.infile

    ########
    # mmseqs2 and hmm decisions
    ########

    # can't have fast and mmseqs2 only
    if args.hmm_only is True and args.mmseqs2_only is True:
        logger.error(
            "You have specified --fast or --hmm_only and --mmseqs2_only. This is impossible. Please choose one or the other."
        )

    # by default
    mmseqs_flag = True
    hmm_flag = True

    if args.hmm_only is True:
        logger.info("You have specified --hmm_only. MMseqs2 will not be run.")
        mmseqs_flag = False
    if args.mmseqs2_only is True:
        logger.info("You have specified --mmseqs2_only. Pyhmmer will not be run.")
        hmm_flag = False

    # running mmseqs2 on the 3 databases
    if mmseqs_flag is True:
        logger.info("Starting MMseqs2.")
        run_mmseqs_proteins(
            input_fasta,
            db_dir,
            out_dir,
            args.threads,
            logdir,
            args.evalue,
            db_name="PHROG",
        )
        run_mmseqs_proteins(
            input_fasta,
            db_dir,
            out_dir,
            args.threads,
            logdir,
            args.evalue,
            db_name="CARD",
        )
        run_mmseqs_proteins(
            input_fasta,
            db_dir,
            out_dir,
            args.threads,
            logdir,
            args.evalue,
            db_name="VFDB",
        )

    if hmm_flag is True:
        # runs pyhmmer on PHROGs
        logger.info("Running PyHMMER.")
        best_results_pyhmmer = run_pyhmmer_proteins(
            input_fasta, db_dir, args.threads, args.evalue
        )

    #################################################
    # post processing
    #################################################

    logger.info("Post Processing Output.")

    # instanatiate the class with some of the params
    pharok = Pharok_Prot()
    pharok.out_dir = out_dir
    pharok.db_dir = db_dir
    pharok.prefix = prefix
    pharok.input_fasta = input_fasta
    pharok.mmseqs_flag = mmseqs_flag
    pharok.hmm_flag = hmm_flag
    if pharok.hmm_flag is True:
        pharok.pyhmmer_results_dict = best_results_pyhmmer

    # post process results
    # includes vfdb and card
    # adds the merged df, vfdb and card top hits dfs to the class objec
    # no need to specify params as they are in the class :)
    pharok.process_dataframes()

    # updates fasta headers
    pharok.update_fasta_headers()

    # cleanup

    remove_directory(os.path.join(out_dir, "mmseqs"))
    remove_file(os.path.join(out_dir, "CARD_results.tsv"))
    remove_file(os.path.join(out_dir, "vfdb_results.tsv"))
    remove_directory(os.path.join(out_dir, "CARD"))
    remove_directory(os.path.join(out_dir, "vfdb"))
    remove_directory(os.path.join(out_dir, "VFDB_dir"))
    remove_directory(os.path.join(out_dir, "tmp_dir"))
    remove_file(os.path.join(out_dir, "mmseqs_results.tsv"))
    remove_file(os.path.join(out_dir, "CARD_results.tsv"))
    remove_file(os.path.join(out_dir, "vfdb_results.tsv"))
    remove_directory(os.path.join(out_dir, "CARD"))
    remove_directory(os.path.join(out_dir, "vfdb"))

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("Pharokka has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")

    logger.info("If you use Pharokka in your research, please cite:")
    logger.info(
        "George Bouras, Roshan Nepal, Ghais Houtak, Alkis James Psaltis, Peter-John Wormald, Sarah Vreugde."
    )
    logger.info("Pharokka: a fast scalable bacteriophage annotation tool.")
    logger.info("Bioinformatics, Volume 39, Issue 1, January 2023, btac776.")
    logger.info("https://doi.org/10.1093/bioinformatics/btac776.")
    logger.info(
        "You should also cite the full list of tools Pharokka uses, which can be found at https://github.com/gbouras13/pharokka#citation."
    )


if __name__ == "__main__":
    main()
