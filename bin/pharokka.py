#!/usr/bin/env python3

import os
import sys
import time
from pathlib import Path

from loguru import logger

from lib.databases import check_db_installation
from lib.hmm import run_pyhmmer
from lib.input_commands import (check_dependencies, get_input,
                                instantiate_dirs, instantiate_split_output,
                                validate_fasta, validate_gene_predictor,
                                validate_meta, validate_terminase,
                                validate_threads)
from lib.post_processing import Pharok, remove_post_processing_files
from lib.processes import (concat_phanotate_meta, concat_trnascan_meta,
                           convert_gff_to_gbk, reorient_terminase, run_aragorn,
                           run_mash_dist, run_mash_sketch, run_minced,
                           run_mmseqs, run_phanotate, run_phanotate_fasta_meta,
                           run_phanotate_txt_meta, run_pyrodigal,
                           run_trna_scan, run_trnascan_meta, split_input_fasta,
                           translate_fastas)
from lib.util import get_version


def main():
    # get the args
    args = get_input()

    if args.citation == True:
        logger.info("If you use pharokka in your research, please cite:")
        logger.info(
            "George Bouras, Roshan Nepal, Ghais Houtak, Alkis James Psaltis, Peter-John Wormald, Sarah Vreugde."
        )
        logger.info("Pharokka: a fast scalable bacteriophage annotation tool.")
        logger.info("Bioinformatics, Volume 39, Issue 1, January 2023, btac776.")
        logger.info("https://doi.org/10.1093/bioinformatics/btac776.")
        logger.info(
            "You should also cite the full list of tools pharokka uses, which can be found at https://github.com/gbouras13/pharokka#citation."
        )
        sys.exit()

    # set the prefix
    if args.prefix == "Default":
        prefix = "pharokka"
    else:
        prefix = args.prefix

    # set the locustag flag for locustag generation
    if args.locustag == "Default":
        locustag = "Random"
    else:
        locustag = args.locustag

    # set the gene_predictor
    gene_predictor = args.gene_predictor

    # instantiate outdir
    out_dir = instantiate_dirs(args.outdir, args.meta, args.force)

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
    log_file = os.path.join(args.outdir, f"{prefix}_pharokka_{start_time}.log")
    # adds log file
    logger.add(log_file)

    # preamble
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Starting pharokka v{get_version()}")
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
            "\nThe database directory was unsuccessfully checked. Please run install_databases.py \n"
        )

    # dependencies
    logger.info("Checking dependencies.")
    check_dependencies(logger)

    # instantiation/checking fasta and gene_predictor
    validate_fasta(args.infile)
    validate_gene_predictor(gene_predictor)
    validate_threads(args.threads)

    # define input - overwrite if terminase reorienter is true
    input_fasta = args.infile

    # terminase reorienting
    if args.terminase == True:
        logger.info(
            "You have chosen to reorient your genome by specifying --terminase. Checking the input."
        )
        validate_terminase(args.infile, args.terminase_strand, args.terminase_start)
        reorient_terminase(
            args.infile,
            out_dir,
            prefix,
            args.terminase_strand,
            args.terminase_start,
            logger,
        )
        input_fasta = os.path.join(
            out_dir, prefix + "_genome_terminase_reoriented.fasta"
        )

    if args.terminase_strand != "nothing" and args.terminase == False:
        logger.info(
            "You have specified a terminase strand using --strand to reorient your genome, but you have not specified --terminase to activate terminase mode. \nContinuing without reorientation."
        )

    if args.terminase_strand != "nothing" and args.terminase_start == False:
        logger.info(
            "You have specified a terminase start coordinate using --terminase_start to reorient your genome, but you have not specified --terminase to activate terminase mode. \nContinuing without reorientation."
        )

    # validates meta mode
    validate_meta(input_fasta, args.meta, args.split, logger)

    # meta mode split input for trnascan and maybe phanotate
    if args.meta == True:
        num_fastas = split_input_fasta(input_fasta, out_dir)
        # will generate split output gffs if meta flag is true
        instantiate_split_output(out_dir, args.meta)

    # CDS predicton
    # phanotate
    if gene_predictor == "phanotate":
        logger.info("Running Phanotate.")
        if args.meta == True:
            logger.info("Applying meta mode.")
            run_phanotate_fasta_meta(
                input_fasta, out_dir, args.threads, num_fastas, logdir
            )
            run_phanotate_txt_meta(
                input_fasta, out_dir, args.threads, num_fastas, logdir
            )
            concat_phanotate_meta(out_dir, num_fastas)
        else:
            run_phanotate(input_fasta, out_dir, logdir)

    if gene_predictor == "prodigal":
        logger.info("Implementing Prodigal using Pyrodigal.")
        run_pyrodigal(input_fasta, out_dir, logger, args.meta, args.coding_table)

    # translate fastas
    logger.info("Translating gene predicted fastas.")
    translate_fastas(out_dir, gene_predictor, args.coding_table)

    # run trna-scan meta mode if required
    if args.meta == True:
        logger.info("Starting tRNA-scanSE. Applying meta mode.")
        run_trnascan_meta(input_fasta, out_dir, args.threads, num_fastas)
        concat_trnascan_meta(out_dir, num_fastas)
    else:
        logger.info("Starting tRNA-scanSE")
        run_trna_scan(input_fasta, args.threads, out_dir, logdir)
    # run minced and aragorn
    run_minced(input_fasta, out_dir, prefix, logdir)
    run_aragorn(input_fasta, out_dir, prefix, logdir)

    # running mmseqs2 on the 3 databases
    logger.info("Starting MMseqs2.")
    run_mmseqs(
        db_dir,
        out_dir,
        args.threads,
        logdir,
        gene_predictor,
        args.evalue,
        db_name="PHROG",
    )
    run_mmseqs(
        db_dir,
        out_dir,
        args.threads,
        logdir,
        gene_predictor,
        args.evalue,
        db_name="CARD",
    )
    run_mmseqs(
        db_dir,
        out_dir,
        args.threads,
        logdir,
        gene_predictor,
        args.evalue,
        db_name="VFDB",
    )

    # runs pyhmmer on PHROGs
    logger.info("Running pyhmmer.")
    best_results_pyhmmer = run_pyhmmer(
        db_dir, out_dir, args.threads, gene_predictor, args.evalue
    )

    print(best_results_pyhmmer)

    #################################################
    # post processing
    #################################################

    logger.info("Post Processing Output.")

    # instanatiate the class with some of the params
    pharok = Pharok()
    pharok.out_dir = out_dir
    pharok.db_dir = db_dir
    pharok.gene_predictor = gene_predictor
    pharok.prefix = prefix
    pharok.locustag = locustag
    pharok.input_fasta = input_fasta
    pharok.meta_mode = args.meta
    pharok.pyhmmer_results_dict = best_results_pyhmmer
    pharok.coding_table = args.coding_table

    # post process results
    # includes vfdb and card
    # adds the merged df, vfdb and card top hits dfs to the class objec
    # no need to specify params as they are in the class :)
    pharok.process_results()

    # gets df of length and gc for each contig
    pharok.get_contig_name_lengths()

    # parse the aragorn output
    # get flag whether there is a tmrna from aragor
    pharok.parse_aragorn()

    # create gff and save locustag to class for table
    pharok.create_gff()

    # create table
    pharok.create_tbl()

    # output single gffs in meta mode
    if args.split == True and args.meta == True:
        # create gffs for each contig
        pharok.create_gff_singles()
        # converts each gff to gbk
        pharok.convert_singles_gff_to_gbk()
        # splits the input fasta into single fastas
        pharok.split_fasta_singles()

    # write vfdb and card tophits
    # needs to be before .create_txt or else won't count properly
    pharok.write_tophits_vfdb_card()

    # write the summary tsv outputs
    pharok.create_txt()

    # convert to genbank
    logger.info("Converting gff to genbank.")
    # not part of the class so from processes.py
    convert_gff_to_gbk(input_fasta, out_dir, out_dir, prefix, args.coding_table)

    # update fasta headers and final output tsv
    pharok.update_fasta_headers()
    pharok.update_final_output()


    # extract terL
    pharok.extract_terl()

    # run mash
    logger.info("Finding the closest match for each contig in INPHARED using mash.")
    # in process.py
    run_mash_sketch(input_fasta, out_dir, logdir)
    run_mash_dist(out_dir, db_dir, logdir)
    # part of the class
    pharok.inphared_top_hits()

    # delete tmp files
    # remove_post_processing_files(out_dir, gene_predictor, args.meta)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("pharokka has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")

    logger.info("If you use pharokka in your research, please cite:")
    logger.info(
        "George Bouras, Roshan Nepal, Ghais Houtak, Alkis James Psaltis, Peter-John Wormald, Sarah Vreugde."
    )
    logger.info("Pharokka: a fast scalable bacteriophage annotation tool.")
    logger.info("Bioinformatics, Volume 39, Issue 1, January 2023, btac776.")
    logger.info("https://doi.org/10.1093/bioinformatics/btac776.")
    logger.info(
        "You should also cite the full list of tools pharokka uses, which can be found at https://github.com/gbouras13/pharokka#citation."
    )


if __name__ == "__main__":
    main()
