#!/usr/bin/env python3

import os
import sys
import time

from loguru import logger

from lib.databases  import check_db_installation
from lib.hmm import run_pyhmmer
from lib.input_commands import (
    get_input,
    check_dependencies,
    instantiate_dirs,
    instantiate_split_output,
    validate_fasta,
    validate_gene_predictor,
    validate_meta,
    validate_terminase,
    validate_threads,
)
from lib.post_processing import (
    create_gff,
    create_tbl,
    get_contig_name_lengths,
    parse_aragorn,
    process_results,
    create_gff_singles,
    convert_singles_gff_to_gbk,
    split_fasta_singles,
    write_tophits_vfdb_card,
    create_txt,
    update_fasta_headers,
    update_final_output,
    extract_terl,
    inphared_top_hits,
    remove_post_processing_files
)
from lib.processes import (
    concat_phanotate_meta,
    concat_trnascan_meta,
    reorient_terminase,
    run_aragorn,
    run_minced,
    run_mmseqs,
    run_mmseqs_card,
    run_mmseqs_vfdb,
    run_phanotate,
    run_phanotate_fasta_meta,
    run_phanotate_txt_meta,
    run_pyrodigal,
    run_trna_scan,
    run_trnascan_meta,
    split_input_fasta,
    translate_fastas,
    convert_gff_to_gbk,
    run_mash_sketch,
    run_mash_dist
)
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
                input_fasta, out_dir, args.threads, num_fastas
            )
            run_phanotate_txt_meta(
                input_fasta, out_dir, args.threads, num_fastas
            )
            concat_phanotate_meta(out_dir, num_fastas)
        else:
            run_phanotate(input_fasta, out_dir, logger)

    if gene_predictor == "prodigal":
        logger.info("Implementing Prodigal using Pyrodigal.")
        run_pyrodigal(
            input_fasta, out_dir, logger, args.meta, args.coding_table
        )

    # translate fastas
    logger.info("Translating gene predicted fastas.")
    translate_fastas(out_dir, gene_predictor, args.coding_table)

    # run trna-scan meta mode if required
    if args.meta == True:
        logger.info("Starting tRNA-scanSE. Applying meta mode.")
        run_trnascan_meta(input_fasta, out_dir, args.threads, num_fastas)
        concat_trnascan_meta(out_dir, num_fastas)
    else:
        print("Running tRNAscan-SE.")
        logger.info("Starting tRNA-scanSE")
        run_trna_scan(input_fasta, args.threads, out_dir, logger)
    # run minced and aragorn
    run_minced(input_fasta, out_dir, prefix, logger)
    run_aragorn(input_fasta, out_dir, prefix, logger)

    # running mmseqs2
    logger.info("Starting MMseqs2.")
    run_mmseqs(
        db_dir, out_dir, args.threads, logger, gene_predictor, args.evalue
    )
    run_mmseqs_card(db_dir, out_dir, args.threads, logger, gene_predictor)
    run_mmseqs_vfdb(db_dir, out_dir, args.threads, logger, gene_predictor)

    # running mmseqs2
    logger.info("Running pyhmmer.")
    run_pyhmmer(db_dir, out_dir, args.threads, logger, gene_predictor, args.evalue)

    # post processing
    logger.info("Post Processing Output.")

    # post process results
    # includes vfdb and card
    # return the merged df, vfdb and card top hits
    (cds_mmseqs_merge_df, vfdb_results, card_results) = process_results(
        db_dir, out_dir, prefix, gene_predictor
    )

    # gets df of length and gc for each contig
    length_df = get_contig_name_lengths(input_fasta)
    tmrna_flag = parse_aragorn(out_dir, length_df, prefix)

    # create gff and return locustag for table
    (locustag, locus_df, gff_df, total_gff) = create_gff(
        cds_mmseqs_merge_df,
        length_df,
        input_fasta,
        out_dir,
        prefix,
        locustag,
        tmrna_flag,
        args.meta,
    )
    create_tbl(
        cds_mmseqs_merge_df,
        length_df,
        out_dir,
        prefix,
        gene_predictor,
        tmrna_flag,
        gff_df,
        args.coding_table,
    )

    # output single gffs in meta mode
    if args.split == True and args.meta == True:
        create_gff_singles(length_df, input_fasta, out_dir, total_gff)
        convert_singles_gff_to_gbk(
            length_df, out_dir, args.coding_table
        )
        split_fasta_singles(input_fasta, out_dir)

    # write vfdb and card tophits
    # needs to be before .create_txt or else won't count properly
    write_tophits_vfdb_card(
        cds_mmseqs_merge_df, vfdb_results, card_results, locus_df, out_dir
    )

    # write the summary tsv outputs
    create_txt(cds_mmseqs_merge_df, length_df, out_dir, prefix)

    # convert to genbank
    logger.info("Converting gff to genbank.")
    convert_gff_to_gbk(
        input_fasta, out_dir, out_dir, prefix, args.coding_table
    )

    # update fasta headers and final output tsv
    update_fasta_headers(locus_df, out_dir, gene_predictor)
    update_final_output(cds_mmseqs_merge_df, locus_df, prefix, out_dir)

    # extract terL
    extract_terl(locus_df, out_dir, gene_predictor, logger)

    # run mash
    logger.info("Finding the closest match for each contig in INPHARED using mash.")
    run_mash_sketch(input_fasta, out_dir, logger)
    run_mash_dist(out_dir, db_dir, logger)
    inphared_top_hits(out_dir, db_dir, length_df, prefix)

    # delete tmp files
    remove_post_processing_files(out_dir, gene_predictor, args.meta)

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

if __name__ == '__main__':
    main()