#!/usr/bin/env python3

import os
import shutil
import sys
import time
import warnings
from pathlib import Path

from Bio import BiopythonDeprecationWarning
from custom_db import run_custom_pyhmmer
from databases import check_db_installation
from hmm import run_pyhmmer
from input_commands import (check_dependencies, get_input, instantiate_dirs,
                            instantiate_split_output,
                            validate_and_extract_genbank, validate_custom_hmm,
                            validate_fasta, validate_gene_predictor,
                            validate_meta, validate_terminase,
                            validate_threads)
from loguru import logger
from post_processing import Pharok, remove_post_processing_files
from processes import (concat_phanotate_meta, concat_trnascan_meta,
                       convert_gff_to_gbk, reorient_terminase, run_aragorn,
                       run_dnaapler, run_mash_dist, run_mash_sketch,
                       run_minced, run_mmseqs, run_phanotate,
                       run_phanotate_fasta_meta, run_phanotate_txt_meta,
                       run_pyrodigal, run_pyrodigal_gv, run_trna_scan,
                       run_trnascan_meta, split_input_fasta, translate_fastas)
from util import count_contigs, get_version

# add this to make sure of deprecation warning with biopython
warnings.simplefilter("default", BiopythonDeprecationWarning)


def main():
    # get the args
    args = get_input()

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
    # set to phanotate by default and prodigal-gv in meta mode
    if gene_predictor == "default":
        if args.meta is True:
            gene_predictor = "prodigal-gv"
        else:
            gene_predictor = "phanotate"

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
    log_file = os.path.join(args.outdir, f"pharokka_{start_time}.log")
    # adds log file
    logger.add(log_file)

    # preamble
    logger.info(f"Starting Pharokka v{get_version()}")
    logger.info("Command executed: {}", args)
    logger.info("Repository homepage is https://github.com/gbouras13/pharokka")
    logger.info("Written by George Bouras: george.bouras@adelaide.edu.au")

    logger.info(f"Checking database installation in {db_dir}.")
    database_installed = check_db_installation(db_dir)
    if database_installed == True:
        logger.info("All databases have been successfully checked.")
    else:
        logger.error(
            "The database directory was unsuccessfully checked. Please run install_databases.py."
        )

    ### custom hmm

    custom_hmm_flag = False
    if args.custom_hmm != "":
        custom_hmm_flag = True

    if custom_hmm_flag is True:
        validate_custom_hmm(args.custom_hmm)

    # dependencies
    logger.info("Checking dependencies.")
    # output versions
    (
        phanotate_version,
        pyrodigal_version,
        pyrodigal_gv_version,
        trna_version,
        aragorn_version,
        minced_version,
    ) = check_dependencies(args.skip_mash)

    # instantiation/checking fasta and gene_predictor
    if args.genbank is True:
        logger.info("You have specified --genbank.")
        logger.info(
            f"Therefore, {args.infile} is a genbank file instead of a FASTA file."
        )
        logger.info("Your custom CDS calls in this genbank file will be preserved.")
        validate_and_extract_genbank(args.infile, out_dir)
        gene_predictor = "genbank"
        input_fasta = f"{out_dir}/genbank.fasta"
    else:
        validate_fasta(args.infile)
        input_fasta = args.infile

    # other validations
    validate_gene_predictor(gene_predictor, args.genbank)
    validate_threads(args.threads)

    ###################
    # define input
    ###################

    # reorient with dnaapler if chosen
    if args.dnaapler is True:
        logger.info(
            "You have chosen to reorient your contigs by specifying --dnaapler. Checking the input."
        )

        # in case both --dnaapler and --terminase is selected
        if args.terminase == True:
            logger.info("Ignoring --terminase. Dnaapler will be run instead.")
            args.terminase = False
            args.terminase_strand = "nothing"
            args.terminase_start = "nothing"

        # count contigs
        contig_count = count_contigs(input_fasta)

        # runs dnaapler
        dnaapler_success = run_dnaapler(
            input_fasta, contig_count, out_dir, args.threads, logdir
        )

        if dnaapler_success == True:
            input_fasta = os.path.join(out_dir, "dnaapler/dnaapler_reoriented.fasta")
            destination_file = os.path.join(
                out_dir, f"{prefix}_dnaapler_reoriented.fasta"
            )
            # copy FASTA to output
            shutil.copy(input_fasta, destination_file)

    # terminase reorienting if chosen and dnaapler is False
    if args.dnaapler is False:
        if args.terminase is True:
            logger.info(
                "You have chosen to reorient your genome by specifying --terminase. Checking the input."
            )
            validate_terminase(input_fasta, args.terminase_strand, args.terminase_start)
            reorient_terminase(
                input_fasta,
                out_dir,
                prefix,
                args.terminase_strand,
                args.terminase_start,
            )
            # overwrite input_fasta if terminase reorienter is true
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

    ########
    # mmseqs2 and hmm decisions
    ########

    # can't have fast and mmseqs2 only
    if args.fast is True and args.mmseqs2_only is True:
        logger.error(
            "You have specified --fast or --hmm_only and --mmseqs2_only. This is impossible. Please choose one or the other."
        )

    # can't have fast and meta_hmm
    if args.fast is True and args.meta_hmm is True:
        logger.error(
            "You have specified --fast or --hmm_only and --meta_hmm. This is impossible. Please choose one or the other."
        )

    # can't have mmseqs2 only and meta_hmm
    if args.mmseqs2_only is True and args.meta_hmm is True:
        logger.error(
            "You have specified --mmseqs2_only and --meta_hmm. This is impossible. Please choose one or the other."
        )

    # by default run both not in meta mode
    mmseqs_flag = True
    hmm_flag = True

    if args.meta is True:  # meta mode default only mmseqs
        if args.meta_hmm is True:  # with --meta_hmm
            logger.info(
                "You have specified --meta_hmm and -m/--meta to run PyHMMER as well as MMseqs2 in meta mode. This may be take a while, please be patient."
            )
            mmseqs_flag = True
            hmm_flag = True
        else:  # just mmseqs2 by default
            mmseqs_flag = True
            hmm_flag = False
    else:  # not in meta mode
        if args.meta_hmm is True:
            logger.warning(
                "You have specified --meta_hmm to run PyHMMER in meta mode, but you have not specified -m to activate meta mode."
            )
            logger.warning("Ignoring --meta_hmm.")

    # overrides if fast/hmm_only is chosen
    if args.fast == True:
        logger.info("You have specified --fast or --hmm_only. MMseqs2 will not be run.")
        if args.meta is True:
            logger.warning(
                "You have specified --fast with -m/--meta. This may be take a while, please be patient."
            )
        mmseqs_flag = False
        hmm_flag = True

    if args.mmseqs2_only == True:
        logger.info("You have specified --mmseqs2_only. PyHMMER will not be run.")
        mmseqs_flag = True
        hmm_flag = False

    # validates meta mode
    validate_meta(input_fasta, args.meta, args.split, args.genbank)

    # meta mode split input for trnascan and maybe phanotate
    if args.meta == True:
        num_fastas = split_input_fasta(input_fasta, out_dir)
        # will generate split output gffs if split flag is true
        instantiate_split_output(out_dir, args.split)

    # CDS predicton
    # phanotate
    if gene_predictor == "phanotate":
        logger.info("Running Phanotate.")
        if args.meta == True:
            logger.info("Applying meta mode.")
            run_phanotate_fasta_meta(input_fasta, out_dir, args.threads, num_fastas)
            run_phanotate_txt_meta(input_fasta, out_dir, args.threads, num_fastas)
            concat_phanotate_meta(out_dir, num_fastas)
        else:
            run_phanotate(input_fasta, out_dir, logdir)
    elif gene_predictor == "prodigal":
        logger.info("Implementing Prodigal using Pyrodigal.")
        run_pyrodigal(
            input_fasta, out_dir, args.meta, args.coding_table, int(args.threads)
        )
    elif gene_predictor == "genbank":
        logger.info("Extracting CDS information from your genbank file.")
    elif gene_predictor == "prodigal-gv":
        logger.info("Implementing Prodigal-gv using Pyrodigal-gv.")
        run_pyrodigal_gv(input_fasta, out_dir, int(args.threads))

    # translate fastas (parse genbank)
    translate_fastas(out_dir, gene_predictor, args.coding_table, args.infile)

    # run trna-scan meta mode if required
    if args.skip_extra_annotations is False:
        if args.meta == True:
            logger.info("Starting tRNA-scanSE. Applying meta mode.")
            run_trnascan_meta(input_fasta, out_dir, args.threads, num_fastas)
            concat_trnascan_meta(out_dir, num_fastas)
        else:
            logger.info("Starting tRNA-scanSE.")
            run_trna_scan(input_fasta, args.threads, out_dir, logdir)
        # run minced and aragorn
        run_minced(input_fasta, out_dir, prefix, args.minced_args, logdir)
        run_aragorn(input_fasta, out_dir, prefix, logdir)

    # running mmseqs2 on the 3 databases
    if mmseqs_flag is True:
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

    if hmm_flag is True:
        # runs pyhmmer on PHROGs
        logger.info("Running PyHMMER on PHROGs.")
        best_results_pyhmmer = run_pyhmmer(
            db_dir, out_dir, args.threads, gene_predictor, args.evalue
        )

    if custom_hmm_flag is True:
        # runs pyhmmer on customer
        logger.info(f"Running PyHMMER on custom HMM database {args.custom_hmm}.")
        best_results_custom_pyhmmer = run_custom_pyhmmer(
            args.custom_hmm, out_dir, args.threads, gene_predictor, args.evalue
        )

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
    pharok.coding_table = args.coding_table
    pharok.mmseqs_flag = mmseqs_flag
    pharok.hmm_flag = hmm_flag
    pharok.custom_hmm_flag = custom_hmm_flag
    pharok.phanotate_version = phanotate_version
    pharok.pyrodigal_version = pyrodigal_version
    pharok.pyrodigal_gv_version = pyrodigal_gv_version
    pharok.trna_version = trna_version
    pharok.aragorn_version = aragorn_version
    pharok.minced_version = minced_version
    pharok.skip_extra_annotations = args.skip_extra_annotations

    if pharok.hmm_flag is True:
        pharok.pyhmmer_results_dict = best_results_pyhmmer
    if pharok.custom_hmm_flag is True:
        pharok.custom_pyhmmer_results_dict = best_results_custom_pyhmmer

    #####################################
    # post processing
    #####################################

    # gets df of length and gc for each contig
    pharok.get_contig_name_lengths()

    # post process results
    # includes vfdb and card
    # adds the merged df, vfdb and card top hits dfs to the class objec
    # no need to specify params as they are in the class :)
    pharok.process_results()

    # parse the aragorn output
    # only if not skipping annots
    if args.skip_extra_annotations is False:
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

    # create and write vfdb and card tophits
    # needs to be before .create_txt or else won't count properly
    # will only write to file if mmseqs2 has been run (in function)

    pharok.write_tophits_vfdb_card()

    # write the summary tsv outputs and functions
    pharok.create_txt()

    # convert to genbank
    logger.info("Converting gff to genbank.")
    # not part of the class so from processes.py
    convert_gff_to_gbk(input_fasta, out_dir, out_dir, prefix, pharok.prot_seq_df)

    # update fasta headers and final output tsv
    pharok.update_fasta_headers()
    pharok.update_final_output()

    # output single gffs in meta mode
    if args.split == True and args.meta == True:
        # splits the faa into single .faa
        pharok.split_faas_singles()

    # extract terL
    pharok.extract_terl()

    # run mash
    if args.skip_mash is False:  # skips mash
        logger.info("Finding the closest match for each contig in INPHARED using mash.")
        # in process.py
        run_mash_sketch(input_fasta, out_dir, logdir)
        run_mash_dist(out_dir, db_dir, args.mash_distance, logdir)
        # part of the class
        pharok.inphared_top_hits()
    else:
        logger.info("You have chosen --skip_mash.")
        logger.info(
            "Skipping finding the closest match for each contig in INPHARED using mash."
        )

    # delete tmp files
    remove_post_processing_files(out_dir, gene_predictor, args.meta)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("Pharokka has finished.")
    logger.info(f"Elapsed time: {elapsed_time} seconds")

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
