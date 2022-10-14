#!/usr/bin/env python3
import sys
import input_commands
import databases
import processes
import post_processing
import os
import logging
import time
import datetime
from version import __version__

if __name__ == "__main__":

    __author__ = 'George Bouras'

    v = __version__

    # get start time
    start_time = time.time()

    # getting time for log file 
    time_for_log = datetime.datetime.now().strftime("%m%d%Y_%H%M%S")

    # get the args
    args = input_commands.get_input()

    # after input so that version is printed only
    print("Starting pharokka v" + v)

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
    out_dir = input_commands.instantiate_dirs(args.outdir, args.force)

    # set the db dir
    if args.database == "Default":
        db_dir = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    else:
        db_dir = args.database

    # start logging
    LOG_FILE = os.path.join(args.outdir, prefix + "_" + str(time_for_log) + ".log")
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO,filename=LOG_FILE,format='%(asctime)s - %(levelname)s - %(message)s')
    logging.captureWarnings(True)
    logger.info("Starting pharokka v" + v)

    # check the database is installed
    print("Checking database installation")
    logger.info("Checking database installation")
    database_installed = databases.check_db_installation(db_dir)
    if database_installed == True:
        print("All databases have been successfully checked")
        logger.info("All databases have been successfully checked")
    else:
        sys.exit("\nThe database directory was unsuccessfully checked. Please run install_databases.py \n") 

    # instantiation/checking fasta and gene_predictor
    input_commands.validate_fasta(args.infile)
    input_commands.validate_gene_predictor(gene_predictor)

    # CDS predicton
    if gene_predictor == "phanotate":
        logger.info("Starting Phanotate")
        processes.run_phanotate(args.infile, out_dir, logger)
    if gene_predictor == "prodigal":
        logger.info("Starting Prodigal")
        processes.run_prodigal(args.infile, out_dir, logger, args.meta, args.coding_table)

    # translate fastas
    logger.info("Translating gene predicted fastas.")
    processes.translate_fastas(out_dir,gene_predictor)

    # run trna-scan and minced
    processes.run_trna_scan(args.infile, out_dir, logger)
    processes.run_minced(args.infile, out_dir, prefix, logger)
    processes.run_aragorn(args.infile, out_dir, prefix, logger)

    # running mmseqs2
    logger.info("Starting mmseqs2.")
    processes.run_mmseqs(db_dir, out_dir, args.threads, logger, gene_predictor, args.evalue)
    processes.run_mmseqs_card(db_dir, out_dir, args.threads, logger, gene_predictor)
    processes.run_mmseqs_vfdb(db_dir, out_dir, args.threads, logger, gene_predictor)

    # post processing
    logger.info("Post Processing Output.")
    print("Post Processing Output.")

    # post process results
    # includes vfdb and card
    cds_mmseqs_merge_df = post_processing.process_results(db_dir, out_dir, prefix, gene_predictor)

    # gets df of length and gc for each contig
    length_df = post_processing.get_contig_name_lengths(args.infile)
    tmrna_flag = post_processing.parse_aragorn(out_dir,length_df, prefix)

    # create gff and return locustag for table
    locustag = post_processing.create_gff(cds_mmseqs_merge_df, length_df, args.infile, out_dir, prefix, locustag, tmrna_flag)
    post_processing.create_tbl(cds_mmseqs_merge_df, length_df, out_dir, prefix, gene_predictor, tmrna_flag, locustag)
    # write the summary tsv outputs
    post_processing.create_txt(cds_mmseqs_merge_df, length_df,out_dir, prefix)
    
    # convert to genbank
    logger.info("Converting gff to genbank.")
    print("Converting gff to genbank.")
    processes.convert_gff_to_gbk(args.infile, out_dir, prefix)
    
    # delete tmp files
    #post_processing.remove_post_processing_files(out_dir, gene_predictor)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("pharokka has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("pharokka has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")

    




