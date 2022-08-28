#!/usr/bin/env python3
import sys
import input_commands
import processes
import post_processing
import os
import subprocess as sp
import logging
import time
import datetime
from version import __version__

if __name__ == "__main__":

    v = __version__

    # get start time
    start_time = time.time()

    # getting time for log file 

    time_for_log = datetime.datetime.now().strftime("%m%d%Y_%H%M%S")

    args = input_commands.get_input()
    # after input so that version is printed only
    print("Starting pharokka.")

        # set the prefix
    if args.prefix == "Default":
        prefix = "pharokka"
    else:
        prefix = args.prefix
    
    if args.locustag == "Default":
        locustag = "Random"
    else:
        locustag = args.locustag

    out_dir = input_commands.instantiate_dirs(args.outdir, args.force) # incase there is already an outdir
    gene_predictor = args.gene_predictor

    LOG_FILE = os.path.join(args.outdir, prefix + "_" + str(time_for_log) + ".log")
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO,filename=LOG_FILE,format='%(asctime)s - %(levelname)s - %(message)s')
    logging.captureWarnings(True)
    logger.info("Starting pharokka v " + v)

    # instantiation/checking fasta and gene_predictor
    input_commands.validate_fasta(args.infile)
    input_commands.validate_gene_predictor(gene_predictor)

    # gene predictor
    if gene_predictor == "phanotate":
        logger.info("Starting Phanotate")
        processes.run_phanotate(args.infile, out_dir, logger, args.TAG, args.TGA, args.TAA)
    if gene_predictor == "prodigal":
        logger.info("Starting Prodigal")
        processes.run_prodigal(args.infile, out_dir, logger, args.meta, args.coding_table)

    # translate fastas
    logger.info("Translating gene predicted fastas.")
    processes.translate_fastas(out_dir,gene_predictor)

    # trna and minced
    processes.run_trna_scan(args.infile, out_dir, logger)
    processes.run_minced(args.infile, out_dir, prefix, logger)
    processes.run_aragorn(args.infile, out_dir, prefix, logger)

    # set the db dir
    if args.database == "Default":
        DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    else:
        DBDIR = args.database


    # running mmseqs2
    logger.info("Starting mmseqs2")
    processes.run_mmseqs(DBDIR, out_dir, args.threads, logger, gene_predictor, args.evalue)

    # post processing
    phan_mmseq_merge_df = post_processing.process_results(DBDIR, out_dir, prefix, gene_predictor)
    logger.info("Post Processing Output.")
    print("Post Processing Output.")
    length_df = post_processing.get_contig_name_lengths(args.infile)
    tmrna_flag = post_processing.parse_aragorn(out_dir,length_df, prefix)
    # return locus_tag for table
    locustag = post_processing.create_gff(phan_mmseq_merge_df, length_df, args.infile, out_dir, prefix, locustag, tmrna_flag)
    post_processing.create_tbl(phan_mmseq_merge_df, length_df, out_dir, prefix, gene_predictor, tmrna_flag, locustag)
    post_processing.create_txt(phan_mmseq_merge_df, length_df,out_dir, prefix, tmrna_flag)
    
    # convert to genbank
    logger.info("Converting gff to genbank.")
    print("Converting gff to genbank.")
    processes.convert_gff_to_gbk(args.infile, out_dir, prefix, logger)
    
    # delete tmp files
    post_processing.remove_post_processing_files(out_dir, gene_predictor)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("pharokka has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("pharokka has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")

    




