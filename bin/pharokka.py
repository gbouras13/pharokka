#!/usr/bin/env python3
import sys
from modules import input_commands
from modules import processes
from modules import post_processing
import os
import subprocess as sp
import logging
import time
import datetime

if __name__ == "__main__":



    print("Starting pharokka.")

    # get start time
    start_time = time.time()

    # getting time for log file 

    time_for_log = datetime.datetime.now().strftime("%m%d%Y_%H%M%S")

    args = input_commands.get_input()
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
    logger.info("Starting pharokka")

    # instantiation/checking fasta and gene_predictor
    input_commands.validate_fasta(args.infile)
    input_commands.validate_gene_predictor(gene_predictor)

    # gene predictor
    if args.gene_predictor == "phanotate":
        logger.info("Starting Phanotate")
        processes.run_phanotate(args.infile, out_dir, logger)
    if gene_predictor == "prodigal":
        logger.info("Starting Prodigal")
        processes.run_prodigal(args.infile, out_dir, logger)

    logger.info("Translating gene predicted fastas.")
    processes.translate_fastas(out_dir,gene_predictor)
    logger.info("Starting tRNA-scanSE")
    processes.run_trna_scan(args.infile, out_dir, logger)

    # set the db dir
    if args.database == "Default":
        DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    else:
        DBDIR = args.database

    processes.remove_delim_fastas(out_dir,gene_predictor)

    # runnin mmseqs2
    logger.info("Starting mmseqs2")
    processes.run_mmseqs(DBDIR, out_dir, args.threads, logger, gene_predictor)
    logger.info("Starting hhsuite")
    processes.run_hmmsuite(DBDIR, out_dir, args.threads, logger, args.gene_predictor)

    # post processing
    phan_mmseq_merge_df = post_processing.process_results(DBDIR, out_dir, prefix, gene_predictor)
    logger.info("Post Processing Data")
    length_df = post_processing.get_contig_name_lengths(args.infile, out_dir, prefix)
    post_processing.create_gff(phan_mmseq_merge_df, length_df, args.infile, out_dir, prefix, locustag)
    post_processing.create_tbl(phan_mmseq_merge_df, length_df, out_dir, prefix)
    post_processing.create_txt(phan_mmseq_merge_df, length_df,out_dir, prefix)
    logger.info("Converting gff to genbank using seqret")
    processes.convert_gff_to_gbk(args.infile, out_dir, prefix, logger)
    
    # delete tmp files
    sp.run(["rm", "-rf", os.path.join(out_dir, "target_dir") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "tmp_dir/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "cleaned_" + gene_predictor + ".tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "input_fasta_delim.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs_results.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_hhsuite.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_mmseqs.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "hhsuite_target_dir") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "phanotate_out.txt") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "trnascan_out.gff") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_aas_tmp.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_out_tmp.fasta") ])

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Show elapsed time for the process
    logger.info("pharokka has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("pharokka has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")

    




