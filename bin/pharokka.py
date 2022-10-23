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
    out_dir = input_commands.instantiate_dirs(args.outdir, args.meta, gene_predictor, args.force)

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

    logging.info("Input args: %r", args)

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

    # validates meta mode 
    input_commands.validata_meta(args.infile, args.meta)

    # meta mode split input for trnascan and maybe phanotate 
    if args.meta == True:
        num_fastas = processes.split_input_fasta(args.infile, out_dir)


    # CDS predicton
    # phanotate 
    if gene_predictor == "phanotate":
        print("Running Phanotate.")
        logger.info("Running Phanotate.")
        if args.meta == True:
            print("Applying meta mode.")
            logger.info("Applying meta mode.")
            processes.run_phanotate_fasta_meta(args.infile, out_dir, args.threads, num_fastas)
            processes.run_phanotate_txt_meta(args.infile, out_dir, args.threads, num_fastas)
            processes.concat_phanotate_meta(out_dir, num_fastas)
        else:
            processes.run_phanotate(args.infile, out_dir, logger)
    if gene_predictor == "prodigal":
        logger.info("Starting Prodigal")
        processes.run_prodigal(args.infile, out_dir, logger, args.meta, args.coding_table)

    # translate fastas
    logger.info("Translating gene predicted fastas.")
    processes.translate_fastas(out_dir,gene_predictor)


    # run trna-scan meta mode if required
    if args.meta == True:
        print("Running tRNAscan-SE. Applying meta mode.")
        logger.info("Starting tRNA-scanSE. Applying meta mode.")
        processes.run_trnascan_meta(args.infile, out_dir, args.threads, num_fastas)
        processes.concat_trnascan_meta(out_dir, num_fastas)
    else:
        print("Running tRNAscan-SE.")
        logger.info("Starting tRNA-scanSE")
        processes.run_trna_scan(args.infile,args.threads, out_dir, logger)
    # run minced and aragorn 
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
    # return the merged df, vfdb and card top hits
    (cds_mmseqs_merge_df,vfdb_results, card_results) = post_processing.process_results(db_dir, out_dir, prefix, gene_predictor)

    # gets df of length and gc for each contig
    length_df = post_processing.get_contig_name_lengths(args.infile)
    tmrna_flag = post_processing.parse_aragorn(out_dir,length_df, prefix)

    # create gff and return locustag for table
    (locustag, locus_df) = post_processing.create_gff(cds_mmseqs_merge_df, length_df, args.infile, out_dir, prefix, locustag, tmrna_flag)
    post_processing.create_tbl(cds_mmseqs_merge_df, length_df, out_dir, prefix, gene_predictor, tmrna_flag, locustag)
    # write the summary tsv outputs
    post_processing.create_txt(cds_mmseqs_merge_df, length_df,out_dir, prefix)
    
    # convert to genbank
    logger.info("Converting gff to genbank.")
    print("Converting gff to genbank.")
    processes.convert_gff_to_gbk(args.infile, out_dir, prefix)

    # update fasta headers and final output tsv
    post_processing.update_fasta_headers(locus_df, out_dir, gene_predictor )
    post_processing.update_final_output(cds_mmseqs_merge_df, vfdb_results, card_results, locus_df, prefix, out_dir )
    # extract terL
    post_processing.extract_terl(locus_df, out_dir, gene_predictor, logger )
    
    # delete tmp files
    post_processing.remove_post_processing_files(out_dir, gene_predictor, args.meta)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("pharokka has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("pharokka has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")

    




