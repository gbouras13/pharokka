#!/usr/bin/env python3
import sys
from modules import input_commands
from modules import processes
from modules import post_processing
import os
import subprocess as sp
import logging

if __name__ == "__main__":


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

    LOG_FILE = os.path.join(args.outdir, prefix + ".log")
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO,filename=LOG_FILE,format='%(asctime)s - %(levelname)s - %(message)s')
    
    input_commands.validate_fasta(args.infile)
    processes.run_phanotate(args.infile, out_dir, logger)
    processes.translate_fastas(out_dir)
    processes.run_trna_scan(args.infile, out_dir, logger)

    # set the db dir
    if args.database == "Default":
        DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    else:
        DBDIR = args.database

    processes.remove_delim_fastas(out_dir)
    processes.run_mmseqs(DBDIR, out_dir, args.threads, logger)
    processes.run_hmmsuite(DBDIR, out_dir, args.threads, logger)
    phan_mmseq_merge_df = post_processing.process_results(DBDIR, out_dir, prefix)
    length_df = post_processing.get_contig_name_lengths(args.infile, out_dir, prefix)
    post_processing.create_gff(phan_mmseq_merge_df, length_df, args.infile, out_dir, prefix, locustag)
    post_processing.create_tbl(phan_mmseq_merge_df, length_df, out_dir, prefix)
    post_processing.create_txt(phan_mmseq_merge_df, length_df,out_dir, prefix)
    post_processing.convert_gff_to_gbk(args.infile, out_dir, prefix)
    
    # delete tmp files
    sp.run(["rm", "-rf", os.path.join(out_dir, "target_dir") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "tmp_dir/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "cleaned_phanotate.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "input_fasta_delim.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs_results.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_hhsuite.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_mmseqs.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "hhsuite_target_dir") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "phanotate_out.txt") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "trnascan_out.gff") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "phanotate_aas_tmp.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "phanotate_out_tmp.fasta") ])
    print("pharokka has finished")

    




