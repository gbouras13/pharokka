#!/usr/bin/env python3
import sys
from modules import input_commands
from modules import processes
from modules import post_processing
import os
import subprocess as sp

DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  

if __name__ == "__main__":
    args = input_commands.get_input()
    input_commands.validate_fasta(args.infile)
    out_dir = input_commands.instantiate_dirs(args.outdir) # incase there is already an outdir
    #processes.run_phanotate(args.infile, out_dir)
    #processes.translate_fastas(out_dir)
    #processes.run_trna_scan(args.infile, out_dir)
    #processes.run_mmseqs(DBDIR, out_dir)
    #processes.run_hmmsuite(DBDIR, out_dir)
    phan_mmseq_merge_df = post_processing.process_results(DBDIR, out_dir)
    length_df = post_processing.get_contig_name_lengths(args.infile)
    post_processing.create_gff(phan_mmseq_merge_df, length_df, args.infile, out_dir)
    post_processing.create_tbl(phan_mmseq_merge_df, length_df, out_dir)
    post_processing.create_txt(phan_mmseq_merge_df, length_df, out_dir)

    # delete tmp
    sp.call(["rm", "-r", os.path.join(os.path.dirname(__file__), "../", "tmp/") ])
     
    sys.exit("phrokka has finished")  



