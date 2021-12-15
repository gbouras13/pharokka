#!/usr/bin/env python3
import sys
from modules import input_commands
from modules import processes
from modules import post_processing
import os

DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  

if __name__ == "__main__":
    args = input_commands.get_input()
    input_commands.validate_fasta(args.infile)
    input_commands.instantiate_dirs(args.outdir)
    processes.run_phanotate(args.infile, args.outdir)
    processes.translate_fastas(args.outdir)
    processes.run_trna_scan(args.infile, args.outdir)
    processes.run_mmseqs(DBDIR, args.outdir)
    phan_mmseq_merge_df = post_processing.process_mmseqs_results(DBDIR, args.outdir)
    length_df = post_processing.get_contig_name_lengths(args.infile)
    post_processing.create_gff(phan_mmseq_merge_df, length_df, args.infile, args.outdir)
    post_processing.create_tbl(phan_mmseq_merge_df, length_df, args.outdir)
    sys.exit("phrokka has finished")  



