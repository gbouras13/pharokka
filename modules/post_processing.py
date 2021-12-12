import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

def process_mmseqs_results():
    input_dir = "output/"
    mmseqs_file = input_dir + "mmseqs_results.tsv"
    print("Processing mmseqs output")
    col_list = ["phrog", "gene", "alnScore", "seqIdentity", "eVal", "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen"] 
    mmseqs_df = pd.read_csv(mmseqs_file, delimiter= '\t', index_col=False , names=col_list) 
    genes = mmseqs_df.gene.unique()
    print(genes)

    tophits = []

    for gene in genes:
        tmp_df = mmseqs_df.loc[mmseqs_df['gene'] == gene].sort_values('eVal').reset_index(drop=True).loc[0]
        tophits.append([tmp_df.phrog, tmp_df.gene, tmp_df.alnScore, tmp_df.seqIdentity, tmp_df.eVal])
        print(tophits)

    tophits_df = pd.DataFrame(tophits, columns=['phrog', 'gene', 'alnScore', 'seqIdentity', 'eVal'])
    tophits_df.to_csv(input_dir + "top_hits.tsv", sep="\t")
    # left join 
    input_dir = "output/"
    phan_file = input_dir + "cleaned_phanotate.tsv"
    col_list = ["ind" ,"start", "stop", "frame", "contig", "score", "gene"] 
    phan_df = pd.read_csv(phan_file, sep="\t", index_col=False, names=col_list)
    phan_df['gene']=phan_df['gene'].astype(str)
    tophits_df['gene']=tophits_df['gene'].astype(str)
    merged_df = phan_df.merge(tophits_df, on='gene', how='left')
    print(merged_df)




