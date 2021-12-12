import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

def process_mmseqs_results():
    mmseqs_file =  "output/mmseqs_results.tsv"
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
    tophits_df.to_csv("output/top_hits.tsv", sep="\t", index=False)
    # left join mmseqs top hits to phanotate
    phan_file = "output/cleaned_phanotate.tsv"
    # automatically picks up the names
    #col_list = ["start", "stop", "frame", "contig", "score", "gene"] 
    phan_df = pd.read_csv(phan_file, sep="\t", index_col=False )
    phan_df['gene']=phan_df['gene'].astype(str)
    tophits_df['gene']=tophits_df['gene'].astype(str)
    # merge top hit
    merged_df = phan_df.merge(tophits_df, on='gene', how='left')
    merged_df[['phrog','top_hit']] = merged_df['phrog'].str.split(' ## ',expand=True)
    merged_df["phrog"] = merged_df["phrog"].str.replace("phrog_", "")
    # get phrog annotaion file
    phrog_annot_df = pd.read_csv('databases/phrog_annot_v3.tsv', sep="\t", index_col=False )
    # merge phrog
    phrog_annot_df['phrog']=phrog_annot_df['phrog'].astype(str)
    merged_df = merged_df.merge(phrog_annot_df, on='phrog', how='left')

    merged_df.to_csv("output/final_merged_output.tsv", sep="\t", index=False)
    print(merged_df)






