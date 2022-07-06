import os
from re import T
import subprocess as sp
from Bio import SeqIO
import random
import string
from Bio.SeqUtils import GC
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None

def process_results(db_dir,out_dir, prefix):

    ##mmseqs

    mmseqs_file =  os.path.join(out_dir, "mmseqs_results.tsv")
    print("Processing mmseqs output")
    col_list = ["phrog", "gene", "alnScore", "seqIdentity", "eVal", "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen"] 
    mmseqs_df = pd.read_csv(mmseqs_file, delimiter= '\t', index_col=False , names=col_list) 
    genes = mmseqs_df.gene.unique()

    tophits = []

    for gene in genes:
        tmp_df = mmseqs_df.loc[mmseqs_df['gene'] == gene].sort_values('eVal').reset_index(drop=True).loc[0]
        tophits.append([tmp_df.phrog, tmp_df.gene, tmp_df.alnScore, tmp_df.seqIdentity, tmp_df.eVal])

    tophits_df = pd.DataFrame(tophits, columns=['phrog', 'gene', 'alnScore', 'seqIdentity', 'eVal'])
    tophits_df.to_csv(os.path.join(out_dir, "top_hits_mmseqs.tsv"), sep="\t", index=False)
    # left join mmseqs top hits to phanotate
    phan_file = os.path.join(out_dir, "cleaned_phanotate.tsv") 
    # automatically picks up the names
    phan_df = pd.read_csv(phan_file, sep="\t", index_col=False )
    phan_df['gene']=phan_df['gene'].astype(str)
    tophits_df['gene']=tophits_df['gene'].astype(str)
    # merge top hit
    merged_df = phan_df.merge(tophits_df, on='gene', how='left')
    merged_df[['phrog','top_hit']] = merged_df['phrog'].str.split(' ## ',expand=True)
    merged_df["phrog"] = merged_df["phrog"].str.replace("phrog_", "")
    
    # get phrog annotaion file
    phrog_annot_df = pd.read_csv( os.path.join(db_dir, "phrog_annot_v3.tsv"), sep="\t", index_col=False )
    # merge phrog
    phrog_annot_df['phrog']=phrog_annot_df['phrog'].astype(str)
    merged_df = merged_df.merge(phrog_annot_df, on='phrog', how='left')
    merged_df = merged_df.replace(np.nan, 'No_PHROG', regex=True)
    merged_df['annot'] = merged_df["annot"].str.replace("No_PHROG", "hypothetical protein")
    merged_df['category'] = merged_df["category"].str.replace("No_PHROG", "unknown function")

    # add columns
    merged_df['Method'] = "PHANOTATE"
    merged_df['Region'] = "CDS"

    ############################
    ########## hhsuite

    hhs_dir = out_dir + "/hhsuite_target_dir/"

    hhsuite_file =  os.path.join(hhs_dir, "results_tsv_file.ffdata")
    print("Processing hhsuite output")
    col_list = ["gene_hmm", "phrog_hmm", "seqIdentity_hmm", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "eVal_hmm", "alnScore_hmm"] 
    hhsuite_df = pd.read_csv(hhsuite_file, delimiter= '\t', index_col=False , names=col_list) 
    genes = hhsuite_df.gene_hmm.unique()
    # remove nan
    genes = [x for x in genes if str(x) != 'nan']
    tophits = []

    for gene in genes:
        tmp_df = hhsuite_df.loc[hhsuite_df['gene_hmm'] == gene].sort_values('eVal_hmm').reset_index(drop=True).iloc[0]
        tophits.append([tmp_df.phrog_hmm, tmp_df.gene_hmm, tmp_df.alnScore_hmm, tmp_df.seqIdentity_hmm, tmp_df.eVal_hmm])

    tophits_hmm__df = pd.DataFrame(tophits, columns=['phrog_hmm', 'gene_hmm', 'alnScore_hmm', 'seqIdentity_hmm', 'eVal_hmm'])

    # filter from 0 to end for savings
    tophits_hmm__df[['spl','ind']] = tophits_hmm__df['gene_hmm'].str.split('delim',expand=True)
    tophits_hmm__df[['ind']] = tophits_hmm__df[['ind']].astype(int)
    tophits_hmm__df = tophits_hmm__df.sort_values(by=['ind']).drop(columns = ['spl', 'ind'])
    tophits_hmm__df.to_csv(os.path.join(out_dir, "top_hits_hhsuite.tsv"), sep="\t", index=False)
    
    
    ################
    ### merge in hmm

    # add match type
    merged_df['match_type'] = np.where(merged_df['phrog'] == "No_PHROG", 'hmm', 'mmseqs')

    merged_df[['gene_hmm','loca']] = merged_df['gene'].str.split(' ',expand=True)
    merged_df = merged_df.merge(tophits_hmm__df, on='gene_hmm', how='left')

    # replace with hmm if nothing found for mmseqs
    merged_df.loc[merged_df['phrog'] == 'No_PHROG', 'phrog'] = merged_df['phrog_hmm']
    merged_df.loc[merged_df['alnScore'] == 'No_PHROG', 'alnScore'] = merged_df['alnScore_hmm']
    merged_df.loc[merged_df['seqIdentity'] == 'No_PHROG', 'seqIdentity'] = merged_df['seqIdentity_hmm']
    merged_df.loc[merged_df['eVal'] == 'No_PHROG', 'eVal'] = merged_df['eVal_hmm']
    merged_df.loc[merged_df['top_hit'] == 'No_PHROG', 'top_hit'] = 'NA'
    merged_df.loc[merged_df['color'] == 'No_PHROG', 'color'] = 'NA'
    
    # get phrog
    merged_df["phrog"] = merged_df["phrog"].str.replace("phrog_", "")
    merged_df['phrog']=merged_df['phrog'].astype(str)
    # drop existing color annot category cols
    merged_df = merged_df.drop(columns = ['color', 'annot', 'category'])
    merged_df = merged_df.merge(phrog_annot_df, on='phrog', how='left')
    merged_df["annot"] = merged_df["annot"].replace(np.nan, 'hypothetical protein', regex=True)

    # get rid of "delimiter"
    merged_df["contig"] = merged_df["contig"].str.replace("delim", "")
    merged_df["gene"] = merged_df["gene"].str.replace("delim", "_")
    merged_df["gene_hmm"] = merged_df["gene_hmm"].str.replace("delim", "_")

    merged_df.to_csv( os.path.join(out_dir, prefix + "_final_merged_output.tsv"), sep="\t", index=False)
    
    return merged_df

def get_contig_name_lengths(fasta_input, out_dir, prefix):
    fasta_sequences = SeqIO.parse(open(fasta_input),'fasta')
    contig_names = []
    lengths = []
    gc = []
    for fasta in fasta_sequences:
        contig_names.append(fasta.id)
        lengths.append(len(fasta.seq))
        gc.append(GC(fasta.seq))
    length_df = pd.DataFrame(
    {'contig': contig_names,
     'length': lengths,
     'gc_perc': gc,
    })
    length_df.to_csv(os.path.join(out_dir, prefix + "_length_gc.tsv"), sep="\t", index=False)
    return(length_df)

def create_txt(phanotate_mmseqs_df, length_df, out_dir, prefix):
    contig_count = len(length_df)
    # with open( os.path.join(out_dir, "pharokka_summary.txt"), 'w') as f:
    #         f.write('Total Number of Contigs: ' + str(contig_count) + '\n')
    #         f.write('------------------------------------\n\n')
    #         f.close()   

    contigs = length_df["contig"]
    description_list = []

    for contig in contigs:
        phanotate_mmseqs_df_cont = phanotate_mmseqs_df[phanotate_mmseqs_df['contig'] == contig]
        cds_count = len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['Region'] == 'CDS'])
        trna_count = len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['Region'] == 'tRNA'])
        # get function
        phanotate_mmseqs_df_cont[['attributes2']] = phanotate_mmseqs_df_cont[['attributes']]
        phanotate_mmseqs_df_cont[['attributes2','function']] = phanotate_mmseqs_df_cont['attributes2'].str.split(';function=',expand=True)
        phanotate_mmseqs_df_cont = phanotate_mmseqs_df_cont.drop(columns=['attributes2'])
        phanotate_mmseqs_df_cont[['function','product']] = phanotate_mmseqs_df_cont['function'].str.split(';product=',expand=True)
        phanotate_mmseqs_df_cont = phanotate_mmseqs_df_cont.drop(columns=['product'])
        # get counts of functions and cds 
        regions = phanotate_mmseqs_df_cont['Region'].value_counts()
        functions = phanotate_mmseqs_df_cont['function'].value_counts()
        # reset index gets the names, then drop drops the first row (a blank index)
        description_df = pd.concat([regions, functions]).to_frame(name = "Test").reset_index()
        description_df.columns = ['Description', 'Count']
        description_df['contig'] = contig
        description_list.append(description_df)
        # remove pharokka summary for now until debug
        # with open( os.path.join(out_dir, "pharokka_summary.txt"), 'a') as f:
        #     f.write('Contig: ' + str(contig) + '\n')
        #     f.write('------------------------------------\n')
        #     f.write('CDS: ' + str(cds_count) + '\n')
        #     f.write('tRNA: ' + str(trna_count) + '\n\n')
        #     f.write('CDS Function Summary\n')
        #     f.write('------------------------------------\n')
        #     f.write('head and packaging: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("head and packaging")])) + '\n')
        #     f.write('connector: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("connector")])) + '\n')
        #     f.write('tail: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("tail")])) + '\n')
        #     f.write('DNA, RNA and nucleotide metabolism: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("DNA, RNA and nucleotide metabolism")])) + '\n')
        #     f.write('integration and excision: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("integration and excision")])) + '\n')
        #     f.write('lysis: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("lysis")])) + '\n')
        #     f.write('transcription regulation: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("transcription regulation")])) + '\n')
        #     f.write('moron, auxiliary metabolic gene and host takeover: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("moron, auxiliary metabolic gene and host takeover")])) + '\n')
        #     f.write('other: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("other")])) + '\n')
        #     f.write('unknown function: ' + str(len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['attributes'].str.contains("unknown function")])) + '\n')
        #     f.write('=====================================\n\n')
        #     f.close()

    #description_df = description_df.drop(index=description_df.index[0], axis=0, inplace=True)
    description_total_df = pd.concat(description_list)
    # add tRNA count
    trna_row = pd.DataFrame({ 'Description':['tRNAs'], 'Count':[trna_count], 'contig':[contig] })
    description_total_df = pd.concat([description_total_df, trna_row])
    #description_total_df = description_total_df.append({'Description':'tRNAs', 'Count':trna_count, 'contig':contig}, ignore_index=True)
    description_total_df.to_csv(os.path.join(out_dir, prefix + "_cds_functions.tsv"), sep="\t", index=False)





    # save as tsv

  
def create_gff(phanotate_mmseqs_df, length_df, fasta_input, out_dir, prefix, locustag):
    # write the headers of the gff file
    with open(os.path.join(out_dir, prefix + ".gff"), 'w') as f:
        f.write('##gff-version 3\n')
        for index, row in length_df.iterrows():
            f.write('##sequence-region ' + row['contig'] + ' 1 ' + str(row['length']) +'\n')
  
    # rearrange start and stop so that start is always less than stop for gff
    cols = ["start","stop"]
    #indices where start is greater than stop
    ixs = phanotate_mmseqs_df['frame'] == '-'
    # Where ixs is True, values are swapped
    phanotate_mmseqs_df.loc[ixs,cols] = phanotate_mmseqs_df.loc[ixs, cols].reindex(columns=cols[::-1]).values

    if locustag == "Random":
        # locus tag header 8 random letters
        locustag = ''.join(random.choice(string.ascii_uppercase) for _ in range(8))
        
    
    phanotate_mmseqs_df['phase'] = 0
    phanotate_mmseqs_df['attributes'] = "ID=" + locustag + "_" + phanotate_mmseqs_df.index.astype(str)  + ";" + "phrog=" + phanotate_mmseqs_df["phrog"] + ";" + "top_hit=" + phanotate_mmseqs_df["top_hit"] + ";" + "locus_tag=" + locustag + "_" + phanotate_mmseqs_df.index.astype(str) + ";" + "function=" + phanotate_mmseqs_df["category"] + ";"  + "product=" + phanotate_mmseqs_df["annot"]

    # get gff dataframe in correct order 
    gff_df = phanotate_mmseqs_df[["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]]

    # change start and stop to int 
    gff_df["start"] = gff_df["start"].astype('int')
    gff_df["stop"] = gff_df["stop"].astype('int')

    with open(os.path.join(out_dir, prefix + ".gff"), 'a') as f:
        gff_df.to_csv(f, sep="\t", index=False, header=False)
      
    ### trnas

    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
    trna_df = pd.read_csv(os.path.join(out_dir,"trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
    # keep only trnas
    trna_df = trna_df[(trna_df['Region'] == 'tRNA') | (trna_df['Region'] == 'pseudogene')]
    trna_df.start = trna_df.start.astype(int)
    trna_df.stop = trna_df.stop.astype(int)
    with open(os.path.join(out_dir, prefix + ".gff"), 'a') as f:
        trna_df.to_csv(f, sep="\t", index=False, header=False)


    # write fasta on the end 

    ##FASTA
    with open(os.path.join(out_dir, prefix + ".gff"), 'a') as f:
        f.write('##FASTA\n')
        fasta_sequences = SeqIO.parse(open(fasta_input),'fasta')
        SeqIO.write(fasta_sequences, f, "fasta")

def create_tbl(phanotate_mmseqs_df, length_df, out_dir, prefix):

    ### readtrnas

    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
     # check if no trnas
    empty = False
    if os.stat(os.path.join(out_dir, "trnascan_out.gff")).st_size == 0:
        empty = True
    if empty == False:    
        trna_df = pd.read_csv(os.path.join(out_dir, "trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
        # keep only trnas and pseudogenes 
        trna_df = trna_df[(trna_df['Region'] == 'tRNA') | (trna_df['Region'] == 'pseudogene')]
        trna_df.start = trna_df.start.astype(int)
        trna_df.stop = trna_df.stop.astype(int)

        trna_df[['attributes','isotypes']] = trna_df['attributes'].str.split(';isotype=',expand=True)
        trna_df[['isotypes','anticodon']] = trna_df['isotypes'].str.split(';anticodon=',expand=True)
        trna_df[['anticodon','rest']] = trna_df['anticodon'].str.split(';gene_biotype',expand=True)
        trna_df['trna_product']='tRNA-'+trna_df['isotypes']+"("+trna_df['anticodon']+")"

    with open( os.path.join(out_dir, prefix + ".tbl"), 'w') as f:
        for index, row in length_df.iterrows():
            contig = row['contig']
            f.write('>' + contig + '\n')
            subset_df = phanotate_mmseqs_df[phanotate_mmseqs_df['contig'] == contig]
            for index, row in subset_df.iterrows():
                f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + row['Region'] + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "PHANOTATE\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "phrog=" + str(row['phrog']) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['annot']) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if empty == False:
                subset_trna_df = trna_df[trna_df['contig'] == contig]
                for index, row in subset_trna_df.iterrows():
                    f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + row['Region'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "tRNAscan-SE")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['trna_product']) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")







