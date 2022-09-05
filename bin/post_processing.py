from cmath import nan
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

def process_results(db_dir,out_dir, prefix, gene_predictor):

    ##mmseqs

    mmseqs_file =  os.path.join(out_dir, "mmseqs_results.tsv")
    print("Processing mmseqs2 output.")
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
    phan_file = os.path.join(out_dir, "cleaned_" +  gene_predictor +  ".tsv") 
    # automatically picks up the names
    phan_df = pd.read_csv(phan_file, sep="\t", index_col=False )
    phan_df['gene']=phan_df['gene'].astype(str)
    tophits_df['gene']=tophits_df['gene'].astype(str)
    phan_df = phan_df[phan_df['start'].notna()]
    phan_df = phan_df.dropna()
    # merge top hits into the phanotate df
    merged_df = phan_df.merge(tophits_df, on='gene', how='left')
    # add test if empty
    if len(tophits_df['phrog']) == 0:
        merged_df['top_hit'] = 'No_PHROG'
    else:
        merged_df[['phrog','top_hit']] = merged_df['phrog'].str.split(' ## ',expand=True)
    merged_df["phrog"] = merged_df["phrog"].str.replace("phrog_", "")
    
    # get phrog annotaion file
    phrog_annot_df = pd.read_csv( os.path.join(db_dir, "phrog_annot_v4.tsv"), sep="\t", index_col=False )
    # merge phrog
    phrog_annot_df['phrog']=phrog_annot_df['phrog'].astype(str)
    merged_df = merged_df.merge(phrog_annot_df, on='phrog', how='left')
    merged_df = merged_df.replace(np.nan, 'No_PHROG', regex=True)
    merged_df['annot'] = merged_df["annot"].str.replace("No_PHROG", "hypothetical protein")
    merged_df['category'] = merged_df["category"].str.replace("No_PHROG", "unknown function")

    # add columns
    if gene_predictor == "phanotate":
        merged_df['Method'] = "PHANOTATE"
    if gene_predictor == "prodigal":
        merged_df['Method'] = "PRODIGAL"
    merged_df['Region'] = "CDS"

    # # replace with NA if nothing found for mmseqs
    merged_df.loc[merged_df['phrog'] == 'No_PHROG', 'phrog'] = 'No_PHROG'
    merged_df.loc[merged_df['alnScore'] == 'No_PHROG', 'alnScore'] = 'No_PHROG'
    merged_df.loc[merged_df['seqIdentity'] == 'No_PHROG', 'seqIdentity'] = 'No_PHROG'
    merged_df.loc[merged_df['eVal'] == 'No_PHROG', 'eVal'] = 'No_PHROG'
    merged_df.loc[merged_df['top_hit'] == 'No_PHROG', 'top_hit'] = 'No_PHROG'
    merged_df.loc[merged_df['color'] == 'No_PHROG', 'color'] = 'No_PHROG'
    
    # get phrog
    merged_df["phrog"] = merged_df["phrog"].str.replace("phrog_", "")
    merged_df['phrog']=merged_df['phrog'].astype(str)
    # drop existing color annot category cols
    merged_df = merged_df.drop(columns = ['color', 'annot', 'category'])
    merged_df = merged_df.merge(phrog_annot_df, on='phrog', how='left')
    merged_df["annot"] = merged_df["annot"].replace(nan, 'hypothetical protein', regex=True)
    merged_df["category"] = merged_df["category"].replace(nan, 'unknown function', regex=True)
    merged_df["color"] = merged_df["color"].replace(nan, 'none', regex=True)

    merged_df.to_csv( os.path.join(out_dir, prefix + "_final_merged_output.tsv"), sep="\t", index=False)
    
    return merged_df

def get_contig_name_lengths(fasta_input):
    fasta_sequences = SeqIO.parse(open(fasta_input),'fasta')
    contig_names = []
    lengths = []
    gc = []
    for fasta in fasta_sequences:
        contig_names.append(fasta.id)
        lengths.append(len(fasta.seq))
        gc.append(round(GC(fasta.seq),2))
    length_df = pd.DataFrame(
    {'contig': contig_names,
     'length': lengths,
     'gc_perc': gc,
    })
    return(length_df)

def create_txt(phanotate_mmseqs_df, length_df, out_dir, prefix, tmrna_flag):

    contigs = length_df["contig"]
    # instantiate the length_df['cds_coding_density']
    length_df['cds_coding_density'] = 0.0
    description_list = []

    # read in trnascan

    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
    trna_df = pd.read_csv(os.path.join(out_dir,"trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
    # keep only trnas
    trna_df = trna_df[(trna_df['Region'] == 'tRNA') | (trna_df['Region'] == 'pseudogene')]
    # get crispr count
    crispr_df = pd.read_csv(os.path.join(out_dir, prefix + "_minced.gff"), delimiter= '\t', index_col=False, names=col_list, comment = '#' ) 
    tmrna_df = pd.read_csv(os.path.join(out_dir, prefix + "_aragorn.gff"), delimiter= '\t', index_col=False, names=col_list ) 
    # get vfdb 
    # gene has the contig 
    vfdb_df = pd.read_csv(os.path.join(out_dir, "top_hits_vfdb.tsv"), delimiter= '\t', index_col=False ) 
    vfdb_df[['gene','coordinate']] = vfdb_df['gene'].str.split(' ',expand=True)
    # strip off the index with the last period delimiter
    vfdb_df['contig'] = vfdb_df['gene'].str.rsplit("\\.",1)

    for contig in contigs:
        phanotate_mmseqs_df_cont = phanotate_mmseqs_df[phanotate_mmseqs_df['contig'] == contig]
        # counts of the cds trnas
        cds_count = len(phanotate_mmseqs_df_cont[phanotate_mmseqs_df_cont['Region'] == 'CDS'])
        trna_count = len(trna_df[trna_df['contig'] == contig])
        tmrna_count = len(tmrna_df[tmrna_df['contig'] == contig])
        crispr_count = len(crispr_df[crispr_df['contig'] == contig])
        vfdb_count = len(vfdb_df[vfdb_df['contig'] == contig])
        

        vfdb_count = vfdb_df['gene'].str.contains(contig).sum()
        # get the total length of the contig
        contig_length = length_df[length_df["contig"] == contig]['length']
        if cds_count > 0:
            # gets the total cds coding length
            cds_lengths = abs(phanotate_mmseqs_df_cont['start'] - phanotate_mmseqs_df_cont['stop']).sum()
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
        else:
            description_df = pd.DataFrame({'Description': ["CDS"], 'Count': [0], 'contig': contig})
            cds_lengths = 0
        # add trna count
        trna_row = pd.DataFrame({ 'Description':['tRNAs'], 'Count':[trna_count], 'contig':[contig] })
        crispr_row = pd.DataFrame({ 'Description':['CRISPRs'], 'Count':[crispr_count], 'contig':[contig] })
        tmrna_row = pd.DataFrame({ 'Description':['tmRNAs'], 'Count':[tmrna_count], 'contig':[contig] })
        vfdb_row = pd.DataFrame({ 'Description':['VFDB_Virulence_Factors'], 'Count':[vfdb_count], 'contig':[contig] })
        # calculate the cds coding density and add to length_df
        cds_coding_density = cds_lengths * 100 / contig_length
        cds_coding_density = round(cds_coding_density, 2)
        length_df.loc[length_df['contig'] == contig, 'cds_coding_density'] = cds_coding_density
        # append it all
        description_list.append(description_df)
        description_list.append(trna_row)
        description_list.append(crispr_row)
        description_list.append(tmrna_row)
        description_list.append(vfdb_row)

    # save the output
    description_total_df = pd.concat(description_list)
    #description_total_df = description_total_df.append({'Description':'tRNAs', 'Count':trna_count, 'contig':contig}, ignore_index=True)
    description_total_df.to_csv(os.path.join(out_dir, prefix + "_cds_functions.tsv"), sep="\t", index=False)
    # save the length_gc.tsv also
    length_df.to_csv(os.path.join(out_dir, prefix + "_length_gc_cds_density.tsv"), sep="\t", index=False)


  
def create_gff(phanotate_mmseqs_df, length_df, fasta_input, out_dir, prefix, locustag, tmrna_flag, vfdb_flag):
  
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
    phanotate_mmseqs_df['attributes'] = "ID=" + locustag + "_CDS_" + phanotate_mmseqs_df.index.astype(str)  + ";" + "phrog=" + phanotate_mmseqs_df["phrog"] + ";" + "top_hit=" + phanotate_mmseqs_df["top_hit"] + ";" + "locus_tag=" + locustag + "_" + phanotate_mmseqs_df.index.astype(str) + ";" + "function=" + phanotate_mmseqs_df["category"] + ";"  + "product=" + phanotate_mmseqs_df["annot"]

    # if vfdb is true
    if vfdb_flag == True:
        phanotate_mmseqs_df.loc[phanotate_mmseqs_df['vfdb_short_name'] != "None", 'attributes'] = phanotate_mmseqs_df['attributes']  + ";"  + "vfdb_short_name=" + phanotate_mmseqs_df['vfdb_short_name'] + ";"  + "vfdb_description="  +  phanotate_mmseqs_df['vfdb_description'] + ";" + "vfdb_species=" + phanotate_mmseqs_df['vfdb_species']

    # get gff dataframe in correct order 
    gff_df = phanotate_mmseqs_df[["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]]
    

    # change start and stop to int 
    gff_df["start"] = gff_df["start"].astype('int')
    gff_df["stop"] = gff_df["stop"].astype('int')

    # write to tmp gff
    with open(os.path.join(out_dir, "phrokka_tmp.gff"), 'w') as f:
        gff_df.to_csv(f, sep="\t", index=False, header=False)
      
    ### trnas
    # check if no trnas
    trna_empty = is_trna_empty(out_dir)
    if trna_empty == False:   
        col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
        trna_df = pd.read_csv(os.path.join(out_dir,"trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
        # keep only trnas
        trna_df = trna_df[(trna_df['Region'] == 'tRNA') | (trna_df['Region'] == 'pseudogene')]
        trna_df = trna_df.reset_index(drop=True)
        trna_df.start = trna_df.start.astype(int)
        trna_df.stop = trna_df.stop.astype(int)
        trna_df[['attributes','isotypes']] = trna_df['attributes'].str.split(';isotype=',expand=True)
        trna_df[['isotypes','anticodon']] = trna_df['isotypes'].str.split(';anticodon=',expand=True)
        trna_df[['anticodon','rest']] = trna_df['anticodon'].str.split(';gene_biotype',expand=True)
        trna_df['trna_product']='tRNA-'+trna_df['isotypes']+"("+trna_df['anticodon']+")"
        trna_df = trna_df.drop(columns=['attributes'])
        trna_df['attributes'] = "ID=" + locustag + "_tRNA_" + trna_df.index.astype(str)  + ";" + "trna=" + trna_df["trna_product"] + ";" + "isotype=" + trna_df["isotypes"] + ";" + "anticodon=" + trna_df["anticodon"] + ";" + "locus_tag=" + locustag + "_tRNA_" + trna_df.index.astype(str)
        trna_df = trna_df.drop(columns=['isotypes', 'anticodon', 'rest', 'trna_product'])
        with open(os.path.join(out_dir, "phrokka_tmp.gff"), 'a') as f:
            trna_df.to_csv(f, sep="\t", index=False, header=False)

    ### crisprs
    crispr_count = get_crispr_count(out_dir, prefix)
    # add to gff if yes
    if crispr_count > 0:
        col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
        minced_df = pd.read_csv(os.path.join(out_dir, prefix + "_minced.gff"), delimiter= '\t', index_col=False, names=col_list, comment='#' ) 
        minced_df.start = minced_df.start.astype(int)
        minced_df.stop = minced_df.stop.astype(int)
        minced_df[['attributes','rpt_unit_seq']] = minced_df['attributes'].str.split(';rpt_unit_seq=',expand=True)
        minced_df[['attributes','rpt_family']] = minced_df['attributes'].str.split(';rpt_family=',expand=True)
        minced_df[['attributes','rpt_type']] = minced_df['attributes'].str.split(';rpt_type=',expand=True)
        minced_df = minced_df.drop(columns=['attributes'])
        minced_df['attributes'] = "ID=" + locustag + "_CRISPR_" + minced_df.index.astype(str)  + ";" + "rpt_type=" + minced_df["rpt_type"] + ";" + "rpt_family=" + minced_df["rpt_family"] + ";" + "rpt_unit_seq=" + minced_df["rpt_unit_seq"] + ";" + "locus_tag=" + locustag + "_CRISPR_" + minced_df.index.astype(str)
        minced_df = minced_df.drop(columns=['rpt_unit_seq', 'rpt_family', 'rpt_type'])
        with open(os.path.join(out_dir, "phrokka_tmp.gff"), 'a') as f:
            minced_df.to_csv(f, sep="\t", index=False, header=False)

    ### tmrna
    if tmrna_flag == True:
        col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
        tmrna_df = pd.read_csv(os.path.join(out_dir, prefix + "_aragorn.gff"), delimiter= '\t', index_col=False, names=col_list ) 
        tmrna_df.start = tmrna_df.start.astype(int)
        tmrna_df.stop = tmrna_df.stop.astype(int)
        tmrna_df['attributes'] = "ID=" + locustag + "_tmRNA_" + tmrna_df.index.astype(str)  + ";" + tmrna_df['attributes'] + ';locus_tag=' + locustag + "_tmRNA_" + tmrna_df.index.astype(str)
        with open(os.path.join(out_dir, "phrokka_tmp.gff"), 'a') as f:
            tmrna_df.to_csv(f, sep="\t", index=False, header=False)

    # write header
    with open(os.path.join(out_dir, prefix + ".gff"), 'w') as f:
        f.write('##gff-version 3\n')
        for index, row in length_df.iterrows():
            f.write('##sequence-region ' + row['contig'] + ' 1 ' + str(row['length']) +'\n')
    
    # sort gff
    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
    total_gff = pd.read_csv(os.path.join(out_dir, "phrokka_tmp.gff"), delimiter= '\t', index_col=False, names=col_list ) 
    total_gff.start = total_gff.start.astype(int)
    total_gff.stop = total_gff.stop.astype(int)
    total_gff = total_gff.sort_values(['contig', 'start'])

    with open(os.path.join(out_dir, prefix + ".gff"), 'a') as f:
        total_gff.to_csv(f, sep="\t", index=False, header=False)

     # write fasta on the end 
     
    ##FASTA
    with open(os.path.join(out_dir, prefix + ".gff"), 'a') as f:
        f.write('##FASTA\n')
        fasta_sequences = SeqIO.parse(open(fasta_input),'fasta')
        SeqIO.write(fasta_sequences, f, "fasta")

    return locustag

def create_tbl(phanotate_mmseqs_df, length_df, out_dir, prefix, gene_predictor, tmrna_flag, locustag):

    ### readtrnas

    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]

     # check if no trnas
    trna_empty = is_trna_empty(out_dir)
    if trna_empty == False:    
        trna_df = pd.read_csv(os.path.join(out_dir, "trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
        # keep only trnas and pseudogenes 
        trna_df = trna_df[(trna_df['Region'] == 'tRNA') | (trna_df['Region'] == 'pseudogene')]
        trna_df.start = trna_df.start.astype(int)
        trna_df.stop = trna_df.stop.astype(int)

        trna_df[['attributes','isotypes']] = trna_df['attributes'].str.split(';isotype=',expand=True)
        trna_df[['isotypes','anticodon']] = trna_df['isotypes'].str.split(';anticodon=',expand=True)
        trna_df[['anticodon','rest']] = trna_df['anticodon'].str.split(';gene_biotype',expand=True)
        trna_df['trna_product']='tRNA-'+trna_df['isotypes']+"("+trna_df['anticodon']+")"

    # check if no crisprs 
    # check if the file has more than 1 line (not empty)
    crispr_count = get_crispr_count(out_dir, prefix)
    if crispr_count > 0:
        crispr_df = pd.read_csv(os.path.join(out_dir, prefix + "_minced.gff"), delimiter= '\t', index_col=False, names=col_list, comment = "#"  ) 
        crispr_df.start = crispr_df.start.astype(int)
        crispr_df.stop = crispr_df.stop.astype(int)
        crispr_df[['attributes','rpt_unit_seq']] = crispr_df['attributes'].str.split(';rpt_unit_seq=',expand=True)

    if tmrna_flag == True:
        col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
        tmrna_df = pd.read_csv(os.path.join(out_dir, prefix + "_aragorn.gff"), delimiter= '\t', index_col=False, names=col_list) 
        tmrna_df.start = tmrna_df.start.astype(int)
        tmrna_df.stop = tmrna_df.stop.astype(int)

    if gene_predictor == "phanotate":
        inf = "PHANOTATE"
    else:
        inf = "PRODIGAL"
    with open( os.path.join(out_dir, prefix + ".tbl"), 'w') as f:
        for index, row in length_df.iterrows():
            contig = row['contig']
            f.write('>' + contig + '\n')
            subset_df = phanotate_mmseqs_df[phanotate_mmseqs_df['contig'] == contig]
            for index, row in subset_df.iterrows():
                f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + row['Region'] + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ inf + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "phrog=" + str(row['phrog']) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['annot']) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ locustag + "_CDS_" + str(index) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if trna_empty == False:
                subset_trna_df = trna_df[trna_df['contig'] == contig]
                for index, row in subset_trna_df.iterrows():
                    f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + row['Region'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "tRNAscan-SE" + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['trna_product']) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ locustag + "_tRNA_" + str(index) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if crispr_count > 0:
                subset_crispr_df = crispr_df[crispr_df['contig'] == contig]
                for index, row in subset_crispr_df.iterrows():
                    f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + row['Region'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "MinCED" +"\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ locustag + "_CRISPR_" + str(index) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['rpt_unit_seq']) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if tmrna_flag == True:
                subset_tmrna_df = tmrna_df[tmrna_df['contig'] == contig]
                for index, row in subset_tmrna_df.iterrows():
                    f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + 'tmRNA' + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "Aragorn" + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ locustag + "_tmRNA_" + str(index) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ 'transfer-messenger RNA, SsrA' + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")


def remove_post_processing_files(out_dir, gene_predictor):
    sp.run(["rm", "-rf", os.path.join(out_dir, "target_dir") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "tmp_dir/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "vfdb_tmp_dir/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "vfdb/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "vfdb_results.tsv") ])  
    sp.run(["rm", "-rf", os.path.join(out_dir, "cleaned_" + gene_predictor + ".tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "input_fasta_delim.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs_results.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_mmseqs.tsv") ])
    # leave in tophits
    sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_aas_tmp.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_out_tmp.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "phrokka_tmp.gff") ])  
    if gene_predictor == "phanotate":
        sp.run(["rm", "-rf", os.path.join(out_dir, "phanotate_out.txt") ])
    if gene_predictor == "prodigal":
        sp.run(["rm", "-rf", os.path.join(out_dir, "prodigal_out.gff") ])

# check if the crispr file has more than 1 line (not empty)
def get_crispr_count(out_dir, prefix):
    crispr_file = os.path.join(out_dir, prefix + "_minced.gff")
    with open(crispr_file) as file:
        lines = file.readlines()
    crispr_count = 0
    for line in lines:
        if line[0] != "#":
            crispr_count +=1
    return crispr_count

# check if the trna file has more than 1 line (not empty)
def is_trna_empty(out_dir):
    trna_empty = False
    if os.stat(os.path.join(out_dir, "trnascan_out.gff")).st_size == 0:
        trna_empty = True
    return trna_empty


# # check if the trna file has more than 1 line (not empty)
# def is_tmrna_empty(out_dir):
#     tmrna_empty = False
#     if os.stat(os.path.join(out_dir, prefix + "_aragorn.gff")).st_size == 0:
#         tmrna_empty = True
#     return tmrna_empty

def parse_aragorn(out_dir,length_df, prefix):
    aragorn_file = os.path.join(out_dir, prefix + "_aragorn.txt")
    f=open(aragorn_file)
    lines=f.readlines()
    contig_count = len(length_df["contig"])
    tmrna_flag = False
    contig_names = []
    methods = []
    regions = []
    starts = []
    stops = []
    scores = []
    frames = []
    phases = []
    attributes = []
    # if there is only one contig
    if contig_count == 1:
        # if no trnas
        if int(lines[1][0]) == 0:
            tmrna_df = pd.DataFrame(
            {'contig': '',
            'Method': '',
            'Region': '',
            'start': '',
            'stop': '',
            'score': '',
            'frame': '',
            'phase': '',
            'attributes': '',
            }, index=[0])
        else:
            tmrna_flag = True
            # get all lines with tmrnas
            tmrna_lines = lines[2:]
            for line in tmrna_lines:
                split = line.split()
                start_stops = split[2].replace("[", "").replace("]", "").split(',')
                contig = length_df["contig"][0]
                method = "Aragorn"
                region = "tmRNA"
                start = start_stops[0].replace("c", "")  # tmrna output is [start,stop] or c[start, stop] so need to remove c also for some phages
                stop = start_stops[1]
                score = "."
                frame = "."
                phase = "."
                tag_peptide = split[3]
                tag_peptide_seq = split[4]
                attribute = "product=transfer-messenger RNA SsrA;tag_peptide=" + tag_peptide + ";tag_peptide_sequence=" + tag_peptide_seq
                contig_names.append(contig)
                methods.append(method)
                regions.append(region)
                starts.append(start)
                stops.append(stop)
                scores.append(score)
                frames.append(frame)
                phases.append(phase)
                attributes.append(attribute)
            tmrna_df = pd.DataFrame(
            {'contig': contig_names,
            'Method': methods,
            'Region': regions,
            'start': starts,
            'stop': stops,
            'score': scores,
            'frame': frames,
            'phase': phases,
            'attributes': attributes,
            })
    # two or more contigs
    else:
        i = 0 # line counter
        j = 0 # contig counter
        for line in lines:
            # contig meets these reqs
            if line[0] == ">" and line[1:4] != "end":
                if lines[i+1][0] != 0:
                    tmrna_flag = True
                    # number of trnas for this contig
                    tmrna_count = int(lines[i+1][0])
                    # iterate over them
                    for k in range(tmrna_count):
                        tmrna_line = lines[i+2+k]
                        split = tmrna_line.split()
                        start_stops = split[2].replace("[", "").replace("]", "").split(',')
                        contig = length_df["contig"][j]
                        method = "Aragorn"
                        region = "tmRNA"
                        start = start_stops[0].replace("c", "")  # tmrna output is [start,stop] or c[start, stop] so need to remove c also for some phages
                        stop = start_stops[1]
                        score = "."
                        frame = "."
                        phase = "."
                        tag_peptide = split[3]
                        tag_peptide_seq = split[4]
                        attribute = "product=transfer-messenger RNA SsrA;tag_peptide=" + tag_peptide + ";tag_peptide_sequence=" + tag_peptide_seq
                        contig_names.append(contig)
                        methods.append(method)
                        regions.append(region)
                        starts.append(start)
                        stops.append(stop)
                        scores.append(score)
                        frames.append(frame)
                        phases.append(phase)
                        attributes.append(attribute)
                j +=1 # iterate contig
            # iterate line
            i += 1 
        # write out the df
        tmrna_df = pd.DataFrame(
        {'contig': contig_names,
        'Method': methods,
        'Region': regions,
        'start': starts,
        'stop': stops,
        'score': scores,
        'frame': frames,
        'phase': phases,
        'attributes': attributes,
        })
    tmrna_df.to_csv(os.path.join(out_dir, prefix + "_aragorn.gff"), sep="\t", index=False, header=False)
    return tmrna_flag


#### process vfdb files
def process_vfdb_results(out_dir, prefix):
    ##vfdb
    vfdb_file = os.path.join(out_dir, "vfdb_results.tsv")
    print("Processing vfdb output.")
    col_list = ["vfdb_hit", "gene", "vfdb_alnScore", "vfdb_seqIdentity", "vfdb_eVal", "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen"] 
    vfdb_df = pd.read_csv(vfdb_file, delimiter= '\t', index_col=False , names=col_list) 
    genes = vfdb_df.gene.unique()

    # get top hit 
    tophits = []

    for gene in genes:
        tmp_df = vfdb_df.loc[vfdb_df['gene'] == gene].sort_values('vfdb_eVal').reset_index(drop=True).loc[0]
        tophits.append([tmp_df.vfdb_hit, tmp_df.gene, tmp_df.vfdb_alnScore, tmp_df.vfdb_seqIdentity, tmp_df.vfdb_eVal])

    tophits_df = pd.DataFrame(tophits, columns=['vfdb_hit', 'gene', 'vfdb_alnScore', 'vfdb_seqIdentity', 'vfdb_eVal'])
    tophits_df.to_csv(os.path.join(out_dir, "top_hits_vfdb.tsv"), sep="\t", index=False)

    # left join vfdb to final_output_file
    final_output_file = os.path.join(out_dir, prefix + "_final_merged_output.tsv")
    # automatically picks up the names
    final_df = pd.read_csv(final_output_file, sep="\t", index_col=False)
    final_df['gene']=final_df['gene'].astype(str)
    tophits_df['gene']=tophits_df['gene'].astype(str)

    # merge top hits into the merged df
    merged_df = final_df.merge(tophits_df, on='gene', how='left')
    merged_df["vfdb_hit"] = merged_df["vfdb_hit"].replace(nan, 'None', regex=True)
    merged_df["vfdb_alnScore"] = merged_df["vfdb_alnScore"].replace(nan, 'None', regex=True)
    merged_df["vfdb_seqIdentity"] = merged_df["vfdb_seqIdentity"].replace(nan, 'None', regex=True)
    merged_df["vfdb_eVal"] = merged_df["vfdb_eVal"].replace(nan, 'None', regex=True)
    # move around
    merged_df[['genbank','desc_tmp', 'vfdb_species']] = merged_df['vfdb_hit'].str.split('[',expand=True)
    merged_df['vfdb_species'] = merged_df['vfdb_species'].str.replace("]", "")
    merged_df['vfdb_species'] = merged_df['vfdb_species'].str.strip()
    merged_df[['genbank','vfdb_short_name', 'vfdb_description']] = merged_df['genbank'].str.split(')',expand=True)
    merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].str.replace("(", "")
    merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].str.strip()
    merged_df["vfdb_description"] = merged_df["vfdb_description"].str.strip()
    merged_df = merged_df.drop(columns = ['genbank', 'desc_tmp'])
    merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].replace(nan, 'None', regex=True)
    merged_df["vfdb_description"] = merged_df["vfdb_description"].replace(nan, 'None', regex=True)
    merged_df["vfdb_species"] = merged_df["vfdb_species"].replace(nan, 'None', regex=True)
    merged_df.to_csv( os.path.join(out_dir, prefix + "_final_merged_output.tsv"), sep="\t", index=False)

    return merged_df
    


