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

def create_mmseqs_tophits(out_dir):
    """
    """

    ##mmseqs
    mmseqs_file =  os.path.join(out_dir, "mmseqs_results.tsv")
    print("Processing mmseqs2 output.")
    col_list = ["phrog", "gene", "alnScore", "seqIdentity", "eVal", "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen"] 
    mmseqs_df = pd.read_csv(mmseqs_file, delimiter= '\t', index_col=False , names=col_list) 
    # get list of genes
    genes = mmseqs_df.gene.unique()

    # instantiate tophits list
    tophits = []

    for gene in genes:
        tmp_df = mmseqs_df.loc[mmseqs_df['gene'] == gene].sort_values('eVal').reset_index(drop=True).loc[0]
        tophits.append([tmp_df.phrog, tmp_df.gene, tmp_df.alnScore, tmp_df.seqIdentity, tmp_df.eVal])

    # create tophits df
    tophits_df = pd.DataFrame(tophits, columns=['phrog', 'gene', 'alnScore', 'seqIdentity', 'eVal'])
    tophits_df.to_csv(os.path.join(out_dir, "top_hits_mmseqs.tsv"), sep="\t", index=False)
    return tophits_df


def process_results(db_dir,out_dir, prefix, gene_predictor):
    """
    Processes and combines PHROGS, CARD and VFDB mmseqs output
    :param db_dir: database directory path
    :param out_dir: output directory path 
    :param prefix: output prefix
    :param gene_predictor: CDS predictor (phanotate or prodigal)
    :return: merged_df a pandas dataframe of the final merged output
    """

    tophits_df = create_mmseqs_tophits(out_dir)

    # left join mmseqs top hits to cds df
    cds_file = os.path.join(out_dir, "cleaned_" +  gene_predictor +  ".tsv") 
    cds_df = pd.read_csv(cds_file, sep="\t", index_col=False )
    # convert the gene to string for the merge
    cds_df['gene']=cds_df['gene'].astype(str)
    tophits_df['gene']=tophits_df['gene'].astype(str)
    cds_df = cds_df[cds_df['start'].notna()]
    cds_df = cds_df.dropna()

    # merge top hits into the cds df
    merged_df = cds_df.merge(tophits_df, on='gene', how='left')

    # get phrog from top hit
    # add test if empty - crashes if no gene call hits
    if len(tophits_df['phrog']) == 0:
        merged_df['top_hit'] = 'No_PHROG'
    else:
        merged_df[['phrog','top_hit']] = merged_df['phrog'].str.split(' ## ',expand=True)
    # strip off phrog_
    merged_df["phrog"] = merged_df["phrog"].str.replace("phrog_", "")
    
    # read in phrog annotaion file
    phrog_annot_df = pd.read_csv( os.path.join(db_dir, "phrog_annot_v4.tsv"), sep="\t", index_col=False )
    phrog_annot_df['phrog']=phrog_annot_df['phrog'].astype(str)

    # merge phrog
    merged_df = merged_df.merge(phrog_annot_df, on='phrog', how='left')
    merged_df = merged_df.replace(np.nan, 'No_PHROG', regex=True)
    # convert no phrog to hyp protein
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
    merged_df['phrog'] = merged_df['phrog'].astype(str)

    # drop existing color annot category cols
    merged_df = merged_df.drop(columns = ['color', 'annot', 'category'])
    merged_df = merged_df.merge(phrog_annot_df, on='phrog', how='left')
    merged_df["annot"] = merged_df["annot"].replace(nan, 'hypothetical protein', regex=True)
    merged_df["category"] = merged_df["category"].replace(nan, 'unknown function', regex=True)
    merged_df["color"] = merged_df["color"].replace(nan, 'None', regex=True)

    # process vfdb results
    (merged_df, vfdb_results) = process_vfdb_results(out_dir, merged_df)
    # process CARD results
    (merged_df, card_results) = process_card_results(out_dir, merged_df, db_dir)

    return (merged_df,vfdb_results, card_results)

def get_contig_name_lengths(fasta_input):
    """
    Gets contig name and length in the input fasta file and calculates gc
    :param fasta_input: input fasta file
    :return: length_df a pandas dataframe 
    """
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

def create_txt(cds_mmseqs_merge_df, length_df, out_dir, prefix):
    """
    Creates the _cds_functions.tsv and _length_gc_cds_density.tsv outputs
    :param cds_mmseqs_merge_df: a pandas df as output from process_results()
    :param length_df: a pandas df as output from get_contig_name_lengths()
    :param out_dir: output directory path
    :param prefix: output prefix
    :return:
    """

    # get contigs - convert to string to make the match work for integer contigs
    contigs = length_df["contig"].astype("string")
    # convert to string to make the match work for integer contigs
    cds_mmseqs_merge_df['contig'] = cds_mmseqs_merge_df['contig'].astype("string")

    # instantiate the length_df['cds_coding_density']
    length_df['cds_coding_density'] = 0.0

    # list of all dataframes with functions (for later)
    combo_list = []

    #### trna scan 
    # read in trnascan
    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
    trna_df = pd.read_csv(os.path.join(out_dir,"trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
    # keep only trnas and pseudogenes
    trna_df = trna_df[(trna_df['Region'] == 'tRNA') | (trna_df['Region'] == 'pseudogene')]

    #### crispr
    crispr_df = pd.read_csv(os.path.join(out_dir, prefix + "_minced.gff"), delimiter= '\t', index_col=False, names=col_list, comment = '#' ) 

    #### tmrna
    tmrna_df = pd.read_csv(os.path.join(out_dir, prefix + "_aragorn.gff"), delimiter= '\t', index_col=False, names=col_list ) 

    #### vfdb
    # gene has the contig 
    try:
        vfdb_df = pd.read_csv(os.path.join(out_dir, "top_hits_vfdb.tsv"), delimiter= '\t', index_col=False ) 
        vfdb_df[['gene','coordinate']] = vfdb_df['gene'].str.split(' ',expand=True)
        # strip off the index with the last period delimiter
        vfdb_df['contig'] = vfdb_df['gene'].str.rsplit("\\.",1)
    except:
        # instantiate empty dataframe if there are no hits
        vfdb_df = pd.DataFrame(columns=['vfdb_hit', 'gene', 'coordinate', 'contig', 'vfdb_alnScore', 'vfdb_seqIdentity', 'vfdb_eVal'])

    #### card
    # gene has the contig 
    try:
        card_df = pd.read_csv(os.path.join(out_dir, "top_hits_card.tsv"), delimiter= '\t', index_col=False ) 
        card_df[['gene','coordinate']] = card_df['gene'].str.split(' ',expand=True)
        # strip off the index with the last period delimiter
        card_df['contig'] = card_df['gene'].str.rsplit("\\.",1)
    except:
        # instantiate empty dataframe if there are no hits
        card_df = pd.DataFrame(columns=['CARD_hit', 'gene', 'coordinate', 'contig', 'CARD_alnScore', 'CARD_seqIdentity', 'CARD_eVal'])
    

    # write descriptions for each contig
    for contig in contigs:
        # get cds's in the contig
        cds_mmseqs_merge_cont_df = cds_mmseqs_merge_df[cds_mmseqs_merge_df['contig'] == contig]
        # counts of the cds trnas
        cds_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['Region'] == 'CDS'])
        trna_count = len(trna_df[trna_df['contig'] == contig])
        tmrna_count = len(tmrna_df[tmrna_df['contig'] == contig])
        crispr_count = len(crispr_df[crispr_df['contig'] == contig])
        if len(vfdb_df['contig']) != 0:
            vfdb_count = len(vfdb_df[vfdb_df['contig'] == contig])
            vfdb_count = vfdb_df['gene'].str.contains(contig).sum()
        else:
            vfdb_count = 0
        if len(card_df['contig']) != 0:
            CARD_count = len(card_df[card_df['contig'] == contig])
            CARD_count = card_df['gene'].str.contains(contig).sum()
        else:
            CARD_count = 0
        # get the total length of the contig
        contig_length = length_df[length_df["contig"] == contig]['length']
        if cds_count > 0:
            # gets the total cds coding length
            cds_lengths = abs(cds_mmseqs_merge_cont_df['start'] - cds_mmseqs_merge_cont_df['stop']).sum()
            # get function
            cds_mmseqs_merge_cont_df[['attributes2']] = cds_mmseqs_merge_cont_df[['attributes']]
            cds_mmseqs_merge_cont_df[['attributes2','function']] = cds_mmseqs_merge_cont_df['attributes2'].str.split(';function=',expand=True)
            cds_mmseqs_merge_cont_df = cds_mmseqs_merge_cont_df.drop(columns=['attributes2'])
            cds_mmseqs_merge_cont_df[['function','product']] = cds_mmseqs_merge_cont_df['function'].str.split(';product=',expand=True)
            cds_mmseqs_merge_cont_df = cds_mmseqs_merge_cont_df.drop(columns=['product'])
            
            # get counts of functions and cds 
            # all 10 PHROGs categories
            # modified for v1.0.0 in case some have 0s (to make consistent for downstream)
            connector_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'connector'])
            metabolism_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'DNA, RNA and nucleotide metabolism'])
            head_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'head and packaging'])
            integration_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'integration and excision'])
            lysis_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'lysis'])
            moron_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'moron, auxiliary metabolic gene and host takeover'])
            other_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'other'])
            tail_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'tail'])
            transcription_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'transcription regulation'])
            unknown_count = len(cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df['function'] == 'unknown function'])
            # create count list  for the dataframe 
            count_list = [cds_count, connector_count, metabolism_count, head_count, integration_count, lysis_count,
            moron_count, other_count, tail_count, transcription_count, unknown_count]
        else:
            cds_lengths = 0

            count_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        description_list = ['CDS', 'connector','DNA, RNA and nucleotide metabolism',  'head and packaging', 'integration and excision','lysis',
        'moron, auxiliary metabolic gene and host takeover','other', 'tail', 'transcription regulation', 'unknown function']
        contig_list = [contig, contig, contig, contig, contig, contig, contig, contig, contig, contig, contig]
        # cds df
        cds_df = pd.DataFrame({'Description': description_list, 'Count': count_list, 'contig': contig_list})

        # add other features 
        trna_row = pd.DataFrame({ 'Description':['tRNAs'], 'Count':[trna_count], 'contig':[contig] })
        crispr_row = pd.DataFrame({ 'Description':['CRISPRs'], 'Count':[crispr_count], 'contig':[contig] })
        tmrna_row = pd.DataFrame({ 'Description':['tmRNAs'], 'Count':[tmrna_count], 'contig':[contig] })
        vfdb_row = pd.DataFrame({ 'Description':['VFDB_Virulence_Factors'], 'Count':[vfdb_count], 'contig':[contig] })
        CARD_row = pd.DataFrame({ 'Description':['CARD_AMR_Genes'], 'Count':[CARD_count], 'contig':[contig] })
        # calculate the cds coding density and add to length_df
        cds_coding_density = cds_lengths * 100 / contig_length
        cds_coding_density = round(cds_coding_density, 2)
        length_df.loc[length_df['contig'] == contig, 'cds_coding_density'] = cds_coding_density
        # eappend it all to combo_list
        combo_list.append(cds_df)
        combo_list.append(trna_row)
        combo_list.append(crispr_row)
        combo_list.append(tmrna_row)
        combo_list.append(vfdb_row)
        combo_list.append(CARD_row)

    # combine all contigs into one final df
    description_total_df = pd.concat(combo_list)
    # save the output
    description_total_df.to_csv(os.path.join(out_dir, prefix + "_cds_functions.tsv"), sep="\t", index=False)
    # save the length_gc.tsv also
    length_df.to_csv(os.path.join(out_dir, prefix + "_length_gc_cds_density.tsv"), sep="\t", index=False)


  
def create_gff(cds_mmseqs_df, length_df, fasta_input, out_dir, prefix, locustag, tmrna_flag, meta_mode):
    """
    Creates the pharokka.gff file 
    :param cds_mmseqs_merge_df: a pandas df as output from process_results()
    :param length_df: a pandas df as output from get_contig_name_lengths()
    :param fasta_input: input fasta file
    :param out_dir: output directory path
    :param prefix: output prefix
    :param locustag: whether or not to create a random locustag - will be Random is so. Otherwise it is parsed
    :tmrna_flag boolean whether there are tmRNAs or not
    :return: locustag for the creation of the .tbl file, locus_df as df with consistent locus_tags
    """

    # create locus tag and get temp df

    # locustag creation
    if locustag == "Random":
        # locus tag header 8 random letters
        locustag = ''.join(random.choice(string.ascii_uppercase) for _ in range(8))
    
    contigs = length_df["contig"].astype("string")

    ############ locus tag #########
    # write df for locus tag parsing

    # zfill - makes the CDS 4 digits trailing zeroes for vcontact
    # in meta mode, 

    locus_df = cds_mmseqs_df
    subset_dfs = []

    # if meta mode is true
    if meta_mode == True:
        for contig in contigs:
            subset_df = locus_df[locus_df['contig'] == contig].reset_index()
            subset_df['count'] = subset_df.index
            # so not 0 indexed
            subset_df['count'] = subset_df['count'] + 1 
            # z fill to make the locus tag 4
            subset_df['count'] = subset_df['count'].astype(str).str.zfill(4)
            subset_dfs.append(subset_df)
        locus_df = pd.concat(subset_dfs, axis=0, ignore_index=True)
    else:
        locus_df['count'] = locus_df.index
        # so not 0 indexed
        locus_df['count'] = locus_df['count'] + 1 
        # z fill to make the locus tag 4
        locus_df['count'] = locus_df['count'].astype(str).str.zfill(4)

    # get the 

    if meta_mode == False:
        locus_df['locus_tag'] = locustag + "_CDS_" + locus_df['count']
    else:
        locus_df['locus_tag'] = locus_df.contig + "_CDS_" + locus_df['count']

    #################################

    #########
    # rearrange start and stop so that start is always less than stop for gff
    #########

    cols = ["start","stop"]
    #indices where start is greater than stop
    ixs = cds_mmseqs_df['frame'] == '-'
    # Where ixs is True, values are swapped
    cds_mmseqs_df.loc[ixs,cols] = cds_mmseqs_df.loc[ixs, cols].reindex(columns=cols[::-1]).values

    # set phase to be 0
    cds_mmseqs_df['phase'] = 0
    # create attributes
    cds_mmseqs_df['attributes'] = "ID=" + locus_df['locus_tag'].astype(str)  + ";" + "phrog=" + cds_mmseqs_df["phrog"].astype(str) + ";" + "top_hit=" + cds_mmseqs_df["top_hit"].astype(str) + ";" + "locus_tag=" + locus_df['locus_tag'].astype(str) + ";" + "function=" + cds_mmseqs_df["category"].astype(str) + ";"  + "product=" + cds_mmseqs_df["annot"].astype(str)
    # adds VFDB
    cds_mmseqs_df.loc[cds_mmseqs_df['vfdb_short_name'] != "None", 'attributes'] = cds_mmseqs_df['attributes'].astype(str)  + ";"  + "vfdb_short_name=" + cds_mmseqs_df['vfdb_short_name'].astype(str) + ";"  + "vfdb_description="  +  cds_mmseqs_df['vfdb_description'].astype(str) + ";" + "vfdb_species=" + cds_mmseqs_df['vfdb_species'].astype(str)
    # adds CARD
    cds_mmseqs_df.loc[cds_mmseqs_df['CARD_short_name'] != "None", 'attributes'] = cds_mmseqs_df['attributes'].astype(str)  + ";"  + "CARD_short_name=" + cds_mmseqs_df['CARD_short_name'].astype(str) + ";"  + "AMR_Gene_Family="  +  cds_mmseqs_df['AMR_Gene_Family'].astype(str) + ";" + "CARD_species=" + cds_mmseqs_df['CARD_species'].astype(str)

    # save back 

    # get gff dataframe in correct order 
    gff_df = cds_mmseqs_df[["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]]
    
    # change start and stop to int 
    gff_df["start"] = gff_df["start"].astype('int')
    gff_df["stop"] = gff_df["stop"].astype('int')

    with open(os.path.join(out_dir, "pharokka_tmp.gff"), 'a') as f:
        gff_df.to_csv(f, sep="\t", index=False, header=False)

    ### trnas
    # check if no trnas
    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
    trna_empty = is_trna_empty(out_dir)
    if trna_empty == False:   
        trna_df = pd.read_csv(os.path.join(out_dir,"trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
        # index hack if meta mode
        if meta_mode == True:
            subset_dfs = []
            for contig in contigs:
                subset_df = trna_df[trna_df['contig'] == contig].reset_index()
                # keep only trnas before indexing
                subset_df = subset_df[(subset_df['Region'] == 'tRNA') | (subset_df['Region'] == 'pseudogene')]
                subset_df = subset_df.reset_index(drop=True)
                subset_df['count'] = subset_df.index
                # so not 0 indexed
                subset_df['count'] = subset_df['count'] + 1 
                # z fill to make the locus tag 4
                subset_df['locus_tag'] = contig + "_tRNA_" + subset_df['count'].astype(str).str.zfill(4)
                subset_df = subset_df.drop(columns=['count'])
                subset_dfs.append(subset_df)
            trna_df = pd.concat(subset_dfs, axis=0, ignore_index=True) 
            trna_df = trna_df.drop(columns=['index'])
        else:
            # keep only trnas
            trna_df = trna_df[(trna_df['Region'] == 'tRNA') | (trna_df['Region'] == 'pseudogene')]
            trna_df = trna_df.reset_index(drop=True)
            trna_df['count'] = trna_df.index 
            trna_df['count'] = trna_df['count'] + 1 
            trna_df['locus_tag'] = locustag + "_tRNA_" + trna_df['count'].astype(str).str.zfill(4)    
            trna_df = trna_df.drop(columns=['count'])

        trna_df.start = trna_df.start.astype(int)
        trna_df.stop = trna_df.stop.astype(int)
        trna_df[['attributes','isotypes']] = trna_df['attributes'].str.split(';isotype=',expand=True)
        trna_df[['isotypes','anticodon']] = trna_df['isotypes'].str.split(';anticodon=',expand=True)
        trna_df[['anticodon','rest']] = trna_df['anticodon'].str.split(';gene_biotype',expand=True)
        trna_df['trna_product']='tRNA-'+trna_df['isotypes']+"("+trna_df['anticodon']+")"
        trna_df = trna_df.drop(columns=['attributes'])
        trna_df['attributes'] = "ID=" + trna_df['locus_tag'] + ";" + "trna=" + trna_df["trna_product"].astype(str) + ";" + "isotype=" + trna_df["isotypes"].astype(str) + ";" + "anticodon=" + trna_df["anticodon"].astype(str) + ";" + "locus_tag="  + trna_df['locus_tag']
        trna_df = trna_df.drop(columns=['isotypes', 'anticodon', 'rest', 'trna_product', 'locus_tag'])
        # append to the end of the gff
        with open(os.path.join(out_dir, "pharokka_tmp.gff"), 'a') as f:
            trna_df.to_csv(f, sep="\t", index=False, header=False)

    ### crisprs
    crispr_count = get_crispr_count(out_dir, prefix)
    # add to gff if > 0 
    if crispr_count > 0:
        minced_df = pd.read_csv(os.path.join(out_dir, prefix + "_minced.gff"), delimiter= '\t', index_col=False, names=col_list, comment='#' ) 
        minced_df.start = minced_df.start.astype(int)
        minced_df.stop = minced_df.stop.astype(int)
        minced_df[['attributes','rpt_unit_seq']] = minced_df['attributes'].str.split(';rpt_unit_seq=',expand=True)
        minced_df[['attributes','rpt_family']] = minced_df['attributes'].str.split(';rpt_family=',expand=True)
        minced_df[['attributes','rpt_type']] = minced_df['attributes'].str.split(';rpt_type=',expand=True)
        minced_df = minced_df.drop(columns=['attributes'])
        # index hack if meta mode
        subset_dfs = []
        if meta_mode == True:
            for contig in contigs:
                subset_df = minced_df[minced_df['contig'] == contig].reset_index()
                subset_df['count'] = subset_df.index
                # so not 0 indexed
                subset_df['count'] = subset_df['count'] + 1 
                # z fill to make the locus tag 4
                subset_df['count'] = subset_df['count'].astype(str).str.zfill(4)
                subset_df['locus_tag'] = contig + "_CRISPR_" + subset_df['count'].astype(str).str.zfill(4)
                subset_df = subset_df.drop(columns=['count'])
                subset_dfs.append(subset_df)
            minced_df = pd.concat(subset_dfs, axis=0, ignore_index=True)
            minced_df = minced_df.drop(columns=['index'])
        else:
            minced_df['count'] = minced_df.index 
            minced_df['count'] = minced_df['count'] + 1 
            minced_df['locus_tag'] = locustag + "_CRISPR_" + minced_df['count'].astype(str).str.zfill(4)    
            minced_df = minced_df.drop(columns=['count'])

        minced_df['attributes'] = "ID=" + minced_df['locus_tag']   + ";" + "rpt_type=" + minced_df["rpt_type"].astype(str) + ";" + "rpt_family=" + minced_df["rpt_family"].astype(str) + ";" + "rpt_unit_seq=" + minced_df["rpt_unit_seq"].astype(str) + ";" + "locus_tag=" + minced_df['locus_tag'] 
        minced_df = minced_df.drop(columns=['rpt_unit_seq', 'rpt_family', 'rpt_type', 'locus_tag'])
        # append to the end of the gff
        with open(os.path.join(out_dir, "pharokka_tmp.gff"), 'a') as f:
            minced_df.to_csv(f, sep="\t", index=False, header=False)

    ### tmrna
     # add to gff there is a tmrna
    if tmrna_flag == True:
        tmrna_df = pd.read_csv(os.path.join(out_dir, prefix + "_aragorn.gff"), delimiter= '\t', index_col=False, names=col_list ) 
        tmrna_df.start = tmrna_df.start.astype(int)
        tmrna_df.stop = tmrna_df.stop.astype(int)

        # index hack if meta mode
        subset_dfs = []
        if meta_mode == True:
            for contig in contigs:
                subset_df = tmrna_df[tmrna_df['contig'] == contig].reset_index()
                subset_df['count'] = subset_df.index
                # so not 0 indexed
                subset_df['count'] = subset_df['count'] + 1 
                # z fill to make the locus tag 4
                subset_df['count'] = subset_df['count'].astype(str).str.zfill(4)
                subset_df['locus_tag'] = contig + "_tmRNA_" + subset_df['count'].astype(str).str.zfill(4)
                subset_df = subset_df.drop(columns=['count'])
                subset_dfs.append(subset_df)
            tmrna_df = pd.concat(subset_dfs, axis=0, ignore_index=True)
            tmrna_df = tmrna_df.drop(columns=['index'])
        else:  
            tmrna_df['count'] = tmrna_df.index 
            tmrna_df['count'] = tmrna_df['count'] + 1 
            tmrna_df['locus_tag'] = locustag + "_tmRNA_" + tmrna_df['count'].astype(str).str.zfill(4)    
            tmrna_df = tmrna_df.drop(columns=['count'])


        tmrna_df['attributes'] = "ID=" + tmrna_df['locus_tag']   + ";" + tmrna_df['attributes'].astype(str) + ';locus_tag=' + tmrna_df['locus_tag']  
        tmrna_df = tmrna_df.drop(columns=['locus_tag'])
        # append to the end of the gff
        with open(os.path.join(out_dir, "pharokka_tmp.gff"), 'a') as f:
            tmrna_df.to_csv(f, sep="\t", index=False, header=False)

    # write header of final gff files 
    with open(os.path.join(out_dir, prefix + ".gff"), 'w') as f:
        f.write('##gff-version 3\n')
        for index, row in length_df.iterrows():
            f.write('##sequence-region ' + row['contig'] + ' 1 ' + str(row['length']) +'\n')
    
    # sort gff by start
    total_gff = pd.read_csv(os.path.join(out_dir, "pharokka_tmp.gff"), delimiter= '\t', index_col=False, names=col_list, low_memory=False ) 
    total_gff.start = total_gff.start.astype(int)
    total_gff.stop = total_gff.stop.astype(int)
    # sorts all features
    total_gff = total_gff.groupby(['contig'], sort=False, as_index=False).apply(pd.DataFrame.sort_values, 'start', ascending =True)

    # write final gff to file
    with open(os.path.join(out_dir, prefix + ".gff"), 'a') as f:
        total_gff.to_csv(f, sep="\t", index=False, header=False)

    # write fasta on the end
    with open(os.path.join(out_dir, prefix + ".gff"), 'a') as f:
        f.write('##FASTA\n')
        fasta_sequences = SeqIO.parse(open(fasta_input),'fasta')
        SeqIO.write(fasta_sequences, f, "fasta")

    return (locustag, locus_df, gff_df)


def update_fasta_headers(locus_df, out_dir, gene_predictor ):
    """
    Updates the fasta output headers to have a consistent locus tag & gene description for downstrea use
    :param locus_df a pandas df as output from create_gff()
    :param out_dir: output directory path
    :gene_predictor: string 'phanotate' or 'prodigal' with the gene predictor used
    """

    #define outputs
    fasta_input_nts_tmp = gene_predictor + "_out_tmp.fasta"
    fasta_input_aas_tmp = gene_predictor + "_aas_tmp.fasta"
    fasta_output_nts_gd = gene_predictor + ".ffn"
    fasta_output_aas_gd = gene_predictor + ".faa"

    # nucleotides

    with open(os.path.join(out_dir, fasta_output_nts_gd), 'w') as nt_fa:   
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_nts_tmp), 'fasta'): 
            dna_record.id = str(locus_df['locus_tag'].iloc[i]) 
            dna_record.description = str(locus_df['annot'].iloc[i])
            SeqIO.write(dna_record, nt_fa, 'fasta')
            i += 1

    # amino acids

    with open(os.path.join(out_dir, fasta_output_aas_gd), 'w') as aa_fa:   
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_aas_tmp), 'fasta'): 
            dna_record.id = str(locus_df['locus_tag'].iloc[i]) 
            dna_record.description = str(locus_df['annot'].iloc[i])
            SeqIO.write(dna_record, aa_fa, 'fasta')
            i += 1


def extract_terl(locus_df, out_dir, gene_predictor, logger ):
    """
    Extract large terminase subunit
    :param locus_df a pandas df as output from create_gff()
    :param out_dir: output directory path
    :gene_predictor: string 'phanotate' or 'prodigal' with the gene predictor used
    :logger: logger
    """

    #phanotate
    fasta_input_nts_tmp = gene_predictor + "_out_tmp.fasta"
    fasta_input_aas_tmp = gene_predictor + "_aas_tmp.fasta"

    # nucleotide

    with open(os.path.join(out_dir, "terL.ffn"), 'w') as aa_fa:   
        # i loops over the rows of the dataframe
        # j counts the number of terLs
        i = 0 
        j = 0
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_nts_tmp), 'fasta'): 
            dna_record.id = str(locus_df['locus_tag'].iloc[i]) 
            dna_record.description = str(locus_df['annot'].iloc[i])
            if locus_df['annot'].iloc[i] == "terminase large subunit":
                SeqIO.write(dna_record, aa_fa, 'fasta')
                # report terL found the first time
                if j < 1:
                    print('Terminase large subunit found.')
                    logger.info("Terminase large subunit found.")
                    # report multiple found the second time
                if j == 1:
                    print('More than one CDS annotated as terminase large subunit found. \nSaving all.') 
                    logger.info('More than one CDS annotated as terminase large subunit found. \nSaving all.') 
                j += 1
            i += 1


    # amino acid no need to print

    with open(os.path.join(out_dir, "terL.faa"), 'w') as aa_fa:   
        # i loops over the rows of the dataframe
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_aas_tmp), 'fasta'): 
            dna_record.id = str(locus_df['locus_tag'].iloc[i]) 
            dna_record.description = str(locus_df['annot'].iloc[i])
            if locus_df['annot'].iloc[i] == "terminase large subunit":
                SeqIO.write(dna_record, aa_fa, 'fasta')
            i += 1


def update_final_output(cds_mmseqs_merge_df, vfdb_results, card_results, locus_df, prefix, out_dir ):
    """
    Updates the fasta output headers to have a consistent locus tag & gene description for downstrea use
    :param cds_mmseqs_merge_df: a pandas df as output from process_results()
    :param locus_df: a pandas df as output from create_gff()
    :param out_dir: output directory path
    :param prefix: output prefix
    :return: 
    """

    # needed for vfdb card matching later
    genes_for_vfdb_card = cds_mmseqs_merge_df['gene']

    # return back the cds_mmseqs_merge_df but with the locus tag instead of gene
    # rename gene with locus_tag
    locus_tag_series = locus_df['locus_tag']
    cds_mmseqs_merge_df['gene'] = locus_tag_series


    #########
    # rearrange start and stop for neg strang
    #########

    st_cols = ["start","stop"]
    #indices where start is greater than stop
    ixs = cds_mmseqs_merge_df['frame'] == '-'
    # Where ixs is True, values are swapped
    cds_mmseqs_merge_df.loc[ixs,st_cols] = cds_mmseqs_merge_df.loc[ixs, st_cols].reindex(columns=st_cols[::-1]).values


    # get a list of columns
    cols = list(cds_mmseqs_merge_df)
    # move the column to head of list using index, pop and insert
    cols.insert(0, cols.pop(cols.index('gene')))
    cds_mmseqs_merge_df = cds_mmseqs_merge_df.loc[:, cols]
    # drop last 2 cols 
    cds_mmseqs_merge_df = cds_mmseqs_merge_df.drop(columns=['phase', 'attributes'])

    # write output
    final_output_path = os.path.join(out_dir, prefix + "_cds_final_merged_output.tsv")
    cds_mmseqs_merge_df.to_csv( final_output_path, sep="\t", index=False)

    ######################################
    ##### update vfdb with locus tag #####
    ###### merge locus into the vfdb ####
    locus_df['gene'] = genes_for_vfdb_card
    vfdb_results = vfdb_results.merge(locus_df, how='left', on='gene')
    # get a list of columns
    cols = list(vfdb_results)
    # move the column to head of list using index, pop and insert
    cols.insert(0, cols.pop(cols.index('locus_tag')))
    vfdb_results = vfdb_results.loc[:, cols]
    # keep only desired columns  and save
    vfdb_results = vfdb_results[['locus_tag', 'vfdb_hit_x', 'vfdb_alnScore_x', 'vfdb_seqIdentity_x', 'start', 'stop', 'frame']]
    vfdb_results.columns = ['gene', 'vfdb_hit', 'vfdb_alnScore', 'vfdb_seqIdentity', 'start', 'stop', 'frame']
    vfdb_results = vfdb_results.sort_values(by=['start'])
    vfdb_results.to_csv(os.path.join(out_dir, "top_hits_vfdb.tsv"), sep="\t", index=False)

    ######################################
    ##### update card with locus tag #####
    ###### merge locus into the card ####
    card_results = card_results.merge(locus_df, how='left', on='gene')
    # get a list of columns
    cols = list(card_results)
    # move the column to head of list using index, pop and insert
    cols.insert(0, cols.pop(cols.index('locus_tag')))
    card_results = card_results.loc[:, cols]

    # keep only desired columns   sand save
    card_results = card_results[['locus_tag', 'CARD_hit_x', 'CARD_alnScore_x', 'CARD_seqIdentity_x', 'start', 'stop', 'frame']]
    card_results.columns = ['gene', 'card_hit', 'card_alnScore', 'card_seqIdentity', 'start', 'stop', 'frame']
    card_results = card_results.sort_values(by=['start'])
    card_results.to_csv(os.path.join(out_dir, "top_hits_card.tsv"), sep="\t", index=False)



def create_tbl(cds_mmseqs_df, length_df, out_dir, prefix, gene_predictor, tmrna_flag, gff_df):
    """
    Creates the pharokka.tbl file 
    :param cds_mmseqs_df: a pandas df as output from process_results()
    :param length_df: a pandas df as output from get_contig_name_lengths()
    :param out_dir: output directory path
    :param prefix: output prefix
    :param gene_predictor: phanotate or prodigal
    :param locustag: output from create_gff()
    :tmrna_flag boolean whether there are tmRNAs or not
    :return: 
    """

    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]

    # read in the gff
    total_gff = pd.read_csv(os.path.join(out_dir, "pharokka_tmp.gff"), delimiter= '\t', index_col=False, names=col_list, low_memory=False ) 

    if gene_predictor == "phanotate":
        cds_df = total_gff[total_gff['Method'] == 'PHANOTATE']
    else:
        cds_df = total_gff[total_gff['Method'] == 'PRODIGAL']

    cds_df[['attributes','locus_tag']] = cds_df['attributes'].str.split(';locus_tag=',expand=True)
    cds_df[['locus_tag','rest']] = cds_df['locus_tag'].str.split(';function=',expand=True)
    cds_mmseqs_df['locus_tag'] = cds_df['locus_tag']


    ### trnas
    # check if no trnas
    trna_empty = is_trna_empty(out_dir)
    if trna_empty == False:    
        trna_df = total_gff[total_gff['Method'] == 'tRNAscan-SE']
        # keep only trnas and pseudogenes 
        trna_df.start = trna_df.start.astype(int)
        trna_df.stop = trna_df.stop.astype(int)
        trna_df[['attributes','locus_tag']] = trna_df['attributes'].str.split(';locus_tag=',expand=True)
        trna_df[['attributes','isotypes']] = trna_df['attributes'].str.split(';isotype=',expand=True)
        trna_df[['isotypes','anticodon']] = trna_df['isotypes'].str.split(';anticodon=',expand=True)
        trna_df['trna_product']='tRNA-'+trna_df['isotypes']+"("+trna_df['anticodon']+")"

    #### CRISPRs
    # check if the file has more than 1 line (not empty)
    crispr_count = get_crispr_count(out_dir, prefix)
    if crispr_count > 0:
        crispr_df = total_gff[total_gff['Region'] == 'repeat_region']
        crispr_df.start = crispr_df.start.astype(int)
        crispr_df.stop = crispr_df.stop.astype(int)
        crispr_df[['attributes','locus_tag']] = crispr_df['attributes'].str.split(';locus_tag=',expand=True)
        crispr_df[['attributes','rpt_unit_seq']] = crispr_df['attributes'].str.split(';rpt_unit_seq=',expand=True)

    ### TMRNAs
    if tmrna_flag == True:
        tmrna_df = total_gff[total_gff['Region'] == 'tmRNA']
        tmrna_df.start = tmrna_df.start.astype(int)
        tmrna_df.stop = tmrna_df.stop.astype(int)
        tmrna_df[['attributes','locus_tag']] = tmrna_df['attributes'].str.split(';locus_tag=',expand=True)
    
    # set inference
    if gene_predictor == "phanotate":
        inf = "PHANOTATE"
    else:
        inf = "PRODIGAL"
    with open( os.path.join(out_dir, prefix + ".tbl"), 'w') as f:
        for index, row in length_df.iterrows():
            contig = row['contig']
            f.write('>Feature ' + contig + '\n')
            subset_df = cds_mmseqs_df[cds_mmseqs_df['contig'] == contig]
            for index, row in subset_df.iterrows():
                start = str(row['start'])
                stop = str(row['stop'])
                if row['frame'] == '-':
                    start = str(row['stop'])
                    stop = str(row['start'])
                f.write(start + "\t" + stop + "\t" + row['Region'] + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['annot']) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ row['locus_tag'] + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if trna_empty == False:
                subset_trna_df = trna_df[trna_df['contig'] == contig]
                for index, row in subset_trna_df.iterrows():
                    start = str(row['start'])
                    stop = str(row['stop'])
                    if row['frame'] == '-':
                        start = str(row['stop'])
                        stop = str(row['start'])
                    f.write(start + "\t" + stop + "\t" + row['Region'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['trna_product']) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ row['locus_tag'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if crispr_count > 0:
                subset_crispr_df = crispr_df[crispr_df['contig'] == contig]
                for index, row in subset_crispr_df.iterrows():
                    start = str(row['start'])
                    stop = str(row['stop'])
                    if row['frame'] == '-':
                        start = str(row['stop'])
                        stop = str(row['start'])
                    f.write(start + "\t" + stop + "\t" + row['Region'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ row['locus_tag'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['rpt_unit_seq']) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if tmrna_flag == True:
                subset_tmrna_df = tmrna_df[tmrna_df['contig'] == contig]
                for index, row in subset_tmrna_df.iterrows():
                    start = str(row['start'])
                    stop = str(row['stop'])
                    if row['frame'] == '-':
                        start = str(row['stop'])
                        stop = str(row['start'])
                    f.write(start + "\t" + stop + "\t" + 'tmRNA' + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"locus_tag" + "\t"+ row['locus_tag'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ 'transfer-messenger RNA, SsrA' + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")


def remove_post_processing_files(out_dir, gene_predictor, meta):
    """
    Cleans temporary files up
    :param out_dir: output directory path
    :param gene_predictor: phanotate or prodigal
    :return: 
    """
    sp.run(["rm", "-rf", os.path.join(out_dir, "target_dir") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "tmp_dir/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "vfdb_tmp_dir/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "vfdb/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "vfdb_results.tsv") ])  
    sp.run(["rm", "-rf", os.path.join(out_dir, "CARD_tmp_dir/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "CARD/") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "CARD_results.tsv") ])  
    sp.run(["rm", "-rf", os.path.join(out_dir, "cleaned_" + gene_predictor + ".tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "input_fasta_delim.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs_results.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_mmseqs.tsv") ])
    # leave in tophits
    sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_aas_tmp.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_out_tmp.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "pharokka_tmp.gff") ])  
    sp.run(["rm", "-rf", os.path.join(out_dir,"mash_out.tsv") ])
    sp.run(["rm", "-rf", os.path.join(out_dir, "input_mash_sketch.msh") ])

    
    if gene_predictor == "phanotate":
        sp.run(["rm", "-rf", os.path.join(out_dir, "phanotate_out.txt") ])
    if gene_predictor == "prodigal":
        sp.run(["rm", "-rf", os.path.join(out_dir, "prodigal_out.gff") ])
    # delete the tmp meta files
    if meta == True:
        sp.run(["rm", "-rf", os.path.join(out_dir, "input_split_tmp/") ])


def get_crispr_count(out_dir, prefix):
    """
    Gets number of crisprs
    :param out_dir: output directory path
    :param prefix: prefix
    :return: crispr_count integer
    """
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
    """
    Determines if trna output file is empty
    :param out_dir: output directory path
    :return: trna_empty Boolean
    """
    trna_empty = False
    if os.stat(os.path.join(out_dir, "trnascan_out.gff")).st_size == 0:
        trna_empty = True
    return trna_empty


def parse_aragorn(out_dir,length_df, prefix):
    """
    Parses aragorn output file
    :param out_dir: output directory path
    :param prefix: prefix
    :param length_df: a pandas df as output from get_contig_name_lengths()
    :return: tmrna_flag Boolean whethere there is tmrna or not
    """
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
def process_vfdb_results(out_dir, merged_df):
    """
    Processes VFDB results
    :param out_dir: output directory path
    :param merged_df: merged_df in process_results
    :return: merged_df merged_df updated with VFDB results
    """
    ##vfdb
    vfdb_file = os.path.join(out_dir, "vfdb_results.tsv")
    col_list = ["vfdb_hit", "gene", "vfdb_alnScore", "vfdb_seqIdentity", "vfdb_eVal", "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen"] 
    vfdb_df = pd.read_csv(vfdb_file, delimiter= '\t', index_col=False , names=col_list) 
    genes = vfdb_df.gene.unique()

    # get top hit 
    tophits = []

    for gene in genes:
        tmp_df = vfdb_df.loc[vfdb_df['gene'] == gene].sort_values('vfdb_eVal').reset_index(drop=True).loc[0]
        tophits.append([tmp_df.vfdb_hit, tmp_df.gene, tmp_df.vfdb_alnScore, tmp_df.vfdb_seqIdentity, tmp_df.vfdb_eVal])

    tophits_df = pd.DataFrame(tophits, columns=['vfdb_hit', 'gene', 'vfdb_alnScore', 'vfdb_seqIdentity', 'vfdb_eVal'])

    # left join vfdb to merged_df
    tophits_df['gene']=tophits_df['gene'].astype(str)

    # merge top hits into the merged df
    merged_df = merged_df.merge(tophits_df, on='gene', how='left')
    merged_df["vfdb_hit"] = merged_df["vfdb_hit"].replace(nan, 'None', regex=True)
    merged_df["vfdb_alnScore"] = merged_df["vfdb_alnScore"].replace(nan, 'None', regex=True)
    merged_df["vfdb_seqIdentity"] = merged_df["vfdb_seqIdentity"].replace(nan, 'None', regex=True)
    merged_df["vfdb_eVal"] = merged_df["vfdb_eVal"].replace(nan, 'None', regex=True)


    # if there is a hit extract information about it
    if len(tophits_df['vfdb_hit']) > 0:
        number_vfs = len(tophits_df['vfdb_hit'])
        print( str(number_vfs) + ' VFDB virulence factors identified.')
        merged_df[['genbank','desc_tmp', 'vfdb_species']] = merged_df['vfdb_hit'].str.split('[',expand=True)
        merged_df['vfdb_species'] = merged_df['vfdb_species'].str.replace("]", "")
        merged_df['vfdb_species'] = merged_df['vfdb_species'].str.strip()
        # genbank has the info 
        merged_df['vfdb_short_name'] = merged_df['genbank'].str.split(')', 1).str[1]
        merged_df['vfdb_description'] = merged_df['genbank'].str.split(')', 2).str[2]
        merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].str.replace("(", "")
        merged_df['vfdb_short_name'] = merged_df['vfdb_short_name'].str.split(')', 1).str[0]
        merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].str.strip()
        merged_df["vfdb_description"] = merged_df["vfdb_description"].str.strip()
        # remove and add None
        merged_df = merged_df.drop(columns = ['genbank', 'desc_tmp'])
        merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].replace(nan, 'None', regex=True)
        merged_df["vfdb_description"] = merged_df["vfdb_description"].replace(nan, 'None', regex=True)
        merged_df["vfdb_species"] = merged_df["vfdb_species"].replace(nan, 'None', regex=True)
    else:
        print('0 VFDB virulence factors identified.')
        merged_df["vfdb_short_name"] = 'None'
        merged_df["vfdb_description"] = 'None'
        merged_df["vfdb_species"] = 'None'
    return (merged_df, tophits_df)


    

#### process CARD files
def process_card_results(out_dir, merged_df, db_dir):
    """
    Processes card results
    :param out_dir: output directory path
    :param merged_df: merged_df in process_results
    :return: merged_df merged_df updated with card results
    """
    ##card
    card_file = os.path.join(out_dir, "CARD_results.tsv")
    print("Processing CARD output.")
    col_list = ["CARD_hit", "gene", "CARD_alnScore", "CARD_seqIdentity", "CARD_eVal", "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen"] 
    card_df = pd.read_csv(card_file, delimiter= '\t', index_col=False , names=col_list) 
    genes = card_df.gene.unique()

    # get top hit 
    tophits = []

    for gene in genes:
        tmp_df = card_df.loc[card_df['gene'] == gene].sort_values('CARD_eVal').reset_index(drop=True).loc[0]
        tophits.append([tmp_df.CARD_hit, tmp_df.gene, tmp_df.CARD_alnScore, tmp_df.CARD_seqIdentity, tmp_df.CARD_eVal])

    tophits_df = pd.DataFrame(tophits, columns=['CARD_hit', 'gene', 'CARD_alnScore', 'CARD_seqIdentity', 'CARD_eVal'])
    
    # left join tophits_df to merged_df
    tophits_df['gene']=tophits_df['gene'].astype(str)

    # merge top hits into the merged df
    merged_df = merged_df.merge(tophits_df, on='gene', how='left')
    merged_df["CARD_hit"] = merged_df["CARD_hit"].replace(nan, 'None', regex=True)
    merged_df["CARD_alnScore"] = merged_df["CARD_alnScore"].replace(nan, 'None', regex=True)
    merged_df["CARD_seqIdentity"] = merged_df["CARD_seqIdentity"].replace(nan, 'None', regex=True)
    merged_df["CARD_eVal"] = merged_df["CARD_eVal"].replace(nan, 'None', regex=True)
    # if there is a hit extract info
    if len(tophits_df['CARD_hit']) > 0:
        number_cards = len(tophits_df['CARD_hit'])
        print( str(number_cards) + ' CARD AMR genes identified.')
        merged_df[['genbank', 'CARD_species']] = merged_df['CARD_hit'].str.split('[',expand=True)
        merged_df['CARD_species'] = merged_df['CARD_species'].str.replace("]", "")
        merged_df['CARD_species'] = merged_df['CARD_species'].str.strip()
        merged_df[['gb','genbank', 'ARO_Accession', 'CARD_short_name']] = merged_df['genbank'].str.split('|',expand=True)
        merged_df["CARD_short_name"] = merged_df["CARD_short_name"].str.strip()
    # read in aro_index 
        CARD_index_file = os.path.join(db_dir, "aro_index.tsv")
        col_list = ["ARO_Accession", "CVTERM_ID", "Model_Sequence_ID", "Model_ID", "Model_Name", "ARO_Name", "Protein_Accession", "DNA_Accession", "AMR_Gene_Family", "Drug_Class", "Resistance_Mechanism", "CARD_Short_Name"] 
        card_index_df = pd.read_csv(CARD_index_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
        card_index_df = card_index_df.drop(columns = ['CVTERM_ID', 'Model_Sequence_ID', 'Model_ID', 'Model_Name', 'ARO_Name', 'CARD_Short_Name'])
        merged_df = merged_df.merge(card_index_df, on='ARO_Accession', how='left')
        merged_df = merged_df.drop(columns = ['gb', 'genbank'])
        merged_df["CARD_species"] = merged_df["CARD_species"].replace(nan, 'None', regex=True)
        merged_df["ARO_Accession"] = merged_df["ARO_Accession"].replace(nan, 'None', regex=True)
        merged_df["CARD_short_name"] = merged_df["CARD_short_name"].replace(nan, 'None', regex=True)
        merged_df["Protein_Accession"] = merged_df["Protein_Accession"].replace(nan, 'None', regex=True)
        merged_df["DNA_Accession"] = merged_df["DNA_Accession"].replace(nan, 'None', regex=True)
        merged_df["AMR_Gene_Family"] = merged_df["AMR_Gene_Family"].replace(nan, 'None', regex=True)
        merged_df["Drug_Class"] = merged_df["Drug_Class"].replace(nan, 'None', regex=True)
        merged_df["Resistance_Mechanism"] = merged_df["Resistance_Mechanism"].replace(nan, 'None', regex=True)
    # if no hits just Nones
    else: 
        print('0 CARD AMR genes identified.')
        merged_df["CARD_species"] = 'None'
        merged_df["ARO_Accession"] = 'None'
        merged_df["CARD_short_name"] = 'None'
        merged_df["Protein_Accession"] = 'None'
        merged_df["DNA_Accession"] = 'None'
        merged_df["AMR_Gene_Family"] = 'None'
        merged_df["Drug_Class"] = 'None'
        merged_df["Resistance_Mechanism"] = 'None'

    return (merged_df, tophits_df)

# check if a file is empt
def is_file_empty(file):
    """
    Determines if file is empty
    :param file: file path
    :return: empty Boolean
    """
    empty = False
    if os.stat(file).st_size == 0:
        empty = True
    return empty

def inphared_top_hits(out_dir, db_dir, length_df, prefix):
    """
    Process mash output to get inphared top hits
    :param out_dir: output directory
    :param length_df: a pandas df as output from get_contig_name_lengths()
    """
    
    mash_tsv = os.path.join(out_dir,"mash_out.tsv")
    col_list = ["contig", "Accession", "mash_distance", "mash_pval", "mash_matching_hashes"]


    # get contigs - convert to string to make the match work for integer contigs
    contigs = length_df["contig"].astype("string")

    # instantiate tophits list
    tophits_mash_df = []


    mash_df = pd.read_csv(mash_tsv, delimiter= '\t', index_col=False, names=col_list ) 

    # instantiate tophits list
    tophits = []

    for contig in contigs:
        hit_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True)
        hits = len(hit_df['mash_distance'])
        # add only if there is a hit
        if hits > 0:
            # top hit 
            top_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True).loc[0]
            tophits.append([top_df.contig, top_df.Accession, top_df.mash_distance, top_df.mash_pval, top_df.mash_matching_hashes])
        else:
            tophits.append([contig, "no_inphared_mash_hit", "no_inphared_mash_hit", "no_inphared_mash_hit", "no_inphared_mash_hit"])
        # create tophits df
    tophits_mash_df = pd.DataFrame(tophits, columns=["contig", "Accession", "mash_distance", "mash_pval", "mash_matching_hashes"])

    # read in the plasdb tsv 
    inphared_tsv_file = os.path.join(db_dir, "5Jan2023_data.tsv")
    # with open(plsdb_tsv_file, 'rb') as f:
    #     result = chardet.detect(f.readline())
    #     print(result)
    cols = ["Accession","Description","Classification","Genome_Length_(bp)","Jumbophage","molGC_(%)","Molecule","Modification_Date","Number_CDS","Positive_Strand_(%)","Negative_Strand_(%)","Coding_Capacity_(%)","Low_Coding_Capacity_Warning","tRNAs","Host","Lowest_Taxa","Genus","Sub-family","Family","Order","Class","Phylum","Kingdom","Realm","Baltimore_Group","Genbank_Division","Isolation_Host_(beware_inconsistent_and_nonsense_values)"]
    inphared_tsv_file = pd.read_csv(inphared_tsv_file, delimiter= '\t', index_col=False, names=cols, skiprows=1,low_memory=False) 

    combined_df = tophits_mash_df.merge(inphared_tsv_file, on='Accession', how='left')

    combined_df.to_csv(os.path.join(out_dir, prefix + "_top_hits_mash_inphared.tsv"), sep="\t", index=False)        

 



        



