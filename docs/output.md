`pharokka` creates a number of output files in different formats.

Main Output
----------
The main output is a gff3 file that is suitable for use downstream pangenomic pipelines such as [Roary](https://sanger-pathogens.github.io/Roary/)) to generate pangenomes.

* The 'phrog=' section shows the closest matching PHROG, or "No_PHROG" if there are no matching PHROGs below the E-value threshold. The 'top_hit=' section shows the closest matching protein in the PHROGs database.

* The 6th column of the gff is the PHANOTATE w_orf score (the more negative, the more likely the CDS is a gene) or the prodigal coding score (the larger, the more likely the CDS is a gene). Please see the PHANOTATE and Prodigal papers for more details.

Other Files
------
* A .gbk genbank formatted file, which is converted from the gff. 

* A .log file, which holds the output from tRNA-scanSE, mmseqs2 and post-processing. It is time-stamped in the "%m%d%Y_%H%M%S" format.

* A .tbl file, which is a flat-file table suitable to be uploaded to the NCBI's Bankit.

* A `_cds_functions.tsv` file, which includes counts of CDSs, tRNAs, CRISPRs and tmRNAs and functions assigned to CDSs according to the PHROGs database.

* A `_length_gc_cds_density.tsv` file, which outputs the phage's length, GC percentage and CDS coding density.

* `phanotate.ffn` or `prodigal.ffn` which will hold all nucleotide sequences of predicted CDSs.

* `phanotate.faa` or `prodigal.faa` which will hold all amino acid sequences of predicted CDSs.

* `_aragorn.txt` and `_aragorn.gff` files, which hold the raw and parsed output from Aragorn, respectively.

* `_minced_spacers.txt` and `_minced.gff`, which hold the output from MinCED.

* A `_trnascan.gff` which holds the output from tRNAscan-SE 2.

* A `_cds_final_merged_output.tsv`, which gives the parsed output from mmseqs2. In general, using the default E-value threshold of 1E-05, mmseqs2 should identify a PHROG for most CDSs, while small (80-200bp) hypothetical proteins often will have no matching PHROG. This may also be the case for phages from uncommon sources where few phage have been isolates (such as environment samples).  Pharokka should be used as a rough guide only in these cases. It is also recommended that pharokka be re-run with a less restrictive e-value threshold in these cases (e.g. -e 0.1).

* A `top_hits_card.tsv` file, which contains any CARD database hits.

* A `top_hits_vfdb.tsv` file, which contains any VFDB database hits.

* A `terL.ffn` file, which contains the nulceotide sequences of all identified large terminase subunit (terL) CDSs.

* A `terL.faa` file, which contains the amino acid sequences of all identified large terminase subunit (terL) CDSs.

* Further, the 'score' column contains the PHANOTATE score for each CDS. In general, the closer the score to 0, the smaller the CDS and the more likely that a PHROG will not be identified by mmseqs2.

* A `_top_hits_mash_inphared.tsv` file which from v1.2.0 holds the top hits of the INPHARED search.

* Optionally, if you reorient the input contig, a `_genome_terminase_reoriented.fasta` with the reoriented genome FASTA.

* Optionally, if you specify `-s` or split mode, folders called `single_gbks`, `single_gffs` and `single_fastas` will be created and contain genbank, gff and FASTA files for each respective input contig, named by the contig header. 

* For more information about PHROGs please consult the website [https://phrogs.lmge.uca.fr](https://phrogs.lmge.uca.fr) and paper [https://doi.org/10.1093/nargab/lqab067](https://doi.org/10.1093/nargab/lqab067).
