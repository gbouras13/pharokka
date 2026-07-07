# Output Files

`pharokka` creates a number of output files in different formats.

## Main Output

The main output is a `.gff` GFF3 file suitable for use in downstream pangenomic pipelines such as [Roary](https://sanger-pathogens.github.io/Roary/) or [panaroo](https://doi.org/10.1186/s13059-020-02090-4).

* The `phrog=` attribute shows the closest matching PHROG, or `No_PHROG` if there are no matching PHROGs below the E-value threshold. The `top_hit=` attribute shows the closest matching protein in the PHROGs database.

* The 6th column of the `.gff` file is the PHANOTATE w_orf score (the more negative, the more likely the CDS is a gene) or the prodigal coding score (the larger, the more likely the CDS is a gene). Please see the PHANOTATE and Prodigal papers for more details.

## Other Files

* A `.gbk` genbank formatted file, converted from the gff.

* A `.log` file, which holds all logging output. It is time-stamped in the `%m%d%Y_%H%M%S` format.

* A `.tbl` file, which is a flat-file table suitable for upload to NCBI's BankIt.

* A `_cds_functions.tsv` file, which includes counts of CDSs, tRNAs, CRISPRs and tmRNAs, and functions assigned to CDSs according to the PHROGs database.

* A `_length_gc_cds_density.tsv` file, which outputs the phage's length, GC percentage, translation table and CDS coding density.

* `phanotate.ffn` or `prodigal.ffn` which will hold all nucleotide sequences of predicted CDSs.

* `phanotate.faa` or `prodigal.faa` which will hold all amino acid sequences of predicted CDSs.

* `_aragorn.txt` and `_aragorn.gff` files, which hold the raw and parsed output from Aragorn, respectively.

* `_minced_spacers.txt` and `_minced.gff`, which hold the output from MinCED.

* A `_trnascan.gff` which holds the output from tRNAscan-SE 2.

* A `trnascan_out.sec` file with tRNAscan-SE 2 secondary structures.

* A `_cds_final_merged_output.tsv`, which gives the parsed output from MMseqs2 and PyHMMER. The key columns are:

  | Column | Description |
  |--------|-------------|
  | `gene` | CDS locus tag |
  | `start`, `stop`, `strand` | Genomic coordinates |
  | `score` | PHANOTATE w_orf score or Prodigal coding score |
  | `mmseqs_phrog` | Best PHROGs hit from MMseqs2 |
  | `mmseqs_seqIdentity`, `mmseqs_eVal` | MMseqs2 hit identity and E-value |
  | `pyhmmer_phrog` | Best PHROGs hit from PyHMMER |
  | `pyhmmer_bitscore`, `pyhmmer_evalue` | PyHMMER hit scores |
  | `phrog` | Final PHROG assignment (MMseqs2 preferred over PyHMMER when both find a hit) |
  | `Method` | Which method assigned the PHROG (`mmseqs2` or `pyhmmer`) |
  | `annot` | Functional annotation string |
  | `category` | PHROGs functional category |
  | `vfdb_hit`, `vfdb_seqIdentity`, `vfdb_eVal` | Top VFDB hit and scores |
  | `CARD_hit`, `CARD_seqIdentity`, `CARD_eVal` | Top CARD hit and scores |
  | `custom_hmm_id`, `custom_hmm_bitscore`, `custom_hmm_evalue` | Custom HMM hit (if `--custom_hmm` was used) |
  | `transl_table` | Translation table used for this CDS |

  In general, using the default E-value threshold of 1E-05:
  * MMseqs2 or PyHMMER should identify a PHROG for most CDSs if the genome is a phage. Small (80–200 bp) hypothetical proteins often have no matching PHROG.
  * The MMseqs2 phrog is preferred over the PyHMMER phrog in the rare case of disagreement between the two methods (as of v1.4.0).
  * The `score` column contains the PHANOTATE score for each CDS. In general, the closer the score to 0, the smaller the CDS and the more likely that a PHROG will not be identified by MMseqs2.

* A `top_hits_card.tsv` file, which contains any CARD database hits.

* A `top_hits_vfdb.tsv` file, which contains any VFDB database hits.

* A `terL.ffn` file, which contains the nucleotide sequences of all identified large terminase subunit (terL) CDSs.

* A `terL.faa` file, which contains the amino acid sequences of all identified large terminase subunit (terL) CDSs.

* A `_top_hits_mash_inphared.tsv` file which (from v1.2.0) holds the top hits of the INPHARED search.

* Optionally, if you reorient the input contig using `--terminase`, a `_genome_terminase_reoriented.fasta` with the reoriented genome FASTA.

* Optionally, if you specify `-s` or split mode, folders called `single_gbks`, `single_gffs` and `single_fastas` will be created and contain genbank, gff and FASTA files for each respective input contig, named by the contig header.

* Optionally, if you specify `--dnaapler`, a `_dnaapler_reoriented.fasta` file containing the reoriented FASTA format genome, and a `dnaapler` directory containing the output from Dnaapler.

* For more information about PHROGs please consult the website [https://phrogs.lmge.uca.fr](https://phrogs.lmge.uca.fr) and paper [https://doi.org/10.1093/nargab/lqab067](https://doi.org/10.1093/nargab/lqab067).
