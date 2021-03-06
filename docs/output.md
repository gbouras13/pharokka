pharokka creates a number of output files in different formats.

Main Output
----------
The main output is a gff3 file that is suitable for use downstream pangenomic pipelines such as Roary ([https://sanger-pathogens.github.io/Roary/](https://sanger-pathogens.github.io/Roary/)) to generate pangenomes.

* The 'phrog=' section shows the closest matching PHROG, or "No_PHROG" if there are no matching PHROGs below the E-value threshold. The 'top_hit=' section shows the closest matching protein in the PHROGs database.

Other Files
------
* A .genbank file, which is converted from the gff using seqret [http://emboss.open-bio.org](http://emboss.open-bio.org).

* A .log file, which holds the output from tRNA-scanSE, mmseqs2 and hh-suite. It is time-stamped in the "%m%d%Y_%H%M%S" format.

* A .tbl file, which is a flat-file table suitable to be uploaded to the NCBI's Bankit.

* A `cds_functions.tsv` file, which includes counts of CDSs, tRNAs, CRISPRs and tmRNAs and functions assigned to CDSs according to the PHROGs database.

* A `_length_gc_cds_density.tsv` file, which outputs the phage's length, GC percentage and CDS coding density.

* A `_nts.fasta` which will hold all nucleotide sequences of predicted CDSs.

* A `_aas.fasta` which will hold all Amino Acid sequences of predicted CDSs.

* `_aragorn.txt` and `_aragorn.gff` files, which hold the raw and parsed output from Aragorn, respectively.

* `_minced_spacers.txt` and `_minced.gff`, which hold the output from MinCED.

* A `_trnascan.gff` which holds the output from tRNAscan-SE 2.

* A `_final_merged_output.tsv`, which gives the parsed output from mmseqs2. In general, using the default E-value threshold, mmseqs2 should identify a PHROG for most CDSs, while small (80-200bp) hypothetical proteins often will have no matching PHROG. This may also be the case for phages from uncommon sources where few phage have been isolates (such as environment samples).  pharokka should be used as a rough guide only in these cases. It is also recommended that pharokka be re-run with a less restrictive e-value threshold in these cases (e.g. -e 0.1).

* Further, the 'score' column contains the PHANOTATE score for each CDS. In general, the closer the score to 0, the smaller the CDS and the more likely that a PHROG will not be identified by mmseqs2.

* For more information about PHROGs please consult the website [https://phrogs.lmge.uca.fr](https://phrogs.lmge.uca.fr) and paper [https://doi.org/10.1093/nargab/lqab067](https://doi.org/10.1093/nargab/lqab067).
