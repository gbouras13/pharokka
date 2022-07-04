pharokka creates a number of output files in different formats.

Main Output
----------
The main output is a gff3 file that is suitable for use downstream pangenomic pipelines such as Roary ([https://sanger-pathogens.github.io/Roary/](https://sanger-pathogens.github.io/Roary/)) to generate pangenomes.

* The 'phrog=' section shows the closest matching PHROG. The 'top_hit=' section shows the closest matching protein in the PHROGs database.

Other Files
------
* A .genbank file, which is converted from the gff using seqret [http://emboss.open-bio.org](http://emboss.open-bio.org).

* A .tbl file, which is a flat-file table suitable to be uploaded to the NCBI's Bankit.

* A `cds_functions.tsv` file, which includes counts of CDSs, tRNAs, and functions assigned to CDSs according to the PHROGs database.

* A `_length_gc.tsv` file, which outputs the phage's length and GC percentage.

* A `_final_merged_output.tsv`, which gives the full output from mmseqs2 and hh-suite. In particular, the 'match_type' column shows whether mmseqs2 or hh-suite was used to conduct the match. In general, mmseqs2 identifies most CDSs, while small (80-200bp) hypothetical proteins are identified using hh-suite matching the closest PHROG. pharokka should be used as a rough guide only in these cases.

* Further, the 'score' column contains the PHANOTATE score for each CDS. In general, the closer the score to 0, the smaller the CDS and the more likely that a PHROG will not be identified by mmseqs2.

* For more information about PHROGs please consult the website [https://phrogs.lmge.uca.fr](https://phrogs.lmge.uca.fr) and paper [https://doi.org/10.1093/nargab/lqab067](https://doi.org/10.1093/nargab/lqab067).
