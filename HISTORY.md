History
=======

1.7.4 (2024-11-23)
------------------

* Adds `--trna_scan_model`  parameter with two accepted options: `--trna_scan_model general` (this will be run by default - what Pharokka has always been running) and `--trna_scan_model bacterial`. See the [tRNAscan-SE paper](https://doi.org/10.1093/nar/gkab688) for more information.
* Bumps the `dnaapler` dependency to v1.0.1 due to a breaking dependency change in `dnaapler`.


1.7.3 (2024-07-10)
------------------

* Fixes issue with genbank creation if certain CARD hits are found [issue #339](https://github.com/gbouras13/pharokka/issues/339)).
    * Due to some semi-colons in the CARD metadata, new qualifier keys were being made in error.
    * Solved by removing semicolons in the updated CARD metadata sheet

1.7.2 (2024-05-27)
------------------

* Identifies issue if `#` is input contig header - Pharokka will error if your contig headers contain this character and prompt you to remove them

1.7.1 (2024-03-13)
------------------

* Adds Google Colab notebook that can run pharokka and [phold](https://github.com/gbouras13/phold). 
* The notebook is  [https://colab.research.google.com/github/gbouras13/pharokka/blob/master/run_pharokka_and_phold.ipynb](https://colab.research.google.com/github/gbouras13/pharokka/blob/master/run_pharokka_and_phold.ipynb)
* Fixes #334 issues with contig ids if they were in scientific notation or lead with 0s.
* Fixes issues with `pharokka_proteins.py` not outputting PHROG annotations.


1.7.0 (2024-03-04)
------------------

* Adds `pharokka_multiplotter.py` to plot multiple phage contigs at once
* Adds separate contig FASTA files if `-s -m` is run (in `single_fastas`)

1.6.1 (2024-01-17)
------------------

* Fixes a bug that was removing tRNAs from the `.tbl` output format #323.

1.6.0 (2024-01-11)
------------------

* Fixes a variety of bugs (#300 `pharokka_proteins.py` crashing if it found VFDB hits, #303 errors in the `.tbl` format, #316 errors with types and where custom HMM dbs had identical scored hits, #317 types and #320 deprecated GC function)
* Adds `--mash_distance` and `--minced_args` as parameters (#299 thanks @iferres).


1.5.1 (2023-10-26)
------------------

* Fixes `dnaapler` version to `>=0.4.0` with new changes to dnaapler
* Adds `.svg` format output with `pharokka_plotter.py`

1.5.0 (2023-09-20)
------------------

* Adds support for `pyrodigal-gv` implementing `prodigal-gv` as a gene predictor ([pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) and [prodigal-gv](https://github.com/apcamargo/prodigal-gv)). This can be specified with `-g prodigal-gv`.
* Thanks to @[althonos](https://github.com/althonos) and @[apcamargo](https://github.com/apcamargo) for making this possible, and to @[asierFernandezP](https://github.com/asierFernandezP) for raising this as an issue in the first place [here](https://github.com/gbouras13/pharokka/issues/290) in #290.
* Adds checks to determine if your input FASTA has duplicated contig headers from #293 [here](https://github.com/gbouras13/pharokka/issues/293). Thanks @[thauptfeld](https://github.com/thauptfeld) for raising this.
* `-g prodigal` and `-g prodigal-gv` should be much faster thanks to multithread support added by the inimitable @althonos.
* Genbank format output will be designated with PHG not VRL (following this issue https://github.com/RyanCook94/inphared/issues/22).
* The `_length_gc_cds_density.tsv` and `_cds_final_merged_output.tsv` files now contain the translation table/genetic code for each contig (usually 11 but now not always if you use `pyrodigal-gv`). 
* `--skip_mash` flag added to skip finding the closest match for each contig in INPHARED using mash.
* `--skip_extra_annotations` flag added to skip running tRNA-scanSE, MinCED and Aragorn in case you only want CDS predictions and functional annotations.

1.4.1 (2023-09-04)
------------------

Pharokka v1.4.1 is a small patch fix release fixing #286, where if you specified `--dnaapler` and `-m`, pharokka would not find the correct output file from `dnaapler` and would crash. 


1.4.0 (2023-08-27)
------------------

Pharokka v1.4.0 is a large update implementing:
 
* More sensitive search for PHROGs using Hidden Markov Models (HMMs) using the amazing [PyHMMER](https://github.com/althonos/pyhmmer). Thanks to @althonos for some of the best written and well documented software I have ever used.
* By default, `pharokka` will now run searches using both MMseqs2 (PHROGs, CARD and VFDB) and HMMs (PHROGs). MMseqs2 was kept for PHROGs as it provides more information than the HMM results (e.g. sequence alignment identities & top hit PHROG protein) if it finds a hit.
* `--fast` or `--hmm_only` parameter, which only runs PyHMMER on PHROGs. It will not run MMseqs2 at all on PHROGs, CARD or VFDB. For phage isolates, this will be much faster than v1.3.2, but you will not get CARD or VFDB annotations. For metagenomes, this will be (much) slower though!
* Updated databases as of 23 August 2023. You will need to download the new `pharokka v1.4.0` databases because these now contain PHROG HMM profiles. The VFDB database is now clustered at 50% sequence identity (which speeds up runtime).
* Other changes in the codebase should make `pharokka v1.4.0` run somewhat faster than v1.3.2, even if PyHMMER is not used i.e. `--mmseqs2_only` is specified.
* The print screen and log files are neater and more information rich using loguru. There is also a new `logs` directory containing separate log files for each tool in the pipeline. This is thanks to taking and modifying some code from @mbhall88 [tbpore](https://github.com/mbhall88/tbpore).
* `install_databases.py` has been modified to be more robust and somewhat faster. This is thanks to taking ideas and modifying some code from @oschwengers [bakta](https://github.com/oschwengers/bakta).
* `--mmseqs2_only` which will essentially run `pharokka` as it was v1.3.2. It is default in meta mode `-m` or `--meta`.
* `pharokka_proteins.py`, which takes an input file of amino acid proteins in FASTA format and runs MMseqs2 (PHROGs, CARD, VFDB) and PyHMMER (PHROGs). See the [proteins documentation](docs/proteins.md) for more details. Thanks to Brady Cress for the idea.
* `--custom_hmm` parameter, which allows for custom HMM profile databases to be used with `pharokka`. Thanks to @pck00 for the idea.
* `create_custom_hmm.py` which facilitates  the creation of a HMM profile database from multiple sequence alignments.  See the [documentation](docs/custom.md) for more details about how to create a compatible HMM profile database.
* `--dnaapler` flag, which automatically detects and reorients your phage to start with the large terminase subunit. For more information, see [dnaapler](https://github.com/gbouras13/dnaapler).
* `--genbank` flag, which allows for genbank format input with `-i`. This will take all (customised) CDS calls in genbank file and PHANOTATE/pyrodigal will not be run. So if you have done manual custom gene curation and want to functionally annotate your customised CDS, this option is recommended. Thanks to @pck00 for the idea.
* Fixes to `-c`, which should now work properly with `-g prodigal` (thanks @alegione for the fixes).

1.3.2 (2023-04-26)
------------------

* Fixes bug with pharokka_plotter.py, which would crash if the phage had tmRMAs or CRISPRs.
* Fixes bug where integration & excision fwd strand CDS would not be plotted in the correct colour
* Adds tmRNAs and CRISPRs to pharokka_plotter.py.

1.3.1 (2023-04-20)
------------------

* Adds tRNAs to pharokka_plotter.py.
* Adds the -s split mode option with metagenome mode, this will output separate single fastas, gff and genbank files along with -m. It is ideally used for situations where you have bulk phage isolates you want to annotate in one go. 

1.3.0 (2023-04-11)
------------------

* Adds pharokka_plotter.py to create plots with pyCirclize.
* Fixes issue with VFDB and CARD counts in _cds_functions.tsv being 0 even is a virulence factor or AMR gene is detected.
* Adds better error checking for --threads.

1.2.1 (2023-02-20)
------------------

* Minor update to fix Biopython version <=1.80, due to a breaking change with 1.81 and bcbio-gff this [issue](https://github.com/chapmanb/bcbb/issues/136).

1.2.0 (2023-01-17)
------------------

* Adds the functionality of mapping each contig against the INPHARED database using mash (https://github.com/RyanCook94/inphared). The top hit for each contig (under a maximum mash distance threshold of 0.2) is kept.
* New database adding INPHARED.
* Replaced prodigal with pyrodigal as it is being actively maintained and used by bakta.
* Adds --citation.
* Adds checks for dependencies.
* Adds --terminase terminase mode to re-orient a single contig phage to begin with a certain orientation and coordinate (most commonly, the large terminase subunit). With this, you must also specify --terminase_strand the strand of the terL gene and --terminase_start the start coordinate.
* All locus tags end with 4 digits (trailing zeros) in order to play nice with vConTACT2 and start with 1 not 0.
* In meta mode, the locus tags now begin with the contig header, not a random string (or chosen prefix).
* Cleans up the .tbl so it should automatically be accepted by NCBI Bankit.


1.1.0 (2022-10-20)
------------------

* Renames the CDS output files to *.faa for amino acids and *.ffn for nulceotide sequences
* Implementation of consistent CDS name (equal to the locus_tag) across all output files
* Creates terL.faa and terL.ffn, which contain the sequences of any identified terminase large subunit CDSs
* Passes multithreading to PHANOTATE and tRNAscan-SE in meta mode indicated by flag -m, which provides approximately a t-fold improvement in run-time for large metavirome datasets, where t is the number of threads. 

1.0.1 (2022-10-10)
------------------

* Minor release to fix a string-parsing bug where pharokka v1.0.0 would crash when certain VFDB virulence factors were detected.

1.0.0 (2022-09-16)
------------------

* Removes errors (with post_processing functions not being parsed as strings) to improve robustness.
* Codebase more reliable and consistent
* Overhaul of install_databases.py
* Adds pre-existing Pharokka Database available at https://zenodo.org/record/7081772

0.1.11 (2022-09-13)
------------------

*  Adds CARD and VFDB databases.

0.1.10 (2022-08-31)
------------------

* Fixes issues with Genbank output files.


0.1.9 (2022-08-04)
------------------

* Fixes bug with parsing Aragorn output for some phages.
* Fixes bug with phages with no mmseqs2 PHROGs hits (for some very small phages).

0.1.8 (2022-07-24)
------------------

* Minor release.
* Fixes bug with install_databases.py.
* Adds locustag to .tbl output.

0.1.7 (2022-07-24)
------------------

* Removes hh-suite dependency to reduce run-time (redundant with e-value option).
* Adds e-value option for passing into mmseqs2.
* Adds CRISPR detection using MinCED.
* Adds tmRNA detection using Aragorn.
* Adds CDS coding density to _length_gc_cds_density.tsv.

0.1.6 (2022-07-16)
------------------

* Version update for bioconda release

0.1.5 (2022-07-09)
------------------

* Adds Prodigal option for gene prediction -g prodigal
* Handles metagenome assembled virus contigs (that may have no genes) with -m flag
* Fixes tRNA count bug in cds_functions.tsv

0.1.4 (2022-07-08)
------------------

* Updates PHROGs annotation to version 4.

0.1.3 (2022-07-06)
------------------

* Adds .log file
* Adds -l option

0.1.2 (2022-07-04)
------------------

* Adds .genbank file

0.1.1 (2022-06-28)
------------------

* Second release

0.1.0 (2021-12-09)
------------------

* First release
