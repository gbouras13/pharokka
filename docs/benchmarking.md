# Benchmarking v1.1.0 (from the Manuscript)

pharokka (v1.1.0) has been benchmarked on an Intel Xeon CPU E5-4610 v2 @ 2.30 specifying 16 threads. Below is benchamarking comparing Pharokka run with PHANOTATE and Prodigal against Prokka v1.14.6 run with PHROGs HMM profiles, as modified by Andrew Millard (https://millardlab.org/2021/11/21/phage-annotation-with-phrogs/).

Benchmarking was conducted on Enterbacteria Phage Lambda (Genbank accession J02459) Staphylococcus Phage SAOMS1 (Genbank Accession MW460250) and 673 crAss-like phage genomes in one multiFASTA input taken from Yutin, N., Benler, S., Shmakov, S.A. et al. Analysis of metagenome-assembled viral genomes from the human gut reveals diverse putative CrAss-like phages with unique genomic features. Nat Commun 12, 1044 (2021) https://doi.org/10.1038/s41467-021-21350-w.

For the crAss-like phage genomes, Pharokka meta mode `-m` was enabled.

| Phage Lambda            | pharokka PHANOTATE | pharokka Prodigal | Prokka with PHROGs |
|------------------------|--------------------|-------------------|--------------------|
| Time (min)             | 4.19               | 3.88              | 0.27               |
| CDS                    | 88                 | 61                | 62                 |
| Coding Density (%)     | 94.55              | 83.69             | 84.96              |
| Annotated Function CDS | 43                 | 37                | 45                 |  
| Unknown Function CDS   | 45                 | 24                | 17                 |

| Phage SAOMS1           | pharokka PHANOTATE | pharokka Prodigal | Prokka with PHROGs |   
|------------------------|--------------------|-------------------|--------------------|
| Time (min)             | 4.26               | 3.89              | 0.93               |
| CDS                    | 246                | 212               | 212                |
| Coding Density (%)     | 92.27              | 89.69             | 89.31              |  
| Annotated Function CDS | 92                 | 93                | 92                 |
| Unknown Function CDS   | 154                | 119               | 120                |  

| 673 crAss-like genomes from Yutin et al., 2021 | pharokka PHANOTATE Meta Mode | pharokka Prodigal Meta Mode  | Prokka with PHROGs |
|------------------------------------------------|------------------------------|------------------------------|--------------------|
| Time (min)                                     | 106.55                       | 11.88                        | 252.33             |
| Time Gene Prediction (min)                     | 96.21                        | 3.4                          | 5.12               |
| Time tRNA Prediction (min)                     | 1.25                         | 1.08                         | 0.3                |
| Time Database Searches (min)                   | 6.75                         | 5.58                         | 238.77             |
| CDS                                            | 138628                       | 90497                        | 89802              |
| Contig Min Coding Density (%)                  | 66.01                        | 46.18                        | 46.13              |
| Contig Max Coding Density (%)                  | 98.86                        | 97.85                        | 97.07              |
| Annotated Function CDS                         | 9341                         | 9228                         | 14461              |
| Unknown Function CDS                           | 129287                       | 81269                        | 75341              |

`pharokka` scales well for large metavirome datasets due to the speed of MMseqs2. In fact, as the size of the input file increases, the extra time taken is required for running gene prediction (particularly PHANOTATE) and tRNA-scan SE2 - the time taken to conduct MMseqs2 searches remain small due to its many vs many approach.

If you require  fast annotations of extremely large datasets (i.e. thousands of input contigs), running `pharokka` with Prodigal is recommended.

# Benchmarking v1.4.0 

`pharokka` v1.4.0 has also been run on phage SAOMS1 and also the same 673 crAss phage dataset to showcase:

1. The improved sensitivity of gene annotation with PyHMMER and a demonstration of how `--fast` is slower for metagenomes. 
    * If you can deal with the compute cost (especially for large metagenomes), I highly recommend `--fast` or  `--meta_hmm` for metagenomes given how much more sensitive HMM search is.
2. The large speed-up over v1.3.2 with `--fast` for phage isolates - with the proviso that no virulence factors or AMR genes will be detected. 
3. The slight speed-up over v1.3.2 with `--mmseqs2_only`.

All benchmarking was conducted on a Intel® Core™ i7-10700K CPU @ 3.80GHz on a machine running Ubuntu 20.04.6 LTS with 16 threads (`-t 16`). 

SAOMS1 was run with Phanotate

| Phage SAOMS1           | pharokka v1.4.0 | pharokka v1.4.0 `--fast` | pharokka v1.3.2 |   
|------------------------|-----------------|--------------------------|-----------------|
| Time (min)             | 3.73            | 0.70                     | 5.08            | 
| CDS                    | 246             | 246                      | 246             | 
| Annotated Function CDS | 93              | 93                       | 92              | 
| Unknown Function CDS   | 153             | 153                      | 154             |  

The 673 crAss-like genomes were run with `-m` (defaults to `--mmseqs2_only` in v 1.4.0) and with `-g prodigal` (i.e. pyrodigal v2.3.0).

| 673 crAss-like genomes | pharokka v1.4.0 `--fast`  | pharokka v1.4.0 `--mmseqs2_only` | pharokka v1.3.2 |
|------------------------|---------------------------|----------------------------------|-----------------|
| Time (min)             | 35.62                     | 11.05                            | 13.27           |
| CDS                    | 91999                     | 91999                            | 91999           |
| Annotated Function CDS | **16713**                 | 9150                             | 9150            |
| Unknown Function CDS   | 75286                     | 82849                            | 82849           |


# Benchmarking v1.5.0

`pharokka v1.5.0` was run on the 673 crAss phage dataset to showcase the improved CDS prediction of `-g prodigal-gv` for metagenomic datasets where some phages likely have alternative genetic codes. 

All benchmarking was conducted on a Intel® Core™ i7-10700K CPU @ 3.80GHz on a machine running Ubuntu 20.04.6 LTS with 8 threads (`-t 8`). `pyrodigal-gv v0.1.0` and `pyrodigal v3.0.0` were used respectively with `--fast`.

| 673 crAss-like genomes | `pharokka` v1.5.0 `-g prodigal-gv`  | `pharokka` v1.5.0 `-g prodigal` | 
|------------------------|------------------------------------|----------------------------------|
| Total CDS              | 81730                              | 91999                            | 
| Annotated Function CDS | **20344**                          | 17458                            | 
| Unknown Function CDS   | 61386                              | 74541                            |
| Contigs with genetic code 15 | 229                          | NA                               | 
| Contigs with genetic code 4 | 38                            | NA                               | 
| Contigs with genetic code 11 | 406                          | 673                              | 

Fewer larger CDS were predicted more accurately, leading to an increase in the number of coding sequences with annotated functions. Approximately 40% of contigs in this dataset were predicted to use non-standard genetic codes according to `pyrodigal-gv`.