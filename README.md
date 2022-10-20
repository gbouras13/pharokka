
pharokka
===============

>pharokka is a fast standardized annotation pipeline for bacteriophage genomes.

<p align="center">
  <img src="img/pharokka_workflow.png" alt="pharokka Workflow" height=600>
</p>

pharokka uses [PHANOTATE](https://github.com/deprekate/PHANOTATE) as a default program for gene prediction and later assigns functional annotation by aligning prediction to the [PHROGs](https://phrogs.lmge.uca.fr), [CARD](https://card.mcmaster.ca) and [VFDB](http://www.mgc.ac.cn/VFs/main.htm) databases thorough [mmseqs2](https://github.com/soedinglab/MMseqs2). pharokka's main output is a GFF file suitable for using in downstream pangenomic pipelines like [Roary](https://sanger-pathogens.github.io/Roary/). Moreover, pharokka generates a `cds_functions.tsv` file, which includes counts of CDSs, tRNAs, tmRNAs, CRISPRs and functions assigned to CDSs according to the PHROGs database. See the full [usage](#usage). Check the full documentation [here](https://pharokka.readthedocs.io). If looking for rapid standardized annotation of bacteria, please use [prokka](https://github.com/tseemann/prokka) or [bakta](https://github.com/oschwengers/bakta).

# Installation

```
conda install -c bioconda pharokka
```
For different installation details see the [install](docs/install.md) section.

# Usage

First the install all databases with:

```
install_databases.py -d
```
Or if you would like to specify a different database directory (recommended), that can be achieved as follows: `install_databases.py -o <path/to/databse_dir>`

Then simply run

```
pharokka.py -i <fasta file> -o <output folder> -t <threads>
```

For more advanced commands see [usage](docs/run.md).

# Database advanced installation

If normal installation does not work, you an alternatively download the databases from [Zenodo]( https://zenodo.org/record/7081772/files/pharokka_database_v1.0.0.tar.gz) and untar the directory in a location of your choice.

If you prefer to use the command line:

```
wget "https://zenodo.org/record/7081772/files/pharokka_database_v1.0.0.tar.gz"
```
and
```
tar -xzf pharokka_database_v1.0.0.tar.gz
```

which will create a directory called "pharokka_database_v1.0.0" containing all databases.

# Version Log

A brief description of what is new in each update of pharokka can be found in the [HISTORY.md](HISTORY.md) file.

# System

pharokka has been tested on Linux and MacOS (M1 and Intel).

# Time

On a standard 16GB RAM laptop specifying 8 threads, pharokka should take between 3-10 minutes to run for a single phage, depending on the genome size.

# Benchmarking

pharokka (v1.1.0) has been benchmarked on an Intel Xeon CPU E5-4610 v2 @ 2.30 specifying 16 threads. Below is benchamarking comparing pharokka run with PHANOTATE and Prodigal against Prokka v1.14.6 run with PHROGs HMM profiles, as modified by Andrew Millard (https://millardlab.org/2021/11/21/phage-annotation-with-phrogs/).

Benchmarking was conducted on Enterbacteria Phage Lambda (Genbank accession J02459) Staphylococcus Phage SAOMS1 (Genbank Accession MW460250) and 673 crAss-like phage genomes in one multiFASTA input taken from Yutin, N., Benler, S., Shmakov, S.A. et al. Analysis of metagenome-assembled viral genomes from the human gut reveals diverse putative CrAss-like phages with unique genomic features. Nat Commun 12, 1044 (2021) https://doi.org/10.1038/s41467-021-21350-w.

For the crAss-like phage genomes, pharokka meta mode `-m` was enabled.

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

pharokka scales well for large metavirome datasets due to the speed of mmseqs2. In fact, as the size of the input file increases, the extra time taken is required for running gene prediction (particularly PHANOTATE) and tRNA-scan SE2 - the time taken to conduct mmseqs2 searches remain small due to its many vs many approach. 

If you require  fast annotations of extremely large datasets (i.e. thousands of input contigs), running pharokka with Prodigal is recommended.
 

# Bugs and Suggestions

If you come across bugs with pharokka, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au

# Citation

If you use pharokka, please also cite:

* McNair K., Zhou C., Dinsdale E.A., Souza B., Edwards R.A. (2019) "PHANOTATE: a novel approach to gene identification in phage genomes", Bioinformatics, https://doi.org/10.1093/bioinformatics/btz26.
* Chan, P.P., Lin, B.Y., Mak, A.J. and Lowe, T.M. (2021) "tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes", Nucleic Acids Res., https://doi.org/10.1093/nar/gkab688.
* Steinegger M. and Soeding J. (2017), "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets", Nature Biotechnology https://doi.org/10.1038/nbt.3988.
* Terzian P., Olo Ndela E., Galiez C., Lossouarn J., Pérez Bucio R.E., Mom R., Toussaint A., Petit M.A., Enault F., "PHROG : families of prokaryotic virus proteins clustered using remote homology", NAR Genomics and Bioinformatics, (2021), https://doi.org/10.1093/nargab/lqab067.
* Bland C., Ramsey L., Sabree F., Lowe M., Brown K., Kyrpides N.C., Hugenholtz P. , "CRISPR Recognition Tool (CRT): a tool for automatic detection of clustered regularly interspaced palindromic repeats", BMC Bioinformatics, (2007), https://doi.org/10.1186/1471-2105-8-209.
* Laslett D., Canback B., "ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences.", Nucleic Acids Research (2004) https://doi.org/10.1093/nar/gkh152.
* Chen L., Yang J., Yao Z., Sun L., Shen Y., Jin Q., "VFDB: a reference database for bacterial virulence factors", Nucleic Acids Research (2005) https://doi.org/10.1093/nar/gki008.
* Alcock et al, "CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database." Nucleic Acids Research (2020) https:doi.org/10.1093/nar/gkz935.
