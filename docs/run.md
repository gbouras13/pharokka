Once `install_databases.py` has been run, pharokka requires an input fasta file. An output directory can optionally be specified using -o. Otherwise, an `output/` directory will be created in your current working directory.

You can run pharokka using the following command:

`pharokka.py -i <fasta file> -o <output folder> -t <threads> -p <prefix>`

A prefix is not required - by default it is `pharokka` e.g. `pharokka.gff`

To specify a different database directory (recommended) from the default:

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>`

To overwrite an existing output directory, use -f or `--force`:

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>  -f `

To specify a custom locus tag in the gff/genbank file, use -l:

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>  -f -l <locus_tag>`

In versions from v0.1.7 and afterwards, the ability to specify an E-value threshold for CDS functional assignment using MMseqs2 was added using the -e flag. It defaults to 1E-5.

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -e <E-value>`

To use Prodigal gene predictions instead of PHANOTATE use `-g prodigal`

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal`

If you are annotating more than 1 contig, it is recommended that you run pharokka in meta mode using the `-m` flag, which will enable pharokka to finish faster by making full use of all available threads when running PHANOTATE and tRNAscan-SE 2.

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -m`

There is also experimental support for alternative genetic codes if pharokka is run with prodigal as a gene predictor using the `-c` flag. See Prodigal's [documentation](https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#alternate-genetic-codes), along with [Yutin et al. 2021](https://doi.org/10.1038/s41467-022-32979-6) and [Peters et al. 2022](https://doi.org/10.1038/s41467-022-32979-6) for more information and background on phage stop codon reassignment.


`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal`

If you need to reorient your phage to begin at a certain point (for example, to reorient the phage to begin with the large terminase subunit), you can now acheive this using the --terminase command along with specifying the   --terminase_strand denoting the strandedness of the large terminase subunit (must be 'pos' or 'neg') and --terminase_start denoting the start co-ordinate. This requires the input to be a single contig, or else pharokka will fail and warn you. For example, if you annotate your phage using pharokka, and then see that your terminase large subunit begins at coordinate 9503 on the positive strand, you could run the following to reorient it to begin at the large terminase subunit:

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --terminase --terminase_start 9503 --terminase_strand pos `

Of course, you can also use this functionality to reorient your phage however you wish!

pharokka defaults to 1 thread.

```
usage: pharokka.py [-h] [-i INFILE] [-o OUTDIR] [-d DATABASE] [-t THREADS] [-f] [-p PREFIX] [-l LOCUSTAG] [-g GENE_PREDICTOR] [-m]
                   [-c CODING_TABLE] [-e EVALUE] [--terminase] [--terminase_strand TERMINASE_STRAND] [--terminase_start TERMINASE_START]
                   [-V] [--citation]

pharokka: fast phage annotation program

options:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input genome file in fasta format.
  -o OUTDIR, --outdir OUTDIR
                        Directory to write the output to.
  -d DATABASE, --database DATABASE
                        Database directory. If the databases have been installed in the default directory, this is not required. Otherwise specify the path.
  -t THREADS, --threads THREADS
                        Number of threads for mmseqs and hhsuite. Defaults to 1.
  -f, --force           Overwrites the output directory.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. This is not required.
  -l LOCUSTAG, --locustag LOCUSTAG
                        User specified locus tag for the gff/gbk files. This is not required. A random locus tag will be generated instead.
  -g GENE_PREDICTOR, --gene_predictor GENE_PREDICTOR
                        User specified gene predictor. Use "-g phanotate" or "-g prodigal". Defaults to phanotate (not required unless prodigal is desired).
  -m, --meta            meta mode for metavirome input samples
  -c CODING_TABLE, --coding_table CODING_TABLE
                        translation table for prodigal. Defaults to 11. Experimental only.
  -e EVALUE, --evalue EVALUE
                        E-value threshold for mmseqs2 PHROGs database search. Defaults to 1E-05.
  --terminase           Runs "terminase large subunit" re-orientation mode. Single genome input only and requires --terminase_strand and --terminase_start to be specified.
  --terminase_strand TERMINASE_STRAND
                        Strand of terminase large subunit. Must be "pos" or "neg".
  --terminase_start TERMINASE_START
                        Start coordinate of the terminase large subunit.
  -V, --version         Print pharokka Version
  --citation            Print pharokka Citation
  ```