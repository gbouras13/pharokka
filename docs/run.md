Once `install_databases.py` has been run, pharokka requires an input fasta file. An output directory can optionally be specified using -o. Otherwise, an `output/` directory will be created in your current working directory.


If you have already downloaded databases for earlier versions of pharokka, these will need to be re-downloaded.

Once the databases have finished downloading, to run pharokka



To specify a prefix for the output files:

```
pharokka.py -i <fasta file> -o <output folder> -t <threads> -p <prefix>
```

To specify a different database directory (recommended):

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -p <prefix>
```

To overwrite an existing output directory, use -f

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -p <prefix> -f
```

To use Prodigal instead of PHANOTATE use `-g prodigal`

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal
```

The `-m` option indicates meta mode designed for metavirome input. pharokka should work with metagenome assembled viral contigs with PHANOTATE automatically. With prodigal, please add the `-m` flag. As of v1.1.0, `-m` has added multi-threaded support for tRNAscan-SE2 and PHANOTATE, speeding their runtime considerably. 

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -m
```

In v0.1.7, the ability to specify an E-value threshold for PHROGs CDS functional assignment using mmseqs2 was added using the -e flag. It defaults to 1E-5.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -e <E-value>
```

pharokka defaults to 1 thread.

```
usage: pharokka.py [-h] -i INFILE [-o OUTDIR] [-d DATABASE] [-t THREADS] [-f] [-p PREFIX] [-l LOCUSTAG]
                   [-g GENE_PREDICTOR] [-m] [-c CODING_TABLE] [-e EVALUE] [-V]

pharokka: fast phage annotation program

optional arguments:
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
  -V, --version         Version
  ```