# Running `pharokka`

Once `install_databases.py` has been run, pharokka requires an input fasta file. An output directory can be specified using -o. Otherwise, an `output/` directory will be created in your current working directory.

## Simple Parameters

You can run `pharokka` using the following command:

```
pharokka.py -i <fasta file> -o <output folder> -t <threads> -p <prefix>
```

`pharokka` defaults to 1 thread.

As of v1.4.0, you can speed up `pharokka` for small inputs (e.g. phage isolates or a small number of contigs) using `-fast`. Note: this will run only PyHMMER against PHROGs (so no VFDB or CARD annotations will be done).

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --fast
```

A prefix is not required - by default it is `pharokka` e.g. `pharokka.gff`

To specify a different database directory (recommended) from the default:

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>
```

To overwrite an existing output directory, use -f or `--force`:

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>  -f
```

To specify a custom locus tag in the gff/genbank file, use -l:

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>  -f -l <locus_tag>
```

In versions from v0.1.7 and afterwards, the ability to specify an E-value threshold for CDS functional assignment using MMseqs2 was added using the -e flag. It defaults to 1E-5.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -e <E-value>
```

To use Prodigal (pyrodigal) gene predictions instead of PHANOTATE use `-g prodigal`

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal
```

To use Prodigal-gv (pyrodigal-gv) gene predictions instead of PHANOTATE use `-g prodigal-gv`. This is recommended for metagenomic datasets where some phages likely have [alternate genetic codes](https://github.com/apcamargo/prodigal-gv).

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal-gv
```

If you are annotating more than 1 contig, it is recommended that you run pharokka in meta mode using the `-m` flag, which will enable pharokka to finish faster by making full use of all available threads when running tRNAscan-SE 2.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -m  -g prodigal-gv
```

As of v1.5.0, you can skip running mash to find the closest match for each contig in INPHARED using `skip_mash`.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --skip_mash
```

As of v1.5.0, you can skip running tRNAscan-SE 2, MinCED and Aragorn fusingg `skip_extra_annotations`.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --skip_extra_annotations
```



## Advanced Parameters

As of v1.4.0,  `pharokka` will automatically run MMseqs2 (PHROGs, CARD, VFDB) and PyHMMER (PHROGs). To turn off PyHMMER, use `--mmseqs2_only`.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --mmseqs2_only
```

As of v1.4.0,  `pharokka` will automatically run `--mmseqs2_only` if you specify `-m` or `--meta` for meta mode. To force `pharokka` to use PyHMMER on PHROGs in meta mode as well, use `--meta_hmm`. Note: this might take a long time if you have a large metagenome!

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -m --meta_hmm
```

You can also force `pharokka` to only use PyHMMER on PHROGs in meta mode using `--fast`.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -m --fast
```

If you are running meta mode, you can optionally specify split mode using the `-s` flag. This will add output directories with separated FASTA, genbank and gff files for each input contig called `single_fastas`, `single_gbks` and `single_gffs`.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -m -s
```

There is also support for alternative genetic codes if `pharokka` is run with prodigal as a gene predictor using the `-c` flag. See Prodigal's [documentation](https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#alternate-genetic-codes), along with [Yutin et al. 2021](https://doi.org/10.1038/s41467-022-32979-6) and [Peters et al. 2022](https://doi.org/10.1038/s41467-022-32979-6) for more information and background on phage stop codon reassignment.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal -c 4
```

As of v1.4.0 you can specify a genbank input file with `--genbank`. This will take all (custom) CDS calls in genbank file and PHANOTATE/pyrodigal will not be run. So if you have done manual gene curation, this option is recommended.

```
pharokka.py -i <GENBANK file> -o <output folder> -d <path/to/database_dir> -t <threads>  --genbank
```

As of v1.4.0 you want to automatically reorient your phage to begin with the large terminase subunit, you can use `--dnaapler` which will run [dnaapler](https://github.com/gbouras13/dnaapler) to automatically reorient your phage. The reoriented FASTA will be named `<prefix>_dnaapler_reoriented.fasta` in the output directory.

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --dnaapler
```

`--dnaapler` is recommended for reorientation as of v1.4.0. However, if you need to reorient your phage to begin at a certain point (for example, to reorient the phage to begin with the large terminase subunit) or if `dnaapler` isn't working for you, you can acheive this using the `--terminase` command along with specifying the  `--terminase_strand` denoting the strandedness of the large terminase subunit (must be 'pos' or 'neg') and `--terminase_start` denoting the start co-ordinate. This requires the input to be a single contig, or else `pharokka` will fail and warn you. For example, if you annotate your phage using `pharokka`, and then see that your terminase large subunit begins at coordinate 9503 on the positive strand, you could run the following to reorient it to begin at the large terminase subunit:

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --terminase --terminase_start 9503 --terminase_strand pos 
```

Additionally, it should be noted that if you want to reorient to begin with the large terminase subunit that is on the negative strand (`--terminase_strand neg`), you will  need to specify a `--terminase_start` coordinate that is actually the end (right) coordinate from the gff when --terminase_strand neg is specified. For example (thanks to Albert Vill) if your large terminase subunit is:

`ON631220.1    PHANOTATE    CDS    80073    82226    -1.721062366005452e+16    -    ...`

Then you need to specify the following to reorient the phage to begin with the large terminase subunit as follows:

```
pharokka.py \
  -i <fasta file> \
  -o <output folder> \
  -d <path/to/database_dir> \
  -t <threads> \
  --terminase \
  --terminase_start 82226 \
  --terminase_strand neg
```

Of course, you can also use this functionality to reorient your phage however you wish!

As of v1.6.0 you can specify a custom mash distance threshold vs the INPHARED database with `--mash_distance` 

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --mash_distance 0.05
```

You can also specify custom arguments to be passed to [MINced](https://doi.org/10.1186/1471-2105-8-209) with `--minced_args`

Note 2 things: (1) that you need to leave off the leading hyphen (i.e. `"minNR 2"` not `"-minNR 2"`) and (2) you will need to use quotation marks (as you extra arguments will contain spaces)

```
pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --minced_args "minNR 2 -minRL 21"
```


```
usage: pharokka.py [-h] [-i INFILE] [-o OUTDIR] [-d DATABASE] [-t THREADS] [-f] [-p PREFIX] [-l LOCUSTAG] [-g GENE_PREDICTOR] [-m] [-s]
                   [-c CODING_TABLE] [-e EVALUE] [--fast] [--mmseqs2_only] [--meta_hmm] [--dnaapler] [--custom_hmm CUSTOM_HMM] [--genbank]
                   [--terminase] [--terminase_strand TERMINASE_STRAND] [--terminase_start TERMINASE_START] [--skip_extra_annotations]
                   [--skip_mash] [--minced_args MINCED_ARGS] [--mash_distance MASH_DISTANCE] [-V] [--citation]

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
                        Number of threads. Defaults to 1.
  -f, --force           Overwrites the output directory.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. This is not required.
  -l LOCUSTAG, --locustag LOCUSTAG
                        User specified locus tag for the gff/gbk files. This is not required. A random locus tag will be generated instead.
  -g GENE_PREDICTOR, --gene_predictor GENE_PREDICTOR
                        User specified gene predictor. Use "-g phanotate" or "-g prodigal" or "-g prodigal-gv" or "-g genbank". 
                        Defaults to phanotate (not required unless prodigal is desired).
  -m, --meta            meta mode for metavirome input samples
  -s, --split           split mode for metavirome samples. -m must also be specified. 
                        Will output separate split FASTA, gff and genbank files for each input contig.
  -c CODING_TABLE, --coding_table CODING_TABLE
                        translation table for prodigal. Defaults to 11.
  -e EVALUE, --evalue EVALUE
                        E-value threshold for MMseqs2 database PHROGs, VFDB and CARD and PyHMMER PHROGs database search. Defaults to 1E-05.
  --fast, --hmm_only    Runs PyHMMER (HMMs) with PHROGs only, not MMseqs2 with PHROGs, CARD or VFDB. 
                        Designed for phage isolates, will not likely be faster for large metagenomes.
  --mmseqs2_only        Runs MMseqs2 with PHROGs, CARD and VFDB only (same as Pharokka v1.3.2 and prior). Default in meta mode.
  --meta_hmm            Overrides --mmseqs2_only in meta mode. Will run both MMseqs2 and PyHMMER.
  --dnaapler            Runs dnaapler to automatically re-orient all contigs to begin with terminase large subunit if found. 
                        Recommended over using '--terminase'.
  --custom_hmm CUSTOM_HMM
                        Run pharokka with a custom HMM profile database suffixed .h3m. 
                        Please use create this with the create_custom_hmm.py script.
  --genbank             Flag denoting that -i/--input is a genbank file instead of the usual FASTA file. 
                         The CDS calls in this file will be preserved and re-annotated.
  --terminase           Runs terminase large subunit re-orientation mode. 
                        Single genome input only and requires --terminase_strand and --terminase_start to be specified.
  --terminase_strand TERMINASE_STRAND
                        Strand of terminase large subunit. Must be "pos" or "neg".
  --terminase_start TERMINASE_START
                        Start coordinate of the terminase large subunit.
  --skip_extra_annotations
                        Skips tRNAscan-se, MINced and Aragorn.
  --skip_mash           Skips running mash to find the closest match for each contig in INPHARED.
  --minced_args MINCED_ARGS
                        extra commands to pass to MINced (please omit the leading hyphen for the first argument). You will need to use quotation marks e.g. --minced_args "minNR 2 -minRL 21"
  --mash_distance MASH_DISTANCE
                        mash distance for the search against INPHARED. Defaults to 0.2.
  -V, --version         Print pharokka Version
  --citation            Print pharokka Citation
  ```