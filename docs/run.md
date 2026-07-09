# Protocols

We have recently published a [protocols paper](https://doi.org/10.1002/cpz1.70405) in _Current Protocols_ outlining how to run Pharokka, [Phold](https://github.com/gbouras13/phold), [Phynteny](https://github.com/susiegriggo/Phynteny_transformer) for annotation and visualisation via our [Phold Plot Wasm application](https://gbouras13.github.io/phold-plot-wasm-app/).

We highly recommend reading and following this protocol for users new to phage annotation.

If you use this protocol, please cite

> Bouras G., Grigson S.R., Durr L., Papudeshi B., Vreugde S.,
> Mallawaarachchi V., Vreugde S., Edwards R.A. 
>  
> *Decoding Viral Dark Matter: Metagenomic Prokaryotic Virus Characterization With Pharokka, Phold, and Phynteny*  
> **Current Protocols**, Volume 6, Number 7, 6 July 2026  
> [https://doi.org/10.1002/cpz1.70405](https://doi.org/10.1002/cpz1.70405)

# Running `pharokka`

Once `pharokka install` has been run, `pharokka run` requires an input FASTA file. An output directory can be specified using `-o`. Otherwise, an `output/` directory will be created in your current working directory.

## Quick-start Examples

**Single phage isolate** (most common case):

```bash
pharokka run -i phage.fasta -o output_dir -d /path/to/database -t 8
```

**Metavirome / multiple contigs** (recommended for >1 contig):

```bash
pharokka run -i contigs.fasta -o output_dir -d /path/to/database -t 8 -m
```

In meta mode (`-m`), `pharokka` defaults to `prodigal-gv` for gene prediction and `--mmseqs2_only` for annotation, which is the fastest and most appropriate combination for metagenomic datasets.

---

## Simple Parameters

`pharokka` defaults to 1 thread.

A prefix is not required — by default it is `pharokka` (e.g. `pharokka.gff`).

To specify a different database directory (recommended) from the default:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -p <prefix>
```

To overwrite an existing output directory, use `-f` or `--force`:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -p <prefix> -f
```

To specify a custom locus tag in the gff/genbank file, use `-l`:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -p <prefix> -f -l <locus_tag>
```

To speed up annotation of phage isolates (small inputs), use `--fast`. This runs only PyHMMER against PHROGs and skips VFDB and CARD:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --fast
```

To specify an E-value threshold for CDS functional assignment (MMseqs2 and PyHMMER), use `-e`. Defaults to 1E-5:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -e <E-value>
```

---

## Gene Predictors

`pharokka` defaults to PHANOTATE for single-phage mode and `prodigal-gv` for meta mode.

To use Prodigal (pyrodigal):

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -g prodigal
```

To use Prodigal-gv (pyrodigal-gv) — recommended for metagenomic datasets where phages may use [alternate genetic codes](https://github.com/apcamargo/prodigal-gv):

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -g prodigal-gv
```

To use Pyrodigal-rv — recommended for RNA phages (and potentially other RNA viruses):

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -g pyrodigal-rv
```

---

## Meta Mode

If you are annotating more than 1 contig, it is recommended that you run `pharokka` in meta mode using the `-m` flag. This enables `pharokka` to finish faster by making full use of all available threads when running tRNAscan-SE 2, and defaults to `prodigal-gv` for gene prediction.

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -m
```

If you are running meta mode, you can optionally specify split mode using the `-s` flag. This will add output directories with separated FASTA, genbank and gff files for each input contig called `single_fastas`, `single_gbks` and `single_gffs`.

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -m -s
```

---

## Advanced Parameters

As of v1.4.0, `pharokka` will automatically run MMseqs2 (PHROGs, CARD, VFDB) and PyHMMER (PHROGs). To turn off PyHMMER, use `--mmseqs2_only`:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --mmseqs2_only
```

As of v1.4.0, `pharokka` will automatically run `--mmseqs2_only` if you specify `-m` or `--meta` for meta mode. To force `pharokka` to use PyHMMER on PHROGs in meta mode as well, use `--meta_hmm`. Note: this might take a long time if you have a large metagenome!

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -m --meta_hmm
```

You can also force `pharokka` to only use PyHMMER on PHROGs in meta mode using `--fast`.

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -m --fast
```

As of v1.5.0, you can skip running mash to find the closest match for each contig in INPHARED using `--skip_mash`:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --skip_mash
```

As of v1.5.0, you can skip running tRNAscan-SE 2, MinCED and Aragorn using `--skip_extra_annotations`:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --skip_extra_annotations
```

There is also support for alternative genetic codes if `pharokka` is run with prodigal as a gene predictor using the `-c` flag. See Prodigal's [documentation](https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#alternate-genetic-codes), along with [Yutin et al. 2021](https://doi.org/10.1038/s41467-022-32979-6) and [Peters et al. 2022](https://doi.org/10.1038/s41467-022-32979-6) for more information:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -g prodigal -c 4
```

As of v1.4.0 you can specify a genbank input file with `--genbank`. This will take all (custom) CDS calls in the genbank file and PHANOTATE/pyrodigal will not be run. Useful if you have done manual gene curation:

```bash
pharokka run -i <GENBANK file> -o <output folder> -d <path/to/database_dir> -t <threads> --genbank
```

As of v1.4.0, to automatically reorient your phage to begin with the large terminase subunit, use `--dnaapler` which will run [dnaapler](https://github.com/gbouras13/dnaapler). The reoriented FASTA will be named `<prefix>_dnaapler_reoriented.fasta`:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --dnaapler
```

`--dnaapler` is recommended for reorientation as of v1.4.0. However, if you need to reorient your phage to begin at a certain point (for example, to reorient the phage to begin with the large terminase subunit) or if `dnaapler` isn't working for you, you can achieve this using the `--terminase` command along with specifying `--terminase_strand` (must be 'pos' or 'neg') and `--terminase_start`. This requires the input to be a single contig. For example, if your terminase large subunit begins at coordinate 9503 on the positive strand:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --terminase --terminase_start 9503 --terminase_strand pos
```

Note: if you want to reorient to begin with a terminase on the negative strand (`--terminase_strand neg`), specify the *end* (right) coordinate from the gff as `--terminase_start`. For example (thanks to Albert Vill) if your large terminase subunit is:

`ON631220.1    PHANOTATE    CDS    80073    82226    -1.721062366005452e+16    -    ...`

```bash
pharokka run \
  -i <fasta file> \
  -o <output folder> \
  -d <path/to/database_dir> \
  -t <threads> \
  --terminase \
  --terminase_start 82226 \
  --terminase_strand neg
```

As of v1.6.0 you can specify a custom mash distance threshold vs the INPHARED database with `--mash_distance`:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --mash_distance 0.05
```

You can also specify custom arguments to be passed to [MINced](https://doi.org/10.1186/1471-2105-8-209) with `--minced_args`.

Note 2 things: (1) leave off the leading hyphen for the first argument (i.e. `"minNR 2"` not `"-minNR 2"`), and (2) use quotation marks as the extra arguments will contain spaces:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --minced_args "minNR 2 -minRL 21"
```

As of v1.7.4 you can specify the bacterial `tRNAscan-SE` model using `--trna_scan_model bacterial`. Otherwise `pharokka` uses the general model by default. See the [tRNAscan-SE paper](https://doi.org/10.1093/nar/gkab688) for more information:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --trna_scan_model bacterial
```

As of v1.9.0, for very large datasets you can reverse the MMseqs2 search direction (PHROG profile database as target rather than query) using `--reverse_mmseqs2`. This is only recommended for enormous datasets where the standard orientation is slow:

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --reverse_mmseqs2
```

As of v1.9.0, you can control MMseqs2 profile search sensitivity with `--sensitivity` (default is 8.5):

```bash
pharokka run -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --sensitivity 7.5
```

---

## Full Usage

```
usage: pharokka run [-h] [-i INFILE] [-o OUTDIR] [-d DATABASE] [-t THREADS] [-f] [-p PREFIX] [-l LOCUSTAG] [-g GENE_PREDICTOR] [-m] [-s] [-c CODING_TABLE] [-e EVALUE]
                   [--fast] [--mmseqs2_only] [--meta_hmm] [--dnaapler] [--custom_hmm CUSTOM_HMM] [--genbank] [--terminase] [--terminase_strand TERMINASE_STRAND]
                   [--terminase_start TERMINASE_START] [--skip_extra_annotations] [--skip_mash] [--minced_args MINCED_ARGS] [--mash_distance MASH_DISTANCE]
                   [--trna_scan_model {general,bacterial}] [--keep_raw_prodigal] [--reverse_mmseqs2] [--sensitivity SENSITIVITY] [-V] [--citation]

pharokka: fast phage annotation program

options:
  -h, --help            show this help message and exit
  -i, --infile INFILE   Input genome file in fasta format.
  -o, --outdir OUTDIR   Directory to write the output to.
  -d, --database DATABASE
                        Database directory. If the databases have been installed in the default directory, this is not required. Otherwise specify the path.
  -t, --threads THREADS
                        Number of threads. Defaults to 1.
  -f, --force           Overwrites the output directory.
  -p, --prefix PREFIX   Prefix for output files. This is not required.
  -l, --locustag LOCUSTAG
                        User specified locus tag for the gff/gbk files. This is not required. A random locus tag will be generated instead.
  -g, --gene_predictor GENE_PREDICTOR
                        User specified gene predictor. Use "-g phanotate" or "-g prodigal" or "-g prodigal-gv" or "-g pyrodigal-rv" or "-g genbank". 
                        Defaults to phanotate usually and prodigal-gv in meta mode.
  -m, --meta            meta mode for metavirome input samples
  -s, --split           split mode for metavirome samples. -m must also be specified. 
                        Will output separate split FASTA, gff and genbank files for each input contig.
  -c, --coding_table CODING_TABLE
                        translation table for prodigal. Defaults to 11.
  -e, --evalue EVALUE   E-value threshold for MMseqs2 database PHROGs, VFDB and CARD and PyHMMER PHROGs database search. Defaults to 1E-05.
  --fast, --hmm_only    Runs PyHMMER (HMMs) with PHROGs only, not MMseqs2 with PHROGs, CARD or VFDB. 
                        Designed for phage isolates, will not likely be faster for large metagenomes.
  --mmseqs2_only        Runs MMseqs2 with PHROGs, CARD and VFDB only (same as Pharokka v1.3.2 and prior). Default in meta mode.
  --meta_hmm            Overrides --mmseqs2_only in meta mode. Will run both MMseqs2 and PyHMMER.
  --dnaapler            Runs dnaapler to automatically re-orient all contigs to begin with terminase large subunit if found. 
                        Recommended over using '--terminase'.
  --custom_hmm CUSTOM_HMM
                        Run pharokka with a custom HMM profile database suffixed .h3m. 
                        Please create this with `pharokka create-hmm`.
  --genbank             Flag denoting that -i/--input is a genbank file instead of the usual FASTA file. 
                        The CDS calls in this file will be preserved and re-annotated.
  --terminase           Runs terminase large subunit re-orientation mode. 
                        Single genome input only and requires --terminase_strand and --terminase_start to be specified.
  --terminase_strand TERMINASE_STRAND
                        Strand of terminase large subunit. Must be "pos" or "neg".
  --terminase_start TERMINASE_START
                        Start coordinate of the terminase large subunit.
  --skip_extra_annotations
                        Skips tRNAscan-SE 2, MinCED and Aragorn.
  --skip_mash           Skips running mash to find the closest match for each contig in INPHARED.
  --minced_args MINCED_ARGS
                        Extra commands to pass to MINced (please omit the leading hyphen for the first argument). You will need to use quotation marks e.g. --minced_args "minNR 2 -minRL 21"
  --mash_distance MASH_DISTANCE
                        mash distance for the search against INPHARED. Defaults to 0.2.
  --trna_scan_model {general,bacterial}
                        tRNAscan-SE model. Defaults to general.
  --keep_raw_prodigal   Keeps raw prodigal header information.
  --reverse_mmseqs2     MMseqs2 database as target not query. Only recommended for enormous datasets.
  --sensitivity SENSITIVITY
                        MMseqs2 profile search sensitivity. Defaults to 8.5.
  -V, --version         Print pharokka Version
  --citation            Print pharokka Citation
```
