# Creating and Using a Custom HMM Database

## When to use a custom HMM database

The PHROGs database covers the majority of known phage protein families, but it won't
have entries for newly characterised or niche families. A custom HMM database lets you
annotate those proteins alongside the standard PHROGs search.

Good candidates include:

- A set of structural proteins (e.g. tailspike variants, receptor-binding proteins)
  from a well-characterised phage genus not yet in PHROGs.
- A family of integrases, anti-CRISPR proteins, or any other group you have curated
  multiple sequence alignments (MSAs) for.
- Any group of proteins you want consistently flagged across a large annotation run.

The workflow is: collect FASTA MSAs (one file per family) → run `pharokka create-hmm`
to build the `.h3m` profile database → pass it to `pharokka run` with `--custom_hmm`.
Hits appear in `pharokka_cds_final_merged_output.tsv` under `custom_hmm_id`, and are
included in the `.gff` and `.gbk` outputs as `custom_annotation=...`.

## Building the database

`pharokka` v1.4.0 allows you to input a custom HMM profile database using the `--custom_hmm` parameter.

The easiest way to do this is with `pharokka create-hmm`, which comes with `pharokka`:

```bash
pharokka create-hmm -i <input directory of MSAs> -o <outdir> -p <prefix>
```

The resulting output you need to feed into `pharokka run` with the `--custom_hmm` parameter will be `outdir/prefix.h3m`.

## Input

`pharokka create-hmm` creates HMM profiles for each file in the input directory. It takes the name as the custom annotation for that HMM. It requires FASTA format MSAs (multiple sequence alignments) e.g.

```bash
>gene1
--------------------MERCMYDLSHFSFMTGRIGQLQTLTCIPVVSGDSFSVNFQGVFRLSPLRRNMVIDCMVDLFAFYVPYRHVYG-----DDWENFMLQGTDETVTFGTAD-FTTN-AVEYLGSQ-----R-TTGVAPLWRLASYNQIWNRYFRVPSDTSAILSNTALE---A-GVNKNKWGRYCARPKVIWSTGVDTTTD--AAD--REVTVTANQLDITDLDAVKGRYASEQRRDWFG-QRYNDILKNQFGGGAGTDADERPTLIMRKQHFLSGYDVDGTGDATLGQYSGKSAAVCDLQIPRKWFPEHGSL-WIMALLRFPTIHEEEIHYLFKKSQPTYDEISGDPDIWERKAPIEHQVQDFFTESTDTTSLGFMPYGQWYRYHPNVVHNQFTEVTGFSFVSA-RPT-----SKDTARYIQDDEYTPTFQTTSLGHWQSQARVDVEAVRPIPGPLKSIFAGV------------
>gene2
MSKSSGYATTSKPQEVVTRDIERKPYDLSHWSFKTGHIGRLQTFSVIPVVAGDSIELNMSAVLRLSPLRHFMYLDAVVDLFAFYVPHRHVYG-----SNWTAFLKAGVDEGSTLGTTT-FTGSNRPECLGVY-----FAAGQVLPSWSVHPYTMIWNNYFRDPQVDADAKTDATYLIGLTYPDPILKYGLPVCHLKRSWNTGVASGLS--TAD--YSLALSGGEVDLYQMSQLKGRLQTEQARDWFAFARYRDILKYTWDSEVNIDADQRPELIMRSTAWLSGKDVDGTDDATLGSYTGKSTGIFNLSFPSRFFAEHGSI-WIMGVVRFPPIVEAEAGYLNSKSEPTYKEISGDPSIIKHEPPISLNANEQFQG-ASSVDLGKIPYGQWYRESINMLSADYLEVSGHPFLLKSSIT-----SRATGVYVVPSQYDTVFSDTLMGHWNSQAFVDVRAKRFIPAPEDSIFAGTI-----------
>gene3
----------------MKGTDTRCLYDLDHWSHVSGNIGSLQTLSAIPVVAGDSMELNFTSLFRLSPLRRNLYLDAMVDLFAFYVPYRHVYG-----DTWIDFIKEGYDESQTLGTYT-IAAGEFINCTGAY-----LESQAVIPKWSIASYTRIWNRYFRHPTDSNEKDDDDLI----TASSGNQIYGYNCCHMKNIWSTGVDSELT--DDD--HKVDVTASKVDLLEIIQKKARFKTERQREWFG-QRYTDILDSVWGSNVNIDADERPELIMRTSTWLSGYDVDGTTETNIGTFSGKAFTTARLQFPMKYFNEHGTI-YIVALVRFPTIHAFERHYLFGKSEPTYKEIAGDPNVIRNEPPQAINPTDYIEN-ATATDVGLQPYAQWYRTHPSYTSLGFSNLDGHPFLQK-IIA-----NTDDAVYVDSNDYDEIFQTQQLQQWQSQGYVGLKAKRNIPDPRKSIFAGVK-----------
```

For example, if my input directory was called `my_tail_msa` and contained 3 MSA files named `tail_family_1.fasta`, `tail_family_2.fasta` and `tail_family_3.fasta`, the following command will create HMM profiles for all 3 MSA files, and combine them together into a `h3m` file with the prefix `tail`:

```bash
pharokka create-hmm -i my_tail_msa -o tail_hmms -p tail
```

The resulting file to be used with `--custom_hmm` is `tail_hmms/tail.h3m`.

## Output

So on a genome, you would run `pharokka` as follows:

```bash
pharokka run -i <input fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> --custom_hmm tail_hmms/tail.h3m
```

A gene hit to e.g. the `tail_family_2` HMM profile found by PyHMMER will be indicated in the `custom_hmm_id` column of the `pharokka_cds_final_merged_output.tsv` output file as `tail_family_2`, with information about the `custom_hmm_bitscore` and `custom_hmm_evalue` from PyHMMER. The annotations will supplement any PHROG annotation found and will be included in the `pharokka.gff` and `pharokka.gbk` files as `custom_annotation=tail_family_2` in this example.

```bash
usage: pharokka create-hmm [-h] [-i INDIR] [-o OUTDIR] [-p PREFIX] [-f] [-V]

pharokka create-hmm: Creates HMMs from FASTA formatted MSAs with PyHMMER for use with Pharokka v1.4.0 and higher.

options:
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        Input directory containing FASTA formatted Multiple Sequence Alignments.
  -o OUTDIR, --outdir OUTDIR
                        Output directory to store HMM profiles.
  -p PREFIX, --prefix PREFIX
                        Prefix used to name HMMs. The relevant file will be 'prefix'.h3m
  -f, --force           Overwrites the output directory.
  -V, --version         Print pharokka Version
```
