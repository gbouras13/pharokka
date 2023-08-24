# Creating and Using a Custom HMM Database

`pharokka.py` v1.4.0 allows you to input a custom HMM profile database using the `--custom_hmm` parameter.

The easiest way to do this is with the `create_custom_hmm.py` script that comes with `pharokka` as follows:

`create_custom_hmm.py -i <input directory of MSAs> -o <outdir> -p <prefix>`

The resulting output you need to feed into `pharokka` in that example will be `outdir/prefix.h3m`

## Input

`create_custom_hmm.py` creates HMM profiles for each file in the input directory. It takes the name as the custom annotation for that HMM. It requires FASTA format MSAs (multiple sequence alignments) e.g.

```
>fasta
--------------------MERCMYDLSHFSFMTGRIGQLQTLTCIPVVSGDSFSVNFQGVFRLSPLRRNMVIDCMVDLFAFYVPYRHVYG-----DDWENFMLQGTDETVTFGTAD-FTTN-AVEYLGSQ-----R-TTGVAPLWRLASYNQIWNRYFRVPSDTSAILSNTALE---A-GVNKNKWGRYCARPKVIWSTGVDTTTD--AAD--REVTVTANQLDITDLDAVKGRYASEQRRDWFG-QRYNDILKNQFGGGAGTDADERPTLIMRKQHFLSGYDVDGTGDATLGQYSGKSAAVCDLQIPRKWFPEHGSL-WIMALLRFPTIHEEEIHYLFKKSQPTYDEISGDPDIWERKAPIEHQVQDFFTESTDTTSLGFMPYGQWYRYHPNVVHNQFTEVTGFSFVSA-RPT-----SKDTARYIQDDEYTPTFQTTSLGHWQSQARVDVEAVRPIPGPLKSIFAGV------------
>fasta2
MSKSSGYATTSKPQEVVTRDIERKPYDLSHWSFKTGHIGRLQTFSVIPVVAGDSIELNMSAVLRLSPLRHFMYLDAVVDLFAFYVPHRHVYG-----SNWTAFLKAGVDEGSTLGTTT-FTGSNRPECLGVY-----FAAGQVLPSWSVHPYTMIWNNYFRDPQVDADAKTDATYLIGLTYPDPILKYGLPVCHLKRSWNTGVASGLS--TAD--YSLALSGGEVDLYQMSQLKGRLQTEQARDWFAFARYRDILKYTWDSEVNIDADQRPELIMRSTAWLSGKDVDGTDDATLGSYTGKSTGIFNLSFPSRFFAEHGSI-WIMGVVRFPPIVEAEAGYLNSKSEPTYKEISGDPSIIKHEPPISLNANEQFQG-ASSVDLGKIPYGQWYRESINMLSADYLEVSGHPFLLKSSIT-----SRATGVYVVPSQYDTVFSDTLMGHWNSQAFVDVRAKRFIPAPEDSIFAGTI-----------
>fasta3
----------------MKGTDTRCLYDLDHWSHVSGNIGSLQTLSAIPVVAGDSMELNFTSLFRLSPLRRNLYLDAMVDLFAFYVPYRHVYG-----DTWIDFIKEGYDESQTLGTYT-IAAGEFINCTGAY-----LESQAVIPKWSIASYTRIWNRYFRHPTDSNEKDDDDLI----TASSGNQIYGYNCCHMKNIWSTGVDSELT--DDD--HKVDVTASKVDLLEIIQKKARFKTERQREWFG-QRYTDILDSVWGSNVNIDADERPELIMRTSTWLSGYDVDGTTETNIGTFSGKAFTTARLQFPMKYFNEHGTI-YIVALVRFPTIHAFERHYLFGKSEPTYKEIAGDPNVIRNEPPQAINPTDYIEN-ATATDVGLQPYAQWYRTHPSYTSLGFSNLDGHPFLQK-IIA-----NTDDAVYVDSNDYDEIFQTQQLQQWQSQGYVGLKAKRNIPDPRKSIFAGVK-----------
```

For example, if my input directory was called `my_tail_msa` and contained 3 MSA files names `tail_family_1.fasta`, `tail_family_2.fasta` and `tail_family_3.fasta`, the follow command will create HMM profiles for all 3 MSAs, and combine them together into a `h3m` file:

`create_custom_hmm.py -i my_tail_msa -o tail_hmms -p tail`

The resulting file to be used with `--custom_hmm` is `tail_hmms/tail.h3m`

So on a genome, you would run `pharokka` as follows:

`pharokka.py -i <input fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  --custom_hmm tail_hmms/tail.h3m`

Any hit to the tail_family_1 HMM profile found by PyHMMER will be indicated in the `custom_hmm_id` column of the `pharokka_cds_final_merged_output.tsv` output file, along with `custom_hmm_bitscore`	`custom_hmm_evalue` for bitscore and evalue from PyHMMER. The annotations will supplement any PHROG annotation found and will be included in the `pharokka.gff` and `pharokka.gbk` files as 'custom_annotation'.


```
usage: create_custom_hmm.py [-h] [-i INDIR] [-o OUTDIR] [-p PREFIX] [-f] [-V]

create_custom_hmm.py: Creates HMMs from FASTA formatted MSAs with PyHMMER for use with Pharokka v1.4.0 and higher.

options:
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        Input directory containing FASTA formatted Multiple Sequence Alignments.
  -o OUTDIR, --outdir OUTDIR
                        Output directory to store HMM profiles.
  -p PREFIX, --prefix PREFIX
                        Prefix used to name HMMs. The relevant file be 'prefix'.h3m
  -f, --force           Overwrites the output directory.
  -V, --version         Print pharokka Version
```