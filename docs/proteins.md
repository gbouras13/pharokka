# Proteins

`pharokka` v1.4.0 implements `pharokka proteins`, which allows for the fast annotation of amino acid FASTA formatted input proteins with PHROGs, CARD and VFDB.

You can run `pharokka proteins` in the following way (most basic command below). By default this will run MMSeqs2 on PHROGs, CARD and VFDB, along with PyHMMER on PHROGs

`pharokka proteins -i <input protein FASTA file> -d <database directory> -o <output folder> -t <threads>`

To run only PHROGs with HMMs (PyHMMER), use `--hmm_only`

`pharokka proteins -i <input protein FASTA file> -d <database directory> -o <output folder> -t <threads> --hmm_only `

To run only MMseqs2 on PHROGs, CARD and VFDB, use `--mmseqs2_only`

`pharokka proteins -i <input protein FASTA file> -d <database directory> -o <output folder> -t <threads> --mmseqs2_only `

## Output

`pharokka proteins` generates 3 output files:

* `_proteins.faa` containing the original protein FASTA file with the header updated with the new PHROG annotation
* `_proteins_full_merged_output.tsv` containing the full information for each protein, similar to `_cds_functions.tsv` for `pharokka run` but without the gene prediction information. This will include VFDB and CARD information.
* `_proteins_summary_output.tsv` containing the following 5 columns (where ID is name of each input protein)

`ID	length (AA)	phrog	annot	category`

e.g.

| ID        | length | phrog | annot                  | category           |
|-----------|------- |-------|------------------------|--------------------|
| protein_1 | 456    | 3717  | terminase large subunit| head and packaging |
| protein_2 | 243    | 315   | endolysin              | lysis              |


```bash
usage: pharokka proteins [-h] [-i INFILE] [-o OUTDIR] [-d DATABASE] [-t THREADS] [-f] [-p PREFIX] [-e EVALUE] [--hmm_only] [--mmseqs2_only] [-V] [--citation]

pharokka proteins: fast phage protein annotation script

options:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input proteins file in amino acid fasta format.
  -o OUTDIR, --outdir OUTDIR
                        Directory to write the output to.
  -d DATABASE, --database DATABASE
                        Database directory. If the databases have been installed in the default directory, this is not required. Otherwise specify the path.
  -t THREADS, --threads THREADS
                        Number of threads. Defaults to 1.
  -f, --force           Overwrites the output directory.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. This is not required.
  -e EVALUE, --evalue EVALUE
                        E-value threshold for MMseqs2 database PHROGs, VFDB and CARD and pyhmmer PHROGs database search. Defaults to 1E-05.
  --hmm_only            Runs pyhmmer (HMMs) with PHROGs only, not MMseqs2 with PHROGs, CARD or VFDB.
  --mmseqs2_only        Runs MMseqs2 with PHROGs, CARD and VFDB only (same as Pharokka v1.3.2 and prior).
  -V, --version         Print pharokka Version
  --citation            Print pharokka Citation
```