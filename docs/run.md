Once `install_databases.py` has been run, pharokka requires an input fasta file. An output directory can optionally be specified using -o. Otherwise, an `output/` directory will be created in your current working directory.

You can run pharokka using the following command:

`pharokka.py -i <fasta file> -o <output folder> -t <threads> -p <prefix>`

A prefix is not required - by default it is `pharokka` e.g. `pharokka.gff`

To specify a different database directory from the default:

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>`

To overwrite an existing output directory, use -f

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>  -f `

To specify a custom locus tag in the gff/genbank file, use -l

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>  -f -l <locus_tag>`

In v0.1.7, the ability to specify an E-value threshold for CDS functional assignment using mmseqs2 was added using the -e flag. It defaults to 1E-5.

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -e <E-value>`

To use Prodigal gene predictions instead of PHANOTATE use `-g prodigal`

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal`

pharokka should work with metagenome assembled viral contigs with PHANOTATE automatically. With prodigal, please add the -m flag

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal -m`


pharokka defaults to 1 thread.
