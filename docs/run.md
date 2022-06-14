Once `install_databases.py` has been run, phrokka requires an input fasta file. An output directory can optionally be specified using -o. Otherwise, an `output/` directory will be created in your current working directory.

You can run phrokka using the following command:

`phrokka.py -i <fasta file> -o <output folder> -t <threads>`

To specify a different database directory from the default:

`phrokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads>`

To overwrite an existing output directory, use -f

`phrokka.py -i <fasta file> -o <output folder> -t <threads>  -f`

phrokka defaults to 1 thread.
