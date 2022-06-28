Once `install_databases.py` has been run, pharokka requires an input fasta file. An output directory can optionally be specified using -o. Otherwise, an `output/` directory will be created in your current working directory.

You can run pharokka using the following command:

`pharokka.py -i <fasta file> -o <output folder> -t <threads> -p <prefix>`

To specify a different database directory from the default:

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>`

To overwrite an existing output directory, use -f

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> -p <prefix>  -f `

pharokka defaults to 1 thread.
