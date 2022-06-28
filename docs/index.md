pharokka is a fast phage annotation pipeline.

pharokka uses PHANOTATE (McNair et al 2019 [https://github.com/deprekate/PHANOTATE](PHANOTATE)) to conduct gene prediction and tRNAscan-SE 2 (Chan et al 2021 [https://github.com/UCSC-LoweLab/tRNAscan-SE](tRNAscan-SE)) to call tRNAs.

pharokka then uses the lightweight PHROGS database (Terzian et al 2021 [https://phrogs.lmge.uca.fr](PHROGS)) for functional annotation of all predicted CDSs.

Each predicted CDS is compared against the PHROGS database twice first using mmseqs2 and secondly using hhsuite (which usually assigns PHROGs to short, hypothetical genes that are predicted by PHANOTATE).
