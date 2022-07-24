pharokka is a fast phage annotation pipeline.

pharokka uses PHANOTATE (McNair et al 2019 [https://github.com/deprekate/PHANOTATE](PHANOTATE)) to conduct gene prediction, tRNAscan-SE 2 (Chan et al 2021 [https://github.com/UCSC-LoweLab/tRNAscan-SE](tRNAscan-SE)) to call tRNAs, MinCED (Bland et al 2007 [https://github.com/ctSkennerton/minced](MinCED)) to detect CRISPRs and Aragorn (Laslett et al 2004 [http://www.ansikte.se/ARAGORN/](Aragorn)) to detect tmRNAs .

pharokka then uses the lightweight PHROGS database (Terzian et al 2021 [https://phrogs.lmge.uca.fr](PHROGS)) for functional annotation of all predicted CDSs using mmseqs2 (Steinegger et al 2017 [https://github.com/soedinglab/MMseqs2](mmseqs2)).
