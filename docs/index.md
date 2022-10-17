pharokka is a fast phage annotation pipeline.

pharokka uses PHANOTATE (McNair et al 2019) to conduct gene prediction, tRNAscan-SE 2 (Chan et al 2021) to call tRNAs, MinCED (Bland et al 2007) to detect CRISPRs and Aragorn (Laslett et al 2004) to detect tmRNAs. There is also the option to specify Prodigal (Hyatt et al 2010) instead of PHANOTATE.

pharokka then uses the lightweight PHROGS database (Terzian et al 2021) for functional annotation of all predicted CDSs using mmseqs2 (Steinegger et al 2017). As of v1.0.0, pharokka also runs each predicted CDS against the VFDB and CARD databases to predict virulence factors and antimicrobial resistance, respectively. 
