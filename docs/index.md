`pharokka` is a fast phage annotation pipeline.

![Image](pharokka_logo.png)

`pharokka` uses PHANOTATE (McNair et al 2019) to conduct gene prediction, tRNAscan-SE 2 (Chan et al 2021) to call tRNAs, MinCED (Bland et al 2007) to detect CRISPRs and Aragorn (Laslett et al 2004) to detect tmRNAs. There is also the option to specify Prodigal (Hyatt et al 2010) implemented with Pyrodigal (Larralde, 2022) instead of PHANOTATE.

`pharokka` then uses the lightweight PHROGS database (Terzian et al 2021) for functional annotation of all predicted CDSs using MMseqs2 (Steinegger et al 2017), and as of v1.4.0, PyHMMER (Larralde and Zeller 2023) for more sensitive annotations. `pharokka` also matches each predicted CDS against the VFDB (Chen et al 2005) and CARD (Alcock et al 2020) databases to predict virulence factors and antimicrobial resistance, respectively. 

![Image](pharokka_workflow.png)




