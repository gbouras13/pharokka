pharokka is a fast phage annotation pipeline.

pharokka uses PHANOTATE (McNair et al 2019) to conduct gene prediction, tRNAscan-SE 2 (Chan et al 2021) to call tRNAs, MinCED (Bland et al 2007) to detect CRISPRs and Aragorn (Laslett et al 2004) to detect tmRNAs. There is also the option to specify Prodigal (Hyatt et al 2010) instead of PHANOTATE.

pharokka then uses the lightweight PHROGS database (Terzian et al 2021) for functional annotation of all predicted CDSs using MMseqs2 (Steinegger et al 2017). As of v1.0.0, pharokka also matches each predicted CDS against the VFDB (Chen et al 2005) and CARD (Alcock et al 2020) databases to predict virulence factors and antimicrobial resistance, respectively. 

As of v1.2.0, Pharokka  implements a major new feature. It quickly matches each input contig against the  [INPHARED](https://github.com/RyanCook94/inphared) database ([paper](http://doi.org/10.1089/phage.2021.0007)) using [mash](https://doi.org/10.1186/s13059-016-0997-x) distances, which may be useful if you are annotating novel phages or metagenomic input samples. If you use this feature, please make sure you cite INPHARED. 

v 1.2.0 also adds the ability to re-orient your phage specifying a coordinate and strandedness using the terminase large subunit reorientation mode, then annotate the re-oriented phage. 



