Why Use pharokka?
---------
`pharokka` is primarily intended for 2 use cases:

1. Where you would like to quickly annotate a phage without using a web-server (such as RAST/PATRIC or CPT Galaxy).
2. Where you have many phages or phage contigs that you would like to annotate quickly in a consistent fashion.

Additionally, `pharokka`, because it uses PHANOTATE for gene prediction by default, will likely predict more small phage genes than existing annotation programs that use different gene prediction programs like Prodigal (Hyatt et al 2010), which is designed for prokaryotes. See the PHANOTATE paper for details (McNair et al 2019 [https://doi.org/10.1093/bioinformatics/btz265]( https://doi.org/10.1093/bioinformatics/btz265)). See also [Fremin et al (2022)](https://pubmed.ncbi.nlm.nih.gov/35732113/) for evidence suggesting that unknown small phage genes may be more prevalent than previously thought.

Regarding VFDB and CARD databases,  chose conservative thresholds of 80% identity over 40% coverage were chosen as recommended for AMR gene searches in [Enault et al (2017)](https://doi.org/10.1038/ismej.2016.90).
