phrokka is a fast phage annotation pipeline.

phrokka uses Phanotate (McNair et al (2019)) to conduct gene calling and tRNAscan-SE 2 (Chan et al (2021)) to call tRNAs.

phrokka then uses the lightweight PHROGS database (https://phrogs.lmge.uca.fr Terzian et al (2021)) for annotation of all predicted ORFs.

Specifically, each predicted ORF is compared against the PHROGS database first using mmseqs2 and secondly using hhsuite.
