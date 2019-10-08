% Edited by hand
\docType{data}
\name{ex.cliqueGroups}
\alias{ex.cliqueGroups}
\title{Example m/z processed data}
\format{It is an 'xcmsSet' object of one sample with 126 features. Has been
obtained with parameters ppm = 15, method = "centWave",
peakwidth = c(5,20), snthresh = 10.
then features have been splitted into cliques with getCliques, with default
parameters and filter = T. Before getCliques it was used set.seed(2).}
\usage{
data(ex.cliqueGroups)
}
\description{
This dataset contains a mass sprectrometry data of metabolite standards
MS1 analyses were performed using an UHPLC system 
(1290 series, Agilent Technologies) coupled to
a 6550 ESI-QTOF MS (Agilent Technologies) operated in positive (ESI+)
electrospray ionization mode.

The original mzdata, which can be found at
"CliqueMS: a computational tool for annotating in-source metabolite
ions from LC-MS untargeted metabolomics data based on a coelution
similarity network" Senan et al, 2019 Bioinformatics
https://doi.org/10.1093/bioinformatics/btz207

The raw data was filtered from scan 0 to 700 with Proteowizard
mzconvert in order to have an smaller file

The metabolites in this example set are the following: 
thymine and uracil
}
\keyword{datasets}
