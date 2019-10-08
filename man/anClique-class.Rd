% edited by hand
\docType{class}
\name{anClique-class}
\alias{anClique-class}
\alias{show,anClique-method}
\alias{getPeaklistanClique,anClique-method}
\alias{getPeaklistanClique}
\alias{getNetanClique,anClique-method}
\alias{getNetanClique}
\alias{getIsolistanClique,anClique-method}
\alias{getIsolistanClique}
\alias{getlistofCliques,anClique-method}
\alias{getlistofCliques}
\alias{hasAnnotation,anClique-method}
\alias{hasAnnotation}
\alias{hasCliques,anClique-method}
\alias{hasCliques}
\alias{hasIsotopes,anClique-method}
\alias{hasIsotopes}
\alias{getIsolistanClique<-,anClique-method}
\alias{getIsolistanClique<-}
\alias{getNetanClique<-,anClique-method}
\alias{getNetanClique<-}
\alias{getlistofCliques<-,anClique-method}
\alias{getlistofCliques<-}
\alias{getPeaklistanClique<-,anClique-method}
\alias{getPeaklistanClique<-}
\alias{hasAnnotation<-,anClique-method}
\alias{hasAnnotation<-}
\alias{hasCliques<-,anClique-method}
\alias{hasCliques<-}
\alias{hasIsotopes<-,anClique-method}
\alias{hasIsotopes<-}
\title{'anClique' S4 class for annotating isotopes and adducts}
\usage{
\S4method{show}{anClique}(object)

\S4method{getPeaklistanClique}{anClique}(object)

\S4method{getNetanClique}{anClique}(object)

\S4method{getIsolistanClique}{anClique}(object)

\S4method{getlistofCliques}{anClique}(object)

\S4method{hasAnnotation}{anClique}(object)

\S4method{hasCliques}{anClique}(object)

\S4method{hasIsotopes}{anClique}(object)

\S4method{getIsolistanClique}{anClique}(object) <- value

\S4method{getNetanClique}{anClique}(object) <- value

\S4method{getlistofCliques}{anClique}(object) <- value

\S4method{getPeaklistanClique}{anClique}(object) <- value

\S4method{hasAnnotation}{anClique}(object) <- value

\S4method{hasCliques}{anClique}(object) <- value

\S4method{hasIsotopes}{anClique}(object) <- value
}
\arguments{
\item{object}{'anClique' S4 object.}

\item{value}{Is the new variable which can be a 'peaklist',
a 'network', a 'isotopes' a 'cliques' a 'cliquesFound'
a 'isoFound' or 'anFound' and it's set by the user.}
}
\value{
An 'anClique' object with annotation of isotopes, adducts
and fragments, and information about the annotation process.
}
\description{
S4 Class \code{anClique-class} for annotating isotopes and adducts
in processed m/z data. Features are first
grouped based on a similarity network algorithm and then
annotation of isotopes and adducts is performed in each group.
The class contains the following slots.
}
\section{Methods (by generic)}{
\itemize{
\item \code{show}: show information about the object

\item \code{getPeaklistanClique}: get the list of features with
current annotation

\item \code{getNetanClique}: get the correlation network

\item \code{getIsolistanClique}: get the table of isotopes

\item \code{getlistofCliques}: get the list of the clique groups

\item \code{hasAnnotation}: is 'TRUE' if annotation has been computed

\item \code{hasCliques}: is 'TRUE' if cliques have been computed

\item \code{hasIsotopes}: is 'TRUE' if isotopes have been computed

\item \code{getIsolistanClique<-}: set the table of isotopes

\item \code{getNetanClique<-}: set the network of correlation

\item \code{getlistofCliques<-}: set the list of clique groups

\item \code{getPeaklistanClique<-}: set the list of features

\item \code{hasAnnotation<-}: set if annotation has been computed

\item \code{hasCliques<-}: set if cliques have been computed

\item \code{hasIsotopes<-}: set if isotopes have been computed
}}

\section{Slots}{

\describe{
\item{\code{'peaklist'}}{Is a data.frame with m/z, retention time
and intensity information for each feature. It also contains
adduct and isotope information if annotation has been performed.}

\item{\code{'network'}}{Is an igraph undirected network of similarity 
used to compute groups of features before annotation.}

\item{\code{'cliques'}}{Is a list that contains the groups of features.
Each id corresponds to a row in the peaklist.}

\item{\code{'isotopes'}}{Is a data.frame with the column 'feature' for
feature id,  column 'charge' for the charge, column 'grade' that starts with 0
and it is 1 for the first isotope, 2 for the second and so on and
column 'cluster' which labels each group of features that are
isotopes.}

\item{\code{'cliquesFound'}}{is
TRUE if clique groups have been computed,}

\item{\code{'isoFound'}}{is
TRUE if isotopes have been annotated,}

\item{\code{'anFound'}}{is
TRUE if annotation of adducts have been computed.}
}}

\examples{
mzfile <- system.file("standards.mzXML", package = "cliqueMS")
library(xcms)
mzraw <- readMSData(files = mzfile, mode = "onDisk")
cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
mzData <- findChromPeaks(object = mzraw, param = cpw)
ex.anClique <- createanClique(mzdata = mzData)
show(ex.anClique)
}
\seealso{
\code{\link{createanClique}}
}
