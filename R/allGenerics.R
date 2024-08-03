## C

#' @export
#' @title 'createanClique' generic function to create an object
#' of class 'anClique'.
#'
#' @description
#' \code{createanClique} creates an 'anClique' object from processed m/z data.e
#' @param mzdata An object with processed m/z data. See methods for
#' valid class types.
#' @return
#' An 'anClique' S4 object with all elements to perform clique grouping,
#' isotope annotation and adduct annotation.
#' @examples
#' ## Using a 'XCMSnExp' object
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' library(xcms)
#' mzraw <- MSnbase::readMSData(files = mzfile, mode = "onDisk")
#' cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' mzData <- findChromPeaks(object = mzraw, param = cpw)
#' ex.anClique <- createanClique(mzdata = mzData)
#' show(ex.anClique)
#' 
#' ## Using a 'xcmsSet' object
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' msSet <- xcms::xcmsSet(files = mzfile, method = "centWave",
#' ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' ex.anClique <- createanClique(msSet)
#' @seealso \code{\link{anClique}}
setGeneric("createanClique", function(mzdata) {
    standardGeneric("createanClique")
})

#' @export
#' @title Generic function to create a similarity network
#' from processed m/z data
#'
#' @description
#' This function creates a similarity network with nodes as features
#' and weighted edges as the cosine similarity between those nodes.
#' Edges with weights = 0 are not included in the network. Nodes
#' without edges are not included in the network. This network will
#' be used to define clique groups and find annotation within this
#' groups.
#'
#' @details Signal processing algorithms may output artefact features.
#' Sometimes they produce two artefact features which are almost identical
#' This artefacts may lead to errors in the computation of the clique
#' groups, so it is recommended to set 'filter' = TRUE to drop repeated
#' features.
#' @param mzdata An object of class 'xcmsSet' or 'XCMSnExp'
#' with processed m/z data.
#' @param peaklist Is a data.frame feature info for m/z data.
#' put each feature in a row and a column 'mz' for mass data, 
#' retention time column 'rt' and intensity in column 'maxo'.
#' @param filter If TRUE, filter out very similar features
#' that have a correlation similarity > 0.99 and equal values of m/z,
#' retention time and intensity.
#' @param mzerror Relative error for m/z, if relative error 
#' between two features is below that value that features
#' are considered with similar m/z value.
#' @param rtdiff Relative error for retention time, if 
#' relative error between two features is below that value
#' that features are considered with similar retention time
#' @param intdiff Relative error for intensity, if relative
#' error between two features is below that value that
#' features are considered with similar intensity
#' @return This function returns a list with the similarity
#' network and the filtered peaklist if 'filter' = TRUE. If
#' filter = FALSE the peaklist is returned unmodified.
#' @examples
#' ## Using a 'xcmsSet' object
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' require(xcms)
#' rawMS <- MSnbase::readMSData(files = mzfile, mode = "onDisk")
#' cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' mzData <- findChromPeaks(rawMS, cpw)
#' peaklist = as.data.frame(chromPeaks(mzData))
#' netlist = createNetwork(mzData, peaklist, filter = TRUE)
#' 
#' ## Using a 'XCMSnExp' object
#' require(xcms)
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' rawMS <- MSnbase::readMSData(files = mzfile, mode = "onDisk")
#' cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' mzData <- findChromPeaks(rawMS, cpw)
#' peaklist = as.data.frame(chromPeaks(mzData))
#' netlist = createNetwork(mzData, peaklist, filter = TRUE)
#' @seealso \code{\link{getCliques}}
setGeneric("createNetwork", function(mzdata, peaklist, filter = TRUE,
    mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    standardGeneric("createNetwork")
})

## G

#' @export
#' @describeIn anClique get the list of the clique groups
setGeneric("getlistofCliques", function(object) {
    standardGeneric ("getlistofCliques")
})
#' @describeIn anClique set the list of the clique groups
setGeneric("getlistofCliques<-", function(object,value) {
    standardGeneric ("getlistofCliques<-")
})

#' @export
#' @describeIn anClique get the isotope annotation
setGeneric("getIsolistanClique", function(object) {
    standardGeneric ("getIsolistanClique")
})

#' @export
#' @describeIn anClique set the isotope annotation
setGeneric("getIsolistanClique<-", function(object,value) {
    standardGeneric ("getIsolistanClique<-")
})

#' @export
#' @describeIn get the similarity Network
setGeneric("getNetanClique", function(object) {
    standardGeneric ("getNetanClique")
})
#' @export
#' @describeIn anClique set the similarity Network
setGeneric("getNetanClique<-", function(object, value) {
    standardGeneric ("getNetanClique<-")
})

#' @export
#' @describeIn anClique get the list of features with the current annotation
setGeneric("getPeaklistanClique", function(object) {
    standardGeneric ("getPeaklistanClique")
})
#' @export
#' @describeIn anClique set the list of features and annotation
setGeneric("getPeaklistanClique<-", function(object, value) {
    standardGeneric ("getPeaklistanClique<-")
})


## H

#' @export
setGeneric("hasAnnotation", function(object) {
    standardGeneric("hasAnnotation")
})
#' @export
setGeneric("hasAnnotation<-", function(object,value) {
    standardGeneric("hasAnnotation<-")
})

#' @export
setGeneric("hasCliques", function(object) {
    standardGeneric("hasCliques")
})
#' @export
setGeneric("hasCliques<-", function(object,value) {
    standardGeneric("hasCliques<-")
})

#' @export
setGeneric("hasIsotopes", function(object) {
    standardGeneric("hasIsotopes")
})
#' @export
setGeneric("hasIsotopes<-", function(object,value) {
    standardGeneric("hasIsotopes<-")
})
