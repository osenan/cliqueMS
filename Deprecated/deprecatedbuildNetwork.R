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
#' @param mzData An object with processed m/z data, see 'methods' for
#' valid class types.
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
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' library(xcms)
#' rawMS <- readMSData(files = mzfile, mode = "onDisk")
#' cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' mzData <- findChromPeaks(rawMS, cpw)
#' peaklist = as.data.frame(chromPeaks(mzData))
#' netlist = createNetwork(mzData, peaklist, filter = TRUE)
#' @seealso \code{\link{getCliques}}
createNetwork <- function(mzData, peaklist, filter = TRUE,
    mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    UseMethod("createNetwork") }

#' @export
#' @title Function to create a similarity network from 'xcmsSet'
#' processed m/z data
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
#' CAMERA package has to be installed to use this method.
#' @param mzData A 'xcmsSet' object with processed m/z data
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
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' msSet <- xcms::xcmsSet(files = mzfile, method = "centWave",
#' ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' netlist = createNetwork(msSet, peaks(msSet), filter = TRUE)
#' @seealso \code{\link{getCliques}}
createNetwork.xcmsSet <- function(mzData, peaklist, filter = TRUE,
    mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    ## function to create similarity network from processed ms data
    ## it filters peaks with very high similarity (0.99 >), m/z,
    ## intensity and retention time
    ## get profile matrix from m/z data
    if (!requireNamespace("CAMERA", quietly = TRUE)) {
        stop("Package CAMERA needed for 'xcmsSet' processed data. Please use
            'XCMSnExp' objects or install package CAMERA.",
        call. = FALSE)
    }
    if(is(mzData,"xcmsSet") == FALSE) stop("mzData should be of class xcmsSet")
    xsan <- CAMERA::xsAnnotate(mzData)
    EIC <- CAMERA::getAllPeakEICs(xsan, rep(1,nrow(peaklist)))
    eicmat <- EIC$EIC
    eicmatnoNA <- nato0(eicmat)
    sparseeic <- as(t(eicmatnoNA), "sparseMatrix")
    cosTotal <- qlcMatrix::cosSparse(sparseeic) # compute cosine corr
    if(filter == TRUE) {
        filterOut <- filterFeatures(cosTotal, peaklist,
            mzerror = mzerror,
            rtdiff = rtdiff,
            intdiff = intdiff)
        cosTotal <- filterOut$cosTotal
        peaklist <- filterOut$peaklist
        message(paste("Features filtered:",
                length(filterOut$deleted),
                sep = " "))
    }
    network <- igraph::graph.adjacency(cosTotal, weighted = TRUE,
        diag = FALSE, mode = "undirected")
    igraph::V(network)$id = seq_len(nrow(peaklist))
    ## remove edges that are zero
    nozeroEdges = igraph::E(network)[which(igraph::E(network)$weight != 0)]
    network <- igraph::subgraph.edges(network, nozeroEdges)
    igraph::E(network)$weight <- round(igraph::E(network)$weight,
        digits = 10)
    ## change similarity of 1 to 0.99999999 to non avoid 'nan'
    igraph::E(network)$weight[
        which(igraph::E(network)$weight == 1)] <- 0.99999999999
    return(list(network = network, peaklist = peaklist))
}

#' @export
#' @title Function to create a similarity network from
#' 'XCMSnExp' processed m/z data
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
#' @param mzData A 'XCMSnExp' object with processed m/z data
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
#' require(xcms)
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' rawMS <- readMSData(files = mzfile, mode = "onDisk")
#' cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' mzData <- findChromPeaks(rawMS, cpw)
#' peaklist = as.data.frame(chromPeaks(mzData))
#' netlist = createNetwork(mzData, peaklist, filter = TRUE)
#' @seealso \code{\link{getCliques}}
createNetwork.XCMSnExp <- function(mzData, peaklist, filter = TRUE,
    mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    ## function to create similarity network from processed ms data
    ## it filters peaks with very high similarity (0.99 >),
    ## m/z, intensity and retention time.
    ## get profile matrix from m/z data
    if(is(mzData,"XCMSnExp") == FALSE) {
        stop("mzData should be of class XCMSnExp")
    }
    eicmat <- defineEIC(mzData)
    sparseeic <- as(t(eicmat), "sparseMatrix")
    cosTotal <- qlcMatrix::cosSparse(sparseeic) # compute cosine corr
    if(filter == TRUE) {
        filterOut <- filterFeatures(cosTotal, peaklist,
            mzerror = mzerror,
            rtdiff = rtdiff,
            intdiff = intdiff)
        cosTotal <- filterOut$cosTotal
        peaklist <- filterOut$peaklist
        message(paste("Features filtered:",
                length(filterOut$deleted),
                sep = " "))
    }
    network <- igraph::graph.adjacency(cosTotal, weighted = TRUE,
        diag = FALSE, mode = "undirected")
    igraph::V(network)$id = seq_len(nrow(peaklist))
    ## remove edges that are zero
    nozeroEdges = igraph::E(network)[which(igraph::E(network)$weight != 0)]
    network <- igraph::subgraph.edges(network, nozeroEdges)
    igraph::E(network)$weight <- round(igraph::E(network)$weight,
        digits = 10)
    ## change similarity of 1 to 0.99999999 to non avoid 'nan'
    igraph::E(network)$weight[
        which(igraph::E(network)$weight == 1)] <- 0.99999999999
    return(list(network = network, peaklist = peaklist))
}
