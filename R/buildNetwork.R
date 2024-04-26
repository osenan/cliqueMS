#################################################
#  FUNCTIONS TO FILTER FEATURES IN THE PEAKLIST  #
#       AND BUILD THE NETWORK OF SIMILARITY      #
##################################################

nato0 <- function(mat) {
    newmat = mat
    newmat[is.na(newmat)] = 0
    return(newmat)
}

similarFeatures <- function(cosine, peaklist, mzerror = 0.000005,
    rtdiff = 0.0001, intdiff = 0.0001) {
    ## identify peaks with very similar cosine correlation,
    ## m/z, rt and intensity
    network <- igraph::graph.adjacency(cosine, weighted = TRUE,
        diag = FALSE, mode = "undirected")
    ## identify edges with weight almost 1
    edges0.99 <- igraph::get.edges(
        network,igraph::E(network)[igraph::E(network)$weight > 0.99]) 
    if( nrow(edges0.99) > 0) {
    ## now check if this features have similar values of m/z,
    ## retention time and intensity, if this is true,
    ## filter the repeated feature
        repeated.peaks <- vapply(seq_len(nrow(edges0.99)), function(x) {
            rows <- peaklist[as.numeric(edges0.99[x,]),
                c("mz","rt","maxo")]
            error <- abs(rows[1,] - rows[2,])/rows[1,]
            res <- sum( c(error["mz"] <= mzerror,
                error["rt"] <= rtdiff,
                error["maxo"] <= intdiff) ) == 3
        }, logical(1))
        if( sum(repeated.peaks) == 0 ) {
            nodes.delete = NULL
        } else {
            filtered.edges <- edges0.99[repeated.peaks,]
            if(sum(repeated.peaks) == 1) {
                ## only one peak filtered
                nodes.delete <- min(filtered.edges)
            } else {
                nodes.delete <- matrixStats::rowMins(filtered.edges)
            }
        }
    } else {
        nodes.delete = NULL
    }
    return(nodes.delete)
}

filterFeatures <- function(cosinus, peaklist, mzerror = 5e-6 ,
    rtdiff = 1e-4, intdiff = 1e-4 ) {
    ## function to filter artifacts from signal processing before network
    newpeaklist <- peaklist
    newcosinus <- cosinus
    deleteN <- similarFeatures(cosinus, peaklist, mzerror = mzerror,
        rtdiff = rtdiff, intdiff = intdiff)
    if( !is.null(deleteN) ) {
        newpeaklist <- newpeaklist[-1*deleteN,]
        newcosinus <- newcosinus[-1*deleteN, -1*deleteN]
    }
    return(list(cosTotal = newcosinus, peaklist = newpeaklist,
        deleted = deleteN))
}

defineEIC <- function(xdata) {
    mzs.xdata <- xcms::mz(xdata)
    rts.xdata <- xcms::rtime(xdata)
    its.xdata <- xcms::intensity(xdata)
    peaks <- xcms::chromPeaks(xdata)
    EIC <- matrix(data = 0, nrow = nrow(peaks),
        ncol = length(rts.xdata))
    for( i in seq_len(nrow(peaks)) ){
        peak <- peaks[i,]
        posrtmin <- which(rts.xdata == peak["rtmin"])
        posrtmax <- which(rts.xdata == peak["rtmax"])
        rangepeak <- seq(posrtmin, posrtmax, by = 1)
        peakint <- unlist(lapply(rangepeak,function(y) {
            mzposc <- which(mzs.xdata[[y]] >= peak["mzmin"])
            finalpos <- mzposc[which(mzs.xdata[[y]][mzposc] <= peak["mzmax"])]
            if(length(finalpos) == 0) {
                int <- 0
            } else {
                int <- mean(its.xdata[[y]][finalpos])
            }
            int
        }))
        EIC[i,rangepeak] <- peakint
    }
    return(EIC)
}

#' @export
#' @describeIn createNetwork To use with 'xcmsSet' class
setMethod("createNetwork", "xcmsSet", function(mzdata, peaklist,
    filter = TRUE, mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    if (!requireNamespace("CAMERA", quietly = TRUE)) {
        stop("Package CAMERA needed for 'xcmsSet' processed data. Please use
            'XCMSnExp' objects or install package CAMERA.",
        call. = FALSE)
    }
    xsan <- CAMERA::xsAnnotate(mzdata)
    EIC <- CAMERA::getAllPeakEICs(xsan, rep(1,nrow(peaklist)))
    eicmat <- EIC$EIC
    eicmatnoNA <- nato0(eicmat)
    sparseeic <- slam::as.simple_triplet_matrix(t(eicmatnoNA))
    cosTotal <- coop::cosine(sparseeic) # compute cosine corr
    cosTotal[is.nan(cosTotal)] <- 0
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
    nozeroEdges <- igraph::E(network)[which(igraph::E(network)$weight != 0)]
    network <- igraph::subgraph.edges(network, nozeroEdges)
    igraph::E(network)$weight <- round(igraph::E(network)$weight,
        digits = 10)
    ## change similarity of 1 to 0.99999999 to non avoid 'nan'
    igraph::E(network)$weight[
        which(igraph::E(network)$weight == 1)] <- 0.99999999999
    return(list(network = network, peaklist = peaklist))
})

#' @export
#' @describeIn createNetwork To use with 'XCMSnExp' class
setMethod("createNetwork", "XCMSnExp", function(mzdata, peaklist,
   filter = TRUE, mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {

    eicmat <- defineEIC(mzdata)
    sparseeic <- slam::as.simple_triplet_matrix(t(eicmat))
    cosTotal <- coop::cosine(sparseeic) # compute cosine cor
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
})
