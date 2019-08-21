##################################################
#      FUNCTIONS TO COMPUTE CLIQUE GROUPS THAT   #
#      WILL BE USED TO FIND ISOTOPE AND ADDUCT   #
#                 ANNOTATION                     #
##################################################

updateCliques <- function(anclique, cliques) {
    ## this function is to assign a clique group value to nodes
    ## that do not have links, because they do not appear in
    ## the edgelist
    ## number of nodes without links, not in the clique data.frame
    cliques = cliques[order(cliques[,"node"]),]
    ## assign correct node value to "nodes"
    cliques[,"node"] <- igraph::V(anclique$network)$id
    maxCliqueLabel <- max(cliques[,"clique"]) + 1
    ## search for nodes without clique value
    missingNodes <- seq_len(nrow(anclique$peaklist))[
        is.na(match(seq_len(nrow(anclique$peaklist)),
            cliques[,"node"]))]
    if( length(missingNodes) > 0 ) {
        ## if nodes without clique value exist
        templateC <- rep(1, length(missingNodes))
        templateC[1] <- maxCliqueLabel
        ## assign clique values to this nodes
        newCliques <- cumsum(templateC)
        dfnewCliques <- data.frame(node = missingNodes,
            clique = newCliques)
        dfnewCliques <- as.matrix(dfnewCliques)
        cliques = rbind(cliques, dfnewCliques)
        ## join this nodes and order again the matrix
        cliques = cliques[order(cliques[,"node"]),]
    }
    ## clique label that is not in the returned ones
    return(cliques[,"clique"])
}

#' @export
#' @title Computes clique groups from a similarity network
#'
#' @description This function splits the features in the network
#' in clique groups. The cliques are fully connected components
#' that have high similarity for inner edges and low similarity
#' for edges outside the clique. This function finds the clique
#' groups that better fit this criteria, moving nodes to different
#' groups until we find the groups that have the best log-likelihood.
#' @param anclique This function uses S3 'anClique' object. Gives
#' warning if clique groups have already been computed.
#' @param tol Minimum relative increase in log-likelihood to do
#' a new round of log-likelihood maximisation.
#' @param silent If 'FALSE' print on the console the log-likelihood
#' maximization progress. Default is 'TRUE'.
#' @return It returns an 'anClique' object with the computed
#' clique groups. It adds the column 'cliqueGroup' to the
#' 'peaklist' in the 'anClique' object.
#' @examples
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' msSet <- xcms::xcmsSet(files = mzfile, method = "centWave",
#' ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' ex.anClique <- createanClique(msSet)
#' summary(ex.anClique)
#' netlist <- createNetwork(msSet, peaks(msSet), filter = TRUE)
#' @seealso \code{\link{getCliques}}
computeCliques <- function(anclique, tol = 1e-5, silent = TRUE) {
    ## this function calls the C++ code and the probabilistic model
    ## to find the clique groups with max log-likelihhod.
    ## Then it assigns clique groups to features that have links and
    ## then assigns clique groups to the remaining features
    ## It is recommended first to apply filters to anClique peaklist
    if(anclique$cliquesFound == TRUE) {
        warning("cliques have already been computed\n")
    }
    # create network data.frame
    netdf <- cbind(igraph::as_edgelist(anclique$network),
        igraph::E(anclique$network)$weight)
    netdf <- data.frame(node1 = netdf[,1], node2 = netdf[,2],
        weight = netdf[,3])
    cliquesRaw <- returnCliques(netdf, tol, silent = silent) # computeCliques
    cliquesGood <- updateCliques(anclique, cliquesRaw)
    anclique$peaklist$cliqueGroup = cliquesGood
    anclique$cliquesFound = TRUE
    anclique$cliques <- 
        lapply(unique(anclique$peaklist$cliqueGroup), function(x) {
            which(anclique$peaklist$cliqueGroup == x)
        })
    names(anclique$cliques) <- unique(anclique$peaklist$cliqueGroup)
    return(anclique)
}

#' @export
#' @title Compute clique groups from processed m/z data
#'
#' @description This function splits features in groups to
#' find isotope and adduct annotation within each group. To
#' find them it uses a similarity network.
#' This similarity network has nodes as features and weighted
#' edges as the cosine similarity between features. Once
#' the network is obtained we find clique groups in this network.
#' The clique groups are fully connected components with
#' high similarity in inner edges and lower similarity in
#' edges outside the clique. We move nodes to different groups
#' until we find the groups with the maximum log-likelihood.
#'
#' @details Signal processing algorithms may output artefact features.
#' Sometimes they produce two artefact features which are almost identical
#' This artefacts may lead to errors in the computation of the clique
#' groups, so it is recommended to set 'filter' = TRUE to drop repeated.
#' features.
#' @param mzData An 'object with processed m/z data. Currently supported
#' class types are 'xcmsSet' or 'XCMSnExp.
#' @param filter If TRUE, filter out very similar features
#' that have a correlation similarity > 0.99 and equal values of m/z,
#' retention time and intensity.
#' @param mzerror Relative error for m/z, if relative error 
#' between two features is below that value that features
#' are considered with similar m/z value.
#' @param rtdiff Relative error for retention time, if 
#' relative error between two features is below that value
#' that features are considered with similar retention time.
#' @param intdiff Relative error for intensity, if relative
#' error between two features is below that value that
#' features are considered with similar intensity.
#' @param tol Minimum relative increase in log-likelihood to do
#' a new round of log-likelihood maximisation.
#' @param silent If 'FALSE' print on the console the log-likelihood
#' maximization progress. Default is 'TRUE'.
#' @return It returns an 'anClique' object with the computed
#' clique groups. It adds the column 'cliqueGroup' to the
#' 'peaklist' in the 'anClique' object.
#' @examples
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' msSet <- xcms::xcmsSet(files = mzfile, method = "centWave",
#' ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' ex.cliqueGroups <- getCliques(msSet)
#' @seealso \code{\link{computeCliques}}
#' \code{\link{createNetwork}}
#' \code{\link{anClique}}
getCliques <- function(mzData, filter = TRUE, mzerror = 5e-6,
    intdiff = 1e-4, rtdiff = 1e-4, tol = 1e-5, silent = TRUE) {
    cat("Creating anClique object\n")
    anclique <- createanClique(mzData)
    cat("Creating network\n")
    netlist <- createNetwork(mzData, anclique$peaklist,
        filter = filter,
        mzerror = mzerror,
        intdiff = intdiff,
        rtdiff = rtdiff)
    anclique$peaklist <- netlist$peaklist
    anclique$network <- netlist$network
    cat("Computing cliques\n")
    anclique = computeCliques(anclique, tol, silent)
    return(anclique)                    
}
