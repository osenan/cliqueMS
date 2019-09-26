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
    cliques[,"node"] <- igraph::V(anclique@network)$id
    maxCliqueLabel <- max(cliques[,"clique"]) + 1
    ## search for nodes without clique value
    missingNodes <- seq_len(nrow(anclique@peaklist))[
        is.na(match(seq_len(nrow(anclique@peaklist)),
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
computeCliques <- function(anclique, tol = 1e-5, silent = TRUE) {
    ## this function calls the C++ code and the probabilistic model
    ## to find the clique groups with max log-likelihhod.
    ## Then it assigns clique groups to features that have links and
    ## then assigns clique groups to the remaining features
    ## It is recommended first to apply filters to anClique peaklist
    if(anclique@cliquesFound == TRUE) {
        warning("cliques have already been computed\n")
    }
    # create network data.frame
    netdf <- cbind(igraph::as_edgelist(anclique@network),
        igraph::E(anclique@network)$weight)
    netdf <- data.frame(node1 = netdf[,1], node2 = netdf[,2],
        weight = netdf[,3])
    cliquesRaw <- returnCliques(netdf, tol, silent = silent) # computeCliques
    cliquesGood <- updateCliques(anclique, cliquesRaw)
    anclique@peaklist$cliqueGroup = cliquesGood
    anclique@cliquesFound  <-  TRUE
    cliquesReplace  <-  
        lapply(unique(anclique@peaklist$cliqueGroup), function(x) {
            which(anclique@peaklist$cliqueGroup == x)
        })
    names(cliquesReplace) <- unique(anclique@peaklist$cliqueGroup)
    anclique@cliques <- cliquesReplace
    return(anclique)
}

#' @export
getCliques <- function(mzdata, filter = TRUE, mzerror = 5e-6,
    intdiff = 1e-4, rtdiff = 1e-4, tol = 1e-5, silent = TRUE) {
    message("Creating anClique object")
    anclique <- createanClique(mzdata)
    message("Creating network")
    netlist <- createNetwork(mzdata, anclique@peaklist,
        filter = filter,
        mzerror = mzerror,
        intdiff = intdiff,
        rtdiff = rtdiff)
    anclique@peaklist <- netlist$peaklist
    anclique@network <- netlist$network
    message("Computing cliques")
    anclique = computeCliques(anclique, tol, silent)
    return(anclique)                    
}
