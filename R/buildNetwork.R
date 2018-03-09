##################################################
#        DEFINE CLASSES AND FUNCTIONS            #
##################################################

anClique <- structure(list("peaklist" = data.frame(),
                           "network" = igraph::empty_graph(),
                           "cliques" = list(),
                           "cliquesFound" = FALSE,
                           "isotopes" = list(),
                           "isoFound" = FALSE,
                           "anFound" = FALSE),
                      class = "anClique")

anClique <- function(msSet, filter = T, mzerror, rtdiff, intdiff) {
    if(class(msSet) != "xcmsSet") stop("msSet should be of class msSet")
    peaklist = as.data.frame(msSet@peaks)
    cliques = list()
    isotopes = matrix()
    return(structure(list("peaklist" = peaklist,
                          "network" = igraph::empty_graph(),
                          "cliques" = cliques,
                          "cliquesFound" = FALSE,
                          "isotopes" = list(),
                          "isoFound" = FALSE,
                          "anFound" = FALSE),
                     class = "anClique"))
}

summary.anClique <- function(object, ...)
{
    cat(paste("anClique object with",nrow(object$peaklist),
              "features\n"), sep = " ")
    if(object$cliquesFound) {
        cat(paste("Features have been splitted into",
                  length(object$cliques), "cliques\n", sep = " "))
    } else {
        cat("No cliques found\n")
    }
    if(object$isoFound) {
        cat(paste(nrow(object$isotopes), "Features are isotopes\n", sep = " "))
    } else {
        cat("No isotopes found\n")
    }
    if(object$anFound) {
        cat(paste(sum(!is.na(object$peaklist$an))), "features annotated\n", sep = " "))
    } else {
        cat("Features do not have annotation\n")
    }
}

nato0 <- function(mat) {
    # transform NA to 0 values in EIC matrix
      newmat <- mat
      newmat[is.na(newmat)] <- 0
      return(newmat)
}

similarFeatures <- function(cosine, peaklist, mzerror = 0.000005, rtdiff = 0.0001, intdiff = 0.0001) {
    # identify peaks with very similar cosine correlation, m/z, rt and intensity
    network <- igraph::graph.adjacency(cosine, weighted = T, diag = F, mode = "undirected")
    # identify edges with weight almost 1
    edges0.99 <- igraph::get.edges(
        network,igraph::E(network)[igraph::E(network)$weight > 0.99]) 
    if( nrow(edges0.99) > 0) {
        # now check if this features have similar values of m/z, retention time and intensity, if this is true, filter the repeated feature
    repeated.peaks <- sapply(1:nrow(edges0.99), function(x) {
        rows <- peaklist[as.numeric(edges0.99[x,]),
                         c("mz","rt","maxo")]
        error <- abs(rows[1,] - rows[2,])/rows[1,]
        res <- sum( c(error["mz"] <= mzerror,
                      error["rt"] <= rtdiff,
                      error["maxo"] <= intdiff) ) == 3
    })
    if( sum(repeated.peaks) == 0 ) {
        nodes.delete = NULL } else { filtered.edges <- 
                                         edges0.99[repeated.peaks,]
                                         if(sum(repeated.peaks) == 1) {
                                    #only one peak filtered
                                             nodes.delete <- 
                                                 min(filtered.edges) } else{
                                                                         nodes.delete <- sapply(1:nrow(filtered.edges),
                                                                                                function(x) { min(filtered.edges[x,])
                                                                                                })
                                                                     } 
                            }
    } else { nodes.delete = NULL }
    return(nodes.delete)
}

filterFeatures <- function(cosinus, peaklist, mzerror = 5e-6 , rtdiff = 1e-4, intdiff = 1e-4 ) {
    # function to filter artifacts from signal processing before network
    newpeaklist <- peaklist
    newcosinus <- cosinus
    deleteN <- similarFeatures(cosinus, peaklist, mzerror = mzerror,
                               rtdiff = rtdiff, intdiff = intdiff)
    if( !is.null(deleteN) ) {
        newpeaklist <- newpeaklist[-1*deleteN,]
        newcosinus <- newcosinus[-1*deleteN, -1*deleteN]
    }
    return(list(cosTotal = newcosinus, peaklist = newpeaklist, deleted = deleteN))
}


createNetwork.anClique <- function(msSet, peaklist, filter = T, mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    #function to create similarity network from processed ms data
    # it filters peaks with very high similarity (0.99 >), m/z, intensity and retention time
    if(filter == T) {
        xsCAMERA = CAMERA::xsAnnotate(msSet)
        EIC <- CAMERA::getAllPeakEICs(xsCAMERA, rep(1,nrow(peaklist)))
        EICmat <- EIC$EIC
        EICnoNA <- nato0(EICmat) # transform NA to 0 values
        sparseEIC <- as(t(EICnoNA), "sparseMatrix") # transform normal EIC matrix into a sparse matrix
        cosTotal <- qlcMatrix::cosSparse(sparseEIC) # compute cosine correlation
        filterOut <- filterFeatures(cosTotal, peaklist,
                                    mzerror = mzerror,
                                    rtdiff = rtdiff,
                                    intdiff = intdiff)
        cosTotal <- filterOut$cosTotal
        peaklist <- filterOut$peaklist
        cat(paste("Features filtered:",
                  length(filterOut$deleted),"\n",
                  sep = " "))
    }
    network <- igraph::graph.adjacency(cosTotal, weighted = T, diag = F, mode = "undirected")
    return(list(network = network, peaklist = peaklist))
}

updateCliques <- function(network, cliques) {
    #this function is to put different clique number to extra nodes that do not have links
    #IMPORTANT: network needs all nodes, even if they have no links
    # number of nodes without links, not in the clique data.frame
    extraclique <- length(igraph::V(network)) - nrow(cliques)
    cliques = cliques[order(cliques[,"node"]),]
    maxCliqueLabel <- max(cliques[,"clique"]) + 1
    # clique label that is not in the returned ones
    cliquesRes <- (maxCliqueLabel:(maxCliqueLabel +
                                   (length(igraph::V(network)) -1)))
    for(i in 1:nrow(cliques)) {
        node <- cliques[i,"node"]
        clique <- cliques[i,"clique"]
        cliquesRes[node] <- clique
    }
    return(cliquesRes)
}

computeCliques.anClique <- function(anclique, network) {
                                        #first apply filters to anClique peaklist (if they have been actived)
    if(anclique$cliquesFound == TRUE) {
        warning("cliques have already been computed\n")
    }
    # create network data.frame
    netdf <- cbind(igraph::as_edgelist(network),
                  igraph::E(network)$weight)
    netdf <- data.frame(node1 = netdf[,1], node2 = netdf[,2],
                       weight = netdf[,3])
    cliquesRaw <- returnCliques(netdf) # computeCliques
    cliquesGood <- updateCliques(network, cliquesRaw)
    anclique$peaklist$cliqueGroup = cliquesGood
    anclique$cliquesFound = TRUE
    anclique$cliques <- 
        lapply(unique(anclique$peaklist$cliqueGroup), function(x) {
            which(anclique$peaklist$cliqueGroup == x)
        })
    names(anclique$cliques) <- unique(anclique$peaklist$cliqueGroup)
    return(anclique)
}

getCliques <- function(msSet, filter = T, mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    cat("Creating anClique object\n")
    anclique <- anClique(msSet)
    cat("Creating network\n")
    netlist <- createNetwork.anClique(msSet, anclique$peaklist,
                                      filter = filter,
                                      mzerror = mzerror,
                                      intdiff = intdiff,
                                      rtdiff = rtdiff)
    anclique$peaklist <- netlist$peaklist
    anclique$network <- netlist$network
    cat("Computing Cliques\n")
    anclique = computeCliques.anClique(anclique, netlist$network)
    return(anclique)                    
}
    
