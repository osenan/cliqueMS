##################################################
#  FUNCTIONS TO FILTER FEATURES IN THE PEAKLIST  #
#       AND BUILD THE NETWORK OF SIMILARITY      #
##################################################

similarFeatures <- function(cosine, peaklist, mzerror = 0.000005,
                            rtdiff = 0.0001, intdiff = 0.0001) {
    # identify peaks with very similar cosine correlation, m/z, rt and intensity
    network <- igraph::graph.adjacency(cosine, weighted = T,
                                       diag = F, mode = "undirected")
    # identify edges with weight almost 1
    edges0.99 <- igraph::get.edges(
        network,igraph::E(network)[igraph::E(network)$weight > 0.99]) 
    if( nrow(edges0.99) > 0) {
   # now check if this features have similar values of m/z,
   # retention time and intensity, if this is true, filter the repeated feature
    repeated.peaks <- sapply(1:nrow(edges0.99), function(x) {
        rows <- peaklist[as.numeric(edges0.99[x,]),
                         c("mz","rt","maxo")]
        error <- abs(rows[1,] - rows[2,])/rows[1,]
        res <- sum( c(error["mz"] <= mzerror,
                      error["rt"] <= rtdiff,
                      error["maxo"] <= intdiff) ) == 3
    })
    if( sum(repeated.peaks) == 0 ) {
        nodes.delete = NULL } else {
                                filtered.edges <- 
                                    edges0.99[repeated.peaks,]
                                if(sum(repeated.peaks) == 1) {
                                    # only one peak filtered
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

getProfileMatrix <- function(msSet, peaklist) {
    # function to get the profiles of all features of the m/z data
    xraw = xcms::xcmsRaw(msSet@filepaths, profstep = 0)
    massrange = as.matrix(peaklist[,c("mzmin","mzmax")])
    timerange = as.matrix(peaklist[,c("rtmin", "rtmax")])
    profiles = xcms::getEIC(xraw, mzrange = massrange,
                      rtrange = timerange)
    positions = lapply(1:nrow(peaklist), function(x) {
        rtmin = peaklist[x,"rtmin"]
        rtmax = peaklist[x,"rtmax"]
        posmin = which(xraw@scantime == rtmin)
        posmax = which(xraw@scantime == rtmax)
        allpos = posmin:posmax
        res = data.frame(i = allpos,
                         j = rep(x, length(allpos)),
                         x = profiles@eic$xcmsRaw[[x]][,"intensity"]
                         )
    })
    dataRaw = do.call(rbind, positions)
    eicmat = Matrix::sparseMatrix(i = dataRaw$i,
                                  j = dataRaw$j,
                                  x = dataRaw$x,
                                  dims = c(length(xraw@scantime),
                                           nrow(peaklist))
                                  )
    eicmat = eicmat
    return(eicmat)
}


#' @title Function to create a similarity network from processed m/z data
#'
#' @description
#' This function creates a similarity network with nodes as features
#' and weighted edges as the cosine similarity between those nodes.
#' Edges with weights = 0 are not included in the network. Nodes
#' without edges are not included in the network. This network will
#' be used to define clique groups and find annotation within this
#' groups
#'
#' @details Signal processing algorithms may output artefact features.
#' Sometimes they produce two artefact features which are almost identical
#' This artefacts may lead to errors in the computation of the clique
#' groups, so it is recommended to set 'filter' = TRUE to drop repeated
#' features.
#' @param msSet A 'xcmsSet' object with processed m/z data
#' @param peaklist Is a data.frame feature info for m/z data.
#' put each feature in a rom and a column 'mz' for mass data, 
#' retention time column 'rt' and intensity in column 'maxo'
#' @param filter If TRUE, filter out very similar features
#' that have a correlation similarity > 0.99 and equal values of m/z,
#' retention time and intensity
#' @param mzerror Relative error for m/z, if relative error 
#' between two features is below that value that features
#' are considered with similar m/z value
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
#' library(cliqueMS)
#' netlist = createNetwork(exmsSet, exmsSet@peaks, filter = TRUE)
#' @seealso \code{\link{getCliques}}
createNetwork <- function(msSet, peaklist, filter = TRUE, mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    #function to create similarity network from processed ms data
    # it filters peaks with very high similarity (0.99 >), m/z, intensity and retention time
    # get profile matrix from m/z data
    if(class(msSet) != "xcmsSet") stop("msSet should be of class xcmsSet")
    eicmat <- getProfileMatrix(msSet, peaklist)
    cosTotal <- qlcMatrix::cosSparse(eicmat) # compute cosine corr
    if(filter == T) {
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
    network <- igraph::graph.adjacency(cosTotal, weighted = TRUE,
                                       diag = FALSE, mode = "undirected")
    # remove edges that are zero
    nozeroEdges = igraph::E(network)[which(igraph::E(network)$weight != 0)]
    network <- subgraph.edges(network, nozeroEdges)
    return(list(network = network, peaklist = peaklist))
}
