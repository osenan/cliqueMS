filterCharge <- function(allnodes, df) {
    ## function to filter the isotopes that have inconsistencies
    ## in charge
    badisoCharge <- unlist(lapply(allnodes, function(x) {
        res <- integer()
        posp <- which(df$pfeature == x)
        if( length(posp) > 0 ) {
            posi <- which(df$ifeature == x)
            if( length(posi) > 0 ) {
                if(df[posp,"pcharge"] != df[posi,"icharge"]) {
                    res <- c(posp, posi)
                }
            }
        }
        res
    }))
    return(badisoCharge)
}

filterInlinks <- function(inlinks, df) {
    ## function to drop one parental mass when one isotope
    ## has two parental masses candidates
    badpfeatures <- unlist(lapply(inlinks, function(x) {
        rowpfeatures <- which(df$ifeature == x)
        ## drop the parental feature with less weight
        dropRows <- rowpfeatures[-1*which.max(
            df[rowpfeatures,"weight"])]
    }))
    return(badpfeatures)
}

filterOutlinks <- function(outlinks, df) {
    ## function to drop one isotope when two isotopes
    ## point to the same parental masss
    badifeatures <- unlist(lapply(outlinks, function(x) {
        rowifeatures <- which(df$pfeature == x)
        ## drop the parental feature with less weight
        dropRows <- rowifeatures[-1*which.max(
        df[rowifeatures,"weight"])]
    }))
    return(badifeatures)
}

filterIso <- function(isodf, network) {
    ## Function to filter isotopes data.frame,
    ## and to create a network of isotopes
    netpfeature = vapply(isodf$pfeature, function(x) {
        which(igraph::V(network)$id == x)
    }, numeric(1))
    netifeature = vapply(isodf$ifeature, function(x) {
        which(igraph::V(network)$id == x)
    }, numeric(1))
    isodf$pfeature = netpfeature
    isodf$ifeature = netifeature
    isodfSorted <- isodf[,c("pfeature", "ifeature")]
    isodfSorted <- do.call(rbind,lapply(seq_len(nrow(isodfSorted)),
        function(x) {
            sort(as.matrix(isodfSorted[x,])) }
    ))
    isodf$weight <- igraph::E(network,
        P = as.numeric(t(isodfSorted)))$weight
    ## as.numeric is important to acces correct weight values
    ## First filter isotopes pointing to two different parents
    inlinks <- as.numeric(names(
        which(table(isodf[,"ifeature"]) > 1)))
    badpfeatures <- filterInlinks(inlinks, isodf)
    if( length(badpfeatures) > 0 ) {
        isodf <- isodf[-1*badpfeatures,]
    }
    ## Second filter parents pointed by two diferent isotopes
    outlinks <- as.numeric(names(
        which(table(isodf[,"pfeature"]) > 1)))
    badifeatures <- filterOutlinks(outlinks, isodf)
    if( length(badifeatures) > 0 ) {
        isodf <- isodf[-1*badifeatures,]
    }
    ## Third filter inconsistency in charge
    allnodes <- unique(c(isodf[,"pfeature",],isodf[,"ifeature"]))
    badisoCharge <- filterCharge(allnodes, isodf)
    if( length(badisoCharge) > 0 ) {
        isodf <- isodf[-1*badisoCharge,]
    }
    realpfeature <- vapply(isodf$pfeature, function(x) {
        igraph::V(network)[x]$id }, numeric(1))
    realifeature <- vapply(isodf$ifeature, function(x) {
        igraph::V(network)[x]$id }, numeric(1))
    isodf$pfeature <- realpfeature
    isodf$ifeature <- realifeature
    ## Finally create the filtered isotope network
    isonet <- igraph::graph.data.frame(isodf[,c("ifeature","pfeature")])
    return(list(network = isonet, isodf = isodf))
}

isoGrade <- function(isonet) {
    ## Function to grade and isotope, starting from 0 to the parental isotpe,
    ## 1 the first isotope and further
    grades <- vapply(igraph::V(isonet), function(x) {
        res <- 0
        nei <- igraph::neighbors(isonet, v = x, mode = "out")
        while(length(nei) > 0) {
            res <- res + 1
            nei <- igraph::neighbors(isonet, v = nei, mode = "out")
        }
        res
    }, numeric(1))
    return(grades)
}

correctGrade <- function(isoTable, maxGrade) {
    maxCluster <- max(isoTable$cluster)
    clusters <- unique(isoTable$cluster)
    res <- do.call(rbind,lapply(seq_len(length(clusters)), function(x) {
        cluster <- clusters[x]
        cpos <- isoTable[isoTable$cluster == x,]
        goodP <- cpos[cpos$grade <= maxGrade,]
        badP <- cpos[cpos$grade > maxGrade,]
        if(nrow(badP) > 1) {
            maxrowN <- max(as.numeric(rownames(badP)))
            badP$cluster = maxCluster + maxrowN + 1
            badP$grade = 0:(nrow(badP)-1)
            goodP = rbind(badP, goodP)
        }
        goodP
    }))
    return(res)
}
            

isonetAttributes <- function(isolist, maxGrade) {
    ## Function to set the node attributes for each isotope:
    ## grade, charge, and community
    ## First assign grade (for info look isoGrade function)
    igraph::V(isolist$network)$grade <- isoGrade(isolist$network)
    charge <- vapply(
        as.numeric(igraph::V(isolist$network)$name), function(x) {
            posp <- which(isolist$isodf[,"pfeature"] == x)
            if( length(posp) > 0 ) {
                res <- isolist$isodf[posp,"pcharge"]
            } else {
                posi <- which(isolist$isodf[,"ifeature"] == x)
                res <- isolist$isodf[posi,"icharge"]
            }
            res
        }, numeric(1))
    ## Second assign charge of each feature as isotope
    igraph::V(isolist$network)$charge <- charge
    ## Third label features that belong to the same isotope cluster
    igraph::V(isolist$network)$cluster <- 
        igraph::clusters(isolist$network, "weak")$membership
    ## Final step write a table with all the isotope data
    isoTable <- data.frame(
        feature = as.numeric(igraph::V(isolist$network)$name),
        charge = igraph::V(isolist$network)$charge,
        grade = igraph::V(isolist$network)$grade,
        cluster = igraph::V(isolist$network)$cluster
    )
    ## Fourth, correct the grade
    while(max(isoTable$grade) > maxGrade) {
        isoTable <- correctGrade(isoTable, maxGrade)
    }
    return(isoTable)
}

addIso2peaklist <- function(isoTable, peaklist) {
    peaklist$isotope <- rep("M0", nrow(peaklist))
    peaklist$isotope[isoTable[,"feature"]] <- 
        paste(paste("M",isoTable$grade, sep = ""),
            paste("[",isoTable$cluster, "]",sep = ""),
                sep = " ")
    return(peaklist)
}

computelistofIsoTable <- function(anclique, maxCharge, maxGrade, ppm, isom ) {
    listofisoTable <- lapply(anclique@cliques, function(x) {
        df.clique <- as.data.frame(
            cbind(anclique@peaklist[x, c("mz","maxo")],x)
        )
        colnames(df.clique) <- c("mz","maxo","feature")
        ## Sort df.clique by intensity because isotopes are less
        ## intense than their parental features
        df.clique <- df.clique[order(df.clique$maxo, decreasing = TRUE),]
        ## compute isotopes from clique
        isodf <- returnIsotopes(df.clique, maxCharge = maxCharge,
            ppm = ppm, isom = isom)
        if( nrow(isodf) > 0 ) {
            ## filter the isotope list by charge
            ## and other inconsistencies
            isolist <- filterIso(isodf, anclique@network)
            if( nrow(isolist$isodf) > 0 ) {
                ## write a table with feature, charge, grade and cluster 
                iTable <- isonetAttributes(isolist, maxGrade)
            } else {
                iTable = NULL}
        } else {
            iTable = NULL
        }
        iTable
    })
    return(listofisoTable)
}

#' @export
#' @title Annotate isotopes
#'
#' @description This function annotates features that are carbon
#' isotopes based on m/z and intensity data. The monoisotopic
#' mass has to be more intense than the first isotope, the first
#' isotope more intense than the second isotope and so one so forth.
#' Isotopes are annotated within each clique group.
#' @param anclique An 'anClique' object with clique groups computed
#' @param maxCharge Maximum charge considered when we test two
#' features to see whether they are isotopes
#' @param maxGrade The maximum number of isotopes apart from the
#' monoisotopic mass. A 'maxGrade' = 2 means than we have the
#' monoisotopic mass, first isotope and second isotope
#' @param isom The mass difference of the isotope
#' @param ppm Relative error in ppm to consider that two features
#' have the mass difference of an isotope
#' @return It returns an 'anClique' object with isotope annotation.
#' it adds the column 'isotope' to the peaklist in the anClique object
#' @examples
#' data(ex.cliqueGroups)
#' show(ex.cliqueGroups)
#' ex.isoAn <- getIsotopes(ex.cliqueGroups)
#' show(ex.isoAn)
#' @seealso
#' \code{\link{getCliques}}
getIsotopes <- function(anclique, maxCharge = 3,
    maxGrade = 2, ppm = 10, isom = 1.003355) {
    # Function to get all the isotopes from the m/z data
    # after splitting it into clique groups
    if(anclique@isoFound == TRUE) {
        warning("Isotopes have been already computed for this object")
    }
    if(anclique@cliquesFound == FALSE) {
        warning("Cliques have not been computed for this object.
        This could lead to long computing times
        for isotope annotation")
    }
    message("Computing isotopes")
    listofisoTable <- computelistofIsoTable(anclique, maxCharge,
        maxGrade, ppm, isom )
    ## If there are no isotopes in all dataset
    if( length(listofisoTable) ==
        sum(vapply(listofisoTable, is.null, logical(1))) ) {
            isoTable <- matrix(c(NA,NA,NA,NA), nrow = 1)
            colnames(isoTable) <- c("feature","charge","grade","cluster")
            anclique@peaklist$isotope <- rep("M0", nrow(anclique@peaklist))
    } else {
    ## The cluster label is inconsistent between all isotopes found
    ## let's correct for avoiding confusions
        listofisoTable <- listofisoTable[
            !vapply(listofisoTable, is.null, logical(1))]
        maxC <- max(listofisoTable[[1]]$cluster)
        for(i in 2:length(listofisoTable)) {
            listofisoTable[[i]]$cluster = listofisoTable[[i]]$cluster + maxC + 1
            maxC <- max(listofisoTable[[i]]$cluster)
        }
        isoTable <- do.call(rbind, listofisoTable)
        rownames(isoTable) <- seq_len(nrow(isoTable))
        ## Change the peaklist adding isotope column
        anclique@peaklist <- addIso2peaklist(isoTable, anclique@peaklist)
    }
    message("Updating anClique object")
    ## Now change status of isotopes at anclique object
    anclique@isoFound <- TRUE
    ## Put new isotopes table
    anclique@isotopes <- isoTable
    return(anclique)
}
