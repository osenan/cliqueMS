computelistofIsoTable <- function(anclique, maxCharge, maxGrade, ppm, isom ) {
    listofisoTable <- lapply(anclique$cliques, function(x) {
        df.clique <- as.data.frame(
            cbind(anclique$peaklist[x, c("mz","maxo")],x)
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
            isolist <- filterIso(isodf, anclique$network)
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
#' summary(ex.cliqueGroups)
#' ex.isoAn <- getIsotopes(ex.cliqueGroups)
#' summary(ex.isoAn)
#' @seealso
#' \code{\link{getCliques}}
getIsotopes <- function(anclique, maxCharge = 3,
    maxGrade = 2, ppm = 10, isom = 1.003355) {
    # Function to get all the isotopes from the m/z data
    # after splitting it into clique groups
    if(anclique$isoFound == TRUE) {
        warning("Isotopes have been already computed for this object")
    }
    if(anclique$cliquesFound == FALSE) {
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
            anclique$peaklist$isotope <- rep("M0", nrow(anclique$peaklist))
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
            anclique$peaklist <- addIso2peaklist(isoTable, anclique$peaklist)
    }
    message("Updating anClique object")
    ## Now change status of isotopes at anclique object
    anclique$isoFound <- TRUE
    ## Put new isotopes table
    anclique$isotopes <- isoTable
    return(anclique)
}
