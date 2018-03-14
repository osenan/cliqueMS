isotopeAnnotation <- function(df.annotation, anclique) {
    # Isotopes of grade > 0 are excluded from annotation
    # This function adds the annotation, in case it exists,
    # to all isotopes of grade > 0 as they should have the same
    # annotation than its isotope of grade 0
    isofeatures <- anclique$isotopes$feature[
        which(anclique$isotopes$grade == 0)]
    extraAnList <- lapply(isofeatures, function(x) {
        cluster <- anclique$isotopes$cluster[
            which(anclique$isotopes$feature == x)]
        extraf <- anclique$isotopes$feature[
            which(anclique$isotopes$cluster == cluster)]
        #now drop the isotope of grade0 because it already has annotation
        extraf <- extraf[-1*which(extraf == x)]
        an <- df.annotation[which(df.annotation$feature == x),]
        # annotation for all isotopes should be the same
        extraan <- lapply(extraf, function(x) {
            an$feature = x
            an
        })
        extraan <- do.call(rbind, extraan)
    })
    extraAn <- do.call(rbind, extraAnList)
    return (extraAn)
}

#' @title Annotate adducts and fragments
#'
#' @description This function annotates adducts after
#' isotope annotation. For each clique group,
#' it searches for combinations of two or more features 
#' compatible with the same neutral mass and two or more
#' adducts in 'adinfo' list. For clique groups than have more
#' than one annotation solution, it scores all
#' possibilites and returns the top five solutions.
#' @param anclique Object of class 'anClique' with isotope
#' annotation
#' @param adinfo Ordered data.frame with columns 'adduct' with
#' adduct name, column 'log10freq' with the log10 frequency of each
#' adduct in the list, column 'massdiff' with the adduct mass
#' diff, column 'nummol' with the number of metabolite's
#' molecule necessary for that adduct and column 'charge' with the
#' charge of that adduct.
#' @param topmassf For each feature, we select a number of 
#' 'topmassf' neutral masses that will be considered in the
#' final scoring of the annotation
#' @param topmasstotal For each group we also select a number
#' of 'topmasstotal' neutral masses that will be considered 
#' in the scoring of the annotation apart from the masses in
#' topmassf
#' @param sizeanG If a clique group has a number
#' of monoisotopic features bigger than 'sizeanG', we try
#' to divide the clique group in non-overlapping annotation groups.
#' Each of this subdivisions is annotated independently.
#' @param ppm Relative error in ppm in which we consider two or 
#' more features compatible with a neutral mass and two or more
#' adducts in 'adinfo'.
#' @param filter If two annotations have a relative 
#' mass difference smaller than 'filter' and the same features
#' and adducts, drop the neutral mass with less adducts
#' @param emptyS Score given to non annotated features
#' @details If clique groups have a lot of features, there are
#' many combinations of neutral masses and adducts. This could
#' lead to long running times to score the top annotations.
#' Parameters
#' 'topmassf', 'topmasstotal' and 'sizeanG' are relevant in
#' those cases to drop the less likely neutral masses to
#' speed up the time of computation and still obtain
#' the most plausible annotation. If the clique group is small
#' usually no neutral masses are discarded for the scoring.
#' @return An 'anClique object' with annotation columns added
#' to the peaklist
#' @examples
#' ex.cliqueGroups <- getCliques(exmsSet, filter = T)
#' summary(ex.cliqueGroups)
#' ex.isoAn <- getIsotopes(ex.cliqueGroups)
#' summary(ex.isoAn)
#' ex.adductAn <- getAnnotation(ex.isoAn, positive.adinfo)
#' @seealso \code{\link{getCliques}}
#' \code{\link{getIsotopes}}
getAnnotation <- function(anclique, adinfo, topmassf = 1,
                          topmasstotal = 10, sizeanG = 20,
                          ppm = 10, filter = 1e-4, emptyS = 1e-6) {
    if(anclique$anFound == TRUE) {
        warning("Annotation has already been computed for this object")
    }
    cat("Computing annotation\n")
    # Compute annotation for each clique
    ppm = 1e-6*ppm
    anList <- lapply(anclique$cliques, function(x) {
        df.clique <- as.data.frame(
            cbind(anclique$peaklist[x, c("mz","isotope")],x)
        )
        colnames(df.clique) <- c("mz","isotope","feature")
        df.clique <- df.clique[grep("^M0", df.clique$isotope),]
        df.clique <- df.clique[order(df.clique$mz, decreasing = F),]
        annotation <- returnAnnotation(df.clique, adinfo,
                                       topmassf, topmasstotal,
                                       sizeanG, ppm,
                                       filter, emptyS)
    })
    df.annotation <- do.call(rbind, anList)
    cat("Annotation computed, updating peaklist\n")
    # Now add the annotation for isotopes
    isoAn <- isotopeAnnotation(df.annotation, anclique)
    df.annotation <- rbind(df.annotation, isoAn)
    df.annotation <- df.annotation[order(df.annotation$feature),]
    df.annotation <- df.annotation[,-1]
    anclique$peaklist <- cbind(anclique$peaklist, df.annotation)
    # uptade object information
    anclique$anFound <- TRUE
    return(anclique)
}
