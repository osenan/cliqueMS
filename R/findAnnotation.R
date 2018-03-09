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
        
getAnnotation.anClique <- function(anclique, adinfo, topmassf = 1,
                          topmasstotal = 10, sizeanG = 20,
                          tol = 1e-5, filter = 1e-4, emptyS = 1e-6) {
    if(anclique$anFound == TRUE) {
        warning("Annotation has already been computed for this object")
    }
    cat("Computing annotation\n")
    # Compute annotation for each clique
    anList <- lapply(anclique$cliques, function(x) {
        df.clique <- as.data.frame(
            cbind(anclique$peaklist[x, c("mz","isotope")],x)
        )
        colnames(df.clique) <- c("mz","isotope","feature")
        df.clique <- df.clique[grep("^M0", df.clique$isotope),]
        df.clique <- df.clique[order(df.clique$mz, decreasing = F),]
        annotation <- returnAnnotation(df.clique, adinfo,
                                       topmassf, topmasstotal,
                                       sizeanG, tol,
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
