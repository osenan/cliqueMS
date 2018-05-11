##################################################
#     DEFINITION OF S3 METHODS AND S3 CLASSES    #
##################################################

#' @title 'anClique-class' for annotating isotopes and adducts
#'
#' @aliases anClique-class
#' @description
#' S3 Class \code{anClique-class} for annotating isotopes and adducts
#' in processed m/z data. Features are first
#' grouped based on a similarity network algorithm and then
#' annotation of isotopes and adducts is performed in each group.
#' @details anClique-class contains the following elements:
#' 
#' 'peaklist'
#' Is a data.frame with m/z, retention time
#' and intensity information for each feature. It also contains
#' adduct and isotope information if annotation has been performed.
#'
#' 'network'
#' Is an igraph undirected network of similarity 
#' used to compute groups of features before annotation.
#' 
#' 'cliques'
#' Is a list that contains the groups of features.
#' Each id corresponds to a row in the peaklist.
#'
#' 'isotopes'
#' Is a data.frame with the column 'feature' for feature id, 
#' column 'charge' for the charge, column 'grade' that starts with 0
#' and it is 1 for the first isotope, 2 for the second and so on and
#' column 'cluster' which labels each group of features that are
#' isotopes.
#'
#' 'cliquesFound' is
#' TRUE if clique groups have been computed,
#' 'isoFound' is
#' TRUE if isotopes have been annotated,
#' 'anFound' is
#' TRUE if annotation of adducts have been computed.
#' @seealso \code{\link{anClique}}
anClique <- structure(list("peaklist" = data.frame(),
                           "network" = igraph::empty_graph(),
                           "cliques" = list(),
                           "cliquesFound" = FALSE,
                           "isotopes" = data.frame,
                           "isoFound" = FALSE,
                           "anFound" = FALSE),
                      class = "anClique")

#' @title 'anClique' produces an object of class 'anClique'.
#'
#' @description
#' \code{anClique} creates an 'anClique' object from processed m/z data.
#' @param msSet A 'xcmsSet' object with processed m/z data.
#' @examples
#' \donttest{
#' library(cliqueMS)
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' msSet <- xcms::xcmsSet(files = mzfile, method = "centWave",
#' ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' ex.anClique <- anClique(msSet = msSet)
#' summary(ex.anClique)
#' }
#' @seealso \code{\link{anClique-class}}
anClique <- function(msSet) {
    if(class(msSet) != "xcmsSet") stop("msSet should be of class xcmsSet")
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
        cat("No computed clique groups\n")
    }
    if(object$isoFound) {
        cat(paste(nrow(object$isotopes), "Features are isotopes\n", sep = " "))
    } else {
        cat("No isotope annotation\n")
    }
    if(object$anFound) {
        pos1 = which(!is.na(object$peaklist$an1))
        pos2 = which(!is.na(object$peaklist$an2))
        pos3 = which(!is.na(object$peaklist$an3))
        pos4 = which(!is.na(object$peaklist$an4))
        pos5 = which(!is.na(object$peaklist$an5))
        anFeatures = unique(c(pos1,pos2,pos3,pos4,pos5))
        cat(paste(length(anFeatures), "features annotated\n", sep = " "))
    } else {
        cat("No adduct annotation\n")
    }
}
