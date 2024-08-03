##################################################
#     DEFINITION OF S3 METHODS AND S3 CLASSES    #
##################################################

setOldClass("igraph")

#' @title 'anClique' S4 class for annotating isotopes and adducts
#' @aliases anClique-class
#' @description
#' S4 Class \code{anClique-class} for annotating isotopes and adducts
#' in processed m/z data. Features are first
#' grouped based on a similarity network algorithm and then
#' annotation of isotopes and adducts is performed in each group.
#' The class contains the following slots.
#' @slot 'peaklist'
#' Is a data.frame with m/z, retention time
#' and intensity information for each feature. It also contains
#' adduct and isotope information if annotation has been performed.
#'
#' @slot 'network'
#' Is an igraph undirected network of similarity 
#' used to compute groups of features before annotation.
#' 
#' @slot 'cliques'
#' Is a list that contains the groups of features.
#' Each id corresponds to a row in the peaklist.
#'
#' @slot 'isotopes'
#' Is a data.frame with the column 'feature' for feature id, 
#' column 'charge' for the charge, column 'grade' that starts with 0
#' and it is 1 for the first isotope, 2 for the second and so on and
#' column 'cluster' which labels each group of features that are
#' isotopes.
#'
#' @slot 'cliquesFound' is
#' TRUE if clique groups have been computed,
#' @slot 'isoFound' is
#' TRUE if isotopes have been annotated,
#' @slot 'anFound' is
#' TRUE if annotation of adducts have been computed.
#' @return An 'anClique' object with annotation of isotopes, adducts
#' and fragments, and information about the annotation process.
#' @examples
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' library(xcms)
#' mzraw <- MSnbase::readMSData(files = mzfile, mode = "onDisk")
#' cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' mzData <- findChromPeaks(object = mzraw, param = cpw)
#' ex.anClique <- createanClique(mzdata = mzData)
#' show(ex.anClique)
#' @seealso \code{\link{createanClique}}
setClass(
    Class = "anClique",
    representation(
        peaklist = "data.frame",
        network = "igraph",    
        cliques = "list",
        cliquesFound = "logical",
        isotopes = "data.frame",
        isoFound = "logical",
        anFound = "logical")
)

#' @title 'anClique' class constructor
#' @description
#' S4 Class \code{anClique} for annotating isotopes and adducts
#' in processed m/z data. Features are first
#' grouped based on a similarity network algorithm and then
#' annotation of isotopes and adducts is performed in each group.
#' @param peaklist 'data.frame' with feature and annotation information.
#' @param network 'igraph' undirected network of similarity.
#' @param cliques list with the groups of features
#' @param isotopes 'data.frame' with isotope annotation.
#' @param anFound 'TRUE' if annotation has been computed.
#' @param cliquesFound 'TRUE' if cliques have been computed.
#' @param isoFound 'TRUE' if isotopes have been computed.
#' @details See \code{help("anClique-class")} for information about the
#' slots and methods of the S4 class 'anClique'.
#' @return
#' A new 'anClique' object with variable values set by the user.
#' @examples
#' mzfile <- system.file("standards.mzXML", package = "cliqueMS")
#' library(xcms)
#' mzraw <- MSnbase::readMSData(files = mzfile, mode = "onDisk")
#' cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
#' mzData <- findChromPeaks(object = mzraw, param = cpw)
#' ex.anClique <- createanClique(mzdata = mzData)
#' show(ex.anClique)
#' @seealso \code{\link{createanClique}}
#' @seealso \code{\link{anClique-class}}
anClique <- function(peaklist = data.frame(),
    network = igraph::make_empty_graph(directed = FALSE),
    cliques = list(), cliquesFound = FALSE,
    isotopes = data.frame(), isoFound = FALSE,
    anFound = FALSE) {
    object  <- new("anClique", peaklist = peaklist, network = network,
        cliques = cliques, cliquesFound = cliquesFound,
        isotopes = isotopes, isoFound = isoFound, anFound = anFound)
    return(object)
}

    
#' @export
#' @describeIn anClique show information about the object
setMethod("show", "anClique", function(object) {
    message(paste("anClique S4 object with",nrow(object@peaklist),
    "features"), sep = " ")
    if(object@cliquesFound) {
        message(paste("Features have been splitted into",
        length(object@cliques), "cliques", sep = " "))
    } else {
        message("No computed clique groups")
    }
    if(object@isoFound) {
        if( sum(is.na(unlist(object@isotopes))) ==
            length(unlist(object@isotopes)) ) {
            message("0 Features are isotopes")
        } else {
            message(paste(nrow(object@isotopes), "Features are isotopes",
                sep = " "))
        }
    } else {
        message("No isotope annotation")
    }
    if(object@anFound) {
        pos1 = which(!is.na(object@peaklist$an1))
        pos2 = which(!is.na(object@peaklist$an2))
        pos3 = which(!is.na(object@peaklist$an3))
        pos4 = which(!is.na(object@peaklist$an4))
        pos5 = which(!is.na(object@peaklist$an5))
        anFeatures = unique(c(pos1,pos2,pos3,pos4,pos5))
        message(paste(length(anFeatures), "features annotated", sep = " "))
    } else {
        message("No adduct annotation")
    }
})

#' @export
#' @describeIn createanClique Method for 'xcmsSet' object
setMethod("createanClique", "xcmsSet", function(mzdata) {
    if (!requireNamespace("CAMERA", quietly = TRUE)) {
        stop("Package CAMERA needed for 'xcmsSet' processed data. Please use
'XCMSnExp' objects or install package CAMERA.",
        call. = FALSE)
    }
    peaklist <- as.data.frame(xcms::peaks(mzdata))
    object <- new("anClique",
    peaklist = peaklist,
    network = igraph::make_empty_graph(n = 1),
    cliques = list(),
    cliquesFound = FALSE,
    isotopes = data.frame(),
    isoFound = FALSE,
    anFound = FALSE)
    return(object)
})

#' @export
#' @describeIn createanClique Method for 'XCMSnExp' object
setMethod("createanClique", "XCMSnExp", function(mzdata) {
    peaklist <- as.data.frame(xcms::chromPeaks(mzdata))
    object <- new("anClique",
    peaklist = peaklist,
    network = igraph::make_empty_graph(n = 1),
    cliques = list(),
    cliquesFound = FALSE,
    isotopes = data.frame(),
    isoFound = FALSE,
    anFound = FALSE)
    return(object)
})

## definition of accesors

#' @export
#' @describeIn anClique get the list of features with current annotation
setMethod("getPeaklistanClique", signature = "anClique", function(object) {
    return(object@peaklist)
})

#' @export
#' @describeIn anClique get the correlation network
setMethod("getNetanClique", signature = "anClique", function(object) {
    return(object@network)
})

#' @export
#' @describeIn anClique get the table of isotopes
setMethod("getIsolistanClique", signature = "anClique", function(object) {
    return(object@isotopes)
})

#' @export
#' @describeIn anClique get the list of the clique groups
setMethod("getlistofCliques", signature = "anClique", function(object) {
    return(object@cliques)
})

#' @export
#' @describeIn anClique is 'TRUE' if annotation has been computed
setMethod("hasAnnotation", signature = "anClique", function(object) {
    return(object@anFound)
})

#' @export
#' @describeIn anClique is 'TRUE' if cliques have been computed
setMethod("hasCliques", signature = "anClique", function(object) {
    return(object@cliquesFound)
})

#' @export
#' @describeIn anClique is 'TRUE' if isotopes have been computed
setMethod("hasIsotopes", signature = "anClique", function(object) {
    return(object@isoFound)
})

## definition of setters

#' @export
#' @param object 'anClique' S4 object.
#' @param value Is the new variable which can be a 'peaklist',
#' a 'network', a 'isotopes' a 'cliques' a 'cliquesFound'
#' a 'isoFound' or 'anFound' and it's set by the user.
#' @describeIn anClique set the table of isotopes
setReplaceMethod("getIsolistanClique", signature = "anClique",
    definition = function(object, value) {
        object@isotopes <- value
        return(object)
    }
)

#' @export
#' @describeIn anClique set the network of correlation
setReplaceMethod("getNetanClique", signature = "anClique",
    definition = function(object, value) {
        object@network <- value
        return(object)
    }
)

#' @export
#' @describeIn anClique set the list of clique groups
setReplaceMethod("getlistofCliques", signature = "anClique",
    definition = function(object, value) {
        object@cliques <- value
        return(object)
    }
)

#' @export
#' @describeIn anClique set the list of features
setReplaceMethod(f = "getPeaklistanClique", signature = "anClique",
    definition = function(object, value) {
        object@peaklist <- value
        return(object)
    }
)

#' @export
#' @describeIn anClique set if annotation has been computed
setReplaceMethod("hasAnnotation", signature = "anClique",
    definition = function(object, value) {
        object@anFound <- value
        return(object)
    }
)

#' @export
#' @describeIn anClique set if cliques have been computed
setReplaceMethod("hasCliques", signature = "anClique",
    definition = function(object, value) {
        object@cliquesFound <- value
        return(object)
    }
)


#' @export
#' @describeIn anClique set if isotopes have been computed
setReplaceMethod("hasIsotopes", signature = "anClique",
    definition = function(object, value) {
        object@isoFound <- value
        return(object)
    }
)


