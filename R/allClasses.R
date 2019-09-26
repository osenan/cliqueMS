##################################################
#     DEFINITION OF S3 METHODS AND S3 CLASSES    #
##################################################


setOldClass("igraph")

#' @export
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

#' @export
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
        message(paste(nrow(object@isotopes), "Features are isotopes", sep = " "))
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
setMethod("createanClique", "XCMSnExp", function(mzdata) {
    peaklist <- as.data.frame(xcms::chromPeaks(mzData))
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
setMethod("getPeaklistanClique", signature = "anClique", function(object) {
    return(object@peaklist)
})

#' @export
setMethod("getNetanClique", signature = "anClique", function(object) {
    return(object@network)
})

#' @export
setMethod("getIsolistanClique", signature = "anClique", function(object) {
    return(object@isotopes)
})

#' @export
setMethod("getlistofCliques", signature = "anClique", function(object) {
    return(object@cliques)
})

## definition of setters

#' @export
setReplaceMethod(f = "getPeaklistanClique", signature = "anClique",
    definition = function(object, value) {
        object@peaklist <- value
        return(object)
    }
)

#' @export
setReplaceMethod("getNetanClique", signature = "anClique",
    definition = function(object, value) {
        object@network <- value
        return(object)
    }
)



