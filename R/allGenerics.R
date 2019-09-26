## C

#' @export
setGeneric("createanClique", function(mzdata) {
    standardGeneric("createanClique")
})

#' @export
setGeneric("createNetwork", function(mzdata, peaklist, filter = TRUE,
    mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    standardGeneric("createNetwork")
})

## G

#' @export
setGeneric("getlistofCliques", function(object) {
    standardGeneric ("getlistofCliques")
})

#' @export
setGeneric("getIsolistanClique", function(object) {
    standardGeneric ("getIsolistanClique")
})

#' @export
setGeneric("getNetanClique", function(object) {
    standardGeneric ("getNetanClique")
})
#' @export
setGeneric("getNetanClique<-", function(object, value) {
    standardGeneric ("getNetanClique<-")
})

#' @export
setGeneric("getPeaklistanClique", function(object) {
    standardGeneric ("getPeaklistanClique")
})
#' @export
setGeneric("getPeaklistanClique<-", function(object, value) {
    standardGeneric ("getPeaklistanClique<-")
})
