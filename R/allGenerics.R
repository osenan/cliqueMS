## C

#' @export
setGeneric("createS4anClique", function(mzdata) {
    standardGeneric("createS4anClique")
})

#' @export
setGeneric("createS4Network", function(mzdata, peaklist, filter = TRUE,
    mzerror = 5e-6, intdiff = 1e-4, rtdiff = 1e-4) {
    standardGeneric("createS4Network")
})

## G

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
