#' @title 'cliqueMS' annotates isotopes and adducts in m/z data
#' @description
#' 'cliqueMS' first separates features in the data into different
#' groups. To do this it computes a similarity weighted network
#' from the data, and searches clique groups. This cliques are
#' fully connected components that have higher similarity in
#' inner edges than edges outside cliques.
#'
#' Once clique groups are computed, annotation of isotopes is
#' first performed. After isotope annotation, adducts are
#' annotated within each group
#'
#' @section Functions
#'
#' 'getCliques' for separate features into clique groups
#' 'getIsotopes' for annotate isotopes
#' 'getAnnotation' for annotate adducts
#' @docType package
#' @name cliqueMS
NULL
