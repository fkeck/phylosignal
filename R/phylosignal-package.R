#' phylosignal
#' 
#' @useDynLib phylosignal
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import adephylo
#' @import ape
#' @import RCurl
#' @importFrom phylobase phylo4d extractTree tdata tipLabels
#' @importMethodsFrom phylobase as
#' @importFrom boot boot boot.ci
#' @importFrom igraph graph.adjacency decompose.graph V graph.density clusters degree plot.igraph
#' 
#' @name phylosignal
#' @docType package
NULL
