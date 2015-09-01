#' phylosignal
#' 
#' @useDynLib phylosignal
#' @import RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import adephylo
#' @import ape
#' @importFrom phylobase extractTree tdata tipLabels
#' @importMethodsFrom phylobase as
#' @importFrom boot boot boot.ci
#' @importFrom igraph graph.adjacency decompose.graph V graph.density clusters degree plot.igraph
#' 
#' @name phylosignal
#' @docType package
NULL
