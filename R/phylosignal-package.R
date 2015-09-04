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


#' Phylogeny and pollution sensitivity of diatoms
#' 
#' A phylogenetic tree and the pollution sensitivity of 17 diatoms species from the order Naviculales.
#' 
#' @format a \code{phylo4d} object.
#' @source Keck F., Rimet F., Franc A. & Bouchez A. (In press) Phylogenetic signal in diatom ecology: perspectives for aquatic ecosystems biomonitoring. Ecological Applications.
#' @docType data
#' @keywords datasets
#' @name navic
#' @usage data(navic)
NULL