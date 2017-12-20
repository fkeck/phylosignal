#' phylosignal
#' 
#' @useDynLib phylosignal
#' @import DBI
#' @import grDevices
#' @import graphics
#' @import methods
#' @import stats
#' @import utils
#' @importFrom boot boot boot.ci
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom ape vcv.phylo node.height node.depth.edgelength unrooted.xy .PlotPhyloEnv extract.clade ladderize plot.phylo rTraitCont read.tree read.nexus
#' @importFrom phylobase phylo4d extractTree tdata tipLabels tipData nTips nNodes nodeData<- nodeData
#' @importFrom igraph graph.adjacency decompose.graph V graph.density clusters degree plot.igraph
#' @importFrom adephylo proxTips distTips
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