

#' Phylogenetically constrained  clustering
#'
#' This function extracts clusters of species based on traits values and phylogenetic proximities.
#' 
#' @param p4d a phylo4d object.
#' @param trait the trait(s) in the phylo4d object to use for clustering.
#' Can be a character vector giving the name of the trait(s) or numbers giving the column index
#' in the table of the data slot of the p4d object.
#' @param lim.phylo the maximum phylogenetic distance for edges selection.
#' @param lim.trait the maximum traits based distance for edges selection.
#' @param dist.phylo a character string specifying the method used to compute phylogenetic distances.
#' See Details.
#' @param dist.trait a character string specifying the method used to compute traits distances.
#' See Details.
#' @param select.method a character string specifying the method used to select edges.
#' @param scale.lim logical (default TRUE) indicating if \code{lim.phylo} and \code{lim.trait} are scaled
#' (divided by the max value)
#'
#' @return An object of class \code{graphclust}.
#'
#' @export
graphClust <- function(p4d, trait = names(tdata(p4d)),
                       lim.phylo = 0.2, lim.trait = 0.2, select.method = "ellipse",
                       dist.phylo = "patristic", dist.trait = "euclidean",
                       scale.lim = TRUE){
  
  select.method <- match.arg(select.method, c("line", "rectangle", "ellipse"))
  
  if(is.numeric(trait)){
    trait <- names(tdata(p4d))[trait]
  }
  
  p4 <- extractTree(p4d)
  phy <- as(p4, "phylo")
  new.order <- phy$edge[, 2][!phy$edge[, 2] %in% phy$edge[, 1]]
  tips <- phy$tip.label[new.order]
  n.tips <- length(tips)
  X <- tdata(p4d, type = "tip")
  X <- X[tips, trait]
  X <- scale(X)
  X <- as.data.frame(X)
  colnames(X) <- trait
  rownames(X) <- tips
  n.traits <- ncol(X)
  
  dmat.phylo <- distTips(phy, method = dist.phylo)
  dmat.phylo <- as.matrix(dmat.phylo)
  dmat.phylo <- dmat.phylo[tips, tips]
  
  dmat.trait <- dist(X, method = dist.trait)
  dmat.trait <- as.matrix(dmat.trait)
  dmat.trait <- dmat.trait[tips, tips]
  
  if(scale.lim){
    lim.p <- max(dmat.phylo) * lim.phylo
    lim.t <- max(dmat.trait) * lim.trait
  } else {
    lim.p <- lim.phylo
    lim.t <- lim.trait
  }
  
  if(select.method == "line"){
    adj.mat <- (((-lim.t/lim.p) * dmat.phylo) + lim.t) > dmat.trait
  }
  if(select.method == "rectangle"){
    adj.mat <- (dmat.phylo < lim.p) + (dmat.trait < lim.t) == 2
  }
  if(select.method == "ellipse"){
    adj.mat <- (dmat.phylo/lim.p)^2 + (dmat.trait/lim.t)^2 < 1
  }
  
  gr <- graph.adjacency(adj.mat, mode = "lower", diag = FALSE)
  
  gr.decomp <- decompose.graph(gr)
  gr.decomp.names <- lapply(gr.decomp, function(x) V(x)$name)
  gr.decomp.den <- sapply(gr.decomp, graph.density)
  names(gr.decomp.den) <- 1:length(gr.decomp)
    
  clust.gr <- clusters(gr)
  clust <- clust.gr$membership
  names(clust) <- V(gr)$name
  
  inclusion <- (degree(gr) + 1) / clust.gr$csize[clust.gr$membership]
  
  meta <- list(p4d = p4d, trait = trait, adj.mat = adj.mat,
               dmat.phylo = dmat.phylo, dmat.trait = dmat.trait,
               lim.p = lim.p, lim.t = lim.t,
               select.method = select.method,
               dist.phylo = dist.phylo, dist.trait = dist.trait,
               graph = gr, clust.gr = clust.gr)
  res <- list(clusters = clust,
              clusters.density = gr.decomp.den,
              taxa.inclusion = inclusion,
              meta = meta)
  
  class(res) <- "graphclust"
  return(res)
}



#' Plot phylogenetically constrained  clustering
#'
#' This function produces three plots (selectable by which):
#' a plot of edges selection based on phylogenetic against trait distances of taxa pairs,
#' a plot of the graph produced with the selected edges
#' and a plot of the clustered phylogenetic tree.
#' 
#' @param x a graphclust object as produced by graphClust.
#' @param which a character vector to select plots.
#' Must be one or more of \code{"selection"}, \code{"graph"}, \code{"tree"}.
#' @param logical if TRUE (default), the user is asked before each plot
#' @param colored logical indicating if plots include colors.
#'
#' @export
plot.graphclust <- function(x, which = c("selection", "graph", "tree"), ask = TRUE, colored = TRUE){
  if(class(x) != "graphclust"){
    stop("x must be an object of class 'graphclust'")
  }
  which <- match.arg(which, c("selection", "graph", "tree"), several.ok = TRUE)
  ask0 <- par("ask")
  
  if(x$meta$select.method == "line"){
    gx <- c(0, x$meta$lim.p)
    gy <- c(x$meta$lim.t, 0)
  }
  if(x$meta$select.method == "rectangle"){
    gx <- c(0, x$meta$lim.p, x$meta$lim.p)
    gy <- c(x$meta$lim.t, x$meta$lim.t, 0)
  }
  if(x$meta$select.method == "ellipse"){
    gx <- seq(0, x$meta$lim.p, length.out = 100)
    gy <- sqrt((1-((gx^2)/(x$meta$lim.p^2)))*(x$meta$lim.t^2))
  }
  
  if(ask) par(ask = TRUE)
  if("selection" %in% which){
    if(colored){
      dot.col <- x$meta$adj.mat + 1
      line.col <- 2
    } else {
      dot.col <- 1
      line.col <- 1
    }
  plot(as.vector(x$meta$dmat.phylo), as.vector(x$meta$dmat.trait),
       col = dot.col,
       pch = ifelse(x$meta$adj.mat == 0, 1, 4),
       main = "Selected Edges",
       xlab = paste("Phylogenetic distance (", x$meta$dist.phylo, ")", sep = ""),
       ylab = paste("Trait distance (", x$meta$dist.trait, ")", sep = ""))
  mtext(paste0(sum(x$meta$adj.mat[lower.tri(x$meta$adj.mat, diag = FALSE)]), "/",
               length(x$taxa.inclusion) * length(x$taxa.inclusion) - length(x$taxa.inclusion)),
        line = 0.5, side= 3)
  lines(gx, gy, col = line.col, lwd = 2, lty = "dashed")
  }
  if(colored){
    graph.gcol <- tree.gcol <- evenColors(x$meta$clust.gr$no)[x$clusters]
    names(graph.gcol) <- names(tree.gcol) <- names(x$clusters)
  } else {
    graph.gcol <- "grey"
    tree.gcol <- "grey35"
  }
  
  if("graph" %in% which){
    plot.igraph(x$meta$graph, vertex.color = graph.gcol, vertex.label.color = 1, vertex.label.family = "")
  }
  if("tree" %in% which){
    barplot(x$meta$p4d, x$meta$trait, bar.col = tree.gcol, center = T, scale = T)
  }
  par(ask = ask0)
}

#' @export
print.graphclust <- function(x){
  cat("Phylogenetically constrained clustering of", length(x$taxa.inclusion), "taxa\n")
  cat("\t Included trait(s):", paste(x$meta$trait, collapse=", "), "\n")
  cat("\t Phylogenetic distance:", x$meta$dist.phylo, "\n")
  cat("\t Trait distance:", x$meta$dist.trait, "\n")
  cat("\t Selected edges:", sum(x$meta$adj.mat[lower.tri(x$meta$adj.mat, diag = FALSE)]), "\n")
  cat("\t Number of clusters:", x$meta$clust.gr$no)
}

