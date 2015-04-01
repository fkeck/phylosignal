

#' Computes phylogenetic signal with different methods
#'
#' This function computes phylogenetic signal statistics
#' (Blomberg's K and K*, Abouheif's Cmean, Moran's I, and Pagel's Lambda) for traits in phylo4d objects.
#'
#' @param p4d a phylo4d object (see Details).
#' @param methods a character vector giving the methods to compute phylogenetic signal (see Details).
#' @param reps a numeric value. Number of repetitions for the estimation of p.values with randomization.
#'
#'@details p4d has to be a phylo4d object as defined in phylobase package. phylo4d objects combine
#'phylogenetic tree with data. Each data associated with tips is interpreted as a trait.
#'By default, the \code{methods} argument is set to "all" and all the available methods are used.
#'The user can specify which method(s) to use. Possible values are
#'"\code{I}", "\code{Cmean}", "\code{Lambda}", "\code{K}" and "\code{K.star}".
#'
#'@return A list of two dataframes with the values of statistics and associated
#'p.values for each tested trait and method.
#'
#'@seealso \code{\link{phyloSimSignal}} .
#'
#'@export
phyloSignal <- function(p4d, methods = c("all", "I", "Cmean", "Lambda", "K", "K.star"), reps = 999){
  methods <- match.arg(methods, several.ok = TRUE)
  p4 <- extractTree(p4d)
  phy <- as(p4, "phylo")
  X <- tdata(p4d, type = "tip")
  X <- as.matrix(X)
  res <- list()
  
  if("all" %in% methods){
    methods <- c("I", "Cmean", "Lambda", "K", "K.star")
  }
  
  if("Cmean" %in% methods){
    W <- proxTips(p4d, method = "Abouheif", useC = TRUE)
    tmp <- apply(X, 2, moranTest, Wr = W, reps = reps)
    res$stat$Cmean <- unlist(sapply(tmp, "[", 1))
    res$pvalue$Cmean <- unlist(sapply(tmp, "[", 2))
  }
  
  if("I" %in% methods){
    W <- proxTips(p4d, method = "patristic", useC = TRUE)
    tmp <- apply(X, 2, moranTest, Wr = W, reps = reps)
    res$stat$I <- unlist(sapply(tmp, "[", 1))
    res$pvalue$I <- unlist(sapply(tmp, "[", 2))
  }
  
  if(any(c("K", "K.star", "Lambda") %in% methods)){  
    VCV <- vcv.phylo(phy, model = "Brownian")
    
    if("K" %in% methods){
      tmp <- apply(X, 2, kTest, vcv = VCV, reps = reps)
      res$stat$K <- unlist(sapply(tmp, "[", 1))
      res$pvalue$K <- unlist(sapply(tmp, "[", 2))
    }
    
    if("K.star" %in% methods){
      tmp <- apply(X, 2, kStarTest, vcv = VCV, reps = reps)
      res$stat$K.star <- unlist(sapply(tmp, "[", 1))
      res$pvalue$K.star <- unlist(sapply(tmp, "[", 2))
    }
    
    if("Lambda" %in% methods){
      tmp <- apply(X, 2, lambdaTest, vcv = VCV)
      res$stat$Lambda <- unlist(sapply(tmp, "[", 1))
      res$pvalue$Lambda <- unlist(sapply(tmp, "[", 2))
    }
  }
  
  res$stat <- as.data.frame(res$stat, row.names = colnames(X))
  res$pvalue <- as.data.frame(res$pvalue, row.names = colnames(X))
  return(res)
}





#Force un retour sur R pour l'optimisation... Essayer d'implémenter ça direct en c++ avec la lib GSL.
lambdaTest <- function(x, vcv){
  lambda.max <- max(vcv)/max(vcv[lower.tri(vcv)])
  opt <- suppressWarnings(
    tryCatch(
      optimize(pagelLogLik, interval=c(0, lambda.max), xr = x, vcvr = vcv, maximum = TRUE),
             warning = function(w) return(list(maximum = NA, objective = NA))
      )
    )
  Lambda <- opt$maximum
  if(is.finite(opt$objective)){
    logL0 <- pagelLogLik(0, x, vcv)
    if(is.finite(logL0)){
      pvalue <- as.numeric(pchisq(2 * (opt$objective - logL0), df = 1, lower.tail = FALSE))
    } else {
      pvalue <- NA
    }
  } else {
    pvalue <- NA
  }  
  return(list(Lambda = Lambda, pvalue = round(pvalue, 3)))
}


