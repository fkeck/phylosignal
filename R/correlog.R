

#' Phylogenetic correlograms
#'
#' This function estimate phylogenetic correlograms for one trait
#' (univariate Moran's I correlogram) or several (multivariate Mantel correlogram).
#' 
#' @param p4d a phylo4d object.
#' @param trait the trait(s) in the p4d object to use for correlogram estimation.
#' Can be a character vector giving the name of the trait(s) or a single number giving the column index
#' in the table of the data slot of the p4d object.
#' @param dist.phylo a character string specifying the method used to compute phylogenetic distances.
#' See Details.
#' @param increment a numeric value. Increment parameter for the uniformly distributed distance classes.
#' @param reps a numeric value. Number of repetitions for the estimation of p.values with randomization.
#' 
#' @return an object of class \code{correlog}.
#' @author This is a modified version of the \code{correlog}
#' function in the \pkg{ncf} package written by Ottar N. Bjornstad.
#'@rdname correlogram
#'@export
correlogram.phylo4d <- function (p4d, trait = names(tdata(p4d)),
                                 dist.phylo = "patristic", increment = 2, reps = 1000){
  
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
  n.traits <- ncol(X)
  if(is.numeric(trait)){
    trait <- names(tdata(p4d))[trait]
  }
  
  if(n.traits > 1){
    X <- as.matrix(X)
    moran <- cor2(t(X), circ = FALSE)
    moran <- moran[lower.tri(moran)]
    moran <- moran - mean(moran, na.rm = TRUE)
  } else {
    X <- as.vector(X)
    X <- scale(X)[, 1]
    moran <- t(outer(X, X))
    moran <- moran[lower.tri(moran)]
  }
  
  dmat <- distTips(phy, method = dist.phylo)
  dmat <- as.matrix(dmat)
  dmat <- dmat[tips, tips]
                    
  if (reps != 0) {
    dmat2 <- dmat
    moran2 <- moran
  }
  
  dmat <- dmat[lower.tri(dmat)]
  
  dkl <- ceiling(dmat/increment)        #Affectation  des points à des grps sur la base de leur distance
  nlok <- sapply(split(moran, dkl), length) # Découpage et comptage des n de chaque grp
  dmean <- sapply(split(dmat, dkl), mean, na.rm = TRUE) # Idem : moyenne distance
  moran <- sapply(split(moran, dkl), mean, na.rm = TRUE) #Idem : moyenne moran
  
  p <- NULL
  if (reps != 0) {
    perm <- matrix(NA, ncol = length(moran), nrow = reps)
    for (i in 1:reps) {
      trekk <- sample(1:n.tips)
      dma <- dmat2[trekk, trekk]
      mor <- moran2
      dma <- dma[lower.tri(dma)]
      dkl <- ceiling(dma/increment)
      perm[i, ] <- sapply(split(mor, dkl), mean, na.rm = TRUE)
    } 
    p <- (apply(moran <= t(perm), 1, sum))/(reps + 1)
    p <- apply(cbind(p, 1 - p), 1, min) + 1/(reps + 1)
  }
  
  res <- list(n = nlok, mean.of.class = dmean,
              correlation = moran, p = p)
  
  class(res) <- "correlog"
  return(res)
}
  




#' Phylogenetic spline correlograms
#'
#' This function estimate phylogenetic spline correlograms for one trait
#' (univariate Moran's I spline correlogram) or several (multivariate Mantel spline correlogram).
#' 
#' @param p4d a phylo4d object.
#' @param trait the trait(s) in the p4d object to use for correlogram estimation.
#' Can be a character vector giving the name of the trait(s) or a single number giving the column index
#' in the table of the data slot of the p4d object.
#' @param dist.phylo a character string specifying the method used to compute phylogenetic distances.
#' See Details.
#' @param df numeric. Degrees of freedom for the spline. Default is \code{sqrt(n)}.
#' @param npoints the number of points at which to save the value for the spline function
#' (and confidence envelope).
#' @param max.it a numeric. The maximum iteration for the Newton method used to estimate the intercepts.
#' @param reps a numeric value. Number of repetitions for the bootstraps.
#' 
#' @return an object of class \code{spline.correlog}.
#' @references Bjornstad, O.N. & Falck, W. (2001)
#' Nonparametric spatial covariance functions: estimation and testing.
#' Environmental and Ecological Statistics, 8:53-70.
#' @author This is a modified version of the \code{spline.correlog}
#' function in the \pkg{ncf} package written by Ottar N. Bjornstad.
#' @rdname correlogram
#' @export
spline.correlogram.phylo4d <- function (p4d, trait = names(tdata(p4d)), dist.phylo = "patristic",
                                        df = NULL, reps = 200,  npoints = 100, max.it = 25){
  
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
  n.traits <- ncol(X)
  if(is.numeric(trait)){
    trait <- names(tdata(p4d))[trait]
  }
  
  if(n.traits > 1){
    X <- as.matrix(X)
    X <- (t(scale(X))) / (sqrt((n.tips - 1)/n.tips))
    moran <- cor2(t(X), circ = FALSE)
    moran <- moran - mean(moran[lower.tri(moran)], na.rm = TRUE)
  } else {
    X <- as.vector(X)
    X <- scale(X)[, 1] / (sqrt((n.tips - 1)/n.tips))
    moran <- t(outer(X, X))
  }
  
  if (is.null(df)) {
    df <- sqrt(n.tips)
  }
  
  dmat <- distTips(phy, method = dist.phylo)
  dmat <- as.matrix(dmat)
  dmat <- dmat[tips, tips]
  maxdist <- max(na.omit(dmat))
  
  real <- list(x.intercept = NA, e.intercept = NA, y.intercept = NA, 
               predicted = list(x = matrix(NA, nrow = 1, ncol = npoints), 
                                y = matrix(NA, nrow = 1, ncol = npoints)))
  
  triang <- lower.tri(dmat)
  u <- dmat[triang]
  v <- moran[triang]
  sel <- is.finite(u) & is.finite(v)
  u <- u[sel]
  v <- v[sel]
  sobj <- smooth.spline(u, v, df = df)
  xpoints <- seq(0, maxdist, length = npoints)
  lx <- predict(sobj, x = xpoints)
  real$y.intercept <- lx$y[1]
  real$predicted <- list(x = xpoints, y = lx$y)
  konst <- 1
  if (real$y.intercept < 0) {
    lx$y <- -lx$y
    konst <- -1
  }
  ly <- 1:length(lx$y)
  choise <- ly[lx$y < 0][1]
  pos <- lx$x[choise - 1]
  neg <- lx$x[choise]
  pos <- pos + (neg - pos)/2
  tmp <- smooth.spline(lx)
  for (j in 1:max.it) {
    if (is.na(neg)) {
      pos <- NA
      break
    }
    if (neg == 0) {
      pos <- 0
      break
    }
    neg <- pos - predict(tmp, pos)$y/predict(tmp, pos, deriv = 1)$y
    if (abs(pos - neg) < 1e-06) {
      break
    }
    pos <- neg
  }
  real$x.intercept <- konst * pos
  sobj <- smooth.spline(u, v - 1/exp(1), df = df)
  lx <- predict(sobj, x = xpoints)
  ly <- 1:length(lx$y)
  choise <- ly[lx$y < 0][1]
  pos <- lx$x[choise - 1]
  neg <- lx$x[choise]
  pos <- pos + (neg - pos)/2
  tmp <- smooth.spline(lx)
  for (j in 1:max.it) {
    if (is.na(neg)) {
      pos <- NA
      break
    }
    if (neg == 0) {
      pos <- 0
      break
    }
    neg <- pos - predict(tmp, pos)$y/predict(tmp, pos, deriv = 1)$y
    if (abs(pos - neg) < 1e-06) {
      break
    }
    pos <- neg
  }
  real$e.intercept <- pos
  boot <- list(NULL)
  boot$boot.summary <- list(NULL)
  if (reps != 0) {
    boot$boot.summary$x.intercept <- matrix(NA, nrow = reps, 
                                            ncol = 1)
    boot$boot.summary$e.intercept <- matrix(NA, nrow = reps, 
                                            ncol = 1)
    boot$boot.summary$y.intercept <- matrix(NA, nrow = reps, 
                                            ncol = 1)
    predicted <- list(x = matrix(NA, nrow = 1, ncol = npoints), 
                      y = matrix(NA, nrow = reps, ncol = npoints))
    predicted$x[1, ] <- xpoints
    for (i in 1:reps) {
      trekkx <- sample(1:n.tips, replace = TRUE)
      trekky <- trekkx
      dmatb <- dmat[trekkx, trekkx]
      triang <- lower.tri(dmat)
      dmatb <- dmatb[triang]
      moranb <- moran[trekky, trekky][triang]
      moranb <- moranb[!(dmatb == 0)]
      dmatb <- dmatb[!(dmatb == 0)]
      u <- dmatb
      v <- moranb
      sel <- is.finite(u) & is.finite(v)
      u <- u[sel]
      v <- v[sel]
      sobj <- smooth.spline(u, v, df = df)
      lx <- predict(sobj, x = xpoints)
      boot$boot.summary$y.intercept[i, 1] <- lx$y[1]
      predicted$y[i, ] <- lx$y
      konst <- 1
      if (boot$boot.summary$y.intercept[i, 1] < 0) {
        lx$y <- -lx$y
        konst <- -1
      }
      ly <- 1:length(lx$y)
      choise <- ly[lx$y < 0][1]
      pos <- lx$x[choise - 1]
      neg <- lx$x[choise]
      pos <- pos + (neg - pos)/2
      tmp <- smooth.spline(lx)
      for (j in 1:max.it) {
        if (is.na(neg)) {
          pos <- NA
          break
        }
        if (neg == 0) {
          pos <- 0
          break
        }
        neg <- pos - predict(tmp, pos)$y/predict(tmp, 
                                                 pos, deriv = 1)$y
        if (abs(pos - neg) < 1e-06) {
          break
        }
        pos <- neg
      }
      boot$boot.summary$x.intercept[i, 1] <- konst * pos
      sobj <- smooth.spline(u, v - 1/exp(1), df = df)
      lx <- predict(sobj, x = xpoints)
      ly <- 1:length(lx$y)
      choise <- ly[lx$y < 0][1]
      pos <- lx$x[choise - 1]
      neg <- lx$x[choise]
      pos <- pos + (neg - pos)/2
      tmp <- smooth.spline(lx)
      for (j in 1:max.it) {
        if (is.na(neg)) {
          pos <- NA
          break
        }
        if (neg == 0) {
          pos <- 0
          break
        }
        neg <- pos - predict(tmp, pos)$y/predict(tmp, 
                                                 pos, deriv = 1)$y
        if (abs(pos - neg) < 1e-06) 
          break
        pos <- neg
      }
      boot$boot.summary$e.intercept[i, 1] <- pos
    }
      boot$boot <- NULL

    ty <- apply(predicted$y, 2, quantile,
                probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), 
                na.rm = TRUE)
    dimnames(ty) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 
                           0.75, 0.9, 0.95, 0.975, 1), NULL)
    tx <- predicted$x
    boot$boot.summary$predicted <- list(x = tx, y = ty)
  }
  else {
    boot <- NULL
    boot.summary <- NULL
  }
  res <- list(real = real, boot = boot, max.distance = maxdist)
  class(res) <- "spline.correlog"
  res
}
  
  
  