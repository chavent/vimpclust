#' @title L0 Sparse weighted k-means
#' @export
#' 
#' @description Sparse weighted k-means with L0 penalty. Two modes are proposed :
#' one with the regularization parameter called \code{s} (proposed in Chang & al. 2018),
#' and one with the regularization parameter called \code{lambda}. Here, the data are always 
#' standardized before clustering. Otherwise, variables with high variance are naturally selected.
#' 
#' @param X a dataframe of dimension \code{n} (observations) by \code{p} (variables) with
#'  numerical data. 
#' @param centers an integer representing the number of clusters. 
#' @param mode either "s" or "lambda".
#' @param grid a vector of numerical values (or a single value) providing 
#'  a grid of values for the regularization parameter.
#'  If NULL (by default), a grid of integer values between 1 and \code{maxs} when \code{mode=s},
#'  and a grid of  \code{nlambda} values in [0,\code{maxlambda}] when  \code{mode=lambda}.
#' @param nlambda The size of the grid of values for the regularization parameter (only
#' when  \code{mode=lambda}).
#' @param nstart an integer representing the number of random starts in the k-means algorithm.
#'  By default, \code{nstart=10}. 
#' @param itermaxw an integer indicating the maximum number of iterations for the inside 
#'  loop over the weights \code{w}. By default, \code{itermaxw=20}.
#' @param itermaxkm  an integer representing the maximum number of iterations in the k-means 
#'  algorithm. By default, \code{itermaxkm=10}.
#' @param epsilonw a positive numerical value. It provides the precision of the stopping 
#'  criterion over \code{w}. By default, \code{epsilonw =1e-04}. 
#' @param verbose an integer value. If \code{verbose=0}, the function stays silent, if \code{verbose=1} (default option), it  prints
#'  whether the stopping criterion over the weights \code{w} is satisfied.
#'
#' @returns A list with the following elements:
#'  \describe{
#'  \item{grid}{a numerical vector containing the regularization parameters 
#'  (a grid of values).}
#'  \item{W}{a \code{p} by \code{length(grid)} matrix. It contains the weights (0 or 1) 
#'  associated to each variable.}
#' \item{cluster}{a \code{n} by \code{length(grid)} integer matrix. It contains the 
#' cluster memberships, for each value of the regularization parameter.}
#' \item{bss.per.feature}{a matrix of size \code{p} by \code{length(grid)}. 
#' It contains the between-class variance computed on the \code{p} variables.}
#' \item{Z}{a matrix of size \code{n} by \code{p} containing the scaled data 
#' if \code{scaling=TRUE}, and a copy of \code{X} otherwise.}
#' }

#' @examples
#' # Swiss data:
#' data("swiss")
#' out <- sparsewkm_L0(X = swiss, centers = 3, mode = "s")
#' # grid of regularization parameters
#' out$s
#' # weights of the variables 
#' out$W
#' 
#' plot(out, what = "weights.features")
#' plot(out, what = "expl.var")
#' 
#' # partitioning obtained with s=2 selected variables
#' out$cluster[ , 2]
#' # between-class variance on each variable
#' out$bss.per.feature[, 2]
#' # between-class variance
#' sum(out$bss.per.feature[, 2])
#' 



sparsewkm_L0 <- function(X, centers, mode = "s", grid = NULL, 
                         nlambda = 20,
                         nstart = 10, itermaxw = 20, itermaxkm = 10, 
                         verbose = 1, epsilonw = 1e-04) {

  if (!mode %in% c("s", "lambda")) {
    stop("mode should be 's' or 'lambda'")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)
  X <- scale(X)*sqrt(n/(n-1))
  
  km <- stats::kmeans(X, centers, nstart, itermaxkm)
  clusterini <- km$cluster
  bss.per.featureini <- weightedss(X, clusterini)$bss.per.feature/nrow(X) 
  
  if ((mode == "lambda") & (!is.null(grid))) {
    check_lambda(grid)
    check_maxlambda(max(grid))
  }

   if (is.null(grid)) {
    if (mode == "s")  
      grid <- 1:ncol(X) 
    if (mode == "lambda")  
      grid <- seq(from = 0, to = 0.9*max(bss.per.featureini), length.out = nlambda)
  }
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste("V", 1:ncol(X), sep = "")
  }
  
  
  ls <- sapply(grid, FUN = function(i) {
                 
                 cluster1 <- clusterini
                 bss.per.feature <- bss.per.featureini
                 w1 <- rep(1/sqrt(ncol(X)), ncol(X))  
                 w0 <- abs(stats::rnorm(ncol(X)))
                 W_cl <- list()
                 niter <- 1
                 while (sum(abs(w1 - w0))/sum(abs(w0)) > epsilonw && niter < itermaxw) {
                   if (niter != 1) {
                     Xw <- t(t(X)*sqrt(w1))
                     km <- stats::kmeans(Xw, centers, nstart = nstart, iter.max = itermaxkm)
                     cluster1 <- km$cluster
                     bss.per.feature <- weightedss(X, cluster1)$bss.per.feature/nrow(X)
                   }
                   w0 <- w1
                   cluster0 <- cluster1
                   w1 <- hard.thresholding(bss.per.feature, i, mode = mode)
                   if (norm.vect(w1)==0)
                     niter <- itermaxw
                   niter = niter + 1
                 } 
                 
                 if (verbose == 1) 
                 {
                   if (niter == itermaxw+1) 
                     print(sprintf("all the feature weights are set to 0 for %s = %g", mode, i))
                   if (niter == itermaxw) 
                     print(sprintf("the stopping criterion over w is not satisfied for %s = %g ", mode, i))  
                   if (niter < itermaxw) 
                     print(sprintf("the stopping criterion over w is satisfied for %s = %g", mode, i))  
                 }
                 W_cl[[1]] <- w1
                 if (niter == itermaxw+1)
                 {
                   W_cl[[2]] <- rep(NA, length(cluster0))
                   W_cl[[3]] <- rep(NA, length(bss.per.feature))
                 } else
                 {
                   W_cl[[2]] <- cluster0
                   W_cl[[3]] <- bss.per.feature
                 }
                 W_cl
               }) 
  
  W <- do.call(cbind, (ls[1, ]))
  cluster <- do.call(cbind, (ls[2, ]))
  bss.per.feature <- do.call(cbind, (ls[3, ]))
  
  rownames(W) <- colnames(X)
  colnames(W) <- paste(sprintf("%s=", mode), round(grid, digits = 3), sep = "")
  
  rownames(bss.per.feature) <- colnames(X)
  colnames(bss.per.feature) <- paste(sprintf("%s=", mode), round(grid, digits = 3), sep = "")
  
  rownames(cluster) <- rownames(X)
  colnames(cluster)  <- paste(sprintf("%s=", mode), round(grid, digits = 3), sep = "")
  
  out <- list(type="L0Sparse", grid = grid,  mode = mode, W = W, 
              cluster = cluster,
              Z = X, bss.per.feature = bss.per.feature, index = 1:ncol(X))
  class(out) <- "spwkm"
  return(out)
}

hard.thresholding <- function(b, grid, mode = "s") {
  if (!mode %in% c("s", "lambda")) 
    stop("mode should be 's' or 'lambda'")
  
  if (mode == "s") {
    w <- rep(0, length(b))
    w[order(b, decreasing = TRUE)[1:grid]] <- 1
  }
  
  if (mode == "lambda") {
    w <- rep(0, length(b))
    w[b > grid] <- 1
  }
  return(w)
}


