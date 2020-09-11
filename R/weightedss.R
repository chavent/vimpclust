wss.components <- function(k,x,cl)
{
  res <- x[cl==k,]
  if (is.null(dim(res))==F)
  {
    res <- apply(res, 2, scale, center=T, scale=F)
    res <- apply(res, 2, function(x) sum(x^2)) 
  }
  else 
    res <- rep(0,dim(x)[2])
  return(res)
}


#' @title Weighted sum-of-squares criteria
#' @export
#'
#' @description This function computes various weighted sum-of-squares criteria for a given
#' partition of a dataset described by numerical features. 
#' 
#' @param X a matrice or a dataframe of size \code{n} (observations) by \code{p} (variables) 
#' with numerical features only.
#' @param cl a vector of integers of length \code{n}. It contains the cluster membership of the data.
#' @param w a numerical vector of length \code{p}. It contains the weights to be applied to the features. 
#' By default, \code{w=NULL}, which amounts to setting each weight equal to 1.
#'
#' @return \item{bss.per.feature}{a numerical vector of length \code{p} containing the weighted 
#' between sum-of-squares per feature.}
#' @return \item{wss.per.feature}{a numerical vector of length \code{p} containing the weighted 
#' within sum-of-squares per feature.}
#' @return \item{bss.per.cluster}{a numerical vector of length \code{K} (\code{K} is the number of 
#' clusters) containing the weighted between sum-of-squares per cluster.}
#' @return \item{wss.per.cluster}{a numerical vector of length \code{K} 
#' containing the weighted within sum-of-squares per cluster.}
#' @return \item{bss}{a scalar representing the weighted between sum-of-squares of the partition.
#' It may be computed as the sum over \code{bss.per.feature} or \code{bss.per.cluster}.}
#' @return \item{wss}{a scalar representing the weighted within sum-of-squares of the partition.
#' It may be computed as the sum over \code{wss.per.feature} or \code{wss.per.cluster}.}
#' 
#' @examples
#' data(iris)
#' out <- weightedss(X = iris[,1:4], cl = as.numeric(iris$Species))
#' out$bss.per.feature
#' out$bss.per.cluster
#' out$bss
#' 
#' w <- c(0.3,0.3,0.2,0.2)
#' out <- weightedss(X = iris[,1:4], cl = as.numeric(iris$Species), w=w)
#' out$bss.per.feature
#' out$bss.per.cluster
#' out$bss

weightedss <- function(X, cl, w = NULL) 
{
    check_fun_weightedss(X, cl, w)
    
    X = as.matrix(X)
    sizecolX = ncol(X)
    if (is.null(w)) 
        w = rep(1, sizecolX)
    
    K <- length(unique(cl))
    
    X.scaled=t(t(X)*sqrt(w))
    cluster.centers <- stats::aggregate(X.scaled, by=list(cl), FUN=mean)[,-1]
    cluster.densities <- table(cl)
    matrix.wss <- t(sapply(1:K, FUN=wss.components, x=X.scaled, cl=cl, simplify=T))
    matrix.bss <- (cluster.centers-matrix(apply(X.scaled,2,mean),nrow=K, ncol=sizecolX,byrow=T))^2*cluster.densities
    wss.per.feature <- apply(matrix.wss, 2, sum)
    wss.per.cluster <- apply(matrix.wss, 1, sum)
    bss.per.feature <- apply(matrix.bss, 2, sum)
    bss.per.cluster <- apply(matrix.bss, 1, sum)
    
    bss <- sum(bss.per.cluster)
    wss <- sum(wss.per.cluster)
    
    names(bss.per.feature) <- colnames(X)
    names(wss.per.feature) <- colnames(X)
    names(bss.per.cluster) <- 1:K
    names(wss.per.cluster) <- 1:K
    
    return(list(bss.per.feature = bss.per.feature, wss.per.feature = wss.per.feature, 
                bss.per.cluster = bss.per.cluster, 
                wss.per.cluster = wss.per.cluster, bss = bss, wss = wss))
}





