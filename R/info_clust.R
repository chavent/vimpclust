#' @title Description of a set of partitions
#' @export
#'
#' @description This function computes descriptive statistics of the clustering produced with 
#' group-sparse weighted k-means on numerical data, or with sparse weighted k-means on mixed data.
#' It displays the average of the numerical variables per cluster, and the relative frequencies 
#' of the levels in the categorical variables per cluster. 
#'
#' @param out an object of class \code{spwkm}.
#' @param which.lambda an integer or a vector of integers
#' selecting the clusterings for which summaries are computed.
#' @param X a matrix or a data frame. The initial data set.
#' 
#' @details 
#' The values in \code{which.lambda} must be integers between 1 and \code{length(out$lambda)}. One may
#' thus select the clusterings corresponding to specific regularization parameters, or the whole set
#' of clusterings
#' obtained for the whole grid of \code{out$lambda}.
#' 
#' @return \item{mean.by.clust}{a list of numerical matrices. Each matrix contains the mean values
#' of the numerical variables computed per cluster, for a given value of the regularization parameter.} 
#' @return \item{freq.by.clust}{a list of numerical matrices. Each matrix contains the relative 
#' frequencies of each level associated to categorical variables, 
#' computed per cluster and for a given value of the regularization parameter.}
#' @return \item{lambda}{a scalar or a numerical vector. The selected values of the regularization
#' parameter. }
#'
#' @seealso \code{\link{groupsparsewkm}}, \code{\link{sparsewkm}}
#'
#' @examples
#' data(HDdata)
#' out <- sparsewkm(X = HDdata[,-14], centers = 2)
#' info_clust(out, which.lambda=c(1,10,20), X = HDdata[,-14])


info_clust <- function(out, which.lambda, X) 
{
  varmean_ls <- list()
  freq_ls <- list()
  
  res <- PCAmixdata::splitmix(X)
  X.quanti <- res$X.quanti
  X.quali <- res$X.quali
  rec <- PCAmixdata::recod(X.quanti, X.quali, rename.level = TRUE)
  
  ncl <-  length(unique(as.vector(out$cluster)))
  
  for(i in 1:length(which.lambda)){
    if (!is.null(X.quanti)) 
    {
      varmean_byclust <- stats::aggregate(X.quanti, list(out$cluster[, which.lambda[i]]), mean)
      varmean_byclust <- varmean_byclust[,-1]
      varmean_byclust <- t(varmean_byclust)
      colnames(varmean_byclust) <- paste("cluster=", 1:ncl, sep = "")
      varmean_ls[[i]] <- varmean_byclust
    }
    if (!is.null(X.quali)) 
    {
      freq_byclust <- stats::aggregate(rec$G, list(out$cluster[, which.lambda[i]]), mean)
      freq_byclust <- freq_byclust[,-1]
      freq_byclust <- t(freq_byclust)
      colnames(freq_byclust) <- paste("cluster=", 1:ncl, sep = "")
      freq_ls[[i]] <- freq_byclust
    } 
  }
  
  names(varmean_ls) <- paste("lambda = ", out$lambda[which.lambda], sep = "")
  names(freq_ls) <- paste("lambda = ", out$lambda[which.lambda], sep = "")
  lambda <- out$lambda[which.lambda]
  
  return(list(mean.by.clust = varmean_ls, freq.by.clust = freq_ls, lambda=lambda))
}
