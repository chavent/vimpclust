#' @title Group-sparse weighted k-means
#' @export
#'
#' @description This function performs group-sparse weighted k-means on a set 
#' of observations described by numerical variables organized in groups. 
#' It generalizes the sparse clustering algorithm introduced by 
#' Witten & Tibshirani (2010) to groups. While the algorithm clusters the observations, the groups of variables are supposed priorly known.
#'  The algorithm computes a series of weights associated to the groups of variables, the weights
#' indicating the importance of each group in the clustering process. 
#' 
#' @param X a numerical matrix or a dataframe of dimension \code{n} (observations) by \code{p} 
#' (variables).
#' @param centers an integer representing the number of clusters.
#' @param lambda a vector of numerical values (or a single value) providing 
#' a grid of values for the regularization parameter. If NULL (by default), the function computes its 
#' own lambda sequence of length \code{nlambda} (see details).
#' @param nlambda an integer indicating the number of values for the regularization parameter. 
#' By default, \code{nlambda=20}.
#' @param index a vector of integers of size \code{p} providing the group membership
#'  for each variable. By default, \code{index=1:ncol(X)} i.e. no groups or groups of size 1. 
#' @param sizegroup a boolean. If TRUE, the group sizes (number of variables in each group) are taken into account in the penalty term (see details).
#'  By default, \code{sizegroup=TRUE}.
#' @param nstart an integer representing the number of random starts in the k-means algorithm.
#'  By default, \code{nstart=10}. 
#' @param itermaxw an integer indicating the maximum number of iterations for the inside 
#' loop over the weights \code{w}. By default, \code{itermaxw=20}.
#' @param itermaxkm an integer representing the maximum number of iterations in the k-means 
#' algorithm. By default, \code{itermaxkm=10}.
#' @param scaling a boolean. If TRUE, variables are scaled to zero mean and unit variance. By default, \code{scaling=TRUE}.
#' @param epsilonw a positive numerical value. It provides the precision of the stopping criterion over \code{w}. By default, \code{epsilonw =1e-04}.
#' @param verbose an integer value. If \code{verbose=0}, the function stays silent, if \code{verbose=1} (default option), it  prints
#'  whether the stopping criterion over the weights \code{w} is satisfied.
#' 
#' @return \item{lambda}{a numerical vector containing the regularization parameters (a grid of values).}
#' @return \item{W}{a \code{p} by \code{length(lambda)} numerical matrix. It contains the weights associated to each variable.}
#' @return \item{Wg}{a \code{L} by \code{length(lambda)} numerical matrix, where \code{L} is the number of groups. It contains the weights associated to each group.}
#' @return \item{cluster}{a \code{n} by \code{length(lambda)} integer matrix. It contains the cluster memberships, for each value of the regularization parameter.}
#' @return \item{sel.feat}{a numerical vector of the same length as \code{lambda}, giving the number of selected variables for each value of the regularization parameter.}
#' @return \item{sel.groups}{a numerical vector of the same length as \code{lambda}, giving the number of selected groups of variables for each value of the regularization parameter.}
#' @return \item{Z}{a matrix of size \code{n} by \code{p} containing the scaled data if \code{scaling=TRUE}, and a copy of \code{X} otherwise.}
#' @return \item{bss.per.feature}{a matrix of size \code{p} by \code{length(lambda)}. It contains the between-class variance computed for each variable.}
#' 
#' 
#' @details 
#' Group-sparse weighted k-means performs clustering on data described by numerical variables priorly partitionned into groups, and automatically selects the most discriminant groups by 
#' setting to zero the weights of the non-discriminant ones. 
#' 
#' The algorithm is based on the optimization of a cost function which is the weighted between-class variance penalized by a group L1-norm. The groups must be priorly defined through 
#' expert knowledge. If there is no group structure (each group contains one variable only), the algorithm reduces to the sparse weighted k-means introduced in Witten & Tibshirani (2010).
#' The penalty term may take into account the size of the groups by setting \code{sizegroup=TRUE} (see Chavent et al. (2020) for further details on the mathematical expression of the
#' optimized criterion). The importance of the penalty term may be adjusted through the regularization parameter \code{lambda}. If \code{lambda=0}, there is no penalty applied to the 
#' weighted between-class variance. The larger \code{lambda}, the larger the penalty term and the number of groups with null weights. 
#' 
#' The output of the algorithm is three-folded: one gets a partitioning of the data, a vector of weights associated to each group, and a vector of weights associated to each
#' variable. Weights equal to zero imply that the associated variables or the associated groups do not participate in the clustering process. 
#' 
#' Since it is difficult to chose the regularization parameter \code{lambda} without prior knowledge, the function builds automatically a grid of parameters and finds the partitioning
#' and the vectors of weights associated to each value in the grid. 
#' 
#' Note that when the regularization parameter is equal to 0 (no penalty applied), the output is different from that of a regular k-means, since the optimized criterion is a weighted 
#' between-class variance and not the between-class variance only. 
#'
#' @references Witten, D. M., & Tibshirani, R. (2010). A framework for feature 
#' selection in clustering. Journal of the American Statistical Association, 
#' 105(490), p.713-726.
#' @references Chavent, M. & Lacaille, J. & Mourer, A. & Olteanu, M. (2020). 
#' Sparse k-means for mixed data via group-sparse clustering, ESANN proceedings.
#' 
#' @seealso \code{\link{plot.spwkm}}, \code{\link{info_clust}}
#' 
#'@examples
#' data(iris)
#' # define two groups of variables: 
#' # "Sepal.Length" and "Sepal.Width" in group 1
#' # "Petal.Length" and "Petal.Width"  in group 2
#' index <- c(1, 2, 1, 2)
#' # group-sparse k-means
#' out <- groupsparsewkm(X = iris[,-5], centers = 3, index = index)
#' # grid of regularization parameters
#' out$lambda
#' k <- 10
#' # weights of the variables for the k-th regularization parameter
#' out$W[,k]
#' # weights of the groups for the k-th regularization parameter
#' out$Wg[,k]
#' # partition obtained with for the k-th regularization parameter
#' out$cluster[,k]
#' # between-class variance on each variable
#' out$bss.per.feature[,k]
#' # between-class variance 
#' sum(out$bss.per.feature[,k])/length(index)
#' 
#' # one variable per group (equivalent to sparse k-means)
#' index <- 1:4 # default option in groupsparsewkm
#' # sparse k-means
#' out <- groupsparsewkm(X = iris[,-5], centers = 3, index = index)
#' # or
#' out <- groupsparsewkm(X = iris[,-5], centers = 3)
#' # group weights and variable weights are identical in this case
#' out$Wg 
#' out$W

 

groupsparsewkm <- function(X, centers, lambda = NULL, nlambda = 20, 
                           index = 1:ncol(X), sizegroup = TRUE,  
                           nstart = 10, itermaxw = 20, itermaxkm = 10, 
                           scaling = TRUE, verbose = 1, epsilonw = 1e-04) 
{
    #options(digits = 3)
    call <- match.call()
    check_fun_groupsparsw(X, lambda, nlambda, index, sizegroup, itermaxw, scaling, verbose)
    
    X <- as.matrix(scale(X, center = scaling, scale = scaling))
    
    X <- X[,order(index)]
    index <- sort(index)
    
    km <- stats::kmeans(X, centers, nstart, itermaxkm)
    clusterini <- km$cluster
    bss.per.featureini <- weightedss(X, clusterini)$bss.per.feature/nrow(X)
    maxlam <-  max(stats::aggregate(bss.per.featureini, by=list(index), function(x){norm.vect(x)/sqrt(length(x))})[,2])
    check_maxlambda(lambda, maxlam)
   
    
    if (is.null(lambda)) 
        lambda <- seq(from = 0, to = (maxlam - 0.001), length.out = nlambda) 
    
    # loop over lambda
    ls <- sapply(lambda, FUN = function(i) {
        
        # loop over w
        cluster1 <- clusterini
        bss.per.feature <- bss.per.featureini
        w1 <- rep(1/ncol(X), ncol(X))  # initialization of w1 / w1 gives the weight at the i+1 iterations
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
            # Soft thresholding operator
            w1 <- groupsoft(bss.per.feature, i, index, sizegroup)
            if (norm.vect(w1)==0)
              niter <- itermaxw
            niter = niter + 1
        }  # end loop over w
        if (verbose == 1) 
            if (niter >= itermaxw) 
                print(paste0("the stopping criterion over w is not satisfied for lambda = ", i)) else
                               print(paste0("the stopping criterion over w is satisfied for lambda = ", i))  
        
        
        W_cl[[1]] <- w0
        W_cl[[2]] <- cluster0
        W_cl[[3]] <- bss.per.feature
        W_cl
    })  # end loop over lambda
    
    W <- do.call(cbind, (ls[1, ]))
    cluster <- do.call(cbind, (ls[2, ]))
    bss.per.feature <- do.call(cbind, (ls[3, ]))
    
    rownames(W) <- colnames(X)
    colnames(W) <- paste("Lambda", 1:length(lambda))
    
    Wg <- stats::aggregate(W, by=list(sort(index)), norm.vect)[,-1]
    Wg <- as.matrix(Wg)
    rownames(Wg) <- paste("Group", unique(sort(index)))
    colnames(Wg) <- paste("Lambda", 1:length(lambda))
    
    selected.features <- apply(W,2,function(x){sum(x>0)})
    selected.groups <- apply(Wg,2,function(x){sum(x>0)})
    
    rownames(bss.per.feature) <- colnames(X)
    colnames(bss.per.feature) <- paste("Lambda", 1:length(lambda))
    
    rownames(cluster) <- rownames(X)
    colnames(cluster)  <- paste("Lambda", 1:length(lambda))
    
    out <- list(call = call, type="NumGroupSparse", lambda = lambda, W = W, Wg = Wg, 
                 cluster = cluster,   sel.feat= selected.features, 
                sel.groups=selected.groups,
                 Z = X, index = index, bss.per.feature = bss.per.feature)
    class(out) <- "spwkm"
    return(out)
}
