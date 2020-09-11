#'@title Recoding mixed data
#'@export
#'
#'@description This function transforms and scales a dataset with numerical and/or categorical
#'variables. Numerical variables are scaled to zero mean and unit variance. Categorical variables
#'are first transformed into dummy variables according to their levels, and second centered and 
#'normalized with respect to the square roots of the relative frequencies of the levels. The complete 
#'procedure is described in Chavent et al. (2014).


#'@param X a matrix or a dataframe with numerical and/or categorical variables. Categorical variables 
#'must be given as factors.
#'@param renamelevel a boolean. If TRUE (default value), the levels of the categorical variables
#'are renamed as \code{'variable_name=level_name'}.

#'@return \item{X}{a data frame or a matrix. The input data \code{X} with reordered columns 
#'(numerical first, categorical second).}
#'@return \item{Z}{a data frame. The transformed data matrix with scaled numerical variables
#' and scaled dummy variables coding for the levels.}
#'@return \item{index}{a vector of integers. Contains an implicit partitioning of the transformed
#' variables: each scaled numerical variable represents a group, all scaled dummy variables 
#' summarizing the levels of a categorical variable represent a group. \code{index} allows to 
#' preserve the information on the initial structure of the data, particularly for categorical variables.}
#'
#'@examples
#'head(HDdata)
#'out <- recodmix(HDdata[,-14], renamelevel=TRUE)
#'# reordered data (numerical/categorical)
#'colnames(out$X)
#'# transformed and scaled data
#'colnames(out$Z)
#'# transformed variables partitioning and group membership
#'out$index

#' @references M. Chavent, V. Kuentz-Simonet, A. Labenne and J. Saracco (2014).
#' Multivariate analysis of mixed data: the PCAmixdata R package, arXiv:1411.4911.


recodmix <- function(X, renamelevel = FALSE) {
    split <- PCAmixdata::splitmix(X)
    X.quanti <- split$X.quanti
    X.quali <- split$X.quali
    rec <- PCAmixdata::recod(X.quanti, X.quali, rename.level = renamelevel)
    X <- rec$X
    Z <- rec$Z
    index <- rec$indexj
    return(list(X = X, Z = Z, index = index))
}


