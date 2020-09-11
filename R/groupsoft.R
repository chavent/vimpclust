norm.vect <- function(x)
{
  res <- sqrt(sum(x^2))
  return(res)
}

soft.thresholding <- function(i.index, index, b, lambda, sizegroup = TRUE)
{
  vect.bcv.group <- b[index==i.index]
  if (length(vect.bcv.group)==1)
    {res <- sign(vect.bcv.group)*max(abs(vect.bcv.group)-lambda,0)} else
  {
    if (norm.vect(vect.bcv.group)>0)
    {
      if (sizegroup==T)
        res <- vect.bcv.group*max(norm.vect(vect.bcv.group)-sqrt(length(vect.bcv.group))*lambda,0)/norm.vect(vect.bcv.group)
      else
        res <- vect.bcv.group*max(norm.vect(vect.bcv.group)-lambda,0)/norm.vect(vect.bcv.group)
    }
   else res <- rep(0, length(vect.bcv.group))
  }
  return(unname(res)) 
}

#' @title Group soft-thresholding operator
#' @export
#'
#' @description 
#' This function implements the group soft-thresholding operator for a vector which elements are priorly split into groups. For the complete mathematical 
#' formulation, the reader may refer to the references below.
#'
#' @param b a numerical vector.
#' @param lambda a positive scalar containing the regularization parameter.
#' @param index a vector of integers of size \code{length(b)} containing the group membership for
#'  each element of \code{b}. By default, \code{index=1:length(b)} i.e. each element of \code{b} constitutes its own group. 
#' @param sizegroup a boolean. if TRUE, the size of 
#' the groups is taken into account in the thresholding operation.
#'
#' @return Returns the sparse vector after the group soft-thresholding operation.
#' 
#' @seealso \code{\link{groupsparsewkm}}
#'
#' @examples
#' b <- c(0.1, 0.2, 0.8, 0.1, 0.1, 0.3)
#' index <- c(1,1,2,2,3,3)
#' lambda <- 0.1
#' groupsoft(b=b, lambda=lambda, index=index, sizegroup=TRUE)
#' lambda <- 0.3
#' groupsoft(b=b, lambda=lambda, index=index, sizegroup=TRUE)
#' lambda <- 0.8
#' groupsoft(b=b, lambda=lambda, index=index, sizegroup=TRUE)
#' 
#' @references M. Chavent, J. Lacaille, A. Mourer and M. Olteanu (2020). 
#' Sparse k-means for mixed data via group-sparse clustering, to appear in ESANN proceedings.
#' @references M. Yuan and Y. Lin (2006). Model selection and estimation in regression with grouped variables. J. R. Statist. Soc. B, Vol. 68(1), p. 49-67.

groupsoft <- function(b, lambda, index = 1:length(b), sizegroup = TRUE) 
{
    check_fun_groupsoft(b, lambda, index, sizegroup)
    
    if (lambda==0)
      w <- b/norm.vect(b)
    else
    {
      w <- c(unlist(sapply(unique(sort(index)), soft.thresholding, index=index, b=b, lambda=lambda, 
                         sizegroup = sizegroup)))
      if (norm.vect(w)==0)
      w <- rep(0, length(b))
      else
      w <- w/norm.vect(w)
    }
    return(w)
    }
