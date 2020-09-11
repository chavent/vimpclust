#' @title check_fun_groupsparsw
#' @keywords internal
#'
check_fun_groupsparsw <- function(X, lambda, nlambda, index, sizegroup, 
                                  itermaxw, scaling, verbose) {
    
    check_X(X)
    check_lambda(lambda)
    check_nlambda(lambda, nlambda)
    check_index(index)
    check_sizegroup(sizegroup)
    check_itermaxw(itermaxw)
    check_scaling(scaling)
    check_verbose(verbose)
    # check_nthread(nthread)
}

check_fun_groupsoft <- function(b, s, index, sizegroup) {
    
    check_b(b)
    check_s(s)
    check_index(index)
    check_sizegroup(sizegroup)
}

check_fun_weightedss <- function(X, cl, w) {
    
    check_X(X)
    check_cl(cl)
    check_w(w)
}

check_X <- function(X) {
    # Missing values, Nan, Inf
    if (sum(is.finite(unlist(X))) != length(unlist(X))) 
        stop("Must provide data with only finite values.")
    X <- as.matrix(X)
    if (ncol(X) == 1) {
        stop("With only one variables, you should run standard k-means instead.")
    }
}

check_nlambda <- function(lambda, nlambda) {
    # nlambda
    if (is.null(lambda)) {
        if (!all.equal(nlambda, as.integer(nlambda), check.attributes = F)) 
            stop("nlambda must be an integer.")
        if (nlambda < 1) 
            stop("nlambda must strictly positive")
        if (nlambda > 10000) 
            stop("nlambda is >10000, one should ask why this precision is needed.")
    }
}

check_itermaxw <- function(itermaxw) {
    # itermaxw
    if (!all.equal(itermaxw, as.integer(itermaxw), check.attributes = F)) 
        stop("itermaxw must be an integer.")
    if (itermaxw < 1) 
        stop("itermaxw must strictly positive")
}

check_renamelevel <- function(renamelevel) {
    #  renamelevel
    if (!is.logical(renamelevel)) 
        stop("renamelevel must be a boolean.")
}


check_sizegroup <- function(sizegroup) {
    #  sizegroup
    if (!is.logical(sizegroup)) 
        stop("sizegroup must be a boolean.")
}

check_index <- function(index) {
    # index
    if (!all.equal(index, as.integer(index), check.attributes = F)) 
        stop("index must be a vector of integer.")
    if (max(index) > length(index)) 
        stop("An index of group is bigger than the number of variables. index must be in [1;ncol(X)].")
    if (min(index) < 1) 
        stop("An index of group is less than 1. index must be in [1;ncol(X)].")
}

check_lambda <- function(lambda) {
    # lambda
    if (!is.null(lambda)) 
        if (sum(lambda < 0) > 0) 
            stop("Lambda must be positive.")
}

check_cl <- function(cl) {
    # cl
    if (!all.equal(cl, as.integer(cl), check.attributes = F)) 
        stop("cl must be a vector of integer.")
    if (max(cl) > length(cl)) 
        stop("A number of group is bigger than the length of cl.")
}

check_s <- function(s) {
    # s
    if (!all.equal(s, as.numeric(s), check.attributes = F)) 
        stop("s must be numerical.")
    if (s < 0) 
        stop("s must be positive.")
}

check_b <- function(b) {
    # b
    if (!all.equal(b, as.numeric(b), check.attributes = F)) 
        stop("b must be numerical.")
}

check_maxlambda <- function(lambda, maxlam) {
    # lambda
    if (!is.null(lambda)) 
        if (sum(lambda > maxlam) > 0) 
            stop(paste0("Lambda must be in [0;", maxlam, "]"))
}

check_w <- function(w) {
    # w
    if (!is.null(w)) {
        if (!all.equal(w, as.numeric(w), check.attributes = F)) 
            stop("w must be numerical.")
        if (sum(w < 0) > 0) 
            stop("w must be positive.")
    }
}

check_nthread <- function(nthread) {
    # nthread
    if (!is.null(nthread)) {
        if (!all.equal(nthread, as.integer(nthread), check.attributes = F)) 
            stop("nthread must be an integer.")
        if (nthread < 1) 
            stop("nthread must strictly positive")
    }
}

check_verbose <- function(verbose) {
    # verbose
    if (!all.equal(verbose, as.integer(verbose), check.attributes = F)) 
        stop("verbose must be an integer.")
    if (verbose < 0) 
        stop("verbose must be 0,1,2 or 3.")
    if (verbose > 3) 
        stop("verbose must be 0,1,2 or 3.")
}

check_scaling <- function(scaling) {
    # scaling
    if (!is.logical(scaling)) 
        stop("scaling must be a boolean.")
}
