
#' Cholesky of Random Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according to the Cholesky
#'     decomposition of a Wishart distribution with parameters \code{Sigma} and
#'     \code{df}, \eqn{W_p(Sigma, df)}.
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{(p * p)} "scale" matrix, the matrix parameter of the distribution.
#'
#' @return a numeric array, say R, of dimension \eqn{p * p * n}, where each \code{R[,,i]} is a Cholesky decomposition of a realization of the Wishart distribution \eqn{W_p(Sigma, df)}. Based on a modification of the existing code for the \code{rWishart} function
#'
#' @seealso \code{\link{rWishart}}, \code{\link{rInvCholWishart}}
#' @useDynLib CholWishart
#' @export
#'
#' @examples
#' # How it is parameterized:
#' set.seed(20180211)
#' A <- rCholWishart(1,10,3*diag(5))[,,1]
#'
#' set.seed(20180211)
#' B <- rInvCholWishart(1,10,3*diag(5))[,,1]
#'
#' A %*% B
#'
#' set.seed(20180211)
#' C <- chol(rWishart(1,10,3*diag(5))[,,1])
#'
rCholWishart <- function(n, df, Sigma){
  if (!is.numeric(Sigma))
    stop("Sigma must be numeric.")
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("N must be 1 or larger.")
  if (!is.matrix(Sigma)  || dims[1] != dims[2])
    stop("'Sigma' must be a square, real matrix");
  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0 )
    stop("inconsistent degrees of freedom and dimension")
  .Call("rCholWishart", n, df, Sigma, PACKAGE = "CholWishart")
}



#' Cholesky Factor of Random Inverse Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according
#'    to the Cholesky factor of an inverse Wishart distribution with
#'    parameters \code{Sigma} and \code{df}, \eqn{W_p(Sigma, df)}.
#'    Note there are different ways of parameterizing the Inverse
#'    Wishart distribution, so check which one you need.
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{(p * p)} "scale" matrix, the matrix parameter of the distribution.
#'
#' @return a numeric array, say R, of dimension \eqn{p * p * n}, where each \code{R[,,i]} is a Cholesky decomposition of a realization of the Wishart distribution \eqn{W_p(Sigma, df)}. Based on a modification of the existing code for the \code{rWishart} function
#'
#' @seealso \code{\link{rWishart}} and \code{\link{rCholWishart}}
#' @useDynLib CholWishart
#' @export
#'
#' @examples
#'
#' rInvCholWishart(1,10,diag(5))
#'
rInvCholWishart <- function(n, df, Sigma){
  if (!is.numeric(Sigma))
    stop("Sigma must be numeric.")
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("N must be 1 or larger.")
  if (!is.matrix(Sigma)  || dims[1] != dims[2])
    stop("'Sigma' must be a square, real matrix");
  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0 )
    stop("inconsistent degrees of freedom and dimension")
  .Call("rInvCholWishart", n, df, Sigma, PACKAGE = "CholWishart")
}

#' Random Inverse Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according
#'     to the inverse Wishart distribution with parameters \code{Sigma} and
#'     \code{df}, \eqn{W_p(Sigma, df)}. Note there are different ways
#'     of parameterizing this distribution, so check which one you need.
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{(p * p)} "scale" matrix, the matrix parameter of the distribution.
#'
#' @return a numeric array, say R, of dimension \eqn{p * p * n}, where each \code{R[,,i]} is a realization of the inverse Wishart distribution \eqn{IW_p(Sigma, df)}. Based on a modification of the existing code for the \code{rWishart} function. If \eqn{X \sim IW_p(Sigma, df)} then \eqn{X^{-1} \sim W_p(Sigma^{-1}, df)}
#'
#' @seealso \code{\link{rWishart}} and \code{\link{rCholWishart}}
#' @useDynLib CholWishart
#' @export
#'
#' @examples
#' set.seed(20180221)
#' A<-rInvWishart(1,10,5*diag(5))[,,1]
#' set.seed(20180221)
#' B<-rWishart(1,10,.2*diag(5))[,,1]
#'
#' A %*% B
rInvWishart <- function(n, df, Sigma){
  if (!is.numeric(Sigma))
    stop("Sigma must be numeric.")
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("N must be 1 or larger.")
  if (!is.matrix(Sigma)  || dims[1] != dims[2])
    stop("'Sigma' must be a square, real matrix");
  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0 )
    stop("inconsistent degrees of freedom and dimension")
  .Call("rInvWishart", n, df, Sigma, PACKAGE = "CholWishart")
}


#' Density for Random Wishart Distributed Matrices
#'
#' Compute the density of an observation of a random Wishart distributed matrix (\code{dWishart})
#' or an observation
#' from the inverse Wishart distribution (\code{dInvWishart}).
#'
#' If \eqn{X} is distributed as a \eqn{p * p} Wishart random variable with \eqn{n > p} degrees of
#' freedom and a covariance matrix \eqn{Sigma}, then \eqn{X^{-1} = Y} is distributed as an
#' inverse Wishart with \eqn{n} degrees of freedom and a covarariance matrix \eqn{Sigma^{-1}}. Note there are different ways of parameterizing the
#'    inverse Wishart distribution, check which one you need.
#'
#' @param x positive definite \eqn{p * p} observation for density estimation
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{(p * p} "scale" matrix, the matrix parameter of the distribution.
#' @param log logical, whether to return value on the log scale.
#'
#' @return Density or log of density
#' @export
#'
#' @examples
#' set.seed(20180222)
#' A <- rWishart(1,10,diag(4))[,,1]
#' A
#' dWishart(x = A, df = 10,Sigma = diag(4), log=TRUE)
#' dInvWishart(x = solve(A), df = 10,Sigma = diag(4), log=TRUE)
dWishart <- function(x, df, Sigma, log = TRUE){
  if (!is.numeric(Sigma))
    stop("Sigma must be numeric.")
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (!is.numeric(x))
    stop("x must be numeric.")
  if(dim(x)[1] != dim(x)[2] || dim(x)[1] != dims[1] || dims[1] != dims[2])
    stop("Non-conformable dimensions")
  if(!symm.check(x) || !symm.check(Sigma))
    stop("Non-symmetric input.")
  cholX <- chol(x)
  cholS <- chol(Sigma)
  ldetX <- sum(log(diag(cholX)))*2
  ldetS <- sum(log(diag(cholS)))*2
  ldensity <- .5*(df - dims[1] - 1)*ldetX - .5 * sum(diag(chol2inv(cholS) %*% x)) -
    (df * dims[1]/2 *log(2) ) - .5 * df * ldetS - lmvgamma(df/2, dims[1])
  if(log) return(ldensity)
  else return(exp(ldensity))
}

#' @describeIn dWishart density for the inverse Wishart distribution.
#' @export
dInvWishart <- function(x, df, Sigma, log = TRUE){
  if (!is.numeric(Sigma))
    stop("Sigma must be numeric.")
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (!is.numeric(x))
    stop("x must be numeric.")
  if(dim(x)[1] != dim(x)[2] || dim(x)[1] != dims[1] || dims[1] != dims[2])
    stop("Non-conformable dimensions")
  if(!symm.check(x) || !symm.check(Sigma))
    stop("Non-symmetric input.")

  cholX <- chol(x)
  cholS <- chol(Sigma)
  ldetX <- sum(log(diag(cholX)))*2
  ldetS <- sum(log(diag(cholS)))*2
  ldensity <- -.5*(df + dims[1] + 1)*ldetX + .5 * df * ldetS - .5 * sum(diag(chol2inv(cholX) %*% Sigma)) -
    (df * dims[1]/2 *log(2) ) - lmvgamma(df/2, dims[1])
  if(log) return(ldensity)
  else return(exp(ldensity))

}


#' Multivariate Gamma Function
#'
#' @description A special mathematical function related to the gamma function,
#'     generalized for multivariate gammas.
#'
#' @param x non-negative numeric vector, matrix, or array
#' @param p positive integer, dimension of a square matrix
#'
#' @return log of multivariate gamma for each entry of x. For non-log variant,
#'     use \code{mvgamma}.
#'
#' @seealso \code{\link{gamma}} and \code{\link{lgamma}}
#'
#' @useDynLib CholWishart
#' @export
#'
#' @examples
#' lgamma(1:12)
#' lmvgamma(1:12,1)
#' mvgamma(1:12,1)
#' gamma(1:12)
lmvgamma <- function(x, p) {
  # p only makes sense as an integer but not testing that. x *could* be
  # less than zero - same domain as gamma function making sure that object
  # returned is same shape as object passed
  if (!all(is.numeric(x),is.numeric(p))) stop("Non-numeric input.")
  dims <- if (is.vector(x))
    length(x) else dim(as.array(x))
  if (p < 1)
    stop("p must be greater than or equal to 1. p = ", p)
  if (any(x <= 0))
    stop("x must be greater than 0. x = ", x)

  result <- .Call("lmvgamma", as.numeric(x), as.integer(p), PACKAGE = "CholWishart")

  return(array(result, dim = dims))
}

#' @describeIn lmvgamma Multivariate gamma function.
#' @export
mvgamma <- function(x, p) exp(lmvgamma(x, p))

#' Multivariate Digamma
#'
#' @description A special mathematical function related to the gamma function,
#'     generalized for multivariate distributions.
#'     The digamma is the derivative of the gamma.
#'
#' @param x non-negative numeric vector, matrix, or array
#' @param p positive integer, dimension of a square matrix
#' @return vector of values of multivariate digamma function.
#' @export
#'
#' @examples
#' digamma(1:10)
#' mvdigamma(1:10,1)
mvdigamma <- function(x,p){
  if (!all(is.numeric(x),is.numeric(p))) stop("Non-numeric input.")
  dims <- if (is.vector(x))
    length(x) else dim(as.array(x))
  if (p < 1)
    stop("p must be greater than or equal to 1. p = ", p)
  if (any(x <= 0))
    stop("x must be greater than 0. x = ", x)

  .Call("mvdigamma",as.numeric(x),as.integer(p),PACKAGE = "CholWishart")

}




#' Quick symmetry check
#'
#' Quick check whether matrix input is symmetric -
#' checks sum of absolute differences of transposes
#' as well as dimensions. Not robust, so only an
#' internal function to be used with known safe input.
#'
#' @param A Numeric real matrix. Does not check if real.
#' @param tol tolerance - note that if you have a big matrix
#'    it may need to be specified as it's a sum of entries.
#'
#' @return logical TRUE if symmetric FALSE otherwise.
#' Not as robust as \code{isSymmetric()}.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' A = ARgenerate(5,.9)
#' symm.check(A)
#' A[1,2] = 5
#' symm.check(A)}
symm.check <- function(A, tol = 10 * (.Machine$double.eps)^.5) {
  # if (!is.matrix(A)) return(FALSE)
  # if (!is.numeric(A)) return(FALSE)
  # commented those out because it is always checked before running symm.check.
  dims <- dim(A)
  if (dims[1] != dims[2]) {
    return(FALSE)
  }
  return(sum(abs(A - t(A))) < prod(dims)*tol)
}




.onUnload <- function(libpath) {
  library.dynam.unload("CholWishart", libpath)
}
