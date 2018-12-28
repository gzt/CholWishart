#   wishart.R
#   CholWishart: Sample the Cholesky Factor of the Wishart and Other Functions
#   Copyright (C) 2018  GZ Thompson <gzthompson@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#   along with this program; if not, a copy is available at
#   https://www.R-project.org/Licenses/
#

#' Cholesky Factor of Random Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according
#'     to the Cholesky factorization of a Wishart distribution with
#'     parameters \code{Sigma} and \code{df}, \eqn{W_p(Sigma, df)}
#'     (known as the Bartlett decomposition
#'     in the context of Wishart random matrices).
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{p \times p}{(p * p)} "scale" matrix, the matrix parameter of the
#' distribution.
#'
#' @return a numeric array, say \code{R}, of dimension \eqn{p \times p \times n}{p * p * n},
#'    where each \code{R[,,i]} is a Cholesky decomposition of a sample
#'    from the Wishart distribution \eqn{W_p(Sigma, df)}. Based on a
#'    modification of the existing code for the \code{rWishart} function.
#'
#' @seealso \code{\link{rWishart}}, \code{\link{rInvCholWishart}}
#'
#' @references
#' Anderson, T. W. (2003). \emph{An Introduction to Multivariate Statistical Analysis} (3rd ed.).
#' Hoboken, N. J.: Wiley Interscience.
#'
#' Mardia, K. V., J. T. Kent, and J. M. Bibby (1979) \emph{Multivariate Analysis},
#' London: Academic Press.
#'
#' A. K. Gupta and D. K. Nagar 1999. \emph{Matrix variate distributions}. Chapman and Hall.
#' @useDynLib CholWishart, .registration = TRUE
#' @export
#'
#' @examples
#' # How it is parameterized:
#' set.seed(20180211)
#' A <- rCholWishart(1,10,3*diag(5))[,,1]
#' A
#' set.seed(20180211)
#' B <- rInvCholWishart(1,10,1/3*diag(5))[,,1]
#' B
#' crossprod(A) %*% crossprod(B)
#'
#' set.seed(20180211)
#' C <- chol(stats::rWishart(1,10,3*diag(5))[,,1])
#' C
rCholWishart <- function(n, df, Sigma) {
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("'n' must be 1 or larger.")

  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0)
    stop("inconsistent degrees of freedom and dimension")
  .Call("C_rCholWishart", n, df, Sigma, PACKAGE = "CholWishart")
}



#' Cholesky Factor of Random Inverse Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according
#'    to the Cholesky factor of an inverse Wishart distribution with
#'    parameters \code{Sigma} and \code{df}, \eqn{W_p(Sigma, df)}.
#'
#'    Note there are different ways of parameterizing the Inverse
#'    Wishart distribution, so check which one you need.
#'     Here,  If \eqn{X \sim IW_p(\Sigma, \nu)}{X ~ IW_p(Sigma, df)} then
#'     \eqn{X^{-1} \sim W_p(\Sigma^{-1}, \nu)}{X^{-1} ~ W_p(Sigma^{-1}, df)}.
#'     Dawid (1981) has a different definition: if
#'     \eqn{X \sim W_p(\Sigma^{-1}, \nu)}{X ~ W_p(Sigma^{-1}, df)} and
#'     \eqn{\nu > p - 1}{df > p - 1}, then
#'     \eqn{X^{-1} = Y \sim IW(\Sigma, \delta)}{X^{-1} = Y ~ IW(Sigma, delta)}, where
#'     \eqn{\delta = \nu - p + 1}{delta = df - p + 1}.
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{p \times p}{(p * p)} "scale" matrix, the matrix parameter of
#' the distribution.
#'
#' @return a numeric array, say \code{R}, of dimension \eqn{p \times p \times n}{p * p * n},
#' where each \code{R[,,i]} is a Cholesky decomposition of a realization of the Wishart distribution
#' \eqn{W_p(Sigma, df)}. Based on a modification of the existing code for the \code{rWishart} function
#'
#' @seealso \code{\link{rWishart}} and \code{\link{rCholWishart}}
#' @references
#' Anderson, T. W. (2003). \emph{An Introduction to Multivariate Statistical Analysis} (3rd ed.).
#' Hoboken, N. J.: Wiley Interscience.
#'
#' Dawid, A. (1981). Some Matrix-Variate Distribution Theory: Notational Considerations and a
#' Bayesian Application. \emph{Biometrika}, 68(1), 265-274. \doi{10.2307/2335827}
#'
#' Gupta, A. K.  and D. K. Nagar (1999). \emph{Matrix variate distributions}. Chapman and Hall.
#'
#' Mardia, K. V., J. T. Kent, and J. M. Bibby (1979) \emph{Multivariate Analysis},
#' London: Academic Press.
#' @export
#'
#' @examples
#' # How it is parameterized:
#' set.seed(20180211)
#' A <- rCholWishart(1,10,3*diag(5))[,,1]
#' A
#' set.seed(20180211)
#' B <- rInvCholWishart(1,10,1/3*diag(5))[,,1]
#' B
#' crossprod(A) %*% crossprod(B)
#'
#' set.seed(20180211)
#' C <- chol(stats::rWishart(1,10,3*diag(5))[,,1])
#' C
rInvCholWishart <- function(n, df, Sigma) {
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("'n' must be 1 or larger.")

  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0)
    stop("inconsistent degrees of freedom and dimension")
  .Call("C_rInvCholWishart", n, df, Sigma, PACKAGE = "CholWishart")
}

#' Random Inverse Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according
#'     to the inverse Wishart distribution with parameters \code{Sigma} and
#'     \code{df}, \eqn{W_p(Sigma, df)}.
#'
#'    Note there are different ways of parameterizing the Inverse
#'    Wishart distribution, so check which one you need.
#'     Here,  If \eqn{X \sim IW_p(\Sigma, \nu)}{X ~ IW_p(Sigma, df)} then
#'     \eqn{X^{-1} \sim W_p(\Sigma^{-1}, \nu)}{X^{-1} ~ W_p(Sigma^{-1}, df)}.
#'     Dawid (1981) has a different definition: if
#'     \eqn{X \sim W_p(\Sigma^{-1}, \nu)}{X ~ W_p(Sigma^{-1}, df)} and
#'     \eqn{\nu > p - 1}{df > p - 1}, then
#'     \eqn{X^{-1} = Y \sim IW(\Sigma, \delta)}{X^{-1} = Y ~ IW(Sigma, delta)}, where
#'     \eqn{\delta = \nu - p + 1}{delta = df - p + 1}.
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{p \times p}{(p * p)} "scale" matrix, the matrix parameter of the
#' distribution.
#'
#' @return a numeric array, say \code{R}, of dimension \eqn{p \times p \times n}{p * p * n},
#' where each \code{R[,,i]} is a realization of the inverse Wishart distribution \eqn{IW_p(Sigma, df)}.
#' Based on a modification of the existing code for the \code{rWishart} function.
#'
#' @seealso \code{\link{rWishart}}, \code{\link{rCholWishart}}, and \code{\link{rInvCholWishart}}
#'
#' @references
#' Dawid, A. (1981). Some Matrix-Variate Distribution Theory: Notational Considerations and a
#' Bayesian Application. \emph{Biometrika}, 68(1), 265-274. \doi{10.2307/2335827}
#'
#' Gupta, A. K.  and D. K. Nagar (1999). \emph{Matrix variate distributions}. Chapman and Hall.
#'
#' Mardia, K. V., J. T. Kent, and J. M. Bibby (1979) \emph{Multivariate Analysis},
#' London: Academic Press.

#' @export
#'
#' @examples
#' set.seed(20180221)
#' A<-rInvWishart(1,10,5*diag(5))[,,1]
#' set.seed(20180221)
#' B<-stats::rWishart(1,10,.2*diag(5))[,,1]
#'
#' A %*% B
rInvWishart <- function(n, df, Sigma) {
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("'n' must be 1 or larger.")

  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0)
    stop("inconsistent degrees of freedom and dimension")
  .Call("C_rInvWishart", n, df, Sigma, PACKAGE = "CholWishart")
}


#' Random Pseudo Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according
#'     to the pseudo Wishart distribution with parameters \code{Sigma} and
#'     \code{df}, \eqn{W_p(Sigma, df)}, with sample size \code{df} less than
#'     the dimension \code{p}. This can work with a singular \code{Sigma}.
#'
#'     Let \eqn{X_i}, \eqn{i = 1, 2, ..., df} be \code{df}
#'     observations of a multivariate normal distribution with mean 0 and
#'     covariance \code{Sigma}. Then \eqn{\sum X_i X_i'} is distributed as a pseudo
#'     Wishart \eqn{W_p(Sigma, df)}. Sometimes this is called a singular Wishart
#'     distribution, however, that can be confused with the case where \code{Sigma}
#'     itself is singular. This is able to work with positive semi-definite
#'     \code{Sigma}, meaning this can generate samples from singular pseudo Wishart
#'     distributions.
#'
#' @param n integer sample size.
#' @param df integer parameter, "degrees of freedom", should be less than the
#'    dimension of \code{p}
#' @param Sigma positive semi-definite \eqn{p \times p}{(p * p)} "scale" matrix, the
#' matrix
#'    parameter of the distribution.
#'
#' @return a numeric array, say \code{R}, of dimension \eqn{p \times p \times n}{p * p * n},
#'     where each \code{R[,,i]} is a realization of the pseudo Wishart distribution
#'     \eqn{W_p(Sigma, df)}.
#'
#' @seealso \code{\link{rWishart}}, \code{\link{rCholWishart}}, \code{\link{rInvWishart}}
#'     and \code{\link{rInvCholWishart}}
#'
#' @references
#' Diaz-Garcia, Jose A, Ramon Gutierrez Jaimez, and Kanti V Mardia. 1997.
#' “Wishart and Pseudo-Wishart Distributions and Some Applications to Shape Theory.”
#' Journal of Multivariate Analysis 63 (1): 73–87. \doi{10.1006/jmva.1997.1689}.
#'
#' Uhlig, Harald. On Singular Wishart and Singular Multivariate Beta Distributions.
#' Ann. Statist. 22 (1994), no. 1, 395--405. \doi{10.1214/aos/1176325375}.
#'
#' @export
#'
#' @examples
#' set.seed(20181227)
#' A<-rPseudoWishart(1,4,5.0*diag(5))[,,1]
#' # A should be singular
#' eigen(A)$values
#'
rPseudoWishart <- function(n, df, Sigma) {
  tol = 1e-06
  p <- ncol(Sigma)
  if (!all(dim(Sigma) == c(p, p)))
    stop("Sigma must be square.")
  if (df > p - 1) {
    warning("df > dimension of Sigma - 1, using rWishart.")
    return(stats::rWishart(n, df, Sigma))
  }
  if (!(df == round(df))) stop("df needs to be a whole number.")
  if (df < 1) stop("df needs to be greater than 1.")
  eigens <- eigen(Sigma, symmetric = TRUE)
  ev <- eigens$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive semidefinite.")
  # Note: this can deal with 0 eigenvalues, so it could be singular.
  # It cannot, however, deal with <0 eigenvalues.
  sqrtmatrix <- eigens$vectors %*% diag(sqrt(pmax(ev, 0)), p)
  X <- array(stats::rnorm(p*df*n), dim = c(df, p, n))
  Xresult <- array(0, dim = c(p, p, n))
  for (i in 1:n) Xresult[,,i] = tcrossprod(tcrossprod(sqrtmatrix, X[,,i]))
  return(Xresult)
}


.onUnload <- function(libpath) {
  library.dynam.unload("CholWishart", libpath)
}
