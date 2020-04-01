#   gammas.R
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




#' Multivariate Gamma Function
#'
#' @description A special mathematical function related to the gamma function,
#'     generalized for multivariate gammas. \code{lmvgamma} is the log of the
#'     multivariate gamma, \code{mvgamma}.
#'
#'    The multivariate gamma function for a dimension p is defined as:
#'
#'    \deqn{\Gamma_{p}(a)=\pi^{p(p-1)/4}\prod_{j=1}^{p}
#'     \Gamma [a+(1-j)/2]}{Gamma_p(a)=\pi^{p(p-1)/4}* Prod_{j=1}^{p}
#'     \Gamma[a+(1-j)/2]}
#'    For \eqn{p = 1}, this is the same as the usual gamma function.
#' @param x non-negative numeric vector, matrix, or array
#' @param p positive integer, dimension of a square matrix
#'
#' @return For \code{lmvgamma} log of multivariate gamma of dimension \code{p}
#'     for each entry of \code{x}. For non-log variant,
#'     use \code{mvgamma}.
#'
#' @seealso \code{\link{gamma}} and \code{\link{lgamma}}
#' @references
#' A. K. Gupta and D. K. Nagar 1999. \emph{Matrix variate distributions}.
#' Chapman and Hall.
#'
#' Multivariate gamma function.
#' In \emph{Wikipedia, The Free Encyclopedia},from
#' \url{https://en.wikipedia.org/w/index.php?title=Multivariate_gamma_function}
#'
#' @export
#'
#' @examples
#' lgamma(1:12)
#' lmvgamma(1:12, 1L)
#' mvgamma(1:12, 1L)
#' gamma(1:12)
lmvgamma <- function(x, p) {
  if (!all(is.numeric(x), is.numeric(p))) {
    stop("non-numeric input")
  }

  # making sure that object
  # returned is same shape as object passed
  dims <- if (is.vector(x)) {
    length(x)
  } else {
    dim(as.array(x))
  }

  result <- .Call("C_lmvgamma", as.numeric(x), as.integer(p),
    PACKAGE = "CholWishart"
  )

  return(array(result, dim = dims))
}

#' @describeIn lmvgamma Multivariate gamma function.
#' @export
mvgamma <- function(x, p) {
  exp(lmvgamma(x, p))
}

#' Multivariate Digamma Function
#'
#' @description A special mathematical function related to the gamma function,
#'     generalized for multivariate distributions.
#'     The multivariate digamma function is the derivative of the log of the
#'     multivariate gamma function; for \eqn{p = 1} it is the same as the
#'     univariate digamma function.
#'
#'     \deqn{\psi_{p}(a)=\sum _{i=1}^{p}\psi(a+(1-i)/2)
#'      }{psi_p(a)=\sum psi(a+(1-i)/2)}
#'     where \eqn{\psi}{psi} is the univariate digamma function (the
#'     derivative of the log-gamma function).
#' @param x non-negative numeric vector, matrix, or array
#' @param p positive integer, dimension of a square matrix
#' @return vector of values of multivariate digamma function.
#'
#' @seealso \code{\link{gamma}}, \code{\link{lgamma}},
#'          \code{\link{digamma}}, and \code{\link{mvgamma}}
#' @export
#' @references
#' A. K. Gupta and D. K. Nagar 1999. \emph{Matrix variate distributions}.
#' Chapman and Hall.
#'
#' Multivariate gamma function.
#' In \emph{Wikipedia, The Free Encyclopedia},from
#' \url{https://en.wikipedia.org/w/index.php?title=Multivariate_gamma_function}
#'
#' @examples
#' digamma(1:10)
#' mvdigamma(1:10, 1L)
mvdigamma <- function(x, p) {
  if (!all(is.numeric(x), is.numeric(p))) {
    stop("non-numeric input")
  }
  dims <- if (is.vector(x)) {
    length(x)
  } else {
    dim(as.array(x))
  }

  result <- .Call("C_mvdigamma", as.numeric(x), as.integer(p),
    PACKAGE = "CholWishart"
  )
  return(array(result, dim = dims))
}
