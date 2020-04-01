#   densities.R
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



#' Density for Random Wishart Distributed Matrices
#'
#' Compute the density of an observation of a random Wishart distributed matrix
#' (\code{dWishart}) or an observation
#' from the inverse Wishart distribution (\code{dInvWishart}).
#'
#'    Note there are different ways of parameterizing the Inverse
#'    Wishart distribution, so check which one you need.
#'     Here,  If \eqn{X \sim IW_p(\Sigma, \nu)}{X ~ IW_p(Sigma, df)} then
#'     \eqn{X^{-1} \sim W_p(\Sigma^{-1}, \nu)}{X^{-1} ~ W_p(Sigma^{-1}, df)}.
#'     Dawid (1981) has a different definition: if
#'     \eqn{X \sim W_p(\Sigma^{-1}, \nu)}{X ~ W_p(Sigma^{-1}, df)} and
#'     \eqn{\nu > p - 1}{df > p - 1}, then
#'     \eqn{X^{-1} = Y \sim IW(\Sigma, \delta)}{X^{-1} = Y ~ IW(Sigma, delta)},
#'     where
#'     \eqn{\delta = \nu - p + 1}{delta = df - p + 1}.
#'
#' @param x positive definite \eqn{p \times p}{p * p} observations for density
#'          estimation - either one matrix or a 3-D array.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{p \times p}{p * p} "scale" matrix,
#'              the matrix parameter of the distribution.
#' @param log logical, whether to return value on the log scale.
#'
#' @return Density or log of density
#'
#' @references
#' Dawid, A. (1981). Some Matrix-Variate Distribution Theory:
#' Notational Considerations and a Bayesian Application.
#' \emph{Biometrika}, 68(1), 265-274. \doi{10.2307/2335827}
#'
#' Gupta, A. K.  and D. K. Nagar (1999). \emph{Matrix variate distributions}.
#' Chapman and Hall.
#'
#' Mardia, K. V., J. T. Kent, and J. M. Bibby (1979)
#' \emph{Multivariate Analysis},
#' London: Academic Press.
#' @export
#'
#' @examples
#' set.seed(20180222)
#' A <- rWishart(1, 10, diag(4))[, , 1]
#' A
#' dWishart(x = A, df = 10, Sigma = diag(4L), log = TRUE)
#' dInvWishart(x = solve(A), df = 10, Sigma = diag(4L), log = TRUE)
dWishart <- function(x, df, Sigma, log = TRUE) {
  if (!is.numeric(Sigma)) {
    stop("'Sigma' must be numeric")
  }
  Sigma <- as.matrix(Sigma)
  dims <- dim(Sigma)
  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }
  if (dim(x)[1L] != dim(x)[2L] ||
    dim(x)[1L] != dims[1L] || dims[1L] != dims[2L]) {
    stop("non-conformable dimensions")
  }
  if (!isSymmetric(Sigma)) {
    stop("non-symmetric input")
  }
  if (length(dim(x)) < 3L) x <- array(x, dim = c(dim(x), 1L))
  dimx <- dim(x)
  ldensity <- rep(0, dimx[3L])
  cholS <- chol(Sigma)
  ldetS <- sum(log(diag(cholS))) * 2
  for (i in seq(dimx[3])) {
    if (!isSymmetric(x[, , i])) {
      stop("non-symmetric input")
    }
    cholX <- chol(x[, , i])
    ldetX <- sum(log(diag(cholX))) * 2
    ldensity[i] <-
      .5 * (df - dims[1L] - 1) * ldetX +
      -.5 * sum(diag(chol2inv(cholS) %*% x[, , i])) +
      -(df * dims[1L] / 2 * log(2)) - .5 * df * ldetS +
      -lmvgamma(df / 2, dims[1L])
  }
  if (log) {
    return(ldensity)
  } else {
    return(exp(ldensity))
  }
}

#' @describeIn dWishart density for the inverse Wishart distribution.
#' @export
dInvWishart <- function(x, df, Sigma, log = TRUE) {
  if (!is.numeric(Sigma)) {
    stop("'Sigma' must be numeric")
  }
  Sigma <- as.matrix(Sigma)
  dims <- dim(Sigma)
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (dim(x)[1L] != dim(x)[2L] ||
    dim(x)[1L] != dims[1L] || dims[1L] != dims[2L]) {
    stop("non-conformable dimensions")
  }
  if (!isSymmetric(Sigma)) {
    stop("non-symmetric input")
  }
  if (length(dim(x)) < 3L) x <- array(x, dim = c(dim(x), 1L))
  dimx <- dim(x)
  ldensity <- rep(0, dimx[3L])
  cholS <- chol(Sigma)
  ldetS <- sum(log(diag(cholS))) * 2
  for (i in seq(dimx[3L])) {
    if (!isSymmetric(x[, , i])) {
      stop("non-symmetric input")
    }
    cholX <- chol(x[, , i])
    ldetX <- sum(log(diag(cholX))) * 2
    ldensity[i] <-
      -.5 * (df + dims[1L] + 1) * ldetX + .5 * df * ldetS +
      -.5 * sum(diag(chol2inv(cholX) %*% Sigma)) -
      (df * dims[1L] / 2 * log(2)) - lmvgamma(df / 2, dims[1L])
  }
  if (log) {
    return(ldensity)
  } else {
    return(exp(ldensity))
  }
}
