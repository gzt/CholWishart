<!-- badges: start --> [![Project Status: Active – The project has
reached a stable, usable state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build
Status](https://travis-ci.org/gzt/CholWishart.svg?branch=master)](https://travis-ci.org/gzt/CholWishart)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/gzt/CholWishart?branch=master&svg=true)](https://ci.appveyor.com/project/gzt/CholWishart)
[![codecov](https://codecov.io/gh/gzt/CholWishart/branch/master/graph/badge.svg)](https://codecov.io/gh/gzt/CholWishart)
[![R build
status](https://github.com/gzt/CholWishart/workflows/R-CMD-check/badge.svg)](https://github.com/gzt/CholWishart/actions)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![](https://cranlogs.r-pkg.org/badges/CholWishart)](https://cran.r-project.org/package=CholWishart)

------------------------------------------------------------------------

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.3.2-6666ff.svg)](https://cran.r-project.org/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/CholWishart)](https://cran.r-project.org/package=CholWishart)
[![packageversion](https://img.shields.io/badge/Package%20version-1.1.0-orange.svg?style=flat-square)](https://github.com/gzt/CholWishart/releases)
<!-- badges: end -->

------------------------------------------------------------------------

[![Last-changedate](https://img.shields.io/badge/last%20change-2020--04--01-yellowgreen.svg)](https://github.com/gzt/CholWishart/)

<!-- README.md is generated from README.Rmd. Please edit that file -->

CholWishart
===========

The most common use for this package is likely fast sampling from the
inverse Wishart or use of the multivariate gamma or digamma functions.
This is a package for fast computation of various functions related to
the Wishart distribution, such as sampling from the Cholesky
factorization of the Wishart, sampling from the inverse Wishart,
sampling from the Cholesky factorization of the inverse Wishart,
sampling from the pseudo Wishart, computing densities for the Wishart
and inverse Wishart, and computing a few auxiliary functions such as the
multivariate gamma and digamma functions. Many of these functions are
written in C to maximize efficiency.

The output of the sampling functions is in the same format as the output
of `stats::rWishart()`.

The main idea: sampling for multivariate or matrix variate statistics
often makes use of distributions related to the Wishart. There are
implementations in a few packages but they are often in R and much
slower than the basic `stats::rWishart()` or there is a lot of
associated overhead in the package. Here, then, is a small package with
some of those distributions and related functions. As the name suggests,
the initial purpose was sampling from the Cholesky factorization of a
Wishart distribution.

Usage
=====

The package offers the following distributions:

-   The inverse Wishart:

<!-- -->

    rInvWishart(n = 1,df = 5,Sigma = diag(4)))

-   The Cholesky decomposition of the inverse Wishart:

<!-- -->

    rInvCholWishart(n = 1,df = 5,Sigma = diag(4)))

-   The Cholesky decomposition of the Wishart

<!-- -->

    rCholWishart(n = 1,df = 5,Sigma = diag(4)))

-   The pseudo Wishart and its pseudo-inverse, the generalized inverse
    Wishart (`df` less than the dimension of `Sigma` minus 1):

<!-- -->

    rPseudoWishart(n = 1,df = 3,Sigma = diag(5)))
    rGenInvWishart(n = 1,df = 3,Sigma = diag(5)))

-   The package also offers density functions for the Wishart and
    inverse Wishart distributions:

<!-- -->

    dInvWishart(x = diag(4), df = 5, Sigma = diag(4), log = TRUE)
    dWishart(x = diag(4), df = 5, Sigma = diag(4), log = TRUE)

In addition to this, the package offers the multivariate gamma function
and multivariate digamma function, the multivariate counterparts of the
gamma and digamma (derivative of the log of the gamma) functions.

    mvgamma(x = 5, p = 3)
    mvdigamma(x = 5, p = 3)

I suspect the most useful functions for most people will be the
`mvgamma` and `rInvWishart` functions.

Installation
============

Now available on CRAN, install it at:

    install.packages('CholWishart')

Install the latest development version at:

    # install.packages('devtools')
    devtools::install_github("gzt/CholWishart")

Calling from other packages
===========================

The header `inst/install/CholWishart.h` is available so that the
exported functions can be called in C from other packages.

Contributing
============

Please note that the ‘CholWishart’ project is released with a
[Contributor Code of
Conduct](https://gzt.github.io/CholWishart/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
