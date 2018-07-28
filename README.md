
<!-- README.md is generated from README.Rmd. Please edit that file -->
CholWishart
===========

A package for fast computation of various functions related to the Wishart distribution, such as sampling from the Cholesky factorization of the Wishart, sampling from the inverse Wishart, sampling from the Cholesky factorization of the inverse Wishart, computing densities for the Wishart and inverse Wishart, and computing a few auxiliary functions such as the multivariate gamma and digamma functions. Many of these functions are written in C to maximize efficiency.

The output of the sampling functions is in the same format as the output of `stats::rWishart()`.

The main idea: sampling for multivariate or matrix variate statistics often makes use of distributions related to the Wishart. There are implementations in a few packages but they are often in R and much slower than the basic `stats::rWishart()` or there is a lot of associated overhead in the package. Here, then, is a small package with some of those distributions and related functions. As the name suggests, the initial purpose was sampling from the Cholesky factorization of a Wishart distribution.

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

The package also offers density functions for the Wishart and inverse Wishart distributions:

    dInvWishart(x = diag(4), df = 5, Sigma = diag(4), log = TRUE)
    dWishart(x = diag(4), df = 5, Sigma = diag(4), log = TRUE)

In addition to this, the package offers the multivariate gamma function and multivariate digamma function, the multivariate counterparts of the gamma and digamma (derivative of the log of the gamma) functions.

    mvgamma(x = 5, p = 3)
    mvdigamma(x = 5, p = 3)

I suspect the most useful functions for most people will be the `mvdigamma` and `rInvWishart` functions.

Installation
============

Now available on CRAN, install it at:

    install.packages('CholWishart')

Install the latest development version at:

    # install.packages('devtools')
    devtools::install_github("gzt/CholWishart")
