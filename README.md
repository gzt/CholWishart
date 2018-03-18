
<!-- README.md is generated from README.Rmd. Please edit that file -->
CholWishart
===========

A package for fast computation of various functions related to the Wishart distribution, such as sampling from the Cholesky factorization of the Wishart, sampling from the inverse Wishart, sampling from the Cholesky factorization of the inverse Wishart, computing densities for the Wishart and inverse Wishart, and computing a few auxiliary functions such as the multivariate gamma and digamma functions. Many of these functions are written in C to maximize efficiency.

The output of the sampling functions is in the same format as the output of `stats::rWishart()`.

The main idea: sampling for multivariate or matrix variate statistics often makes use of distributions related to the Wishart. There are implementations in a few packages but they are often in R and much slower than the basic `stats::rWishart()` or there is a lot of associated overhead in the package. Here, then, is a small package with some of those distributions and related functions. As the name suggests, the initial purpose was sampling from the Cholesky factorization of a Wishart distribution.

Now available on CRAN, install it at:

    install.packages('CholWishart')

Install the latest development version at:

    devtools::install_github("gzt/CholWishart")
