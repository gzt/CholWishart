# CholWishart

[![Travis-CI Build Status](https://travis-ci.org/gzt/CholWishart.svg?branch=master)](https://travis-ci.org/gzt/CholWishart) 

A package for fast computation of various functions related to the Wishart distribution, such as sampling from the Cholesky factor of the Wishart, sampling from the inverse Wishart, sampling from the Cholesky factor of the inverse Wishart, computing densities for the Wishart and inverse Wishart, and computing a few auxiliary functions such as the multivariate gamma and digamma functions. Many of these functions are written in C to maximize efficiency. 

The output of the sampling functions is in the same format as the output 
of `stats::rWishart`.

Install it at:

```
devtools::install_github("gzt/CholWishart")
```

