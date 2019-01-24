## Test environments
* local Fedora 29 install, R 3.5.1
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.1, 3.5.2, devel
* OS X Sierra 10.13.3 (on travis-ci), R 3.5.2 
* win-builder (devel and release)


## R CMD check results

1 NOTE: Found the following (possibly) invalid URLs:
URL: http://www.jstor.org/stable/2335827 (moved to https://www.jstor.org/stable/2335827)
    From: inst/doc/wishart.html
    Status: 403
    Message: Forbidden
    
This paper is referred to in the documentation by \doi{10.2307/2335827}, so the problem is not 
with this package, but rather resolution of the DOI reference. It is probably better to leave
it as the DOI reference than to hard-code some other URL.

## Reverse dependencies

This does not affect the 1 reverse dependency.

## Other notes

This is a modification of two functions in R that does not make any changes in usage of functions
from previous versions. As before, information about authorship is included in inst/AUTHORS.
