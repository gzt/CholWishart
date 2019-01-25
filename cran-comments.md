## Test environments
* local Fedora 29 install, R 3.5.1
* ubuntu 14.04.5 LTS (on travis-ci), R 3.4.4, 3.5.2, devel
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

Functions used in reverse dependencies did not change.

## Other notes

This is a minor update of functions and documentation which does not alter the interface (and the version numbering is intended to indicate that the interface is stable now). As before, information about authorship is included in inst/AUTHORS.
