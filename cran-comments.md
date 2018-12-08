## Test environments
* local Fedora 28 install, R 3.5.1
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.1
* OS X Sierra 10.13.3 (on travis-ci), R 3.5.0
* win-builder (devel and release)
* R-hub.io Windows Server 2008 R2 SP1, R-devel, 32/64 bit
           Fedora Linux, R-devel
           Ubuntu Linux 16.04 LTS, R-release
           Debian Linux, R-devel

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

There are no reverse dependencies.

## Other notes

This is a minor set of changes from the previous version that does not 
make any incompatible changes in usage from previous versions. As before,
information about authorship is included in inst/AUTHORS.
