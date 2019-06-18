## Test environments
* local Fedora 30 install, R 3.6.0
* Ubuntu Linux 16.04 LTS R-release, Debian Linux, R-devel, Fedora Linux, 
  R-devel, (r-hub)
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.0, release, old, 3.3
* OS X Sierra macOS High Sierra 10.13.3 (on travis-ci), R 3.6.0 
* win-builder (devel and release)


## R CMD check results

On my own machine I got this NOTE but it doesn't show up elsewhere,
so it's probably just my local configuration

checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_FORTIFY_SOURCE=2’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’


## Reverse dependencies

Functions used in reverse dependencies did not change.

## Other notes

This is a minor update to the documentation and the associated web site of the 
package.
