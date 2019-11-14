## Test environments
* local Fedora 30 install, R 3.6.1
* ubuntu Ubuntu 16.04.6 LTS (on travis-ci), R devel, release, old, 3.3
* OS X Sierra macOS High Sierra 10.13.3 (on travis-ci), R 3.6.1 
* Ubuntu Linux 16.04 LTS, R-release, GCC on rhub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit on rhub
* Fedora Linux, R-devel, clang, gfortran on rhub
* Debian Linux, R-devel, GCC ASAN/UBSAN on rhub


## R CMD check results

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_FORTIFY_SOURCE=2’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’

This is a local configuration issue.


## Reverse dependencies

Functions used in reverse dependencies did not change.

## Other notes

This is a minor update which adds header files so the package
can be linked to from other packages more easily.
