## Resubmission

* Previous submission had bad URIs in the README - this has been corrected.

## Test environments
* local Fedora 30 install, R 3.6.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.0, release, old, 3.3
* OS X Sierra macOS High Sierra 10.13.3 (on travis-ci), R 3.6.0 



## R CMD check results

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_FORTIFY_SOURCE=2’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’

This is a local configuration issue.


## Reverse dependencies

Functions used in reverse dependencies did not change.

## Other notes

This is an update to fix the previous submission using full URLs in the README,
which is why not all the test environments were used.
