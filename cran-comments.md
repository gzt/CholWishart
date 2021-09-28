## Test environments
* local Fedora install, R 4.0.5
* ubuntu Ubuntu 20.04.3), R devel, release, old (on github)
* Microsoft Windows Server 2019, R release (on github)
* Mac OS X 10.15.7, R release (on github)
* win-builder.r-project.org, R release

## R CMD check results

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_FORTIFY_SOURCE=2’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’

This is a local configuration issue.


## Reverse dependencies

Functions used in reverse dependencies did not change.

## Other notes

This is a minor update which fixes the USE_FC_LEN issue.
