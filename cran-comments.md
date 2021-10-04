This is a resubmission. I believe I have fixed all the links apart from one mentioned below
(which is an issue with the reference provided by DOI itself).

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

Some checks provide a NOTE about a JSTOR link but it is the 
canonical one generated from the DOI and the link works.


## Reverse dependencies

I have confirmed that the changes do not impact the reverse dependencies.

## Other notes

This is a minor update which fixes the USE_FC_LEN issue.
