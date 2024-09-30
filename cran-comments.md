This is a minor update fixing the person("john doe") issue.
I also update Calloc to R_Calloc and Free to R_Free.

## Test environments
* local Fedora install, R 4.4.1
* ubuntu (latest), R release, old-rel1 (on github)
* Microsoft Windows latest, R release (on github)
* Mac OS (latest), R release (on github)
* win-builder.r-project.org, R gdevel

## R CMD check results

 checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_GLIBCXX_ASSERTIONS’
    ‘-Wp,-U_FORTIFY_SOURCE,-D_FORTIFY_SOURCE=3’ ‘-march=x86-64’
    ‘-mno-omit-leaf-frame-pointer’


This is a local configuration issue.

## Reverse dependencies

I have confirmed that the changes do not impact the reverse dependencies.

## Other notes

This is a minor update which fixes the person field and some
memory allocation functions that do not impact functionality.
