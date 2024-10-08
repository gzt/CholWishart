# CholWishart 1.1.4

* Update names
* Change PI to M_PI

# CholWishart 1.1.2

* Minor updates to documentation to fix URLs for CRAN submission.

# CholWishart 1.1.1

* Updates to conform with the following recommendation: 
  https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Fortran-character-strings While this should not
  break anything, I believe this may require R>=3.6.2.

# CholWishart 1.1.0

* Latest update to CRAN.

# CholWishart 1.0.3

* Minor updates to documentation


# CholWishart 1.0.2

* Add `R_RegisterCCallable()` interface and header `inst/install/CholWishart.h` 
  so that the exported functions can be called in C from other packages.
  All the exported functions are now also available to be called by 
  external C code in other packages, just pull in the header.

# CholWishart 1.0.1

* Update documentation to use `pkgdown`.

# CholWishart 1.0.0

* Tweak the documentation
* Port pseudo-Wishart to C, gen inv based on pseudo-Wishart

# CholWishart 0.9.4

* Add new functions to the vignette
* Add generalized inverse Wishart (pseudo inverse of the pseudo Wishart)
* Add pseudo-Wishart (Wishart distribution based on fewer observations than the 
dimension of the covariance matrix).
* Add contributor guidelines and code of conduct.

# CholWishart 0.9.3

* Minor update to internal functions

# CholWishart 0.9.2

* Adding possibility of array input to density functions. 

# CholWishart 0.9.1

* Finalize edits to documentation including additional references.

# CholWishart 0.9.0.9002

* Add more documentation, add more references to documentation, clean LaTeX equations in documentation.

# CholWishart 0.9.0.9001

* Add additional tests for `dWishart` and `dInvWishart` functions
* Add references and equations to help files
* Add additional tests for complex entries (should fail) and other erroneous input


# CholWishart 0.9.0

* Feature complete, fully documented, and the math should be correct.

# CholWishart 0.1.0

* Breaking off from `matrixdist`



