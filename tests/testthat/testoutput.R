# https://academic.oup.com/biomet/article-abstract/87/1/99/221352?redirectedFrom=PDF
# https://faculty.mccombs.utexas.edu/carlos.carvalho/CarvalhoMassamWest2007.pdf


context("Checking outputs match")


test_that("Testing helper functions:", {


  expect_equal(lgamma(1:12) - lmvgamma(1:12, 1), array(0, dim = 12),
               tolerance = 1e-7)
  expect_equal(digamma(1:12) - mvdigamma(1:12, 1), array(0, dim = 12),
               tolerance = 1e-7)

  expect_equal(gamma(1:12) - mvgamma(1:12, 1), array(0, dim = 12),
               tolerance = 1e-7)
  p = 2
  expect_equal((p * (p - 1) / 4 * log(pi) + lgamma(5 - 0) + lgamma(5 - .5)),
               as.numeric(lmvgamma(5, 2)))
  expect_equal((3 * (3 - 1) / 4 * log(pi) + lgamma(5 - 0) +
                  lgamma(5 - .5) + lgamma(5 - 1)),
               as.numeric(lmvgamma(5, 3)))
  expect_equal(digamma(1:100),c(mvdigamma(1:100,1)))
})

test_that("Equivalent outputs for different options:", {
  set.seed(20180211)
  A <- rInvCholWishart(1, 10, .5*diag(5))[, , 1]
  set.seed(20180211)
  B <- rCholWishart(1, 10, 2*diag(5))[, , 1]
  set.seed(20180211)
  C <- chol(rWishart(1, 10, 2*diag(5))[, , 1])


  expect_equal(sum(abs(A[lower.tri(A)])), 0)
  expect_equal(sum(abs(B[lower.tri(B)])), 0)
  expect_equal(crossprod(A) %*% crossprod(B), diag(5))
  expect_equal(B, C)

  set.seed(20180221)
  A <- rInvWishart(1,10,5*diag(5))[,,1]
  set.seed(20180221)
  B <- rWishart(1,10,.2*diag(5))[,,1]
  expect_equal(A %*% B, diag(5))

  # this really shouldn't work in general, it only works on diag()
  expect_equal(dWishart(diag(5), 10, 5*diag(5)), dInvWishart(diag(5), 10, .2*diag(5)))
  expect_equal(dWishart(diag(5), 10, 5*diag(5), log = FALSE), dInvWishart(diag(5), 10, .2*diag(5), log = FALSE))

  A <- array(c(diag(3), diag(3)), dim = c(3,3,2))
  B <- dWishart(A, df = 4, Sigma = diag(3))
  C <- dInvWishart(A, df = 4, Sigma = diag(3))
  expect_equal(B[1], B[2])
  expect_equal(C[1], C[2])

})
