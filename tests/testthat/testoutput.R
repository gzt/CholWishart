
context("Checking outputs match")


test_that("Testing helper functions:", {
  expect_equal(lgamma(1:12) - lmvgamma(1:12, 1), array(0, dim = 12),
    tolerance = 1e-7
  )
  expect_equal(digamma(1:12) - mvdigamma(1:12, 1), array(0, dim = 12),
    tolerance = 1e-7
  )

  expect_equal(gamma(1:12) - mvgamma(1:12, 1), array(0, dim = 12),
    tolerance = 1e-7
  )
  p <- 2
  expect_equal(
    (p * (p - 1) / 4 * log(pi) + lgamma(5 - 0) + lgamma(5 - .5)),
    as.numeric(lmvgamma(5, 2))
  )
  expect_equal(
    (3 * (3 - 1) / 4 * log(pi) + lgamma(5 - 0) +
      lgamma(5 - .5) + lgamma(5 - 1)),
    as.numeric(lmvgamma(5, 3))
  )
  expect_equal(digamma(1:100), c(mvdigamma(1:100, 1)))
})

test_that("Equivalent outputs for different options:", {
  set.seed(20180211)
  a_mat <- rInvCholWishart(1, 10, .5 * diag(5))[, , 1]
  set.seed(20180211)
  b_mat <- rCholWishart(1, 10, 2 * diag(5))[, , 1]
  set.seed(20180211)
  c_mat <- chol(rWishart(1, 10, 2 * diag(5))[, , 1])


  expect_equal(sum(abs(a_mat[lower.tri(a_mat)])), 0)
  expect_equal(sum(abs(b_mat[lower.tri(b_mat)])), 0)
  expect_equal(crossprod(a_mat) %*% crossprod(b_mat), diag(5))
  expect_equal(b_mat, c_mat)

  set.seed(20180221)
  a_mat <- rInvWishart(1, 10, 5 * diag(5))[, , 1]
  set.seed(20180221)
  b_mat <- rWishart(1, 10, .2 * diag(5))[, , 1]
  expect_equal(a_mat %*% b_mat, diag(5))

  # this really shouldn't work in general, it only works on diag()
  expect_equal(
    dWishart(diag(5), 10, 5 * diag(5)),
    dInvWishart(diag(5), 10, .2 * diag(5))
  )
  expect_equal(
    dWishart(diag(5), 10, 5 * diag(5), log = FALSE),
    dInvWishart(diag(5), 10, .2 * diag(5), log = FALSE)
  )

  a_mat <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
  b_mat <- dWishart(a_mat, df = 4, Sigma = diag(3))
  c_mat <- dInvWishart(a_mat, df = 4, Sigma = diag(3))
  expect_equal(b_mat[1], b_mat[2])
  expect_equal(c_mat[1], c_mat[2])
  expect_equal(b_mat[1], -7.255196, tolerance = 1e-6)

  set.seed(20180221)
  a_mat <- rGenInvWishart(1, 4, 5 * diag(5))[, , 1]
  set.seed(20180221)
  b_mat <- rPseudoWishart(1, 4, 5 * diag(5))[, , 1]
  expect_equal(a_mat %*% b_mat %*% a_mat, a_mat)
  expect_equal(b_mat %*% a_mat %*% b_mat, b_mat)
})

test_that("Ranks are correct:", {
  rango <- 5
  expect_equal(qr(rGenInvWishart(1, rango, diag(rango + 2))[, , 1])$rank, rango)
  expect_equal(qr(rPseudoWishart(1, rango, diag(rango + 2))[, , 1])$rank, rango)
  expect_equal(
    qr(rInvWishart(1, rango + 2, .5 * diag(rango))[, , 1])$rank,
    rango
  )
})
