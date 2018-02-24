
context("Testing input integrity")


test_that("trying wrong type of input", {


  expect_error(rCholWishart(1, df = 5, Sigma = "A"))
  expect_error(rCholWishart(1, df = "A", Sigma = diag(5)))
  expect_error(rCholWishart("A", df = 5, Sigma = diag(5)))
  expect_error(rCholWishart(.4, df = 5, Sigma = diag(5)))

  expect_error(rInvCholWishart(1, df = 5, Sigma = "A"))
  expect_error(rInvCholWishart(1, df = "A", Sigma = diag(5)))
  expect_error(rInvCholWishart("A", df = 4, Sigma = diag(5)))
  expect_error(rInvCholWishart(.4, df = 4, Sigma = diag(5)))

  expect_error(rInvWishart(1, df = 5, Sigma = "A"))
  expect_error(rInvWishart(1, df = "A", Sigma = diag(5)))
  expect_error(rInvWishart("A", df = 4, Sigma = diag(5)))
  expect_error(rInvWishart(.4, df = 4, Sigma = diag(5)))

  expect_error(dInvWishart(1, df = 5, Sigma = "A"))
  expect_error(dInvWishart(1, df = "A", Sigma = diag(5)))
  expect_error(dInvWishart("A", df = 4, Sigma = diag(5)))
  expect_error(dInvWishart(.4, df = 4, Sigma = diag(5)))

  expect_error(dWishart(1, df = 5, Sigma = "A"))
  expect_error(dWishart(1, df = "A", Sigma = diag(5)))
  expect_error(dWishart("A", df = 4, Sigma = diag(5)))
  expect_error(dWishart(.4, df = 4, Sigma = diag(5)))


  expect_error(rCholWishart(1, 10, matrix(c(3, 1, 1, 1, 1, 3), nrow = 2)))
  expect_error(rInvCholWishart(1, 10, matrix(c(3, 1, 1, 1, 1, 3), nrow =
                                               2)))

}
)

test_that("Out of bounds numeric input: ", {


  expect_error(rCholWishart(1, 10, matrix(c(3, 1, 1, 1, 1, 3), nrow = 2)))
  expect_error(rInvCholWishart(1, 10, matrix(c(3, 1, 1, 1, 1, 3), nrow =
                                               2)))

  expect_error(lmvgamma(-1, 5))
  expect_error(lmvgamma(1, -5))

  expect_error(rCholWishart(1, 10, -diag(5)))
  expect_error(rInvCholWishart(1, 10, -diag(5)))


  expect_error(rCholWishart(1, 4, diag(5)))
  expect_error(rInvCholWishart(1, 4, diag(5)))


  expect_error(rCholWishart(-1, 10, diag(5)))
  expect_error(rInvCholWishart(-1, 10, diag(5)))



  expect_error(rCholWishart(1, 10, matrix(c(1, 1, 1, 1), nrow = 2)))
  expect_error(rInvCholWishart(1, 10, matrix(c(1, 1, 1, 1), nrow = 2)))

  }
)