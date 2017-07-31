library(testthat)
library(zipfextR)

# test_check("zipfextR")
test_that("[moezipf - pmf] The summation of all probabilities must be 1.", {
  expect_equal(1,sum(dmoezipf(1:1000, alpha = 2.5, beta = 0.10)), tolerance = 1*10^(-3))
  expect_equal(1,sum(dmoezipf(1:1000, alpha = 2.5, beta = 1.3)), tolerance = 1*10^(-3))
  expect_equal(1,sum(dmoezipf(1:1000, alpha = 2.5, beta = 3)), tolerance = 1*10^(-3))
})

test_that("[moezipf - pmf] The probability lays on the interval (0, 1).", {
  pvalue <- dmoezipf(x = 9, alpha = 2, beta = 0.10)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)

  pvalue <- dmoezipf(x = 9, alpha = 1.25, beta = 5.3)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)
})

test_that("[moezipf - cdf] The summation of all probabilities must be 1.", {
  expect_equal(1,pmoezipf(1000, alpha = 2.5, beta = 0.10), tolerance = 1*10^(-3))
  expect_equal(1,pmoezipf(1000, alpha = 2.5, beta = 1.3), tolerance = 1*10^(-3))
  expect_equal(1,pmoezipf(1000, alpha = 2.5, beta = 3), tolerance = 1*10^(-3))
})

test_that("[moezipf- cdf] The cumulative probabilities have to be in the interval (0, 1).", {
  pvalue <- pmoezipf(q = 9, alpha = 2, beta = 0.10)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)

  pvalue <- pmoezipf(q = 9, alpha = 1.25, beta = 5.3)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)
})
