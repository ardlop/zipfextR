test_that("[zpss - preconditions] Checking preconditions", {
  expect_error(dzpss(1, "2.5", -1.5), "Incorrect alpha parameter. This parameter should be greater than one.")
  expect_error(dzpss(1, 0.5, 1.5), "Incorrect alpha parameter. This parameter should be greater than one.")
  expect_error(dzpss(1, 2.5, "-1.5"), "Incorrect lambda parameter. You should provide a numeric value.")
  expect_error(dzpss("1", 1.5, 1.5), "The x value is not included into the support of the distribution.")
  expect_error(dzpss(2.5, 2.5, 1.5), "The x value is not included into the support of the distribution.")
  expect_error(dzpss(-1, 2.5, 1.5), "The x value is not included into the support of the distribution.")
  # expect_error(qmoezipf(1.3, 2.5, 1.6), "Wrong values for the p parameter.")
  # expect_error(rmoezipf(1.6, 2.5, 1.5), "The x value is not included into the support of the distribution.")
})

test_that("[zpss - pmf] The summation of all probabilities must be 1.", {
  expect_equal(1,sum(dzpss(1:1000, alpha = 2.5, lambda = 1.3)), tolerance = 1*10^(-3))
  expect_equal(1,sum(dzpss(1:1000, alpha = 2.5, lambda = 3)), tolerance = 1*10^(-3))
})

test_that("[zpss - pmf] The probability lays on the interval (0, 1).", {
  pvalue <- dzpss(x = 9, alpha = 2, lambda = 0.10)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)

  pvalue <- dzpss(x = 9, alpha = 1.25, lambda = 5.3)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)
})

test_that("[zpss - cdf] The summation of all probabilities must be 1.", {
  expect_equal(1,pzpss(1000, alpha = 2.5, lambda = 1.3), tolerance = 1*10^(-3))
  expect_equal(1,pzpss(1000, alpha = 2.5, lambda = 3, isTruncated = T), tolerance = 1*10^(-3))
})

test_that("[zpss- cdf] The cumulative probabilities have to be in the interval (0, 1).", {
  pvalue <- pzpss(q = 9, alpha = 2, lambda = 0.10, isTruncated = T)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)

  pvalue <- pzpss(q = 9, alpha = 1.25, lambda = 5.3)
  expect_lt(pvalue, 1)
  expect_gt(pvalue, 0)
})

