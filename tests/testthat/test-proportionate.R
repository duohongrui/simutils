test_that("tset proportion", {
  ## Save 1 decimal place
  a <- proportionate(number = 355,
                     prop = c(0.2, 0.6, 0.15, 0.36),
                     digits = 1)
  expect_true(sum(a) != 355)
  ## The sum of the proportions is 1
  b <- proportionate(number = 355,
                     prop = c(0.2, 0.6, 0.15, 0.05),
                     prop_sum_strict = 1,
                     digits = 1)
  expect_true(sum(b) == 355)
  ## Save 0 decimal place
  c <- proportionate(number = 355,
                     prop = c(0.2, 0.6, 0.15, 0.05),
                     prop_sum_strict = 1,
                     digits = 0)
  expect_true(sum(c) == 355)
  ## The sum of the results is 355
  d <- proportionate(number = 355,
                     result_sum_strict = 355,
                     prop = c(0.2, 0.6, 0.15, 0.05),
                     prop_sum_strict = 1,
                     digits = 0)
  expect_true(sum(d) == 355)
})
