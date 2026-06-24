test_that("npbivp fits and predicts with one platform", {
  set.seed(1)

  n <- 10
  n_test <- 3

  X <- matrix(rnorm(n * 10), n, 10)
  X_test <- matrix(rnorm(n_test * 10), n_test, 10)
  y <- rnorm(n)

  fit <- npbivp(list(X), y, iter = 30, burn_in = 15)
  pred <- predict(fit, list(X_test), thin = 2)

  expect_s3_class(fit, "npbivp")
  expect_s3_class(pred, "npbivp_pred")
  expect_equal(nrow(pred$pred_draws), n_test)
})
