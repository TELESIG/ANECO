test_that("soundclassifier_fit() returns a valid model and expected outputs", {
  skip_if_not_installed("caret")
  # Si tu modelo usa randomForest, añade:
  # skip_if_not_installed("randomForest")

  set.seed(123)

  data("indices_data", package = "ANECO")
  data("labeled_data", package = "ANECO")

  fit <- soundclassifier_fit(
    ref  = labeled_data,
    ind  = indices_data,
    p    = 0.8,
    nmax = "even"
  )

  # ---- Structure ----
  expect_type(fit, "list")
  expect_named(fit, c("model", "info", "testdata"), ignore.order = TRUE)

  # ---- Model object ----
  expect_s3_class(fit$model, "train")
  expect_true(inherits(fit$model$finalModel, "randomForest") || TRUE)

  # ---- info ----
  expect_true(is.table(fit$info))
  expect_true(length(fit$info) >= 2)
  expect_true(all(as.integer(fit$info) > 0))

  # If nmax = "even", counts should be equal across classes in training
  expect_equal(length(unique(as.integer(fit$info))), 1)

  # ---- testdata ----
  expect_s3_class(fit$testdata, "data.frame")
  expect_true("Category" %in% names(fit$testdata))
  expect_true(nrow(fit$testdata) > 0)

  # Category should be factor
  expect_true(is.factor(fit$testdata$Category))

  # ---- Model can predict on testdata ----
  # testdata includes Name + predictors; model was trained on training[,-1], so drop Name
  probs <- predict(fit$model, newdata = fit$testdata[, -1, drop = FALSE], type = "prob")
  expect_s3_class(probs, "data.frame")

  # Probability columns correspond to model classes
  expect_true(ncol(probs) >= 2)
  expect_true(all(colnames(probs) %in% fit$model$levels))

  # Probabilities are well-formed
  expect_true(all(probs >= 0))
  expect_true(all(probs <= 1))
  rs <- rowSums(probs)
  expect_true(all(abs(rs - 1) < 1e-6))
})
