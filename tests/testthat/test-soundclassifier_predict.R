test_that("soundclassifier_predict() returns classification output", {
  skip_on_cran()  # Optional, to avoid long runs on CRAN

  # Load example data
  data("indices_data", package = "ANECO")
  data("labeled_data", package = "ANECO")
  data("newdata", package = "ANECO")

  # Train a model
  fit <- soundclassifier_fit(
    ref = labeled_data,
    ind = indices_data,
    p = 0.8,
    nmax = "even"
  )

  # Predict using the model
  pred_result <- soundclassifier_predict(data = newdata, mod = fit$model)

  # Check output structure
  expect_s3_class(pred_result, "data.frame")
  expect_true("predicted" %in% names(pred_result))

  # Check dimensions match input data
  expect_equal(nrow(pred_result), nrow(newdata))

  # Check predicted values are factor and match expected levels
  expect_type(pred_result$predicted, "integer")  # factors are integer-based
  expect_true(all(as.character(pred_result$predicted) %in% levels(fit$testdata$Category)))
})
