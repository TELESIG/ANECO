test_that("soundclassifier_predict() returns probabilities and predicted_class", {
  skip_if_not_installed("caret")

  set.seed(123)

  data("indices_data", package = "ANECO")
  data("labeled_data", package = "ANECO")
  data("newdata", package = "ANECO")

  fit <- soundclassifier_fit(
    ref  = labeled_data,
    ind  = indices_data,
    p    = 0.8,
    nmax = "even"
  )

  pred_result <- soundclassifier_predict(data = newdata, mod = fit$model)

  # Output structure
  expect_s3_class(pred_result, "data.frame")
  expect_true("predicted_class" %in% names(pred_result))

  # Prob columns exist
  prob_cols <- grep("^prob_", names(pred_result), value = TRUE)
  expect_true(length(prob_cols) >= 2)

  # Dimensions match input
  expect_equal(nrow(pred_result), nrow(newdata))

  # predicted_class type: character or factor (both acceptable)
  expect_true(is.character(pred_result$predicted_class) || is.factor(pred_result$predicted_class))

  # Predicted labels are valid classes (use prob_ columns as the ground truth)
  classes <- sub("^prob_", "", prob_cols)
  expect_true(all(as.character(pred_result$predicted_class) %in% classes))

  # Probabilities sanity
  expect_true(all(pred_result[, prob_cols, drop = FALSE] >= 0))
  expect_true(all(pred_result[, prob_cols, drop = FALSE] <= 1))
  rs <- rowSums(pred_result[, prob_cols, drop = FALSE])
  expect_true(all(abs(rs - 1) < 1e-6))
})
