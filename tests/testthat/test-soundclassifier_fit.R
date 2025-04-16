test_that("soundclassifier_fit() returns expected structure and model", {
  skip_on_cran()  # Optional: skips on CRAN to avoid long processing or OS-specific issues

  # Load example data provided with the package
  data("indices_data", package = "ANECO")
  data("labeled_data", package = "ANECO")

  # Run the function with a small subset to make the test faster
  fit_result <- soundclassifier_fit(
    ref = labeled_data,
    ind = indices_data,
    p = 0.8,
    nmax = "even"  # Force equal class representation to ensure balanced categories
  )

  # Check that the result is a list with the correct elements
  expect_type(fit_result, "list")
  expect_named(fit_result, c("model", "info", "testdata"))

  # Check that model is a caret train object
  expect_s3_class(fit_result$model, "train")

  # Check that testdata is a data.frame and includes 'Category'
  expect_s3_class(fit_result$testdata, "data.frame")
  expect_true("Category" %in% names(fit_result$testdata))

  # Check that info is a named integer vector with category counts
  expect_type(fit_result$info, "integer")
  expect_true(length(fit_result$info) > 0)
})
