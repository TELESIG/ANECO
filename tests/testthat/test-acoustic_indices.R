test_that("acoustic_indices() returns a valid minute-level table on example wavs", {
  skip_if_not_installed("tuneR")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("progressr")

  set.seed(123)

  data("example_wavs", package = "ANECO")

  # Isolated temp dir
  temp_dir <- withr::local_tempdir()

  file1 <- file.path(temp_dir, "MORN01_20200101_000000.wav")
  file2 <- file.path(temp_dir, "NOON01_20200101_010000.wav")

  tuneR::writeWave(example_wavs$morning, filename = file1)
  tuneR::writeWave(example_wavs$noon, filename = file2)

  # Minimal but representative index set:
  # - one calibrated SPL summary
  # - one entropy
  # - one spectral stat (raw path)
  # - noise feature + band features (efb expands columns)
  sel <- c("msldB_bio", "Hf", "SpecCent", "S2N_chi", "efb")

  res <- acoustic_indices(
    dir = temp_dir,
    calibparam = "SM4",
    sel.ind = sel,
    exclude = "none",
    noise.ind = FALSE,     # keep test small; we explicitly request S2N_chi + efb
    save.file = FALSE,
    ncores = 1,
    parallel = "files"
  )

  # ---- Basic structure ----
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)

  # ---- Required ID fields ----
  id_cols <- c("Name", "Code", "Site", "RecorderID", "Date", "Year", "Month", "Day", "Hour", "Minute", "Duration")
  expect_true(all(id_cols %in% names(res)))

  # ---- Requested indices exist ----
  expect_true(all(c("msldB_bio", "Hf", "SpecCent", "S2N_chi") %in% names(res)))

  # efb should expand to E## columns (at least E01 if sampling rate allows >= 1 kHz)
  efb_cols <- grep("^E\\d+\\d+$", names(res), value = TRUE)
  expect_true(length(efb_cols) >= 1)

  # ---- Sanity checks ----
  expect_true(all(res$Duration > 0))
  expect_true(all(res$Minute >= 0))
})
