test_that("acoustic_indices() works on example data", {
  skip_on_cran()  # opcional: evitar en CRAN por uso de audio/tiempo

  # Cargar el objeto
  data("example_wavs", package = "ANECO")

  # Crear carpeta temporal
  temp_dir <- tempdir()

  # Guardar los dos archivos
  file1 <- file.path(temp_dir, "morning.wav")
  file2 <- file.path(temp_dir, "noon.wav")
  tuneR::writeWave(example_wavs$morning, filename = file1)
  tuneR::writeWave(example_wavs$noon, filename = file2)

  # Ejecutar la función con una selección mínima de índices
  result <- acoustic_indices(
    dir = temp_dir,
    calibparam = "SM4",
    sel.ind = "all",
    noise.ind = TRUE,
    save.file = FALSE,
    ncores = 2,
    whereparallel = "files"
  )

  # Comprobaciones
  expect_s3_class(result, "data.frame")
  expect_true("msldBA_bio" %in% names(result))
  expect_gt(nrow(result), 0)
})
