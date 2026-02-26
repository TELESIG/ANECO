# ============================================================
# Supplementary Material: Practical benchmarking (R)
# Tool comparison: soundecology vs ANECO across CPU cores
# ============================================================

# --------------------------
# 0) Packages
# --------------------------
# install.packages(c("soundecology", "tuneR", "ANECO"))
library(soundecology)
library(ANECO)

set.seed(456)

# --------------------------
# 1) Data source (choose ONE option)
# --------------------------

# Option A (recommended): Download + unpack the Zenodo tar into a temporary folder
tar_url <- "https://zenodo.org/records/18778160/files/Amistosa_soundscape.tar"

tmp_dir <- tempfile("amistosa_audio_")
dir.create(tmp_dir, recursive = TRUE)
tar_file <- file.path(tmp_dir, "Amistosa_soundscape.tar")

options(timeout = 60 * 20)  # increase timeout for large downloads

download.file(tar_url, tar_file, mode = "wb", quiet = FALSE)
untar(tar_file, exdir = tmp_dir)

wavdir <- tmp_dir

# Option B: Use an existing local folder (uncomment and edit)
# wavdir <- "C:/folder/uncompressed_Amistosa_soundscape_files"

# --------------------------
# 2) Benchmark configuration
# --------------------------
cores_to_test <- c(1, 2, 4, 6, 8)

wav_files <- list.files(wavdir, pattern = "\\.wav$", full.names = TRUE)
n_files <- length(wav_files)
stopifnot(n_files > 0)

message("Number of WAV files detected: ", n_files)

# --------------------------
# 3) Helper function: time execution
# --------------------------
time_it <- function(expr) {
  t <- system.time(force(expr))
  data.frame(
    Elapsed_sec = unname(t["elapsed"]),
    User_sec    = unname(t["user.self"]),
    stringsAsFactors = FALSE
  )
}

# --------------------------
# 4) soundecology benchmark
# --------------------------
# Note: multiple_sounds() writes results to a CSV. We use a temporary file.
soundeco_indices <- c(
  "acoustic_complexity",
  "acoustic_diversity",
  "acoustic_evenness",
  "bioacoustic_index",
  "ndsi"
)

bench_soundeco <- lapply(cores_to_test, function(k) {
  
  out_csv <- tempfile(fileext = ".csv")
  
  tdf <- time_it({
    for (idx in soundeco_indices) {
      soundecology::multiple_sounds(
        directory  = wavdir,
        resultfile = out_csv,
        soundindex = idx,
        no_cores   = k
      )
    }
  })
  
  data.frame(
    Tool = "soundecology",
    Cores = k,
    tdf,
    Sec_per_file = tdf$Elapsed_sec / n_files,
    stringsAsFactors = FALSE
  )
})

bench_soundeco <- do.call(rbind, bench_soundeco)

# --------------------------
# 5) ANECO benchmark
# --------------------------
# Keep the same index set used in the manuscript benchmarking section.
aneco_indices <- c("ACI", "ADI", "AEI", "BIO_soundecology", "NDSI_soundecology")

bench_aneco <- lapply(cores_to_test, function(k) {
  
  tdf <- time_it({
    acoustic_indices(
      dir        = wavdir,
      calibparam = "SM4",
      sel.ind    = aneco_indices,
      noise.ind  = FALSE,
      bioband    = c(0, 10000),
      channel    = "left",
      wl         = 512,
      ncores     = k,
      parallel   = "files",
      save.file  = FALSE
    )
  })
  
  data.frame(
    Tool = "ANECO",
    Cores = k,
    tdf,
    Sec_per_file = tdf$Elapsed_sec / n_files,
    stringsAsFactors = FALSE
  )
})

bench_aneco <- do.call(rbind, bench_aneco)

# --------------------------
# 6) Combine results + export
# --------------------------
combined <- rbind(bench_aneco, bench_soundeco)

print(combined)

# Export (edit output path as needed)
out_csv_path <- file.path(getwd(), "benchmark_soundecology_vs_aneco.csv")
write.csv(combined, out_csv_path, row.names = FALSE)

message("Benchmark results saved to: ", out_csv_path)