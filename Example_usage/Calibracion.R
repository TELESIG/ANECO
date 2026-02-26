# ============================================================
# Supplementary Material: Calibration mismatch & residual variability
# Cluster bootstrap by minute (SM4 vs MINI)
#
# NOTE ON DATA AVAILABILITY:
# The raw calibration-test audio datasets used in this analysis are very large and
# are therefore not distributed with the repository. They can be shared upon
# reasonable request by contacting the corresponding author at:
# jimmybarrantesm@gmail.com
# ============================================================

# --------------------------
# 0) Packages
# --------------------------
# install.packages(c("rio", "dplyr", "tidyr", "ggplot2", "purrr"))
library(rio)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

set.seed(456)

# --------------------------
# 1) Load and merge input files
# --------------------------
# Folder containing exported tables (CSV) for two conditions:
# - "no_" prefix: uncalibrated/mismatch setup
# - "ok_" prefix: calibrated/match setup

parent_dir <- "C:/Paisaje acustico UNA/Articulos/Aneco/Ecologial frontiers/Enviado - Revision 1/Suplemental materials/CalibrationTest"

read_and_bind <- function(folder, pattern) {
  files <- list.files(path = folder, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found for pattern: ", pattern)
  
  bind_rows(lapply(files, rio::import))
}

dfraw <- read_and_bind(parent_dir, pattern = "no_")  # uncalibrated (mismatch)
dfcal <- read_and_bind(parent_dir, pattern = "ok_")  # calibrated (match)

# --------------------------
# 2) Add metadata: calibration setting and recorder model
# --------------------------
# IMPORTANT: The original logic is preserved:
# - For dfraw: Calib is SM4 except RecorderID containing "A" which is MINI.
# - Model is MINI except RecorderID containing "A" which is SM4.
# - For dfcal: Calib is MINI except RecorderID containing "A" which is SM4.
# - Model is MINI except RecorderID containing "A" which is SM4.

annotate_dfraw <- function(df) {
  df$Calib <- "SM4"
  df$Calib[grep("A", df$RecorderID)] <- "MINI"
  
  df$Model <- "MINI"
  df$Model[grep("A", df$RecorderID)] <- "SM4"
  df
}

annotate_dfcal <- function(df) {
  df$Calib <- "MINI"
  df$Calib[grep("A", df$RecorderID)] <- "SM4"
  
  df$Model <- "MINI"
  df$Model[grep("A", df$RecorderID)] <- "SM4"
  df
}

dfraw <- annotate_dfraw(dfraw)
dfcal <- annotate_dfcal(dfcal)

df <- bind_rows(dfcal, dfraw) %>%
  mutate(
    Status  = ifelse(Calib == Model, "OK", "NO"),
    minuteID = paste0("m", Hour, Minute)
  )

# --------------------------
# 3) Optional exploratory plots (quick checks)
# --------------------------
# Example: Mamp summary by Code and calibration status
# (kept optional; does not affect results)
# ggplot(df, aes(x = Code, y = Mamp, colour = Status)) +
#   stat_summary() +
#   facet_wrap(~Status) +
#   theme_bw()

# Subsets used in the manuscript logic
dfs  <- df %>% filter(Calib == "SM4")      # treated as "uncal"
dfok <- df %>% filter(Status == "OK")      # calibrated minutes (match)

dfs$state  <- "Uncalibrated"
dfok$state <- "Calibrated"

# --------------------------
# 4) Cluster bootstrap by minute (difference SM4 - MINI)
# --------------------------
cluster_boot_diff <- function(data, indice, nsim = 1000) {
  minutes <- unique(data$minuteID)
  
  sampres <- numeric(nsim)      # mean signed difference
  sampres_abs <- numeric(nsim)  # mean absolute difference
  
  for (s in seq_len(nsim)) {
    
    # Cluster bootstrap at the minute level
    boot_minutes <- sample(minutes, size = length(minutes), replace = TRUE)
    data_b <- dplyr::filter(data, minuteID %in% boot_minutes)
    
    test <- data_b %>%
      group_by(Model, minuteID) %>%
      summarise(
        val = mean(.data[[indice]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_wider(names_from = Model, values_from = val) %>%
      filter(!is.na(SM4) & !is.na(MINI)) %>%
      mutate(
        dif_abs = abs(SM4 - MINI),
        dif     = SM4 - MINI
      )
    
    sampres[s]     <- mean(test$dif, na.rm = TRUE)
    sampres_abs[s] <- mean(test$dif_abs, na.rm = TRUE)
  }
  
  data.frame(
    mean    = mean(sampres, na.rm = TRUE),
    low     = unname(quantile(sampres, probs = 0.025, na.rm = TRUE)),
    upp     = unname(quantile(sampres, probs = 0.975, na.rm = TRUE)),
    absmean = mean(sampres_abs, na.rm = TRUE),
    abslow  = unname(quantile(sampres_abs, probs = 0.025, na.rm = TRUE)),
    absupp  = unname(quantile(sampres_abs, probs = 0.975, na.rm = TRUE))
  )
}

indices_to_test <- c("Mamp", "ACI", "HvSPL", "ADI", "mdBGL", "NDSI")
nsim <- 1000

consolidated <- bind_rows(lapply(indices_to_test, function(ind) {
  
  # Calibrated (Status == OK)
  rok <- cluster_boot_diff(dfok, indice = ind, nsim = nsim) %>%
    mutate(state = "Calibrated")
  
  # Uncalibrated (Calib == SM4 in the original script)
  rs <- cluster_boot_diff(dfs, indice = ind, nsim = nsim) %>%
    mutate(state = "Uncalibrated")
  
  bind_rows(rok, rs) %>%
    mutate(indice = ind)
}))

print(consolidated)

# --------------------------
# 5) Main figure (mean difference with 95% bootstrap CI)
# --------------------------
p <- ggplot(consolidated, aes(x = state, y = mean, ymin = low, ymax = upp)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.05) +
  geom_hline(yintercept = 0) +
  geom_segment(aes(x = state, y = 0, xend = state, yend = mean), linetype = 2) +
  geom_text(
    aes(y = mean, label = round(mean, 2)),
    hjust = -0.2, vjust = -1, size = 3, color = "gray30"
  ) +
  facet_wrap(~ indice, scales = "fixed") +
  ylab("Mean difference between recorder models (SM4 - MINI)") +
  xlab("") +
  theme_bw() +
  ylim(c(-6, 8))

print(p)

# Optional export
# write.csv(consolidated, file = "res_diff.csv", row.names = FALSE)
# ggsave("calibration_bootstrap_differences.png", p, width = 10, height = 6, dpi = 300)