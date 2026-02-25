# ============================================================
# ANECO - Reproducible script
# ============================================================

# --------------------------
# 0) Setup
# --------------------------
# install.packages(c("remotes", "ANECO", "ggplot2", "dplyr", "ggeffects",
#                    "progressr", "caret", "rio", "ggsci"))

# Optional: install ANECO from GitHub (only if needed)
remotes::install_github("TELESIG/ANECO")

library(ANECO)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(progressr)
library(caret)

set.seed(456)

progressr::handlers("txtprogressbar")
progressr::handlers(global = TRUE)

# --------------------------
# 1) Download & unpack data
# --------------------------
tar_url <- "https://zenodo.org/records/18624044/files/Amistosa_soundscape.tar"

tmp_dir <- tempfile("amistosa_audio_")
dir.create(tmp_dir, recursive = TRUE)
tar_file <- file.path(tmp_dir, "Amistosa_soundscape.tar")

options(timeout = 60 * 20)  # 20 min

download.file(tar_url, tar_file, mode = "wb", quiet = FALSE)
untar(tar_file, exdir = tmp_dir)

# quick sanity check
wav_files <- list.files(tmp_dir, pattern = "\\.wav$", full.names = FALSE)
print(head(wav_files, 5))

# --------------------------
# 2) Run acoustic_indices()
# --------------------------
sel_indices <- c("ACI", "ADI", "AAc", "msldB_bio")

ind_df <- acoustic_indices(
  dir = tmp_dir,
  calibparam = "SMmini",
  sel.ind = sel_indices,       # keep same indices used in the manuscript example
  prefix.format = "SSRRR",
  bioband = c(1000, 10000),
  noiseband = c(0, 1000),
  noise.ind = TRUE,            # keep TRUE (needed later for classifier)
  channel = "left",
  parallel = "files",
  ncores = NULL,               # NULL = all cores - 1 
  wl = 512,
  save.file = FALSE            # avoids writing CSVs in temp folder
)

print(head(ind_df))

# --------------------------
# 3) Forest-type comparison (linear models + Figure 1)
# --------------------------
fit_site_models <- function(df) {
  list(
    ACI       = lm(ACI ~ Site, data = df),
    ADI       = lm(ADI ~ Site, data = df),
    AAc       = lm(AAc ~ Site, data = df),
    msldB_bio = lm(msldB_bio ~ Site, data = df)
  )
}

effects_long <- function(model_list, rain_label = NULL) {
  out <- lapply(names(model_list), function(nm) {
    m <- model_list[[nm]]
    x <- as.data.frame(ggeffects::ggeffect(m))
    x$indice <- nm
    if (!is.null(rain_label)) x$rain <- rain_label
    x
  })
  bind_rows(out)
}

models_all <- fit_site_models(ind_df)
df_indices <- effects_long(models_all)

p_fig1 <- ggplot(df_indices, aes(x = Site.x, y = Site.predicted)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Site.conf.low, ymax = Site.conf.high),
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ indice, scales = "free_y") +
  xlab("") + ylab("Index value") +
  theme_bw() +
  theme(legend.position = "none")

print(p_fig1)

# --------------------------
# 4) Rain classifier training (package example data)
# --------------------------
data("indices_data", package = "ANECO")
data("labeled_data", package = "ANECO")

# (Optional) quick class balance inspection
labeled_data %>%
  group_by(Category) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  print()

model_list <- soundclassifier_fit(
  ref = labeled_data,
  ind = indices_data,
  p = 0.80,
  nmax = "even",
  classes = c("Rain", "Clean")
)

print(model_list$info)
print(model_list$model)

# --------------------------
# 5) Internal test performance (hold-out from labeled_data)
# --------------------------
test <- model_list$testdata
pred <- predict(model_list$model, newdata = test[, -c(1, 2)])

cm_internal <- caret::confusionMatrix(data = pred, reference = test$Category, positive = "Rain")
print(cm_internal)

metrics_from_cm <- function(cm) {
  data.frame(
    Accuracy  = unname(cm$overall["Accuracy"]),
    Precision = unname(cm$byClass["Pos Pred Value"]),
    Recall    = unname(cm$byClass["Sensitivity"]),
    F1        = unname(cm$byClass["F1"])
  )
}

perf_internal <- metrics_from_cm(cm_internal)
print(perf_internal)

# --------------------------
# 6) Predict rain over Amistosa example + external evaluation (optional)
# --------------------------
df2 <- soundclassifier_predict(data = ind_df, mod = model_list$model)

# OPTIONAL: If you have a manual label file for the 180 minutes, set the path here.
# (This keeps your manuscript logic but removes hardcoded absolute paths.)
manual_labels_path <- NULL
# manual_labels_path <- "path/to/Labeled_example.xlsx"

if (!is.null(manual_labels_path)) {
  real <- rio::import(manual_labels_path)
  
  cm_external <- caret::confusionMatrix(
    data = df2$predicted_class,
    reference = as.factor(real$Category)
  )
  print(cm_external)
  
  perf_external <- metrics_from_cm(cm_external)
  print(perf_external)
}

# --------------------------
# 7) Filter out predicted rain minutes and refit models (Figure 2)
# --------------------------
# Safer join: explicit key (avoid accidental multi-column joins)
merged_df <- left_join(
  df2,
  ind_df,
  by = c("Name", "Minute", "Site", "RecorderID", "Date", "Hour")
)

filtered_df <- merged_df %>%
  filter(predicted_class != "Rain")

models_no_rain <- fit_site_models(filtered_df)

df_clean <- effects_long(models_no_rain, rain_label = "Without Rain")
df_rain  <- df_indices %>% mutate(rain = "With Rain")

dd <- bind_rows(df_clean, df_rain)

library(ggsci)

p_fig2 <- ggplot(dd, aes(x = rain, y = Site.predicted, colour = Site.x)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Site.conf.low, ymax = Site.conf.high),
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ indice, scales = "free_y") +
  xlab("") + ylab("Index value") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_color_npg()

print(p_fig2)