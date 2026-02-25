#' Fit a Random Forest classifier to detect rain-dominated minutes
#'
#' Trains a supervised classifier (Random Forest via \pkg{caret}) to identify
#' 1-minute segments dominated by rain (\code{"Rain"}) using (1) a table of labeled
#' time intervals and (2) a table of acoustic indices computed with
#' \code{\link[ANECO]{acoustic_indices}} (with \code{noise.ind = TRUE}).
#'
#' Although this function can be used in a multiclass setting (e.g., including
#' \code{"Cicada"} or other labels), the predictor set and intended use are optimized
#' for rain detection and performance on other noise sources may vary.
#'
#' @param ref A data.frame with labeled events. Must contain at least the columns
#'   \code{File}, \code{Category}, \code{Start}, and \code{End}.
#'   \describe{
#'     \item{File}{Audio file name (with extension). If the extension is missing,
#'       \code{.wav} is appended.}
#'     \item{Category}{Class label for the event. At minimum, should include
#'       \code{"Rain"} and at least one non-rain class (e.g., \code{"Clean"}).}
#'     \item{Start}{Event start time in minutes.}
#'     \item{End}{Event end time in minutes. \strong{Exclusive} upper bound.}
#'   }
#'
#' @param ind A data.frame of acoustic indices produced by
#'   \code{\link[ANECO]{acoustic_indices}} with \code{noise.ind = TRUE}. It must contain
#'   at least \code{Name} (file identifier) and \code{Minute} (minute index), plus the
#'   indices used by this function (see \strong{Details}).
#'
#' @param p Numeric in (0, 1]. Proportion of samples used for training. The remainder
#'   is used as a hold-out test set. Default is \code{0.80}.
#'
#' @param nmax Numeric, \code{"even"}, or \code{NULL}. Controls the maximum number of
#'   samples per class used to fit the model:
#'   \describe{
#'     \item{\code{NULL}}{No truncation; all available samples are used.}
#'     \item{\code{"even"}}{Caps all classes to the size of the least frequent class.}
#'     \item{numeric}{Caps each class to at most \code{nmax} randomly sampled rows.}
#'   }
#'
#' @param classes Character. Either \code{"All"} (default) to use all categories in
#'   \code{ref}, or a character vector of class names to keep (e.g.,
#'   \code{c("Clean","Rain")}).
#'
#' @details
#' \strong{Rain definition.} \code{"Rain"} refers to minutes dominated by rain such
#' that other sounds are substantially masked. Because the boundary between light and
#' heavy rain is gradual, labels should reflect whether rain is sufficiently strong
#' to impair biological signal interpretation; when uncertain, labeling as
#' \code{"Rain"} is often preferable for conservative filtering.
#'
#' \strong{Minute expansion rule (exclusive end).} Each labeled interval is expanded
#' to minute-level labels using the half-open interval \code{[Start, End)} at a
#' 1-minute resolution. For example, \code{Start = 2} and \code{End = 5} will label
#' minutes 2, 3, and 4 (minute 5 is not included).
#'
#' \strong{Predictors.} This function uses a fixed subset of indices designed for
#' noise/rain detection:
#' \code{msldB_bio}, \code{msldB_low}, \code{mdBGL}, \code{SpecIQR}, \code{SpecKurt},
#' \code{SpecSkew}, \code{SpecQ3}, \code{SpecCent}, \code{S2N_chi}, \code{dfreq},
#' \code{E01}, \code{E12}, \code{E23}, \code{E34}, \code{E45}, \code{E56}, \code{E67},
#' \code{E78}, \code{E89}, \code{E910}, \code{E1011}.
#'
#' Model fitting uses \code{\link[caret]{train}} with method \code{"rf"} and
#' preprocessing \code{c("center","scale")}. Cross-validation progress is printed
#' during training (\code{verboseIter = TRUE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{model}{A fitted \pkg{caret} model object (Random Forest).}
#'   \item{info}{A table with the number of training samples per class.}
#'   \item{testdata}{Hold-out dataset not used for training.}
#' }
#'
#' @examples
#' \dontrun{
#' library(caret)
#' set.seed(123)
#'
#' data("indices_data", package = "ANECO")
#' data("labeled_data", package = "ANECO")
#'
#' # Binary rain detector (recommended use)
#' mod <- soundclassifier_fit(ref = labeled_data, ind = indices_data,
#'                            p = 0.8, nmax = "even",
#'                            classes = c("Clean","Rain"))
#'
#' mod$model
#' table(mod$testdata$Category)
#'
#' # Multiclass (user-provided additional labels); performance may vary
#' # mod2 <- soundclassifier_fit(ref = labeled_data, ind = indices_data,
#' #                             classes = c("Clean","Rain","Cicada"))
#' }
#'
#' @export
soundclassifier_fit <- function(ref, ind, p = 0.80, nmax = NULL, classes = "All" ){

  if(!all(c("Start", "End", "File", "Category") %in% names(ref))){stop("The labeled dataset must contain the columns File, Category, Start, End")}

  ref$Start <- as.numeric(ref$Start)
  ref$End <- as.numeric(ref$End)

  audio_format <- ".wav"
  if(any(grepl(pattern = ".flac", ind$Name))){
    audio_format <- ".flac"
  }

  for (f in 1:nrow(ref)) {
    if(!grepl(pattern = audio_format, x = ref$File[f])){
      ref$File[f] <- paste0(ref$File[f], ".wav")
    }
  }

  # Renombrar y categorizar cada minuto de la base de datos de etiquetas
  # ------------------------------------

  if(any(classes != "All")){
    ref <- ref[ref$Category %in% classes, ]
  }

  name <- vector()
  minute <- vector()
  cate <- vector()
  cont <- 1
  for (r in 1:nrow(ref)) {
    dur <- (ref[r, "End"] - ref[r, "Start"])
    if(dur < 1) {
      next }
    min <- seq(from = as.integer(ref[r, "Start"]), to = as.integer(ref[r, "End"]), by = 1)
    min <- min[-length(min)]
    for(i in min){
      minute[cont] <- i
      name[cont] <- ref[r, "File"]
      cate[cont] <- as.character(ref[r, "Category"])
      cont <- cont + 1
    }
  }
  datfinal <- data.frame(Name = name, Minute = minute, Category = cate)

  # Unir matriz de etiquetas con matriz de indices acusticos
  #-----------------------------------------------

  sub <- ind[which(ind$Name %in% datfinal$Name), ]

  mrg <- merge(y = sub, x = datfinal, by = c("Name", "Minute"))

  mrg$Category <- as.factor(mrg$Category)

  tabmrg <- table(mrg$Category)

  if(!is.null(nmax)){
    if(nmax == "even"){
      nmax <- min(tabmrg)
    }
    ntable <- data.frame()
    for (c in unique(mrg$Category)) {
      subsetcategory <- mrg[mrg$Category == c, ]
      maxsize <- nrow(subsetcategory)
      if(nmax < maxsize){
        subsetcategory <- subsetcategory[sample(1:nrow(subsetcategory), size = nmax), ]
      }
      ntable <- rbind(ntable, subsetcategory)
    }
  mrg <- ntable
  }

  setdatos <- mrg %>% dplyr::select("Name", "Category", "msldB_bio", "msldB_low", "mdBGL", "SpecIQR", "SpecKurt", "SpecSkew", "SpecQ3","SpecCent", "S2N_chi", "dfreq", "E01", "E12", "E23", "E34", "E45", "E56", "E67", "E78", "E89", "E910", "E1011")

  setdatos$Category <- as.factor(setdatos$Category)
  setdatos <- setdatos %>%
    dplyr::mutate(dplyr::across(msldB_bio:E1011, as.numeric))

  # Control del modelo
  setControl <- caret::trainControl(method = "cv", verboseIter = TRUE)

  # Dividir el set de datos en datos para entrenamiento y datos para prueba
  inTrain = caret::createDataPartition(y = setdatos$Category, p = p, list = FALSE)
  training = setdatos[inTrain,]
  testing = setdatos[-inTrain,]

  mod_rf <- caret::train(Category ~ . , data = training[,-1], method = "rf", preProcess = c("center", "scale"), trControl = setControl)

  trainInfo <- table(training$Category)

  return(list(model = mod_rf, info = trainInfo, testdata = testing))
}

