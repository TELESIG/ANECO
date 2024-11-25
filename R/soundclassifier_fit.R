#' Fitting a model to predict heavy rain and cicada sounds
#'
#'Fit a model to predict rain and cicadas noises based on a dataset of labeled audios used to train the model.
#' @param ref Dataset of labeled acoustic data. It must include the columns File, Category, Start and End.
#' \describe{
#'   \item{File}{Full name of the audio file including the extension .wav at the end. Example: FT02_20190615_110000.wav}
#'   \item{Category}{Category in which the recorded audio segment is classified. Three categories are used: Clean, Rain and Cicada. See section \strong{Details}}
#'   \item{Start}{Start minute of the event of the analyzed category, in a resolution of 0.5 seconds}
#'   \item{End}{Minute of completion of the event of the analyzed category, in a resolution of 0.5 seconds}
#' }
#' @param ind An acoustic index dataset resulting from the function \code{\link[ANECO]{acoustic_indices}} (The \code{noise.ind = TRUE} option must have been included in the function \code{\link[ANECO]{acoustic_indices}})
#' @param p Specifies the proportion of data to be used in training the model. The default value is 0.80, equivalent to 80 percent of the data.
#' @param nmax Numeric. This parameter is used to prevent bias due to an uneven number of elements across categories. For example, if \code{nmax = 1000} is specified, no category will have more than 1000 elements for training. If \code{nmax = "even"} is chosen, the maximum number of elements used in training is set to the number of elements in the least represented category, ensuring that all categories have the same number of elements. If \code{nmax = NULL}, no truncation is applied, and all elements will be used for training. The default is \code{nmax = NULL}.
#' @details The categories that the model allows to classify are Clean, Rain and Cicada.
#' \describe{
#'   \item{Clean}{Refers to audio segments that are not obscured by unwanted noise. In these segments it is possible to distinguish animal sounds (if present) such as bird calls, amphibians or insects}
#'   \item{Rain}{Rain sound that masks other sounds. The difference between heavy rain and light rain (which is considered clean) is gradual, so it must be made at discretion considering whether the noise level is sufficient to drown out other sounds such as bird songs. If in doubt, it is preferable to annotate the audio as heavy rain}
#'   \item{Cicada}{Loud sounds made by cicadas or insects that can generate noise. They generally occupy a small frequency range}}
#'
#' @return Returns a list with three elements. The first element corresponds to the model, the second element is the data that was not used in the construction of the model and can be used for its evaluation. Finally, the third element corresponds to the information on the number of elements of each category used to fit the model.
#'
#' @examples
#'\code{
#' library(caret)
#' # Use the example data
#'head(ejemplo_etiquetas)
#'head(ejemplo_indices)
#' # Para mayor informacion sobre los datos de ejemplo usar:
#' # ?ejemplo_etiquetas
#' # ?ejemplo_indices
#'
#'mod <- soundclassifier_model(ref = ejemplo_etiquetas, ind = ejemplo_indices,  = 2.5, p = 0.9)
#'mod$modelo
#'
#' # Evaluar con los datos incluidos en mod$testdata
#'pred <- predict(mod$modelo, newdata = mod$testdata, type = "raw")
#'confusionMatrix(pred, mod$testdata$Category)
#'}
#' @export
#'
soundclassifier_model <- function(ref, ind, p = 0.80, nmax = NULL){

  if(!all(names(ref) %in% c("Start", "End", "File", "Category"))){stop("The labeled dataset must contain the columns File, Category, Start, End")}

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
  minu <- vector()
  cate <- vector()
  cont <- 1
  for (r in 1:nrow(ref)) {
    dur <- (ref[r, "End"] - ref[r, "Start"])
    if(dur < 1) {
      next }
    min <- seq(from = as.integer(ref[r, "Start"]), to = as.integer(ref[r, "End"]), by = 1)
    min <- min[-length(min)]
    for(i in min){
      minu[cont] <- paste0(ref[r, "File"], i)
      cate[cont] <- as.character(ref[r, "Category"])
      cont <- cont + 1
    }
  }
  datfinal <- data.frame(Name = minu, Category = cate)

  # Unir matriz de etiquetas con matriz de indices acusticos
  #-----------------------------------------------

  sub <- ind[which(ind$Name %in% datfinal$Name), ]

  mrg <- merge(y = sub, x = datfinal, by = "Name")

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

  setdatos <- mrg %>% dplyr::select("Name", "Category", "msldBA_bio", "msldBA_low", "mdBGL", "SpecIQR", "SpecKurt", "SpecSkew", "SpecQ3","SpecCent", "S2N_chi", "dfreq", "E01", "E12", "E23", "E34", "E45", "E56", "E67", "E78", "E89", "E910", "E1011")

  setdatos$Category <- as.factor(setdatos$Category)
  setdatos <- setdatos %>%
    dplyr::mutate(across(msldBA_bio:E1011, as.numeric))

  # Control del modelo
  setControl <- caret::trainControl(method = "cv", verboseIter = TRUE)

  # Dividir el set de datos en datos para entrenamiento y datos para prueba
  inTrain = caret::createDataPartition(y = setdatos$Category, p = p, list = FALSE)
  training = setdatos[inTrain,]
  testing = setdatos[-inTrain,]

  mod_rf <- caret::train(Category ~ . , data = training[,-1], method = "rf", preProcess = c("center", "scale"), trControl = setControl)

  trainInfo <- table(training$Category)

  return(list(modelo = mod_rf, info = trainInfo, testdata = testing))
}

