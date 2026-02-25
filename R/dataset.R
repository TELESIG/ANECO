#' Example of labeled data.
#'
#' A dataset indicating the periods with occurrences of rain and cicada sounds.
#'
#' @format A dataframe with the columns File, Category, Start, and End. Where:
#' \describe{
#'   \item{File}{Full name of the audio file, including the .wav extension. For example: FT02_20190615_110000.wav}
#'   \item{Category}{Category of the audio segment. Three categories are available: Clean, Rain, Cicada}
#'   \item{Start}{Start time of the event in the analyzed category, with a resolution of 0.5 seconds}
#'   \item{End}{End time of the event in the analyzed category, with a resolution of 0.5 seconds}
#' }
"labeled_data"


#' Example dataset of acoustic indices
#'
#' A dataset of acoustic indices resulting from the function \code{\link[ANECO]{acoustic_indices}}, with the parameter \code{"noise.ind = TRUE"}
#'
"indices_data"
