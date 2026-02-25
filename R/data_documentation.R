#' Example labeled data for sound classification
#'
#' A dataset containing manually labeled audio segments used
#' to train the random forest classifier in ANECO.
#'
#' @format A data frame with 4 columns:
#' \describe{
#'   \item{File}{Audio file name}
#'   \item{Category}{Label (Clean, Rain, Cicada)}
#'   \item{Start}{Start minute}
#'   \item{End}{End minute}
#' }
#' @source Example dataset included in ANECO
"labeled_data"


#' Example acoustic indices dataset
#'
#' Acoustic indices calculated using acoustic_indices()
#'
#' @format A data frame where rows represent minutes and
#' columns represent acoustic indices and metadata.
#' @source Example dataset included in ANECO
"indices_data"
