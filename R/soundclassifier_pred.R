#' Predict rain-dominated minutes from acoustic indices
#'
#' Uses a previously fitted classifier (typically a Random Forest model trained with \code{\link[ANECO]{soundclassifier_fit}}) to predict class probabilities and the final class label for each 1-minute segment in an acoustic index dataset.
#'
#' The workflow is optimized for rain detection (i.e., identifying minutes dominated by \code{"Rain"}). However, if the fitted model includes additional classes (e.g., user-defined noise sources), probabilities are returned for all classes learned by the model.
#'
#' @param data A data.frame produced by \code{\link[ANECO]{acoustic_indices}} with \code{noise.ind = TRUE}. It must include:
#'   \itemize{
#'     \item identification columns: \code{Name}, \code{Code}, \code{Site},
#'       \code{RecorderID}, \code{Date}, \code{Year}, \code{Month}, \code{Day}, \code{Hour}, \code{Minute}
#'     \item predictor columns required by the classifier (see \strong{Details})
#'   }
#'
#' @param mod A fitted classification model compatible with
#'   \code{predict(mod, newdata, type = "prob")} and
#'   \code{predict(mod, newdata, type = "raw")}. In the ANECO workflow this is the \code{$model} element returned by
#'   \code{\link[ANECO]{soundclassifier_fit}} (a \pkg{caret} \code{train} object).
#'
#' @details
#' \strong{Predictors.} The function uses a fixed subset of noise-related indices:
#' \code{msldB_bio}, \code{msldB_low}, \code{mdBGL}, \code{SpecIQR}, \code{SpecKurt},
#' \code{SpecSkew}, \code{SpecQ3}, \code{SpecCent}, \code{S2N_chi}, \code{dfreq},
#' \code{E01}, \code{E12}, \code{E23}, \code{E34}, \code{E45}, \code{E56}, \code{E67},
#' \code{E78}, \code{E89}, \code{E910}, \code{E1011}.
#'
#' \strong{Output probabilities.} Class probabilities are returned in columns prefixed
#' with \code{prob_} (e.g., \code{prob_Rain}). The number and names of probability columns depend on the classes present in the fitted model.
#'
#' If columns named \code{prob_*} or \code{predicted_class} already exist in \code{data},
#' they will be overwritten in the returned object.
#'
#' @return A data.frame with one row per minute, containing:
#' \itemize{
#'   \item the identification columns from \code{data}
#'   \item one probability column per class, prefixed with \code{prob_}
#'   \item \code{predicted_class}: the predicted class label (character)
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' data("indices_data", package = "ANECO")
#' data("labeled_data", package = "ANECO")
#'
#' # Train a binary rain detector (recommended use)
#' fit <- soundclassifier_fit(ref = labeled_data,
#'                            ind = indices_data,
#'                            classes = c("Clean","Rain"),
#'                            nmax = "even")
#'
#' # Predict on a full index dataset
#' pred <- soundclassifier_predict(data = indices_data,
#'                                  mod = fit$model)
#'
#' head(pred)
#' }
#'
#' @export
soundclassifier_predict <- function(data = NULL, mod = NULL){

  if (is.null(mod)) {
    stop("Argument 'mod' must be a fitted classification model (e.g., fit$model from soundclassifier_fit).")
  }

  # check columns
  indnames <- c(
    "msldB_bio", "msldB_low", "mdBGL", "SpecIQR", "SpecKurt", "SpecSkew", "SpecQ3", "SpecCent", "S2N_chi", "dfreq", "E01", "E12", "E23", "E34", "E45", "E56", "E67", "E78", "E89", "E910", "E1011"
  )

  if(!all(indnames %in% colnames(data))){
    stop("Incomplete database. Must contain acoustic indices for noise detection (noise.ind = TRUE) in the acoustic_indices function: 'msldB_bio', 'msldB_low', 'mdBGL', 'SpecIQR', 'SpecKurt', 'SpecSkew', 'SpecQ3','SpecCent', 'S2N_chi', 'dfreq', 'E01', 'E12', 'E23', 'E34', 'E45', 'E56', 'E67', 'E78', 'E89', 'E910', 'E1011'")
  }

  clmnsID <- c("Name", "Code", "Site", "RecorderID", "Date", "Year","Month", "Day", "Hour", "Minute")

  if(!all(clmnsID %in% colnames(data))){
    stop("Incomplete database. It must contain the columns for identification that are included in the result of the acoustic_indices function: 'Name', 'Code', 'Site', 'RecorderID', 'Date', 'Year', 'Month', 'Day', 'Hour', 'Minute' ")
  }

  # Select acoustic indices
  datos <- data[, which(colnames(data) %in% indnames), drop = FALSE]

  # select ID columns
  id <- data[, which(colnames(data) %in% clmnsID), drop = FALSE]

  # predict class probabilities
  rf_prob <- stats::predict(mod, newdata = datos, type = "prob")

  rf_raw <-  stats::predict(mod, newdata = datos, type = "raw")

  colnames(rf_prob) <- paste0("prob_", colnames(rf_prob))
  # save results
  result <- cbind(id, rf_prob)
  result$predicted_class <- rf_raw
  rownames(result) <- NULL
  return(result)
}
