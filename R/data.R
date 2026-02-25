#' Catalog of acoustic indices available in ANECO
#'
#' A lookup table describing all acoustic indices implemented in ANECO.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{code}{Index code used in `sel.ind`.}
#'   \item{label}{Human-readable name.}
#'   \item{description}{Short explanation of what the index measures.}
#'   \item{reference}{Primary bibliographic reference (short form).}
#'   \item{family}{Functional family of the index (e.g., level_spl, entropy, complexity).}
#'   \item{calibrated}{Logical. TRUE if computed from calibrated SPL spectra.}
#'   \item{source}{Implementation origin (ANECO / seewave / soundecology).}
#'   \item{implementation}{Short note describing how it is computed.}
#' }
#'
#' @usage data(indices_catalog)
"indices_catalog"


#' Index codes available in ANECO
#'
#' Character vector of valid index codes that can be passed to
#' `acoustic_indices()` through the `sel.ind` argument.
#'
#' @format Character vector.
#' @usage data(indices_codes)
"indices_codes"


#' List available acoustic indices
#'
#' Returns the index catalog.
#'
#' @export
aneco_list_indices <- function() {
  indices_catalog
}
