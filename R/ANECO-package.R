#'ANECO: an R package for easy calculation of acoustic indices and identification of heavy rain and cicadas sounds
#'
#' The ANECO package is a tool developed by the TELESIG laboratory from the Instituto Internacional en Conservación y Manejo de Vida Silvestre, Universidad Nacional de Costa Rica. It is designed to facilitate ecoacoustic analyses with minimal programming expertise. ANECO integrates core functionalities for calculating acoustic indices, calibrating sound pressure levels, and identifying audio segments with heavy rain. The package provides three primary functions, enabling users to compute multiple indices, build classification models, and apply these models to new data sets seamlessly.
#'Its main function are:
#'\itemize{
#'\item{Simplification of the calculation of a set of acoustic indices on audio files.}{}
#'\item{Obtaining calibrated acoustic indices with absolute sound pressure level values}{}
#'\item{Provides measurements of dominance and energy in frequency bands at 1 kHz resolution}{}
#'\item{Automatic detection and classification of audio minutes saturated with heavy rain noise}{}}

#'The ANECO package uses the calibration process described in Merchant et al. (2015), which provides absolute sound pressure level measurements. The code for the calibrated indices is based on the work of Buxton et al. (2018), while the traditional (uncalibrated) indices are based on functions from the seewave (Sueur et al., 2016), soundecology (Villanueva-Rivera & Pijanowski, 2016) and Sinax (https://github.com/osoramirez/Sinax4) packages.
#'
#'@section Functions:
#'\itemize{
#'\item{\code{\link[ANECO]{acoustic_indices}}}
#'\item{\code{\link[ANECO]{soundclassifier_predict}}}
#'\item{\code{\link[ANECO]{soundclassifier_model}}}
#'}
#'
#' @name ANECO
#' @references \itemize{
#' \item{Buxton, R., McKenna, M. F., Clapp, M., Meyer, E., Stabenau, E., Angeloni, L. M., ... & Wittemyer, G. 2018. Efficacy of extracting indices from large‐scale acoustic recordings to monitor biodiversity. Conservation Biology.}{}
#'
#' \item{Merchant, N. D., Fristrup, K. M., Johnson, M. P., Tyack, P. L., Witt, M. J., Blondel, P., & Parks, S. E. (2015). Measuring acoustic habitats. Methods in Ecology and Evolution, 6(3), 257-265.}{}
#'
#' \item{Sueur, J., Simonis, C., Brown, E., Depraetere, M., Desjonqueres, C., Fabianek, F., Gasc, A., LaZerte, S., Lees, J., Marchal, J., Pavoine, S., Stotz, A., Villanueva-Rivera, L., Ross, Z., Witthoft, C., & Zhivomirov, H. (2016). Paquete para R: Seewave, Sound Analysis and Synthesis. Version 2.0.5}{}
#'
#' \item{Villanueva-Rivera, L. J. & Pijanowski, B. C. 2016. soundecology: Soundscape Ecology. R package version 1.3.2.}{} }
NULL
#> NULL
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr
#' @import seewave
#' @importFrom caret confusionMatrix
#' @importFrom caret createDataPartition
#' @importFrom caret train
#' @importFrom caret trainControl
#' @importFrom dplyr slice_sample
#' @importFrom future multicore
#' @importFrom future multisession
#' @importFrom future plan
#' @importFrom future sequential
#' @importFrom future.apply future_lapply
#' @importFrom ineq Gini
#' @importFrom moments kurtosis
#' @importFrom moments skewness
#' @importFrom NbClust NbClust
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @importFrom progressr handlers
#' @importFrom progressr make_progression_handler
#' @importFrom progressr progressor
#' @importFrom randomForest randomForest
#' @importFrom rio import
#' @importFrom soundecology acoustic_complexity
#' @importFrom soundecology acoustic_diversity
#' @importFrom soundecology acoustic_evenness
#' @importFrom soundecology bioacoustic_index
#' @importFrom soundecology ndsi
#' @importFrom stats median
#' @importFrom stats mvfft
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats var
#' @importFrom tuneR readWave
#' @importFrom tuneR writeWave
#' @importFrom utils capture.output
#' @importFrom vegan diversity
## usethis namespace: end
NULL
