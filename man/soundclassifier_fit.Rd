% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/soundclassifier_fit.R
\name{soundclassifier_fit}
\alias{soundclassifier_fit}
\title{Fitting a model to predict heavy rain and cicada sounds}
\usage{
soundclassifier_fit(ref, ind, p = 0.8, nmax = NULL)
}
\arguments{
\item{ref}{Dataset of labeled acoustic data. It must include the columns File, Category, Start and End.
\describe{
  \item{File}{Full name of the audio file including the extension .wav at the end. Example: FT02_20190615_110000.wav}
  \item{Category}{Category in which the recorded audio segment is classified. Three categories are used: Clean, Rain and Cicada. See section \strong{Details}}
  \item{Start}{Start minute of the event of the analyzed category, in a resolution of 0.5 seconds}
  \item{End}{Minute of completion of the event of the analyzed category, in a resolution of 0.5 seconds}
}}

\item{ind}{An acoustic index dataset resulting from the function \code{\link[ANECO]{acoustic_indices}} (The \code{noise.ind = TRUE} option must have been included in the function \code{\link[ANECO]{acoustic_indices}})}

\item{p}{Specifies the proportion of data to be used in training the model. The default value is 0.80, equivalent to 80 percent of the data.}

\item{nmax}{Numeric. This parameter is used to prevent bias due to an uneven number of elements across categories. For example, if \code{nmax = 1000} is specified, no category will have more than 1000 elements for training. If \code{nmax = "even"} is chosen, the maximum number of elements used in training is set to the number of elements in the least represented category, ensuring that all categories have the same number of elements. If \code{nmax = NULL}, no truncation is applied, and all elements will be used for training. The default is \code{nmax = NULL}.}
}
\value{
Returns a list with three elements. The first element corresponds to the model, the second element is the data that was not used in the construction of the model and can be used for its evaluation. Finally, the third element corresponds to the information on the number of elements of each category used to fit the model.
}
\description{
Fit a model to predict rain and cicadas noises based on a dataset of labeled audios used to train the model.
}
\details{
The categories that the model allows to classify are Clean, Rain and Cicada.
\describe{
  \item{Clean}{Refers to audio segments that are not obscured by unwanted noise. In these segments it is possible to distinguish animal sounds (if present) such as bird calls, amphibians or insects}
  \item{Rain}{Rain sound that masks other sounds. The difference between heavy rain and light rain (which is considered clean) is gradual, so it must be made at discretion considering whether the noise level is sufficient to drown out other sounds such as bird songs. If in doubt, it is preferable to annotate the audio as heavy rain}
  \item{Cicada}{Loud sounds made by cicadas or insects that can generate noise. They generally occupy a small frequency range}}
}
\examples{
library(caret)
# Use the example data
data("indices_data", package = "ANECO")
data("labeled_data", package = "ANECO")

mod <- soundclassifier_fit(ref = labeled_data, ind = indices_data, p = 0.8, nmax = "even")

mod$model

}
