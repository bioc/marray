\name{maTwoSamples}
\alias{maTwoSamples}
\alias{widget.TwoSamples}

\title{Changing signs for two sample analysis}
\description{
  Taking target file information and flip the dye swaps experiments.
}
\usage{
maTwoSamples(targetfile, normdata, Trt, Ctl, targetID = "TargetName", slidesID = "Slides", dyesID = "Dyes", RedID = 5, path = ".", output = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{targetfile}{A data.frame containing target samples information.}
  \item{normdata}{A R object of class 'marrayNorm'}
  \item{Trt}{A character string representing "treatment" sample.}
  \item{Ctl}{A character string representing "controls" sample.}
  \item{targetID}{A character string representing the column name
    in 'targetfile' containing target samples information.}
  \item{slidesID}{A character string representing the column name
    in 'targetfile' containing the slide label.}
  \item{dyesID}{A character string representing the column name
    in 'targetfile' containing dye labeled information.}
  \item{RedID}{The character use to represent the Cy5 dye.}
  \item{path}{A character string representing the data directory. By
    default this is set to the current working directory ("."). }
  \item{output}{Save and tab delimited file}
}

\value{
  An objects of 'marrayNorm' with the dye assignment adjusted.
}

\author{Yee Hwa (Jean) Yang}

\keyword{manip}% at least one, from doc/KEYWORDS
