\name{maDotPlots}
\alias{maDotPlots}

\title{Produce a dot plot for replicated controls genes}
\description{
  A dot plot showing the values of replicated control genes
}
\usage{
maDotPlots(data, x = list("maA"), id = "ID", pch, col, nrep = 3, ...)
}

\arguments{
  \item{data}{Microarray object of class `"marrayRaw"'. By default, plots
    are done for the first array in the batch.}
  \item{x}{Name of accessor method for the spot statistic to be display.
    E.g. maM to display log-ratios and maA to display log-intensities}
  \item{id}{Integer value specifying the index of the array.}
  \item{pch}{point style for different control samples.}
  \item{col}{color code for different control samples.}
  \item{nrep}{Integer value indicating the minimum number of replicated
    spots required.}
  \item{\dots}{Optional graphical parameters, see `par'.}
}

\value{
  Produce a dot plot showing the values of replicated control genes.
}

\examples{
  data(swirl)
  maControls(swirl) <- maGenControls(maGnames(swirl), id="Name")
  maDotPlots(swirl[,1])
}

\author{Yee Hwa (Jean) Yang}

\keyword{hplot}
