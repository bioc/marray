\name{maDiagnPlots}
\alias{maDiagnPlots}

\title{Basic diagnostic plots}
\description{
  These three functions produce sets of basic diagnostic plots for pre-
  and post-normalization cDNA microarray data. 
}
\usage{
maDiagnPlots(mraw, mNorm = NULL, save = TRUE, fname = NULL,
dev= "png", pch, col, DEBUG = FALSE, ...)
}

\arguments{
  \item{mraw}{Microarray object of class `"marrayRaw"'.
    By default, plots are done for the first array in the batch.}
  \item{mNorm}{Microarray object of class `"marrayNorm"'. By default,
    this value is NULL and the function maNorm is call to perform normalization.}
  \item{save}{If `TRUE', the figures will be saved to the file named in
          `fname'.}
  \item{fname}{A "character" string naming the output file.}
  \item{dev}{A "character" string naming the graphics device.}
  \item{pch}{point style for different control samples.}
  \item{col}{color code for different control samples.}
  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  \item{\dots}{Optional graphical parameters, see `par'.}
}
\details{
  \code{maDiagnPlots} produces plots of pre- and  post-normalization cDNA
  microarray data: color images of log-ratios (M) and log-intensities
  (A) pre- and post-normalization. MA-plots of pre- and post-normalization
  log-ratios M. Log foreground to background ratios for both Cy5 and Cy3
  channels. Dot plots for replicate control spots. 
}


\author{Yee Hwa (Jean) Yang}

\note{This function may not work properly when "save=FALSE".}

\seealso{\code{\link{maPlot}}, \code{\link{maImage}}, \code{\link{maBoxplot}}}


\examples{
data(swirl)
## maDiagnPlots(swirl, dev="postscript", save=TRUE)
}

\keyword{hplot}

