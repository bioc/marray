\name{marrayTools}
\alias{spotTools}
\alias{gpTools}

\title{Wrapper function}
\description{
  This functions will produce a directory of diagnostic plots, quality
  measures and normalization.
}
\usage{
spotTools(fnames, path = ".", galfile, bg = TRUE, plot = TRUE,
                 quality = TRUE, fill = TURE, raw = FALSE, echo =
                 TRUE, ...)

gpTools(fnames, path = ".", galfile, bg = TRUE, plot = TRUE, quality = TRUE, fill = TRUE, raw = FALSE, echo = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
    \item{fnames}{a vector of character strings containing the file names of
    each spot quantification data file.  These typically end in \code{.spot}
    for the software \code{Spot} or \code{.gpr} for the software \code{GenePix}.}
  \item{path}{a character string representing the data directory. By
    default this is set to the current working directory (".").  In the
    case where \code{fnames} contains the full path name, path should be set
    to NULL.}
  \item{galfile}{a character string representing the Gal file.}
  \item{bg}{logical: if 'TRUE', background correction is performed before normalization}
  \item{plot}{logical: if 'TRUE', maDiagnPlots is call to generate basic
    diagnostic plot per array.}
  \item{quality}{logical: if 'TRUE', maQuality is call to calculate basic quality measures.}
  \item{fill}{logical: if `TRUE', `scan' will implicitly add empty fields
    to any lines with fewer fields than implied by `what'.}
  \item{raw}{logical: if 'TRUE', the microarray object class of
    'marrayRaw' is return rather than 'marrayNorm'.}
  \item{echo}{logical: if 'TRUE' (default), run-time
    information of the array currently being process is printed.}
  \item{\dots}{further arguments to \code{\link{scan}} or \code{\link{maNorm}}.}
}

\value{
  Jpeg files of diagnostic plots, excel files of quality measures, tab
  delimited files of normalized log-ratios and log-intensities values are
  saved to disk.
  In addition, an object of class \code{\link{marrayRaw}} or \code{\link{marrayNorm}} is return.
}

\author{Yee Hwa (Jean) Yang}

\seealso{\code{\link{read.marrayRaw}}}

\examples{
 data(swirl)
 datadir <- system.file("data", package="marray")

## normdata <- spotTools(path=datadir, quality=FALSE)
}

\keyword{file}% at least one, from doc/KEYWORDS

