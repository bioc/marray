\name{maQualityMain}
\alias{maQualityMain}
\alias{maQuality}

\title{Calculate basic array wide statistics}

\description{
  This is based on GenePix output and it calculates summary statistics
  of various signal to noise ratios.
}
\usage{
maQualityMain(mraw, path=".", fname, output = FALSE)
maQuality(mraw, path=".")
}

\arguments{
  \item{mraw}{Microarray object of class `"marrayRaw"'.
    By default, plots are done for the first array in the batch.}
  \item{path}{a character string representing the data directory. By
    default this is set to the current working directory (".").
    In the case where "fname" contains the full path name, path
    should be set to NULL. }
  \item{fname}{a vector of character strings representing the name of
    the output file. This is only used when "output" is TRUE.}
  \item{output}{if TRUE, write the results to a file.}
}

\value{
  Produce a lists of various summary statistics
}

\author{Jean Yee Hwa yang}

\seealso{\code{\link{summary}}, \code{\link{maGenControls}}}

\examples{
##   data(swirl)
##   datadir <- system.file("data", package="marray")
##   maControls(swirl) <- maGenControls(maGnames(swirl), id="Name")
##   maQualityMain(swirl[,1], path=datadir)
}

\keyword{manip}% at least one, from doc/KEYWORDS

