\name{cbind}
\alias{cbind.marrayRaw}
\alias{cbind.marrayNorm}
\alias{rbind.marrayInfo}

\title{Combine marrayRaw,  marrayNorm or marrayInfo Objects}
\description{
  Combine a series of \code{marrayRaw}, \code{marrayNorm} and
  \code{marrayInfo} objects.
}

\usage{
\method{cbind}{marrayRaw}(\dots, deparse.level=1)
\method{cbind}{marrayNorm}(\dots, deparse.level=1)
\method{rbind}{marrayInfo}(\dots, deparse.level=1)
}
\arguments{
  \item{\dots}{\code{marrayRaw} objects or \code{marrayNorm} objects}
  \item{deparse.level}{not currently used, see \code{\link[base]{cbind}} in the base package}
}
\details{
\code{cbind} combines data objects assuming the same gene lists but
different arrays. 
\code{rbind} combines data objects assuming equivalent arrays, i.e., the same RNA targets, but different genes.

For \code{cbind}, the matrices o f expression data from the individual objects are cbinded.
The data.frames of target information, if they exist, are rbinded.
The combined data object will preserve any additional components or attributes found in the first object to be combined.
For \code{rbind}, the matrices of expression data are rbinded while the target information, in any, is unchanged.
}
\value{
}
\author{Jean Yang}
\seealso{
  \code{\link[base]{cbind}} in the base package.
}
\examples{
}
\keyword{manip}
