\name{maCompInd}
\alias{maCompInd}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generate spot indices}
\description{
This function generates spot indices from ranges of rows and columns for
the grid and spot matrices. Spots on the array are numbered consecutively starting
from the top left grid and the top left spot within each grid.

}
\usage{
maCompInd(grows, gcols, srows, scols, L)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{grows}{numeric vector of grid rows.}
  \item{gcols}{numeric vector of grid columns.}
  \item{srows}{numeric vector of spot rows.}
  \item{scols}{numeric vector of spot columns.}
  \item{L}{object of class \code{"\link{marrayLayout}"}.}
}

\value{
a numeric vector of spot indices.
}

\author{Yee Hwa (Jean) Yang, Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayLayout}}, \code{\link{maCoord2Ind}},
\code{\link{maInd2Coord}}, \code{\link{maCompCoord}}.}


\examples{
L <- new("marrayLayout", maNgr=4, maNgc=4, maNsr=22, maNsc=24)
maCompInd(1:2,1,1:4,1:3,L)


}
\keyword{manip}
\keyword{methods}
