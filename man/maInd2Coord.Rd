\name{maInd2Coord}
\alias{maInd2Coord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert spot indices to grid and spot matrix coordinates}
\description{
This functions converts spot indices to grid and spot matrix coordinates (four
coordinates), where spots on the array are numbered consecutively starting
from the top left grid and the top left spot within each grid.
}
\usage{
maInd2Coord(x, L)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of spot indices.}
  \item{L}{an object of class \code{"\link{marrayLayout}"}.}
}
\value{
a matrix of spot four-coordinates, with rows corresponding to
  spots and columns to grid row, grid column, spot row, and spot column coordinates.
}

\author{Yee Hwa (Jean) Yang, Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayLayout}}, \code{\link{maCoord2Ind}},
\code{\link{maCompCoord}}, \code{\link{maCompInd}}.}

\examples{

L <- new("marrayLayout", maNgr=4, maNgc=4, maNsr=22, maNsc=24)
maInd2Coord(c(1:10,529:538), L)


}
\keyword{manip}
\keyword{methods}
