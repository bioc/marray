\name{maCoord2Ind}
\alias{maCoord2Ind}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert grid and spot matrix coordinates to spot indices}
\description{
This functions converts grid and spot matrix coordinates (four
coordinates) to spot
indices, where spots on the array are numbered consecutively starting
from the top left grid and the top left spot within each grid.
}
\usage{
maCoord2Ind(x, L)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a matrix of spot four-coordinates, with rows corresponding to
  spots and columns to grid row, grid column, spot row, and spot column coordinates.}
  \item{L}{an object of class \code{"\link{marrayLayout}"}.}
}

\value{
 a numeric vector of spot indices.
}

\author{Yee Hwa (Jean) Yang, Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayLayout}}, \code{\link{maInd2Coord}},
\code{\link{maCompCoord}}, \code{\link{maCompInd}}.}

\examples{
L <- new("marrayLayout", maNgr=4, maNgc=4, maNsr=22, maNsc=24)
coord<-cbind(rep(2,4),rep(1,4),rep(1,4),1:4)
maCoord2Ind(coord, L)


}
\keyword{manip}
\keyword{methods}
