\name{maCompCoord}
\alias{maCompCoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generate grid and spot matrix coordinates}
\description{
This function generates grid and spot matrix coordinates from ranges of
rows and columns for the grid and spot matrices. Spots on the array are numbered consecutively starting
from the top left grid and the top left spot within each grid.
}
\usage{
maCompCoord(grows, gcols, srows, scols)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{grows}{numeric vector of grid rows.}
  \item{gcols}{numeric vector of grid columns.}
  \item{srows}{numeric vector of spot rows.}
  \item{scols}{numeric vector of spot columns.}
}

\value{

a matrix of spot four-coordinates, with rows corresponding to
  spots and columns to grid row, grid column, spot row, and spot column coordinates.


}

\author{Yee Hwa (Jean) Yang, Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayLayout}}, \code{\link{maCoord2Ind}},
\code{\link{maInd2Coord}}, \code{\link{maCompInd}}.}


\examples{

maCompCoord(1:2,1,1:4,1:3)


}
\keyword{manip}
\keyword{methods}
