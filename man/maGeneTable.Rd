\name{maGeneTable}
\alias{maGeneTable}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Table of spot coordinates and gene names}
\description{
This function produces a table of spot coordinates and gene names for
objects of class \code{"\link{marrayRaw}"} and \code{"\link{marrayNorm}"}}.
}
\usage{
maGeneTable(object)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{microarray object of class \code{"\link{marrayRaw}"} and
    \code{"\link{marrayNorm}"}.}
}

\value{
  an object of class \code{\link{data.frame}}, with rows corresponding
  to spotted probe sequences. The first four columns are the grid matrix
  and spot matrix coordinates, and the remaining columns are the spot
  descriptions stored in the \code{maGnames} slot of the microarray
  object. 
}

\author{Yee Hwa (Jean) Yang}


\seealso{\code{\link{marrayInfo}}, \code{\link{marrayLayout}}, \code{\link{marrayRaw}}, \code{\link{marrayNorm}},  \code{\link{maCompCoord}}.}


\examples{
# Example uses swirl dataset, for description type ? swirl

data(swirl)

tab<-maGeneTable(swirl)
tab[1:10,]

}

\keyword{methods}% at least one, from doc/KEYWORDS

