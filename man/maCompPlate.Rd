\name{maCompPlate}
\alias{maCompPlate}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generate plate IDs}
\description{
This function generates plate IDs from the dimensions of the grid and
spot matrices. Note that this function only applies to arrays with a
regular plate layout, where the number of spots is a multiple of the
number of wells on a plate (usually 96 or 384) and each well contributes 
exactly one spot. It should thus be used with caution.}

\usage{
maCompPlate(x, n=384)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{"\link{marrayLayout}"},
    \code{"\link{marrayRaw}"} and \code{"\link{marrayNorm}"}}
  \item{n}{object of class "numeric", number of wells in each plate,
    usually 384 or 96.} 
}

\details{
  Having plate IDs may be useful for the purpose of
  normalization. Normalization by plate can be done using the function
  \code{\link{maNormMain}}. 
}
\value{
a vector of plate IDs (\code{\link{factor}}).
}

\author{Yee Hwa (Jean) Yang, Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayLayout}}, \code{\link{marrayRaw}},
  \code{\link{marrayNorm}}}

\examples{

L<-new("marrayLayout", maNgr=4, maNgc=4, maNsr=22, maNsc=24)
plate<-maCompPlate(L,384)
table(plate)
maPlate(L)<-plate

}
\keyword{methods}% at least one, from doc/KEYWORDS

