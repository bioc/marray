\name{maCompLayout}
\alias{maCompLayout}

\title{Generate a marrayLayout object}
\description{
  Take a matrix of cooordiates and generate a marrayLayout object.
}
\usage{
maCompLayout(mat, ncolumns = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{a matrix of coordinates, this can either be n by 3 matrix
    with columns (Block, Row, Column) or n by 4 matrix with columns
    (Grid.R, Grid.C, Spot.R, Spot.C)}
  \item{ncolumns}{For n by 3 matrix, the number of meta-grid columns.
    By default, it is set to 4.}
}
\details{

}
\value{
  An object of class \code{"\link{marrayLayout}"}.
}

\author{Jean Yang}

\examples{
X <- cbind(Block = c(1,1,2,2,3,3,4,4), Rows=c(1,2,1,2,1,2,1,2), Columns=rep(1,8))
maCompLayout(X, ncolumns=2)
}
\keyword{manip}% at least one, from doc/KEYWORDS
